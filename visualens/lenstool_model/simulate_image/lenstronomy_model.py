import numpy as np
import math
import os
import matplotlib.pyplot as plt
import pickle
import copy
import corner

from lenstronomy.Plots.model_plot import ModelPlot
from lenstronomy.Plots import chain_plot
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.Data.pixel_grid import PixelGrid  
from lenstronomy.Util import param_util
        
from ...utils.utils_astro.utils_general import world_to_relative, relative_to_world, arcsec_to_kpc
from ...utils.utils_astro.get_cosmology import get_cosmo
from ...utils.utils_Qt.drag_widgets import (
    SelectableCircleROI,
    SelectableEllipseROI,
    attach_sliders,
    _update_dpos_square_overlay,
)
from ...utils.utils_Qt.utils_general import transform_ROI_params_inverse, make_handles


def _apply_kwargs_bounds_to_roi_sliders(roi, lo, hi, eps=1e-8, plot_widget=None):
    """
    Apply ``kwargs_lower`` / ``kwargs_upper`` entries to the ROI's TripleSliders.

    Lenstronomy names are mapped to slider keys (e.g. ``source_amp`` -> ``amp``).
    Expands each slider's ``min_range`` / ``max_range`` if the requested bounds
    lie outside the defaults from ``get_light_model_ranges``.
    """
    if not isinstance(lo, dict) or not isinstance(hi, dict):
        return
    sliders = getattr(roi, 'sliders', None)
    if not sliders:
        return

    def set_triple_bounds(s, lo_v, hi_v):
        lo_v, hi_v = float(lo_v), float(hi_v)
        if lo_v >= hi_v:
            return
        s.min_range = min(s.min_range, lo_v)
        s.max_range = max(s.max_range, hi_v)
        s.vmin = lo_v
        s.vmax = hi_v
        s.vmid = max(min(s.vmid, s.vmax - eps), s.vmin + eps)
        s.update()

    for k_kwargs, k_slider in (
        ('amp', 'amp'),
        ('n_sersic', 'n_sersic'),
        ('sigma', 'sigma'),
        ('R_sersic', 'R_sersic'),
        ('q', 'q'),
        ('phi', 'phi'),
        ('source_amp', 'amp'),
    ):
        if k_kwargs in lo and k_kwargs in hi and k_slider in sliders:
            set_triple_bounds(sliders[k_slider], lo[k_kwargs], hi[k_kwargs])

    if 'dpos' in sliders:
        s = sliders['dpos']
        dpos_val = None
        if all(k in lo and k in hi for k in ('center_x', 'center_y')):
            dx = (float(hi['center_x']) - float(lo['center_x'])) / 2.0
            dy = (float(hi['center_y']) - float(lo['center_y'])) / 2.0
            dpos_val = 0.5 * (dx + dy)
        elif all(k in lo and k in hi for k in ('ra_source', 'dec_source')):
            dx = (float(hi['ra_source']) - float(lo['ra_source'])) / 2.0
            dy = (float(hi['dec_source']) - float(lo['dec_source'])) / 2.0
            dpos_val = 0.5 * (dx + dy)
        if dpos_val is not None and dpos_val > 0:
            s.vmid = max(min(dpos_val, s.max_range - eps), s.min_range + eps)
            s.update()

    # Programmatic dpos changes do not emit TripleSlider.valuesChanged; refresh overlay.
    if plot_widget is not None:
        _update_dpos_square_overlay(plot_widget, roi)


class lenstronomy_model :
    def __init__(self, lm_dict_or_str, imsim, chain_list=None) :
        # lm is a dictionary containing the parameters of a source model, either in local or world coordinates.
        self.imsim = imsim
        self.chain_list = chain_list
        if type(lm_dict_or_str)==str :
            lm_dict_or_str = self.load_lm(path=lm_dict_or_str)
        if 'chain_list' in lm_dict_or_str :
            self.chain_list = lm_dict_or_str['chain_list']
            del lm_dict_or_str['chain_list']
        if 'lens_model_list' in lm_dict_or_str['models'] :
            self.local = lm_dict_or_str
            self.world = self.to_world()
        else :
            self.world = lm_dict_or_str
            self.local = self.to_local()

    def _iter_kwargs_blocks(self, model_dict):
        """
        Yield all sub-dicts that can contain kwargs_source / kwargs_source_ps blocks.

        This includes the standard `model_dict['kwargs']` plus any optimization
        blocks added by `make_lm_dict_opt` (e.g. kwargs_init/sigma/lower/upper/fixed,
        kwargs_opt, kwargs_MCMC, ...).
        """
        if not isinstance(model_dict, dict):
            return
        for k, v in model_dict.items():
            if k == 'models':
                continue
            if isinstance(v, dict) and ('kwargs_source' in v or 'kwargs_source_ps' in v):
                yield v

    def _local_to_world_coords(self, d):
        """Convert a single kwargs dict from local-relative to world coordinates."""
        if not isinstance(d, dict):
            return

        # visualens naming convention: center_x/center_y -> center_ra/center_dec
        if 'center_x' in d and 'center_y' in d:
            ra, dec = relative_to_world(d['center_x'], d['center_y'], self.imsim.center_world)
            d['center_ra'] = ra
            d['center_dec'] = dec
            del d['center_x']
            del d['center_y']

        # lenstronomy naming convention: ra_source/dec_source are local-relative in Lenstronomy
        if 'ra_source' in d and 'dec_source' in d:
            ra, dec = relative_to_world(d['ra_source'], d['dec_source'], self.imsim.center_world)
            d['center_ra'] = ra
            d['center_dec'] = dec
            del d['ra_source']
            del d['dec_source']

    def _world_to_local_coords(self, d):
        """Convert a single kwargs dict from world to local-relative coordinates."""
        if not isinstance(d, dict):
            return

        if 'center_ra' in d and 'center_dec' in d:
            xr, yr = world_to_relative(d['center_ra'], d['center_dec'], self.imsim.center_world)
            d['center_x'] = xr
            d['center_y'] = yr
            del d['center_ra']
            del d['center_dec']

        #if 'center_ra' in d and 'center_dec' in d:
        #    xr, yr = world_to_relative(d['center_ra'], d['center_dec'], self.imsim.center_world)
        #    d['ra_source'] = xr
        #    d['dec_source'] = yr
        #    del d['center_ra']
        #    del d['center_dec']
    
    def plot(self, source_resolution=1000, cmap='grey', axes=None, cutoff=1.) :
        if axes is None :
            fig, axes = plt.subplots(1, 4, figsize=(4, 16), sharex=False, sharey=False)
        else :
            fig = axes[0].figure

        """ Observed image """
        observed_image = self.imsim.image_data.copy()
        observed_image_crop = observed_image.copy()
        observed_image_crop[observed_image>cutoff] = cutoff

        """ Image simulated """
        image_array = self._make_image_array()
        image_array_crop = image_array.copy()
        image_array_crop[image_array>cutoff] = cutoff

        vmin = min(observed_image_crop.min(), image_array_crop.min())
        vmax = max(observed_image_crop.max(), image_array_crop.max())
        axes[0].imshow(observed_image_crop, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
        axes[1].imshow(image_array_crop, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

        """ Residuals """
        residuals = observed_image - image_array
        axes[2].imshow(residuals, origin='lower', cmap='bwr')

        """ Source plane """
        source_array = self._make_source_array(source_resolution=source_resolution)
        #axes[3].imshow(source_array, origin='lower', cmap='inferno')
        source_array_crop = source_array.copy()
        #cutoff_src = 0.008
        #source_array_crop[source_array>cutoff_src] = cutoff_src
        axes[3].imshow(np.log10(source_array), origin='lower', cmap=cmap)

        for ax in axes :
            ax.axis('off')
        fig.tight_layout(pad=0)



        if False :
            kwargs_name = 'kwargs_opt' if 'kwargs_opt' in self.local else 'kwargs'
            kwargs = self.local[kwargs_name].copy()
            kwargs['kwargs_lens'] = self.local['kwargs']['kwargs_lens']
            kwargs['kwargs_lens_light'] = []
            if 'kwargs_source_ps' in kwargs :
                kwargs['kwargs_ps'] = kwargs['kwargs_source_ps']
                del kwargs['kwargs_source_ps']
            else :
                kwargs['kwargs_ps'] = []
                
            ps = []
            kwargs_name = 'kwargs_opt' if 'kwargs_opt' in self.local else 'kwargs'
            for i, src in enumerate(self.local[kwargs_name]['kwargs_source_ps']) :
                ps.append({})
                for param in src :
                    if param=='center_x' :
                        ps[i]['ra_source'] = src[param] - self.imsim.source_center_coordinates[0]
                    elif param=='center_y' :
                        ps[i]['dec_source'] = src[param] - self.imsim.source_center_coordinates[1]
                    else :
                        ps[i][param] = src[param]
            
            kwargs['kwargs_ps'] = ps

            self.modelPlot = ModelPlot([[self.imsim.ImageData_kwargs, self.imsim.PSF_kwargs, self.imsim.kwargs_numerics]], 
                            self.local['models'], kwargs, arrow_size=0.02, cmap_string="gist_heat",
                            linear_solver=False)
                            
            f, axes = plt.subplots(1, 3, figsize=(16, 4), sharex=False, sharey=False)
            self.modelPlot.data_plot(ax=axes[0])
            self.modelPlot.model_plot(ax=axes[1])
            self.modelPlot.normalized_residual_plot(ax=axes[2], v_min=-6, v_max=6)
            #self.modelPlot.source_plot(ax=axes[1, 0], deltaPix_source=0.001, numPix=1000)
            #self.modelPlot.convergence_plot(ax=axes[1, 1], v_max=1)
            #self.modelPlot.magnification_plot(ax=axes[1, 2])
            f.tight_layout()
            f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
            plt.show()
        
            f2, a2 = plt.subplots()
            a2 = self.modelPlot.source_plot(ax=a2, deltaPix_source=0.0005, numPix=1500, cmap='gray') #0.0005, 1500
            source_array = a2.get_children()[0].get_array()
            
            fig, ax = plt.subplots()
            #ax.imshow( source_array-np.min(source_array), origin="lower", cmap='gray' )
            ax.imshow( source_array, vmax=np.log(0.06), origin='lower', cmap='gray' )
            
            f3, a3 = plt.subplots()
            a3.imshow( np.exp(source_array), vmax=0.06, origin="lower", cmap='gray' )
        
        
        #if self.chain_list is not None :
        #    for i in range(len(self.chain_list)):
        #        chain_plot.plot_chain_list(self.chain_list, 0)
        
        #return ax, source_array
        #return source_array
        
    
    def to_world(self) :
        model_world = copy.deepcopy(self.local)
        if 'lens_model_list' in model_world['models'] :
            del model_world['models']['lens_model_list']
        if 'kwargs' in model_world and 'kwargs_lens' in model_world['kwargs'] :
            del model_world['kwargs']['kwargs_lens']

        for block in self._iter_kwargs_blocks(model_world):
            if 'kwargs_source' in block:
                for src in block['kwargs_source']:
                    self._local_to_world_coords(src)
            if 'kwargs_source_ps' in block:
                for ps in block['kwargs_source_ps']:
                    self._local_to_world_coords(ps)
        return model_world


    def to_local(self) :        
        model_local = copy.deepcopy(self.world)
        model_local['models']['lens_model_list'] = self.imsim.LensModel_list
        if 'kwargs' not in model_local:
            model_local['kwargs'] = {}
        model_local['kwargs']['kwargs_lens'] = self.imsim.LensModel_kwargs

        for block in self._iter_kwargs_blocks(model_local):
            if 'kwargs_source' in block:
                for src in block['kwargs_source']:
                    self._world_to_local_coords(src)
            if 'kwargs_source_ps' in block:
                for ps in block['kwargs_source_ps']:
                    self._world_to_local_coords(ps)
        return model_local


    def add(self, model_2) :
        if type(model_2)==str :
            model_2 = self.load_lm(path=model_2)
        if 'lens_model_list' in model_2['models'] :
            self.local['models']['source_light_model_list'] += model_2['models']['source_light_model_list']
            self.local['models']['point_source_model_list'] += model_2['models']['point_source_model_list']
            for name in ['kwargs', 'kwargs_fixed', 'kwargs_lower', 'kwargs_upper', 'kwargs_sigma', 'kwargs_init', 'kwargs_opt', 'kwargs_MCMC'] : #, 'kwargs_PSO'
                self.local[name]['kwargs_source'] += model_2[name]['kwargs_source']
                self.local[name]['kwargs_source_ps'] += model_2[name]['kwargs_source_ps']

            #self.__init__(self.local, self.imsim)
            self.world = self.to_world()
            return self.local
        else :
            self.world['models']['source_light_model_list'] += model_2['models']['source_light_model_list']
            self.world['models']['point_source_model_list'] += model_2['models']['point_source_model_list']
            for name in ['kwargs', 'kwargs_fixed', 'kwargs_lower', 'kwargs_upper', 'kwargs_sigma', 'kwargs_init', 'kwargs_opt', 'kwargs_MCMC'] : #, 'kwargs_PSO'
                self.world[name]['kwargs_source'] += model_2[name]['kwargs_source']
                self.world[name]['kwargs_source_ps'] += model_2[name]['kwargs_source_ps']

            #self.__init__(self.world, self.imsim)
            self.local = self.to_local()
            return self.world


    def save(self, path='./lenstronomy_model.pkl') :
        #if os.path.exists(path) :
        #    with open(path, 'rb') as file :
        #        existing_model = pickle.load(file)
        #else :
        #    existing_model = {'models': {'source_light_model_list': []}, 
        #                      'results': {'kwargs_source': []}}
        
        #full_model = join_lenstronomy_model(existing_model, formatted_model)
        
        path = make_model_path(path)
        
        to_write = self.world.copy()
        
        if self.chain_list is not None :
            to_write['chain_list'] = self.chain_list
        
        with open(path, 'wb') as file :
            pickle.dump(to_write, file)
        print('Model saved at ' + path)


    def load_lm(self, path='./lenstronomy_model.pkl') :
        path = make_model_path(path)
        if os.path.exists(path) :
            print('Source model found at: ' + path)
        
        with open(path, 'rb') as file :
            imported_model_world = pickle.load(file)
        
        # For backward compatibility
        if 'results' in imported_model_world :
            imported_model_world['kwargs'] = copy.deepcopy(imported_model_world['results'])
            del imported_model_world['results']
        
        #self.__init__(self.world, self.imsim)
        #self.local = self.to_local()
        return imported_model_world
    
    
    def plot_mcmc(self, figsize=None, **corner_kwargs):
        """Corner plot of the MCMC parameter samples stored in ``kwargs_MCMC``."""

        if 'kwargs_MCMC' not in self.local:
            raise KeyError("No MCMC kwargs found in the model. Expected key 'kwargs_MCMC'.")

        mcmc_block = self.local['kwargs_MCMC']

        columns, labels = [], []
        for block_key in ('kwargs_source', 'kwargs_source_ps'):
            if block_key not in mcmc_block:
                continue
            prefix = 'src' if block_key == 'kwargs_source' else 'ps'
            for i, kw in enumerate(mcmc_block[block_key]):
                for param, val in kw.items():
                    arr = np.array(val)
                    if arr.size <= 1:
                        continue
                    columns.append(arr)
                    labels.append(f"{prefix}{i}_{param}")
        
        lengths = [len(col) for col in columns]
        if len(np.unique(lengths)) > 1 :
            length_min = np.min(lengths)
            for i in range(len(columns)) :
                columns[i] = columns[i][:length_min]
        
        data = np.column_stack(columns)

        kw = dict(labels=labels, show_titles=True, title_fmt=".3f")
        if figsize is not None:
            kw['fig'] = plt.figure(figsize=figsize)
        kw.update(corner_kwargs)

        fig = corner.corner(data, **kw)
        plt.show()
        return fig

    def remove_source(self, index) :
        del self.local['models']['source_light_model_list'][index]
        for name in ['kwargs', 'kwargs_fixed', 'kwargs_lower', 'kwargs_upper', 'kwargs_sigma', 'kwargs_init', 'kwargs_opt', 'kwargs_MCMC', 'kwargs_PSO'] :
            del self.local[name]['kwargs_source'][index]
            del self.local[name]['kwargs_source_ps'][index]
        self.world = self.to_world()

    def _apply_bounds_if_present(self, roi, block_key, index):
        """Apply ``kwargs_lower``/``kwargs_upper`` for ``block_key`` at ``index`` to ``roi`` sliders."""
        if 'kwargs_lower' not in self.local or 'kwargs_upper' not in self.local:
            return
        kl = self.local['kwargs_lower'].get(block_key)
        ku = self.local['kwargs_upper'].get(block_key)
        if kl is None or ku is None or index >= len(kl) or index >= len(ku):
            return
        _apply_kwargs_bounds_to_roi_sliders(
            roi, kl[index], ku[index], plot_widget=self.imsim.source_plane_widget
        )

    def send_to_imsim(self) :
        PlotWidget = self.imsim.source_plane_widget
        kwargs_name = 'kwargs_opt' if 'kwargs_opt' in self.local else 'kwargs'
        
        for i, model in enumerate(self.local['models']['source_light_model_list']) :
            if model=='SERSIC' :
                kwargs = self.local[kwargs_name]['kwargs_source'][i]
                x = kwargs['center_x'] - self.imsim.source_center_coordinates[0]
                y = kwargs['center_y'] - self.imsim.source_center_coordinates[1]
                radius = kwargs['R_sersic']
                x_corner = x - radius
                y_corner = y - radius
                
                PlotWidget.last_roi = SelectableCircleROI([x_corner, y_corner], radius=radius, pen='purple', invertible=True)
                PlotWidget.last_roi.type = 3
                PlotWidget.last_roi.type_str = 'SERSIC'
                PlotWidget.roi_list.append(PlotWidget.last_roi)
                for handle in PlotWidget.last_roi.handles:
                    PlotWidget.last_roi.removeHandle(handle['item'])
                PlotWidget.addItem(PlotWidget.last_roi)
                
                PlotWidget.last_roi.addScaleHandle([0.5+2**-1.5, 0.5+2**-1.5], [0.5, 0.5])
                
                # Add sliders to control light source parameters
                params = ['amp', 'n_sersic']#, 'R_sersic']
                attach_sliders(PlotWidget, params)
                
                for p, v in kwargs.items() :
                    if p in PlotWidget.last_roi.sliders :
                        PlotWidget.last_roi.sliders[p].vmid = v

                self._apply_bounds_if_present(PlotWidget.last_roi, 'kwargs_source', i)

            if model in ('SERSIC_ELLIPSE_Q_PHI', 'SERSIC_ELLIPSE') :
                kwargs = self.local[kwargs_name]['kwargs_source'][i]
                x_center = kwargs['center_x'] - self.imsim.source_center_coordinates[0]
                y_center = kwargs['center_y'] - self.imsim.source_center_coordinates[1]
                R_sersic = kwargs['R_sersic']
                #e1, e2 = kwargs['e1'], kwargs['e2']
                #e = (e1**2 + e2**2)**0.5
                #q = (1-e) / (1+e)
                #a = 2*R_sersic / (1+q)
                #b = a*q
                #theta = 0.5*np.arctan(e2/e1)
                #x_corner, y_corner, a, b, angle = transform_ROI_params_inverse(x_center, y_center, a, b, theta)
                if model == 'SERSIC_ELLIPSE_Q_PHI' :
                    q = kwargs['q']
                    phi = kwargs['phi']
                else :
                    e1, e2 = kwargs['e1'], kwargs['e2']
                    phi, q = param_util.ellipticity2phi_q(float(e1), float(e2))
                semi_major = R_sersic / np.sqrt(q)
                semi_minor = R_sersic * np.sqrt(q)
                x_corner, y_corner, a, b, angle = transform_ROI_params_inverse(
                    x_center, y_center, semi_major, semi_minor, phi
                )
                
                PlotWidget.last_roi = SelectableEllipseROI([x_corner, y_corner], [a, b], angle=angle, pen='orange', invertible=True)
                PlotWidget.last_roi.type = 1
                PlotWidget.last_roi.type_str = 'SERSIC_ELLIPSE_Q_PHI'
                PlotWidget.roi_list.append(PlotWidget.last_roi)
                for handle in PlotWidget.last_roi.handles:
                    PlotWidget.last_roi.removeHandle(handle['item'])
                PlotWidget.addItem(PlotWidget.last_roi)
                
                make_handles(PlotWidget.last_roi)
                
                # Add sliders to control light source parameters
                params = ['amp', 'n_sersic', 'q', 'phi']
                attach_sliders(PlotWidget, params)
                
                for p, v in kwargs.items() :
                    if p in PlotWidget.last_roi.sliders :
                        PlotWidget.last_roi.sliders[p].vmid = v

                self._apply_bounds_if_present(PlotWidget.last_roi, 'kwargs_source', i)

            if model=='GAUSSIAN' :
                kwargs = self.local[kwargs_name]['kwargs_source'][i]
                x = kwargs['center_x'] - self.imsim.source_center_coordinates[0]
                y = kwargs['center_y'] - self.imsim.source_center_coordinates[1]
                radius = kwargs['sigma'] / 2.0
                x_corner = x - radius
                y_corner = y - radius
                
                PlotWidget.last_roi = SelectableCircleROI([x_corner, y_corner], radius=radius, pen='red', invertible=True)
                PlotWidget.last_roi.type = 2
                PlotWidget.last_roi.type_str = 'GAUSSIAN'
                PlotWidget.roi_list.append(PlotWidget.last_roi)
                for handle in PlotWidget.last_roi.handles:
                    PlotWidget.last_roi.removeHandle(handle['item'])
                PlotWidget.addItem(PlotWidget.last_roi)
                
                PlotWidget.last_roi.addScaleHandle([0.5+2**-1.5, 0.5+2**-1.5], [0.5, 0.5])
                
                # Add sliders to control light source parameters
                params = ['amp']
                attach_sliders(PlotWidget, params)
                
                for p, v in kwargs.items() :
                    if p in PlotWidget.last_roi.sliders :
                        PlotWidget.last_roi.sliders[p].vmid = v

                self._apply_bounds_if_present(PlotWidget.last_roi, 'kwargs_source', i)

        for i, model in enumerate(self.local['models']['point_source_model_list']) :
            if model == 'SOURCE_POSITION' :
                kwargs = self.local[kwargs_name]['kwargs_source_ps'][i]
                x = kwargs['center_x'] - self.imsim.source_center_coordinates[0] #kwargs['ra_source']
                y = kwargs['center_y'] - self.imsim.source_center_coordinates[1] #kwargs['dec_source']
                radius = 1e-3
                x_corner = x - radius
                y_corner = y - radius

                PlotWidget.last_roi = SelectableCircleROI([x_corner, y_corner], radius=radius, pen='cyan', invertible=True)
                PlotWidget.last_roi.type = 4
                PlotWidget.last_roi.type_str = 'SOURCE_POSITION'
                PlotWidget.roi_list.append(PlotWidget.last_roi)
                for handle in PlotWidget.last_roi.handles:
                    PlotWidget.last_roi.removeHandle(handle['item'])
                PlotWidget.addItem(PlotWidget.last_roi)

                params = ['amp']
                attach_sliders(PlotWidget, params)

                amp_val = kwargs.get('source_amp', kwargs.get('amp'))
                if amp_val is not None and 'amp' in PlotWidget.last_roi.sliders :
                    PlotWidget.last_roi.sliders['amp'].vmid = amp_val

                self._apply_bounds_if_present(PlotWidget.last_roi, 'kwargs_source_ps', i)
    
    
    def _make_source_array(self, source_resolution=1000) :
        LensModel_empty = LensModel(lens_model_list=[])
        LightModel_source = LightModel(light_model_list=self.local['models']['source_light_model_list'])
        PSModel_source = PointSource(point_source_type_list=self.local['models']['point_source_model_list'], fixed_magnification_list=[True for src in self.local['models']['point_source_model_list']])

        offset_source = []
        offset_ps = []
        
        kwargs_name = 'kwargs_opt' if 'kwargs_opt' in self.local else 'kwargs'
        
        x_list = [src['center_x'] for src in self.local[kwargs_name]['kwargs_source_ps']] + [src['center_x'] for src in self.local[kwargs_name]['kwargs_source']]
        y_list = [src['center_y'] for src in self.local[kwargs_name]['kwargs_source_ps']] + [src['center_y'] for src in self.local[kwargs_name]['kwargs_source']]
        x_mean, y_mean = np.mean(x_list), np.mean(y_list)
        
        for i, src in enumerate(self.local[kwargs_name]['kwargs_source']) :
            offset_source.append({})
            for param in src :
                if param=='center_x' :
                    offset_source[i][param] = src[param] - x_mean
                elif param=='center_y' :
                    offset_source[i][param] = src[param] - y_mean
                else :
                    offset_source[i][param] = src[param]
        for i, src in enumerate(self.local[kwargs_name]['kwargs_source_ps']) :
            offset_ps.append({})
            for param in src :
                if param=='center_x' :
                    offset_ps[i]['ra_source'] = src[param] - x_mean
                elif param=='center_y' :
                    offset_ps[i]['dec_source'] = src[param] - y_mean
                else :
                    offset_ps[i][param] = src[param]
        
        pix_scale = 0.001
        extent = 0.5
        transform_pix2angle = np.array([[1, 0], [0, 1]]) * pix_scale
        PixelGrid_kwargs = {'nx': round(extent / pix_scale),
                            'ny': round(extent / pix_scale),
                            'ra_at_xy_0': - extent/2,
                            'dec_at_xy_0': - extent/2,
                            'transform_pix2angle': transform_pix2angle} 
        pixel_grid = PixelGrid(**PixelGrid_kwargs)
                    
        image_model_class = ImageModel( data_class=pixel_grid, psf_class=self.imsim.PSF,
                                        lens_model_class=LensModel_empty,
                                        source_model_class=LightModel_source,
                                        point_source_class=PSModel_source,
                                        #lens_light_model_class=,
                                        kwargs_numerics=self.imsim.kwargs_numerics )

        source_array = image_model_class.image(kwargs_source=offset_source,
                                                kwargs_ps=offset_ps,
                                                #kwargs_lens_light=kwargs_light_lens,
                                                kwargs_lens=[{}], unconvolved=False)
        return source_array
        
    def _make_image_array(self) :
        LightModel_source = LightModel(light_model_list=self.local['models']['source_light_model_list'])
        PSModel_source = PointSource(point_source_type_list=self.local['models']['point_source_model_list'], fixed_magnification_list=[True for src in self.local['models']['point_source_model_list']])

        ps = []
        kwargs_name = 'kwargs_opt' if 'kwargs_opt' in self.local else 'kwargs'
        for i, src in enumerate(self.local[kwargs_name]['kwargs_source_ps']) :
            ps.append({})
            for param in src :
                if param=='center_x' :
                    ps[i]['ra_source'] = src[param]
                elif param=='center_y' :
                    ps[i]['dec_source'] = src[param]
                else :
                    ps[i][param] = src[param]

        image_model_class = ImageModel( data_class=self.imsim.PixelGrid, psf_class=self.imsim.PSF,
                                        lens_model_class=self.imsim.LensModel,
                                        source_model_class=LightModel_source,
                                        point_source_class=PSModel_source,
                                        #lens_light_model_class=,
                                        kwargs_numerics=self.imsim.kwargs_numerics )

        image_array = image_model_class.image(kwargs_source=self.local[kwargs_name]['kwargs_source'],
                                                kwargs_ps=ps,
                                                #kwargs_lens_light=kwargs_light_lens,
                                                kwargs_lens=self.local['kwargs']['kwargs_lens'], unconvolved=False)
        return image_array

    def make_latex(self, ref=None, pc=True, uncertainty='hpd', flux_per_amp=None, flux_unit=None):
        """Return (and print) a LaTeX table of all light-model component parameters.

        When ``kwargs_MCMC`` is present, each cell shows the central value with
        asymmetric uncertainties ``$v^{+hi}_{-lo}$``.  The method used depends
        on ``uncertainty``:

        * ``'hpd'`` *(default)* – centre is the best-fit (``kwargs_opt``);
          bounds are the 68.3 % Highest Posterior Density interval of the MCMC
          chain.
        * ``'median'`` – centre is the MCMC median (50th percentile); bounds
          are the 16th and 84th percentiles.  When ``ref`` is an integer index
          the reference position is also taken from the MCMC median of that
          component (not the best-fit).

        Parameters are not MCMC-sampled fall back to the plain point estimate.

        Parameters
        ----------
        ref : None, tuple, or int
            * ``None``  – coordinates are left as-is (arcsec relative to the
              image-cutout centre; no unit conversion is applied).
            * ``(ra, dec)`` – world coordinates (degrees) used as the origin;
              each component's position is expressed relative to this point.
            * ``int``   – index into the combined list of components (source
              light models first, then point-source models); the position of
              that component becomes the origin.
        pc : bool
            When ``ref`` is not ``None``, convert the relative arcsec offsets
            to parsecs using ``imsim.z_source`` and the default cosmology.
            Has no effect when ``ref=None``.
        uncertainty : {'hpd', 'median'}
            Method used to derive uncertainties from the MCMC chain.

        Returns
        -------
        str  – the complete LaTeX ``table`` environment as a string.
        """
        if uncertainty not in ('hpd', 'median'):
            raise ValueError(f"uncertainty must be 'hpd' or 'median', got {uncertainty!r}")
        
        cosmo = get_cosmo()
        z = self.imsim.z_source
        
        if flux_unit is None :
            if 'BUNIT' in self.imsim.individual_filter.header :
                flux_unit = self.imsim.individual_filter.header['BUNIT']
            else :
                flux_unit = 'arbitrary units'
        flux_per_amp = 1. if flux_per_amp is None else flux_per_amp
        
        kwargs_key = 'kwargs_opt' if 'kwargs_opt' in self.local else 'kwargs'
        kw_block = self.local[kwargs_key]
        models = self.local['models']

        # --- collect all components as (kind, local_index, model_name, kwargs) ---
        components = []
        for i, model in enumerate(models.get('source_light_model_list', [])):
            components.append(('src', i, model, dict(kw_block['kwargs_source'][i])))
        for i, model in enumerate(models.get('point_source_model_list', [])):
            components.append(('ps', i, model, dict(kw_block['kwargs_source_ps'][i])))

        # --- MCMC samples (optional) ---
        mcmc_block = self.local.get('kwargs_MCMC')
        mcmc_src   = (mcmc_block or {}).get('kwargs_source',    [])
        mcmc_ps    = (mcmc_block or {}).get('kwargs_source_ps', [])

        def _mcmc_kw(kind, idx):
            """Return the MCMC kwargs dict for this component, or None."""
            if mcmc_block is None:
                return None
            pool = mcmc_src if kind == 'src' else mcmc_ps
            return pool[idx] if idx < len(pool) else None

        def _samples(mcmc_kw, col):
            """Return a 1-D float array of MCMC samples for ``col``, or None."""
            if mcmc_kw is None:
                return None
            keys = ('amp', 'source_amp') if col == 'amp' else (col,)
            for k in keys:
                if k in mcmc_kw:
                    arr = np.asarray(mcmc_kw[k], dtype=float).ravel()
                    if arr.size > 1:
                        return arr
            return None

        # --- resolve reference position in local arcsec ---
        do_offset = ref is not None
        ref_x, ref_y = 0.0, 0.0

        if isinstance(ref, (tuple, list)):
            ref_x, ref_y = world_to_relative(ref[0], ref[1], self.imsim.center_world)
        elif isinstance(ref, int):
            if ref < 0 or ref >= len(components):
                raise IndexError(
                    f"ref={ref} is out of range for {len(components)} component(s)."
                )
            ref_kind, ref_idx, _, kw_ref = components[ref]
            ref_x = float(kw_ref.get('center_x', kw_ref.get('ra_source', 0.0)))
            ref_y = float(kw_ref.get('center_y', kw_ref.get('dec_source', 0.0)))
            # For the median method, override with the MCMC median position
            if uncertainty == 'median' and mcmc_block is not None:
                ref_mk = _mcmc_kw(ref_kind, ref_idx)
                for raw_key, alt_key, attr in (
                    ('center_x', 'ra_source',  'ref_x'),
                    ('center_y', 'dec_source', 'ref_y'),
                ):
                    samp = _samples(ref_mk, raw_key)
                    if samp is None:
                        samp = _samples(ref_mk, alt_key)
                    if samp is not None:
                        if attr == 'ref_x':
                            ref_x = float(np.median(samp))
                        else:
                            ref_y = float(np.median(samp))
        elif ref is not None:
            raise TypeError(f"ref must be None, a (ra, dec) tuple, or an int; got {type(ref)}")

        # arcsec → pc scale factor (used for both coordinate offsets and size params)
        _pc_scale    = float(arcsec_to_kpc(1.0, z).value) * 1000.0 if pc else 1.0
        _coord_scale = _pc_scale if do_offset else 1.0

        coord_unit   = 'pc' if (do_offset and pc) else 'arcsec'
        size_unit    = 'pc' if pc else 'arcsec'
        coord_prefix = r'$\Delta$' if do_offset else ''
        x_header = rf'{coord_prefix}$x$ [{coord_unit}]'
        y_header = rf'{coord_prefix}$y$ [{coord_unit}]'

        # --- parameter display order and LaTeX header names ---
        _latex_name = {
            'flux':     f'Flux density [{flux_unit}]',
            'ab_mag':   r'$m_{\rm AB}$',
            'size':     f'sigma or R_sersic [{size_unit}]',
            'n_sersic': r'$n$',
            'e1':       r'$e_1$',
            'e2':       r'$e_2$',
            'q':        r'$q$',
            'phi':      r'$\phi$',
            'r':        rf'{coord_prefix}$r$ [{coord_unit}]',
        }
        _param_order = ['flux', 'ab_mag', 'size', 'n_sersic', 'e1', 'e2', 'q', 'phi']
        _amp_aliases = {'source_amp'}
        _size_keys   = {'R_sersic', 'sigma'}   # angular sizes converted alongside coords
        _coord_keys  = {'center_x', 'center_y', 'ra_source', 'dec_source'}

        # --- formatters ---
        def _fmt(v):
            if isinstance(v, (float, np.floating)):
                return f'{float(v):.4g}'
            if isinstance(v, (int, np.integer)):
                return str(int(v))
            return str(v).replace('_', r'\_')

        def _hpd(samples, credible_level=0.6827):
            """Return the shortest interval containing ``credible_level`` of the samples."""
            s = np.sort(samples)
            n = len(s)
            n_in = max(1, int(np.floor(credible_level * n)))
            if n_in >= n:
                return s[0], s[-1]
            widths  = s[n_in:] - s[:n - n_in]
            i_min   = int(np.argmin(widths))
            return float(s[i_min]), float(s[i_min + n_in])

        def _fmt_unc(val, lo, hi):
            """Format ``val +hi / -lo`` with precision set by the uncertainty."""
            ref_unc = min(lo, hi)
            if ref_unc <= 0:
                return _fmt(val)
            mag      = math.floor(math.log10(ref_unc))
            decimals = max(0, 1 - mag)          # 2 sig figs on the uncertainty
            fs = f'.{decimals}f'
            return (rf'${val:{fs}}^{{+{hi:{fs}}}}_{{-{lo:{fs}}}}$')

        def _cell_from_samples(samples, point_val):
            """Format a cell using the chosen uncertainty method, or plain point estimate."""
            if samples is None:
                return _fmt(point_val)
            if uncertainty == 'median':
                p16, p50, p84 = np.percentile(samples, [16, 50, 84])
                return _fmt_unc(p50, p50 - p16, p84 - p50)
            else:  # 'hpd'
                v = float(point_val)
                hpd_lo, hpd_hi = _hpd(samples)
                return _fmt_unc(v, v - hpd_lo, hpd_hi - v)

        # --- build table rows ---
        table_rows = []
        extra_cols_seen = []
        extra_cols_set  = set()

        for kind, idx, model, kw in components:
            mk = _mcmc_kw(kind, idx)

            cx = float(kw.get('center_x', kw.get('ra_source', 0.0)))
            cy = float(kw.get('center_y', kw.get('dec_source', 0.0)))

            # Coordinate cells (with possible MCMC)
            coord_data = {}
            for coord_key, raw_val, ref_val, col in (
                ('center_x', cx, ref_x, 'x'),
                ('center_y', cy, ref_y, 'y'),
            ):
                point_out = (raw_val - ref_val) * _coord_scale
                samp_raw  = _samples(mk, coord_key) if mk is not None else None
                if samp_raw is None and mk is not None:
                    # try the alternate coordinate key name
                    alt = 'ra_source' if coord_key == 'center_x' else 'dec_source'
                    samp_raw = _samples(mk, alt)
                samp_out = (samp_raw - ref_val) * _coord_scale if samp_raw is not None else None
                if samp_out is not None:
                    # HPD: centre on best-fit; median: centre on sample median (point_out unused)
                    center = point_out if uncertainty == 'hpd' else None
                    cell = _cell_from_samples(samp_out, center)
                else:
                    cell = _fmt(point_out)
                coord_data[col] = {
                    'cell': cell,
                    'point': float(point_out),
                    'samples': samp_out,
                }

            row = {
                '_type': model.replace('_', r'\_'),
                'x': coord_data['x']['cell'],
                'y': coord_data['y']['cell'],
            }
            if do_offset:
                x_point = coord_data['x']['point']
                y_point = coord_data['y']['point']
                x_samp = coord_data['x']['samples']
                y_samp = coord_data['y']['samples']
                r_point = float(np.hypot(x_point, y_point))
                r_samp = None
                if x_samp is not None and y_samp is not None:
                    n = min(len(x_samp), len(y_samp))
                    if n > 1:
                        r_samp = np.hypot(x_samp[:n], y_samp[:n])
                elif x_samp is not None and len(x_samp) > 1:
                    r_samp = np.hypot(x_samp, y_point)
                elif y_samp is not None and len(y_samp) > 1:
                    r_samp = np.hypot(x_point, y_samp)
                center = r_point if uncertainty == 'hpd' else None
                row['r'] = _cell_from_samples(r_samp, center)

            for k, v in kw.items():
                if k in _coord_keys:
                    continue
                col = 'flux' if k in _amp_aliases or k == 'amp' else k
                if col == 'flux':
                    scaled_v    = float(v) * flux_per_amp
                    samp        = _samples(mk, 'amp')
                    scaled_samp = samp * flux_per_amp if samp is not None else None
                    row[col] = _cell_from_samples(scaled_samp, scaled_v)
                    if flux_unit == 'nJy':
                        if scaled_v > 0:
                            mag_v = 31.4 - 2.5 * np.log10(scaled_v)
                            mag_samp = None
                            if scaled_samp is not None:
                                pos = scaled_samp[scaled_samp > 0]
                                if pos.size > 1:
                                    mag_samp = 31.4 - 2.5 * np.log10(pos)
                            center = mag_v if uncertainty == 'hpd' else None
                            row['ab_mag'] = _cell_from_samples(mag_samp, center)
                        else:
                            row['ab_mag'] = '--'
                elif col in _size_keys:
                    scaled_v    = float(v) * _pc_scale
                    samp        = _samples(mk, col)
                    scaled_samp = samp * _pc_scale if samp is not None else None
                    row['size'] = _cell_from_samples(scaled_samp, scaled_v)
                else:
                    row[col] = _cell_from_samples(_samples(mk, col), v)
                out_col = 'size' if col in _size_keys else col
                if out_col not in extra_cols_set and out_col not in {'x', 'y'}:
                    extra_cols_seen.append(out_col)
                    extra_cols_set.add(out_col)
                if col == 'flux' and flux_unit == 'nJy' and 'ab_mag' not in extra_cols_set:
                    extra_cols_seen.append('ab_mag')
                    extra_cols_set.add('ab_mag')

            table_rows.append(row)

        # --- ordered non-coordinate param columns ---
        param_cols = [k for k in _param_order if k in extra_cols_set]
        leftover   = [k for k in extra_cols_seen if k not in set(_param_order)]
        param_cols += leftover

        all_col_keys  = ['_type', 'x', 'y'] + (['r'] if do_offset else []) + param_cols
        all_col_heads = [
            'Type', x_header, y_header,
        ] + ([_latex_name['r']] if do_offset else []) + [_latex_name.get(k, k.replace('_', r'\_')) for k in param_cols]

        ncols    = len(all_col_keys)
        col_spec = 'l' + 'r' * (ncols - 1)

        # --- assemble LaTeX ---
        lines = [
            r'\begin{table}[h]',
            r'\centering',
            rf'\begin{{tabular}}{{{col_spec}}}',
            r'\hline\hline',
            ' & '.join(all_col_heads) + r' \\',
            r'\hline',
        ]

        for row in table_rows:
            cells = [row.get(k, '--') for k in all_col_keys]
            lines.append(' & '.join(cells) + r' \\')

        lines.append(r'\hline')
        lines.append(r'\end{tabular}')

        # caption / note
        if isinstance(ref, (tuple, list)):
            note = (
                rf'Coordinates relative to '
                rf'$(\alpha,\delta)=({ref[0]:.6f}^\circ,\,{ref[1]:.6f}^\circ)$.'
            )
        elif isinstance(ref, int):
            note = rf'Coordinates relative to component {ref} ({components[ref][2]}).'
        else:
            note = 'Coordinates in arcsec relative to the image-cutout centre.'
        if mcmc_block is not None:
            if uncertainty == 'hpd':
                note += (r' Central values are best-fit (kwargs\_opt); '
                         r'uncertainties are 68.3\,\% HPD intervals.')
            else:
                note += (r' Central values are MCMC medians; '
                         r'uncertainties are 16th/84th percentiles.')

        lines.append(rf'\caption{{{note}}}')
        lines.append(r'\end{table}')

        latex = '\n'.join(lines)
        print(latex)
        return latex



def make_model_path(path) :
    if not os.path.isdir( os.path.dirname(path) ) :
        path = os.path.join('./', path)
        extension = '.pkl'
        if not extension in path :
            path += extension
    return path
