import numpy as np
import os
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt, QObject, QEvent
import pyqtgraph as pg
import pickle
import copy

from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.Sampling.parameters import Param
from lenstronomy.Plots.model_plot import ModelPlot
from lenstronomy.Workflow.fitting_sequence import FittingSequence

from .utils_astro.utils_general import world_to_relative, relative_to_world
from .utils_Qt.utils_general import transform_ROI_params




class SourceFilter(QObject):
    def __init__(self, lm) :
        super().__init__()
        self.lm = lm
        
    def eventFilter(self, obj, event):
        if event.type() == QEvent.KeyPress and event.key() in (Qt.Key_Return, Qt.Key_Enter) :
            self.lm.LightModel_source_list = []
            self.lm.LightModel_source_kwargs = []
            
            for roi in self.lm.source_plane_widget.roi_list :
                if type(roi).__name__ == 'SelectableEllipseROI' :
                    self.lm.LightModel_source_list.append('SERSIC_ELLIPSE')
                    n_sersic = 3
                    x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
                    q = semi_minor / semi_major
                    e = (1-q) / (1+q)
                    e1 = e * np.cos(2*angle)
                    e2 = e * np.sin(2*angle)
                    R_sersic = (semi_major + semi_minor)/2
                    src_xr = x_center + self.lm.source_center_coordinates[0]
                    src_yr = y_center + self.lm.source_center_coordinates[1]
                    to_add = {'amp': 1, 'R_sersic': R_sersic, 'n_sersic': n_sersic, 'e1': e1, 'e2': e2, 'center_x': x_center, 'center_y': y_center}
                    self.lm.LightModel_source_kwargs.append(to_add)
                elif type(roi).__name__ == 'SelectableCircleROI' :
                    self.lm.LightModel_source_list.append('GAUSSIAN')
                    x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
                    src_xr = x_center + self.lm.source_center_coordinates[0]
                    src_yr = y_center + self.lm.source_center_coordinates[1]
                    sigma = abs(roi.size()[0])
                    to_add = {'amp': sigma, 'sigma': self.lm.sigma, 'center_x': src_xr, 'center_y': src_yr}
                    self.lm.LightModel_source_kwargs.append(to_add)
                    
            self.lm.LightModel_source = LightModel(light_model_list=self.lm.LightModel_source_list)
            
            print('calculate image_model')
            self.lm.ImageModel = ImageModel(data_class=self.lm.PixelGrid, psf_class=self.lm.PSF,
                                                        lens_model_class=self.lm.LensModel,
                                                        source_model_class=self.lm.LightModel_source,
                                                        #point_source_class=point_source,
                                                        #lens_light_model_class=,
                                                        kwargs_numerics=self.lm.kwargs_numerics)
            
            print('Start simulating image')
            self.lm.simulated_image = self.lm.ImageModel.image(kwargs_source=self.lm.LightModel_source_kwargs,
                                                                #kwargs_ps=kwargs_ps,
                                                                #kwargs_lens_light=kwargs_light_lens,
                                                                kwargs_lens=self.lm.LensModel_kwargs, unconvolved=False)
            self.lm.image_plane_plot.setImage(self.lm.simulated_image[::-1,:])
            x = (self.lm._lt_curve_coords_relative_broken[0] - self.lm._SquareOfInterest_xr_bottomleft) / (self.lm.fits_image.pix_deg_scale*3600)
            y = self.lm.simulated_image.shape[0] - (self.lm._lt_curve_coords_relative_broken[1] - self.lm._SquareOfInterest_yr_bottomleft) / (self.lm.fits_image.pix_deg_scale*3600)
            self.lm.critical_curve_plot.setData(x, y)
            print('done')
            
            #solver = LensEquationSolver(self.lm.lens_model)
            
            
            
            self.lm.models = {'lens_model_list': self.lm.LensModel_list,
                                            #'lens_light_model_list': self.lm.LensModel_light_list,
                                            'source_light_model_list': self.lm.LightModel_source_list}
            
            
            
            if event.modifiers() == Qt.ShiftModifier :
                print('Starting optimization')
                
                kwargs_source_fixed = []
                kwargs_source_lower = []
                kwargs_source_upper = []
                kwargs_source_sigma = []
                kwargs_source_init = []
                for i, src in enumerate(self.lm.LightModel_source_list) :
                    kwargs = self.lm.LightModel_source_kwargs[i]
                    kwargs_fixed = {}
                    kwargs_lower = {}
                    kwargs_upper = {}
                    kwargs_sigma = {}
                    kwargs_init = {}
                    if src == 'GAUSSIAN' :
                        fixed_params = ['center_x', 'center_y']
                        opt_params = ['amp', 'sigma']
                        for p in fixed_params :
                            #kwargs_fixed[p] = kwargs[p]
                            #kwargs_fixed[p] = kwargs[p]
                            
                            kwargs_lower[p] = kwargs[p] -0.001
                            kwargs_upper[p] = kwargs[p] +0.001
                            kwargs_sigma[p] = 0.0001
                            kwargs_init[p] = kwargs[p]
                        for p in opt_params :
                            kwargs_lower[p] = kwargs[p] /10
                            kwargs_upper[p] = kwargs[p] *10
                            kwargs_sigma[p] = kwargs[p] /10
                            kwargs_init[p] = kwargs[p]
                        kwargs_upper['sigma'] = self.lm.sigma
                        kwargs_init['sigma'] = self.lm.sigma/2
                            
                    if src == 'SERSIC_ELLIPSE' :
                        print(None)
                    
                    kwargs_source_fixed.append(kwargs_fixed)
                    kwargs_source_lower.append(kwargs_lower)
                    kwargs_source_upper.append(kwargs_upper)
                    kwargs_source_sigma.append(kwargs_sigma)
                    kwargs_source_init.append(kwargs_init)

                    
                        
                param = Param(self.lm.models,
                              kwargs_fixed_lens=self.lm.LensModel_kwargs,
                              kwargs_fixed_source=kwargs_source_fixed,#self.lm.fixed_source_kwargs,
                              #kwargs_fixed_lens_light=kwargs_fixed_lens_light,
                              #kwargs_fixed_ps=kwargs_fixed_ps, 
                              kwargs_lower_lens=[{}],
                              kwargs_lower_source=kwargs_source_lower,
                              #kwargs_lower_lens_light=kwargs_lower_lens_light,
                              #kwargs_lower_ps=kwargs_lower_ps,
                              kwargs_upper_lens=[{}],
                              kwargs_upper_source=kwargs_source_upper,
                              #kwargs_upper_lens_light=kwargs_upper_lens_light,
                              #kwargs_upper_ps=kwargs_upper_ps
                              #, kwargs_lens_init=kwargs_lens
                              #, joint_lens_with_light: [[0, 0, ['center_x', 'center_y']]]
                              )
                param.print_setting()
                
                
                kwargs_lens_empty = [{} for _ in self.lm.LensModel_list]
                lens_params = [kwargs_lens_empty, kwargs_lens_empty, self.lm.LensModel_kwargs, kwargs_lens_empty, kwargs_lens_empty]
                source_params = [kwargs_source_init, kwargs_source_sigma, kwargs_source_fixed, kwargs_source_lower, kwargs_source_upper]

                self.lm.optimization_params = {'lens_model': lens_params,
                                                             #'lens_light_model': lens_light_params,
                                                             'source_model': source_params}

                kwargs_likelihood = {'source_marg': False}
                single_band = [[self.lm.ImageData_kwargs, self.lm.PSF_kwargs, self.lm.kwargs_numerics]]
                kwargs_data_joint = {'multi_band_list': single_band, 'multi_band_type': 'multi-linear'}
                
                fitting_seq = FittingSequence(kwargs_data_joint, self.lm.models, {}, kwargs_likelihood, self.lm.optimization_params, mpi=False)

                fitting_kwargs_list = [['PSO', {'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}]]

                chain_list = fitting_seq.fit_sequence(fitting_kwargs_list)
                self.lm.result_kwargs = fitting_seq.best_fit()
                
                make_LENSTRONOMY_plot(self.lm, self.lm.models, self.lm.result_kwargs)
            return True  # Stop propagation
        return False     # Let other events pass through  

            

class ImageFilter(QObject) :
    def __init__(self, lm) :
        super().__init__()
        self.lm = lm
    def eventFilter(self, obj, event) :
        if event.type() == QEvent.KeyPress and event.key() == Qt.Key_Space :
            print("Hello world!")
            return True  # Stop propagation
        if event.type() == QEvent.MouseButtonDblClick :
            if event.button() == Qt.LeftButton :
                im_coords = self.lm.image_plane_plot.getView().mapSceneToView(event.pos())
                x_im_full = self.lm.anchor[0] + im_coords.x()
                y_im_full = self.lm.anchor[1] - im_coords.y()
                ra, dec = self.lm.fits_image.image_to_world(x_im_full, y_im_full) #to do: replace with filter data
                xr, yr = world_to_relative(ra, dec, self.lm.center_world)
                src_xr, src_yr = self.lm.LensModel.ray_shooting(xr, yr, self.lm.LensModel_kwargs)
                src_xr = src_xr - self.lm.source_center_coordinates[0]
                src_yr = src_yr - self.lm.source_center_coordinates[1]
                
                if not hasattr(self.lm, 'interactive_source_plot') :
                    self.lm.interactive_source_plot = pg.ScatterPlotItem()
                    self.lm.source_plane_widget.addItem(self.lm.interactive_source_plot)
                self.lm.interactive_source_plot.setData([src_xr], [src_yr], symbol='x')

                return True  # Stop further processing of this event
        #return super().eventFilter(obj, event)
        return False     # Let other events pass through


def make_LENSTRONOMY_plot(self, models, kwargs) :
    self.modelPlot = ModelPlot([[self.ImageData_kwargs, self.PSF_kwargs, self.kwargs_numerics]], models, kwargs, arrow_size=0.02, cmap_string="gist_heat")
    

    f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

    self.modelPlot.data_plot(ax=axes[0,0])
    self.modelPlot.model_plot(ax=axes[0,1])
    self.modelPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
    self.modelPlot.source_plot(ax=axes[1, 0], deltaPix_source=0.0003, numPix=1200)
    self.modelPlot.convergence_plot(ax=axes[1, 1], v_max=1)
    self.modelPlot.magnification_plot(ax=axes[1, 2])
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    plt.show()


def save_lenstronomy_model(path, models, results, LENSTRONOMY_center_world) :
    if os.path.exists(path) :
        with open(path, 'rb') as file :
            existing_model = pickle.load(file)
    else :
        existing_model = {'models': {'source_light_model_list': []}, 
                          'results': {'kwargs_source': []}}


    results = copy.deepcopy(results)
    for name in ['kwargs_source'] : #, 'kwargs_lens_light'] :
        if name in results :
            for i in range(len(results[name])) :
                if 'center_x' in results[name][i] and 'center_y' in results[name][i] :
                    ra, dec = relative_to_world(results[name][i]['center_x'], results[name][i]['center_y'], LENSTRONOMY_center_world)
                    results[name][i]['center_ra'] = ra
                    results[name][i]['center_dec'] = dec
                    del results[name][i]['center_x']
                    del results[name][i]['center_y']

    full_model = {'models': {}, 'results': {}}
    full_model['models']['source_light_model_list'] = existing_model['models']['source_light_model_list'] + models['source_light_model_list']
    full_model['results']['kwargs_source'] = existing_model['results']['kwargs_source'] + results['kwargs_source']
    
    with open(path, 'wb') as file :
        pickle.dump(full_model, file)
    
    print('Model saved at ' + path)
    return full_model
    

def load_lenstronomy_model(path, LENSTRONOMY_center_world) :
    with open(path, 'rb') as file :
        model = pickle.load(file)
    for i in range(len(model['results']['kwargs_source'])) :
        if 'center_ra' in model['results']['kwargs_source'][i] and 'center_dec' in model['results']['kwargs_source'][i] :
            ra = model['results']['kwargs_source'][i]['center_ra']
            dec = model['results']['kwargs_source'][i]['center_dec']
            xr, yr = world_to_relative(ra, dec, LENSTRONOMY_center_world)
            print('Delta ra: ', LENSTRONOMY_center_world[0] - ra)
            print('Delta dec: ', LENSTRONOMY_center_world[1] - dec)
            model['results']['kwargs_source'][i]['center_x'] = xr
            model['results']['kwargs_source'][i]['center_y'] = yr
            del model['results']['kwargs_source'][i]['center_ra']
            del model['results']['kwargs_source'][i]['center_dec']
            
    return model