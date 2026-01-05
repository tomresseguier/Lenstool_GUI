import numpy as np
from PyQt5.QtCore import Qt, QObject, QEvent
import pyqtgraph as pg
import copy

from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.Sampling.parameters import Param
from lenstronomy.Workflow.fitting_sequence import FittingSequence

from .lenstronomy_model import lenstronomy_model
from ...utils.utils_astro.utils_general import world_to_relative
from ...utils.utils_Qt.utils_general import transform_ROI_params




class SourceFilter(QObject):
    def __init__(self, imsim) :
        super().__init__()
        self.imsim = imsim
        
    def eventFilter(self, obj, event):
        if event.type() == QEvent.KeyPress and event.key() in (Qt.Key_Return, Qt.Key_Enter) :
            self.imsim.LightModel_source_list = []
            self.imsim.LightModel_source_kwargs = []
            
            for roi in self.imsim.source_plane_widget.roi_list :
                if roi.type==1 :
                    self.imsim.LightModel_source_list.append('SERSIC_ELLIPSE')
                    amp = roi.sliders['amp'].vmid
                    n_sersic = roi.sliders['n_sersic'].vmid
                    x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
                    q = semi_minor / semi_major
                    e = (1-q) / (1+q)
                    e1 = e * np.cos(2*angle)
                    e2 = e * np.sin(2*angle)
                    R_sersic = (semi_major + semi_minor)/2
                    src_xr = x_center + self.imsim.source_center_coordinates[0]
                    src_yr = y_center + self.imsim.source_center_coordinates[1]
                    to_add = {'amp': amp, 'R_sersic': R_sersic, 'n_sersic': n_sersic, 'e1': e1, 'e2': e2, 'center_x': src_xr, 'center_y': src_yr}
                    self.imsim.LightModel_source_kwargs.append(to_add)
                elif roi.type==2 :
                    self.imsim.LightModel_source_list.append('GAUSSIAN')
                    amp = roi.sliders['amp'].vmid
                    x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
                    src_xr = x_center + self.imsim.source_center_coordinates[0]
                    src_yr = y_center + self.imsim.source_center_coordinates[1]
                    sigma = abs(roi.size()[0])
                    to_add = {'amp': amp, 'sigma': sigma, 'center_x': src_xr, 'center_y': src_yr}
                    self.imsim.LightModel_source_kwargs.append(to_add)
                if roi.type==3 :
                    self.imsim.LightModel_source_list.append('SERSIC')
                    amp = roi.sliders['amp'].vmid
                    n_sersic = roi.sliders['n_sersic'].vmid
                    x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
                    src_xr = x_center + self.imsim.source_center_coordinates[0]
                    src_yr = y_center + self.imsim.source_center_coordinates[1]
                    R_sersic = abs(roi.size()[0])
                    to_add = {'amp': amp, 'R_sersic': R_sersic, 'n_sersic': n_sersic, 'center_x': src_xr, 'center_y': src_yr}
                    self.imsim.LightModel_source_kwargs.append(to_add)
            
            
            
            self.imsim.LightModel_source = LightModel(light_model_list=self.imsim.LightModel_source_list)
            
            print('calculate image_model')
            self.imsim.ImageModel = ImageModel(data_class=self.imsim.PixelGrid, psf_class=self.imsim.PSF,
                                                        lens_model_class=self.imsim.LensModel,
                                                        source_model_class=self.imsim.LightModel_source,
                                                        #point_source_class=point_source,
                                                        #lens_light_model_class=,
                                                        kwargs_numerics=self.imsim.kwargs_numerics)
            
            print('Start simulating image')
            self.imsim.simulated_image = self.imsim.ImageModel.image(kwargs_source=self.imsim.LightModel_source_kwargs,
                                                                #kwargs_ps=kwargs_ps,
                                                                #kwargs_lens_light=kwargs_light_lens,
                                                                kwargs_lens=self.imsim.LensModel_kwargs, unconvolved=False)
            self.imsim.image_plane_plot.setImage(self.imsim.simulated_image[::-1,:])
            x = (self.imsim._lt_curve_coords_relative_broken[0] - self.imsim._SquareOfInterest_xr_bottomleft) / (self.imsim.fits_image.pix_deg_scale*3600)
            y = self.imsim.simulated_image.shape[0] - (self.imsim._lt_curve_coords_relative_broken[1] - self.imsim._SquareOfInterest_yr_bottomleft) / (self.imsim.fits_image.pix_deg_scale*3600)
            self.imsim.critical_curve_plot.setData(x, y)
            print('done')
            
            #solver = LensEquationSolver(self.imsim.lens_model)
            
            
            
            self.imsim.models = {'lens_model_list': self.imsim.LensModel_list,
                                            #'lens_light_model_list': self.imsim.LensModel_light_list,
                                            'source_light_model_list': self.imsim.LightModel_source_list}
            
            kwargs = {'kwargs_lens': self.imsim.LensModel_kwargs, 'kwargs_source': self.imsim.LightModel_source_kwargs}
            
            lm_dict = format_lm_local(self.imsim.models, kwargs)
            
            self.imsim.lm_current = lenstronomy_model(lm_dict, self.imsim)
            
            
            if event.modifiers() == Qt.ShiftModifier :
                print('Starting optimization')
                
                kwargs_source_fixed = []
                kwargs_source_lower = []
                kwargs_source_upper = []
                kwargs_source_sigma = []
                kwargs_source_init = []
                for i, src in enumerate(self.imsim.LightModel_source_list) :
                    kwargs = self.imsim.LightModel_source_kwargs[i]
                    kwargs_fixed = {}
                    kwargs_lower = {}
                    kwargs_upper = {}
                    kwargs_sigma = {}
                    kwargs_init = {}
                    
                    
                    fixed_params = ['center_x', 'center_y']
                    self.imsim.source_plane_widget.roi_list[i]
                    opt_params = self.imsim.LightModel_source_kwargs.keys()
                    for param in fixed_params :
                        opt_params.remove(param)
                    for param in fixed_params :
                        kwargs_fixed[param] = kwargs[param]
                    for param in opt_params :
                        vmax = self.imsim.source_plane_widget.roi_list[i].sliders[param].vmin
                        vmin = self.imsim.source_plane_widget.roi_list[i].sliders[param].vmax
                        kwargs_lower[param] = vmax
                        kwargs_upper[param] = vmin
                        kwargs_sigma[param] = (vmax + vmin)/10
                        kwargs_init[param] = self.imsim.source_plane_widget.roi_list[i].sliders[param].vmid
                    
                    if False :
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
                            kwargs_upper['sigma'] = self.imsim.sigma
                            kwargs_init['sigma'] = self.imsim.sigma/2
                                
                        if src == 'SERSIC_ELLIPSE' :
                            fixed_params = ['center_x', 'center_y']
                            opt_params = ['amp', 'R_sersic', 'n_sersic', 'e1', 'e2']
                            for p in fixed_params :
                                #kwargs_fixed[p] = kwargs[p]
                                #kwargs_fixed[p] = kwargs[p]
                                
                                kwargs_lower[p] = kwargs[p] -0.001
                                kwargs_upper[p] = kwargs[p] +0.001
                                kwargs_sigma[p] = 0.0001
                                kwargs_init[p] = kwargs[p]
                            for p in opt_params :
                                if p=='n_sersic' :
                                    kwargs_lower[p] = 0.5
                                    kwargs_upper[p] = 10.
                                    kwargs_sigma[p] = 0.2
                                elif p=='e1' or p=='e2' :
                                    kwargs_lower[p] = -1.
                                    kwargs_upper[p] = 1.
                                    kwargs_sigma[p] = 0.1
                                else :
                                    kwargs_lower[p] = kwargs[p] /10
                                    kwargs_upper[p] = kwargs[p] *10
                                    kwargs_sigma[p] = kwargs[p] /10
                                kwargs_init[p] = kwargs[p]                        
                    
                    kwargs_source_fixed.append(kwargs_fixed)
                    kwargs_source_lower.append(kwargs_lower)
                    kwargs_source_upper.append(kwargs_upper)
                    kwargs_source_sigma.append(kwargs_sigma)
                    kwargs_source_init.append(kwargs_init)

                    
                        
                param = Param(self.imsim.models,
                              kwargs_fixed_lens=self.imsim.LensModel_kwargs,
                              kwargs_fixed_source=kwargs_source_fixed,#self.imsim.fixed_source_kwargs,
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
                
                
                kwargs_lens_empty = [{} for _ in self.imsim.LensModel_list]
                lens_params = [kwargs_lens_empty, kwargs_lens_empty, self.imsim.LensModel_kwargs, kwargs_lens_empty, kwargs_lens_empty]
                source_params = [kwargs_source_init, kwargs_source_sigma, kwargs_source_fixed, kwargs_source_lower, kwargs_source_upper]

                self.imsim.optimization_params = {'lens_model': lens_params,
                                                             #'lens_light_model': lens_light_params,
                                                             'source_model': source_params}

                kwargs_likelihood = {'source_marg': False}
                single_band = [[self.imsim.ImageData_kwargs, self.imsim.PSF_kwargs, self.imsim.kwargs_numerics]]
                kwargs_data_joint = {'multi_band_list': single_band, 'multi_band_type': 'multi-linear'}
                
                fitting_seq = FittingSequence(kwargs_data_joint, self.imsim.models, {}, kwargs_likelihood, self.imsim.optimization_params, mpi=False)

                fitting_kwargs_list = [['PSO', {'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}]]

                chain_list = fitting_seq.fit_sequence(fitting_kwargs_list)
                self.imsim.result_kwargs = fitting_seq.best_fit()
                
                lm_dict = format_lm_local(self.imsim.models, self.imsim.result_kwargs)
                
                self.imsim.lm_optimized = lenstronomy_model(lm_dict, self.imsim)
                self.imsim.lm_optimized.plot()
            return True  # Stop propagation
        return False     # Let other events pass through  

            

class ImageFilter(QObject) :
    def __init__(self, imsim) :
        super().__init__()
        self.imsim = imsim
    def eventFilter(self, obj, event) :
        if event.type() == QEvent.KeyPress and event.key() == Qt.Key_Space :
            print("Hello world!")
            return True  # Stop propagation
        if event.type() == QEvent.MouseButtonDblClick :
            if event.button() == Qt.LeftButton :
                im_coords = self.imsim.image_plane_plot.getView().mapSceneToView(event.pos())
                x_im_full = self.imsim.anchor[0] + im_coords.x()
                y_im_full = self.imsim.anchor[1] - im_coords.y()
                ra, dec = self.imsim.fits_image.image_to_world(x_im_full, y_im_full) #to do: replace with filter data
                xr, yr = world_to_relative(ra, dec, self.imsim.center_world)
                src_xr, src_yr = self.imsim.LensModel.ray_shooting(xr, yr, self.imsim.LensModel_kwargs)
                src_xr = src_xr - self.imsim.source_center_coordinates[0]
                src_yr = src_yr - self.imsim.source_center_coordinates[1]
                
                if not hasattr(self.imsim, 'interactive_source_plot') :
                    self.imsim.interactive_source_plot = pg.ScatterPlotItem()
                    self.imsim.source_plane_widget.addItem(self.imsim.interactive_source_plot)
                self.imsim.interactive_source_plot.setData([src_xr], [src_yr], symbol='x')

                return True  # Stop further processing of this event
        #return super().eventFilter(obj, event)
        return False     # Let other events pass through



def join_lm_world(m1, m2) :
    full_model = copy.deepcopy(m1)
    full_model['models']['source_light_model_list'] = m1['models']['source_light_model_list'] + m2['models']['source_light_model_list']
    full_model['kwargs']['kwargs_source'] = m1['kwargs']['kwargs_source'] + m2['kwargs']['kwargs_source']
    return full_model


def join_lm_local(m1, m2) :
    full_model = {'models': {}, 'kwargs': {}}
    full_model['models']['source_light_model_list'] = m1['models']['source_light_model_list'] + m2['models']['source_light_model_list']
    full_model['kwargs']['kwargs_source'] = m1['kwargs']['kwargs_source'] + m2['kwargs']['kwargs_source']
    return full_model


def format_lm_local(models, kwargs) :
    model_local = {'models': {}, 'kwargs': {}}
    model_local['models'] = models
    model_local['kwargs'] = kwargs
    return model_local


def get_light_model_ranges() :
    ranges = {'amp': [0, 100], 'n_sersic': [1, 10]}
    return ranges
    
    


