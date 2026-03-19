import numpy as np
from PyQt5.QtCore import Qt, QObject, QEvent
import pyqtgraph as pg
import corner

from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.Sampling.parameters import Param
from lenstronomy.Workflow.fitting_sequence import FittingSequence
from lenstronomy.Plots import chain_plot

from .lenstronomy_model import lenstronomy_model
from .utils import format_lm_local
from ...utils.utils_astro.utils_general import world_to_relative
from ...utils.utils_Qt.utils_general import transform_ROI_params




class SourceFilter(QObject):
    def __init__(self, imsim) :
        super().__init__()
        self.imsim = imsim
        
    def eventFilter(self, obj, event):
        if event.type() == QEvent.KeyPress and event.key() in (Qt.Key_Return, Qt.Key_Enter) :
            LightModel_source_list = []
            LightModel_source_kwargs = []
            
            for roi in self.imsim.source_plane_widget.roi_list :
                if roi.type==1 :
                    LightModel_source_list.append('SERSIC_ELLIPSE')
                    amp = roi.sliders['amp'].vmid
                    n_sersic = roi.sliders['n_sersic'].vmid
                    x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
                    q = semi_minor / semi_major
                    e = (1-q) / (1+q)
                    e1 = e * np.cos(2*angle)
                    e2 = e * np.sin(2*angle)
                    R_sersic = (semi_major + semi_minor) / 2.0
                    src_xr = x_center + self.imsim.source_center_coordinates[0]
                    src_yr = y_center + self.imsim.source_center_coordinates[1]
                    to_add = {'amp': amp, 'R_sersic': R_sersic, 'n_sersic': n_sersic, 'e1': e1, 'e2': e2, 'center_x': src_xr, 'center_y': src_yr}
                    LightModel_source_kwargs.append(to_add)
                elif roi.type==2 :
                    LightModel_source_list.append('GAUSSIAN')
                    amp = roi.sliders['amp'].vmid
                    x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
                    src_xr = x_center + self.imsim.source_center_coordinates[0]
                    src_yr = y_center + self.imsim.source_center_coordinates[1]
                    sigma = abs(roi.size()[0])
                    to_add = {'amp': amp, 'sigma': sigma, 'center_x': src_xr, 'center_y': src_yr}
                    LightModel_source_kwargs.append(to_add)
                if roi.type==3 :
                    LightModel_source_list.append('SERSIC')
                    amp = roi.sliders['amp'].vmid
                    n_sersic = roi.sliders['n_sersic'].vmid
                    x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
                    src_xr = x_center + self.imsim.source_center_coordinates[0]
                    src_yr = y_center + self.imsim.source_center_coordinates[1]
                    R_sersic = abs(roi.size()[0]) / 2.0
                    to_add = {'amp': amp, 'R_sersic': R_sersic, 'n_sersic': n_sersic, 'center_x': src_xr, 'center_y': src_yr}
                    LightModel_source_kwargs.append(to_add)
            
            
            ######### Create the current lenstronomy model class #########
            lm_dict = {}
            lm_dict['models'] = {'lens_model_list': self.imsim.LensModel_list,
                                 #'lens_light_model_list': self.imsim.LensModel_light_list,
                                 'source_light_model_list': LightModel_source_list}
            lm_dict['kwargs'] = {'kwargs_lens': self.imsim.LensModel_kwargs, 
                                 'kwargs_source': LightModel_source_kwargs}
            self.imsim.lm_current = lenstronomy_model(lm_dict, self.imsim)
            
            
            ######### Plot the simulated image #########
            LightModel_source = LightModel(light_model_list=LightModel_source_list)
            
            self.imsim.ImageModel = ImageModel( data_class=self.imsim.PixelGrid, psf_class=self.imsim.PSF,
                                                lens_model_class=self.imsim.LensModel,
                                                source_model_class=LightModel_source,
                                                #point_source_class=point_source,
                                                #lens_light_model_class=,
                                                kwargs_numerics=self.imsim.kwargs_numerics )
            
            print('Start simulating image')
            self.imsim.simulated_image = self.imsim.ImageModel.image(kwargs_source=LightModel_source_kwargs,
                                                                     #kwargs_ps=kwargs_ps,
                                                                     #kwargs_lens_light=kwargs_light_lens,
                                                                     kwargs_lens=self.imsim.LensModel_kwargs, unconvolved=False)
            self.imsim.image_plane_plot.setImage(self.imsim.simulated_image[::-1,:])
            
            ######### Plot the critical curve #########
            x = (self.imsim._lt_curve_coords_relative_broken[0] - self.imsim._SquareOfInterest_xr_bottomleft) / (self.imsim.fits_image.pix_deg_scale*3600)
            y = self.imsim.simulated_image.shape[0] - (self.imsim._lt_curve_coords_relative_broken[1] - self.imsim._SquareOfInterest_yr_bottomleft) / (self.imsim.fits_image.pix_deg_scale*3600)
            self.imsim.critical_curve_plot.setData(x, y)
            print('done')
            
            #solver = LensEquationSolver(self.imsim.lens_model)
            
            
            if event.modifiers() == Qt.ShiftModifier :
                print('Starting optimization')
                
                N = 40
                
                kwargs_source_fixed = []
                kwargs_source_lower = []
                kwargs_source_upper = []
                kwargs_source_sigma = []
                kwargs_source_init = []
                for i, src in enumerate(LightModel_source_list) :
                    kwargs = LightModel_source_kwargs[i]
                    kwargs_fixed = {}
                    kwargs_lower = {}
                    kwargs_upper = {}
                    kwargs_sigma = {}
                    kwargs_init = {}
                    
                    fixed_params = []#['center_x', 'center_y']                    
                    opt_params = LightModel_source_kwargs[i].keys()
                    
                    for param in fixed_params :
                        opt_params.remove(param)
                    for param in fixed_params :
                        kwargs_fixed[param] = kwargs[param]
                    for param in opt_params :
                        if param in ['center_x', 'center_y'] :
                            dpos = self.imsim.source_plane_widget.roi_list[i].sliders['dpos'].vmid
                            v = kwargs[param]
                            kwargs_lower[param] = v - dpos
                            kwargs_upper[param] = v + dpos
                            kwargs_sigma[param] = 2*dpos/N
                            kwargs_init[param] = v
                        elif param in ['e1', 'e2'] :
                            kwargs_lower[param] = -1.
                            kwargs_upper[param] = 1.
                            kwargs_sigma[param] = 2/N
                            kwargs_init[param] = kwargs[param]
                        else :
                            vmax = self.imsim.source_plane_widget.roi_list[i].sliders[param].vmax
                            vmin = self.imsim.source_plane_widget.roi_list[i].sliders[param].vmin
                            vmid = self.imsim.source_plane_widget.roi_list[i].sliders[param].vmid
                            kwargs_lower[param] = vmin
                            kwargs_upper[param] = vmax
                            kwargs_sigma[param] = (vmax + vmin)/N
                            kwargs_init[param] = vmid
                    
                    
                    kwargs_source_fixed.append(kwargs_fixed)
                    kwargs_source_lower.append(kwargs_lower)
                    kwargs_source_upper.append(kwargs_upper)
                    kwargs_source_sigma.append(kwargs_sigma)
                    kwargs_source_init.append(kwargs_init)

                
                param = Param(lm_dict['models'],
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
                
                fitting_seq = FittingSequence(kwargs_data_joint, lm_dict['models'], {}, kwargs_likelihood, self.imsim.optimization_params, mpi=False)

                fitting_kwargs_list = [['PSO', {'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 200}]]
                
                self.imsim.fitting_seq = fitting_seq
                self.imsim.chain_list = fitting_seq.fit_sequence(fitting_kwargs_list)
                self.imsim.result_kwargs = fitting_seq.best_fit()
                
                
                #for i in range(len(chain_list)):
                chain_plot.plot_chain_list(self.imsim.chain_list, 0)
                
                if False :
                    sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = self.imsim.chain_list[1]
                    print("number of non-linear parameters in the MCMC process: ", len(param_mcmc))
                    print("parameters in order: ", param_mcmc)
                    print("number of evaluations in the MCMC process: ", np.shape(samples_mcmc)[0])
                    n_sample = len(samples_mcmc)
                    print(n_sample)
                    samples_mcmc_cut = samples_mcmc[int(n_sample*1/2.):]
                    n, num_param = np.shape(samples_mcmc_cut)
                    plot = corner.corner(samples_mcmc_cut[:,:], labels=param_mcmc[:], show_titles=True)
                
                
                
                
                lm_dict = format_lm_local(lm_dict['models'], self.imsim.result_kwargs)
                
                self.imsim.lm_optimized = lenstronomy_model(lm_dict, self.imsim, chain_list=self.imsim.chain_list)
                self.imsim.lm_optimized.plot()
                self.imsim.lm_optimized.save('./last_optimized_lm.pkl')
                
            return True  # Stop propagation
        return False     # Let other events pass through  

            

class ImageFilter(QObject) :
    def __init__(self, imsim) :
        super().__init__()
        self.imsim = imsim
    def eventFilter(self, obj, event) :
        if event.type() == QEvent.KeyPress and event.key() == Qt.Key_Space :
            print("Spacebar function not yet implemented.")
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







