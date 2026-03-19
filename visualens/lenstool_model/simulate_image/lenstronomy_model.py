import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
import copy
from lenstronomy.Plots.model_plot import ModelPlot
from lenstronomy.Plots import chain_plot

from ...utils.utils_astro.utils_general import world_to_relative, relative_to_world
from ...utils.utils_Qt.drag_widgets import SelectableCircleROI, SelectableEllipseROI, attach_sliders
from ...utils.utils_Qt.sliders import TripleSlider
from ...utils.utils_Qt.utils_general import transform_ROI_params_inverse, make_handles



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
    
    def plot(self) :        
        self.local['kwargs']['kwargs_lens_light'] = []
        self.local['kwargs']['kwargs_ps'] = []
        
        simulated_image_kwargs = self.imsim.ImageData_kwargs.copy()
        simulated_image_kwargs['image_data'] = self.imsim.simulated_image
        self.modelPlot = ModelPlot([[simulated_image_kwargs, self.imsim.PSF_kwargs, self.imsim.kwargs_numerics]], 
                                   self.local['models'], self.local['kwargs'], arrow_size=0.02, cmap_string="gist_heat",
                                   linear_solver=True)
                                   #linear_solver=False)
        
        if True :
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
        
        
        if self.chain_list is not None :
            for i in range(len(self.chain_list)):
                chain_plot.plot_chain_list(self.chain_list, 0)
        
        #return ax, source_array
        return source_array
        
    
    def to_world(self) :
        model_world = copy.deepcopy(self.local)
        if 'lens_model_list' in model_world['models'] :
            del model_world['models']['lens_model_list']
        if 'kwargs_lens' in model_world['kwargs'] :
            del model_world['kwargs']['kwargs_lens']
        for src in model_world['kwargs']['kwargs_source'] :
            if 'center_x' in src and 'center_y' in src :
                ra, dec = relative_to_world(src['center_x'], src['center_y'], self.imsim.center_world)
                src['center_ra'] = ra
                src['center_dec'] = dec
                del src['center_x']
                del src['center_y']
        return model_world


    def to_local(self) :        
        model_local = copy.deepcopy(self.world)
        model_local['models']['lens_model_list'] = self.imsim.LensModel_list
        model_local['kwargs']['kwargs_lens'] = self.imsim.LensModel_kwargs
        
        for src in model_local['kwargs']['kwargs_source'] :
            if 'center_ra' in src and 'center_dec' in src :
                xr, yr = world_to_relative(src['center_ra'], src['center_dec'], self.imsim.center_world)
                src['center_x'] = xr
                src['center_y'] = yr
                del src['center_ra']
                del src['center_dec']
        return model_local


    def add(self, model_2) :
        if type(model_2)==str :
            model_2 = self.load_lm(path=model_2)
        if 'lens_model_list' in model_2['models'] :
            self.local['models']['source_light_model_list'] += model_2['models']['source_light_model_list']
            self.local['kwargs']['kwargs_source'] += model_2['kwargs']['kwargs_source']
            #self.__init__(self.local, self.imsim)
            self.world = self.to_world()
            return self.local
        else :
            self.world['models']['source_light_model_list'] += model_2['models']['source_light_model_list']
            self.world['kwargs']['kwargs_source'] += model_2['kwargs']['kwargs_source']
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
    
    
    def remove_source(self, index) :
        del self.local['models']['source_light_model_list'][index]
        del self.local['kwargs']['kwargs_source'][index]
        self.world = self.to_world()
    
    
    def send_to_imsim(self) :
        PlotWidget = self.imsim.source_plane_widget
        for i, model in enumerate(self.local['models']['source_light_model_list']) :
            
            if model=='SERSIC' :
                kwargs = self.local['kwargs']['kwargs_source'][i]
                x = kwargs['center_x'] - self.imsim.source_center_coordinates[0]
                y = kwargs['center_y'] - self.imsim.source_center_coordinates[1]
                radius = kwargs['R_sersic']
                x_corner = x - radius
                y_corner = y - radius
                
                PlotWidget.last_roi = SelectableCircleROI([x_corner, y_corner], radius=radius, pen='b', invertible=True)
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
                
            if model=='SERSIC_ELLIPSE' :
                kwargs = self.local['kwargs']['kwargs_source'][i]
                x_center = kwargs['center_x'] - self.imsim.source_center_coordinates[0]
                y_center = kwargs['center_y'] - self.imsim.source_center_coordinates[1]
                R_sersic = kwargs['R_sersic']
                e1, e2 = kwargs['e1'], kwargs['e2']
                
                e = (e1**2 + e2**2)**0.5
                q = (1-e) / (1+e)
                a = 2*R_sersic / (1+q)
                b = a*q
                theta = 0.5*np.arctan(e2/e1)
                
                x_corner, y_corner, a, b, angle = transform_ROI_params_inverse(x_center, y_center, a, b, theta)
                
                PlotWidget.last_roi = SelectableEllipseROI([x_corner, y_corner], [a, b], angle=angle, pen='r', invertible=True)
                PlotWidget.last_roi.type = 1
                PlotWidget.last_roi.type_str = 'SERSIC_ELLIPSE'
                PlotWidget.roi_list.append(PlotWidget.last_roi)
                for handle in PlotWidget.last_roi.handles:
                    PlotWidget.last_roi.removeHandle(handle['item'])
                PlotWidget.addItem(PlotWidget.last_roi)
                
                make_handles(PlotWidget.last_roi)
                
                # Add sliders to control light source parameters
                params = ['amp', 'n_sersic']
                attach_sliders(PlotWidget, params)
                
                for p, v in kwargs.items() :
                    if p in PlotWidget.last_roi.sliders :
                        PlotWidget.last_roi.sliders[p].vmid = v
            
            if model=='GAUSSIAN' :
                kwargs = self.local['kwargs']['kwargs_source'][i]
                x = kwargs['center_x'] - self.imsim.source_center_coordinates[0]
                y = kwargs['center_y'] - self.imsim.source_center_coordinates[1]
                radius = kwargs['sigma'] / 2.0
                x_corner = x - radius
                y_corner = y - radius
                
                PlotWidget.last_roi = SelectableCircleROI([x_corner, y_corner], radius=radius, pen='g', invertible=True)
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




def make_model_path(path) :
    if not os.path.isdir( os.path.dirname(path) ) :
        path = os.path.join('./', path)
        extension = '.pkl'
        if not extension in path :
            path += extension
    return path
