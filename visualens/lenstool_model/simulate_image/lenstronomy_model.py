import numpy as np
import os
import matplotlib.pyplot as plt
import pickle
import copy
from lenstronomy.Plots.model_plot import ModelPlot

from ...utils.utils_astro.utils_general import world_to_relative, relative_to_world
from ...utils.utils_Qt.drag_widgets import SelectableCircleROI, SelectableEllipseROI
from ...utils.utils_Qt.sliders import TripleSlider



class lenstronomy_model :
    def __init__(self, lm_dict_or_str, imsim) :
        # lm is a dictionary containing the parameters of a source model, either in local or world coordinates.
        self.imsim = imsim
        if type(lm_dict_or_str)==str :
            lm_dict_or_str = self.load_lm(path=lm_dict_or_str)
        if 'lens_model_list' in lm_dict_or_str['models'] :
            self.local = lm_dict_or_str
            self.world = self.to_world()
        else :
            self.world = lm_dict_or_str
            self.local = self.to_local()
    
    def plot(self) :        
        self.local['kwargs']['kwargs_lens_light'] = []
        self.local['kwargs']['kwargs_ps'] = []
        
        self.modelPlot = ModelPlot([[self.imsim.ImageData_kwargs, self.imsim.PSF_kwargs, self.imsim.kwargs_numerics]], 
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
            a2 = self.modelPlot.source_plot(ax=a2, deltaPix_source=0.0005, numPix=1500, cmap='gray')
            source_array = a2.get_children()[0].get_array()
            
            fig, ax = plt.subplots()
            #ax.imshow( source_array-np.min(source_array), origin="lower", cmap='gray' )
            ax.imshow( source_array, vmax=np.log(0.06), origin='lower', cmap='gray' )
            
            f3, a3 = plt.subplots()
            a3.imshow( np.exp(source_array), vmax=0.06, origin="lower", cmap='gray' )
            
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
        
        with open(path, 'wb') as file :
            pickle.dump(self.world, file)
        print('Model saved at ' + path)


    def load_lm(self, path='./lenstronomy_model.pkl') :
        with open(path, 'rb') as file :
            imported_model_world = pickle.load(file)
        
        if 'results' in imported_model_world :
            imported_model_world['kwargs'] = copy.deepcopy(imported_model_world['results'])
            del imported_model_world['results']
        
        #self.__init__(self.world, self.imsim)
        #self.local = self.to_local()
        return imported_model_world
    
    
    def remove_source(self, index) :
        del self.local['models']['source_light_model_list'][index]
        del self.local['kwargs']['kwargs_source'][index]
    
    
    def send_to_imsim(self) :
        PlotWidget = self.imsim.source_plane_widget
        for i, model in enumerate(self.local['models']['source_light_model_list']) :
            if model=='GAUSSIAN' :
                kwargs = self.local['kwargs']['kwargs_source'][i]
                x = kwargs['center_x'] - self.imsim.source_center_coordinates[0]
                y = kwargs['center_y'] - self.imsim.source_center_coordinates[1]
                sigma = kwargs['sigma']
                x_corner = x - sigma
                y_corner = y - sigma
                
                PlotWidget.last_roi = SelectableCircleROI([x_corner, y_corner], PlotWidget=PlotWidget, radius=sigma, pen='g', invertible=True)
                PlotWidget.last_roi.type = 2
                PlotWidget.roi_list.append(PlotWidget.last_roi)
                for handle in PlotWidget.last_roi.handles:
                    PlotWidget.last_roi.removeHandle(handle['item'])
                PlotWidget.addItem(PlotWidget.last_roi)
                                
                # Add sliders to control light source parameters
                param = 'amp'
                min_range, max_range = PlotWidget.ranges_dict[param][0], PlotWidget.ranges_dict[param][1]
                PlotWidget.last_roi.sliders = {param : TripleSlider(min_range, max_range, PlotWidget=PlotWidget, label=param, roi=self.last_roi)}
                PlotWidget.last_roi.sliders[param].vmid = kwargs[param]
                PlotWidget.last_roi.sliders[param].setParentItem(PlotWidget.last_roi)
        
            if model=='SERSIC' :
                kwargs = self.local['kwargs']['kwargs_source'][i]
                x = kwargs['center_x'] - self.imsim.source_center_coordinates[0]
                y = kwargs['center_y'] - self.imsim.source_center_coordinates[1]
                R_sersic = kwargs['R_sersic']
                x_corner = x - R_sersic
                y_corner = y - R_sersic
                
                PlotWidget.last_roi = SelectableCircleROI([x_corner, y_corner], PlotWidget=PlotWidget, radius=R_sersic, pen='b', invertible=True)
                PlotWidget.last_roi.type = 3
                PlotWidget.roi_list.append(PlotWidget.last_roi)
                for handle in PlotWidget.last_roi.handles:
                    PlotWidget.last_roi.removeHandle(handle['item'])
                PlotWidget.addItem(PlotWidget.last_roi)
                                
                # Add sliders to control light source parameters
                params = ['amp', 'n_sersic']
                for i, param in enumerate(params) : #, 'R_sersic'
                    min_range, max_range = PlotWidget.ranges_dict[param][0], PlotWidget.ranges_dict[param][1]
                    if not hasattr(PlotWidget.last_roi, 'sliders') :
                        PlotWidget.last_roi.sliders = {param : TripleSlider(min_range, max_range, PlotWidget=PlotWidget, label=param, roi=self.last_roi)}
                    else :
                        offset = PlotWidget.last_roi.sliders[params[i-1]].offset + PlotWidget.last_roi.sliders[params[i-1]].bounding_height + PlotWidget.last_roi.sliders[params[i-1]].offset_text
                        PlotWidget.last_roi.sliders[param] = TripleSlider(min_range, max_range, PlotWidget=PlotWidget, label=param, roi=self.last_roi, offset=offset)
                    PlotWidget.last_roi.sliders[param].vmid = kwargs[param]
                    PlotWidget.last_roi.sliders[param].setParentItem(PlotWidget.last_roi)
                
            if model=='SERSIC_ELLIPSE' :
                kwargs = self.local['kwargs']['kwargs_source'][i]
                x = kwargs['center_x'] - self.imsim.source_center_coordinates[0]
                y = kwargs['center_y'] - self.imsim.source_center_coordinates[1]
                R_sersic = kwargs['R_sersic']
                e1, e2 = kwargs['e1'], kwargs['e2']
                x_corner = x - R_sersic # check if this works with ellipticity...
                y_corner = y - R_sersic
                
                PlotWidget.last_roi = SelectableEllipseROI([x_corner, y_corner], [e1, e2], PlotWidget=PlotWidget, radius=R_sersic, pen='r', invertible=True)
                PlotWidget.last_roi.type = 1
                PlotWidget.roi_list.append(PlotWidget.last_roi)
                for handle in PlotWidget.last_roi.handles:
                    PlotWidget.last_roi.removeHandle(handle['item'])
                PlotWidget.addItem(PlotWidget.last_roi)
                                
                # Add sliders to control light source parameters
                params = ['amp', 'n_sersic']
                for i, param in enumerate(params) : #, 'R_sersic'
                    min_range, max_range = PlotWidget.ranges_dict[param][0], PlotWidget.ranges_dict[param][1]
                    if not hasattr(PlotWidget.last_roi, 'sliders') :
                        PlotWidget.last_roi.sliders = {param : TripleSlider(min_range, max_range, PlotWidget=PlotWidget, label=param, roi=self.last_roi)}
                    else :
                        offset = PlotWidget.last_roi.sliders[params[i-1]].offset + PlotWidget.last_roi.sliders[params[i-1]].bounding_height + PlotWidget.last_roi.sliders[params[i-1]].offset_text
                        PlotWidget.last_roi.sliders[param] = TripleSlider(min_range, max_range, PlotWidget=PlotWidget, label=param, roi=self.last_roi, offset=offset)
                    PlotWidget.last_roi.sliders[param].vmid = kwargs[param]
                    PlotWidget.last_roi.sliders[param].setParentItem(PlotWidget.last_roi)
                    
                    
                    
                    