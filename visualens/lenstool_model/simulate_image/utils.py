import copy
import numpy as np
from ...utils.utils_Qt.utils_general import transform_ROI_params
from ...utils.utils_Qt.drag_widgets import _roi_link_box_slider
import re
from lenstronomy.Sampling.parameters import Param
from lenstronomy.LensModel.lens_model import LensModel




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


def get_ROI_slider_init_dict_light_model() :
    """
    Slider specs: ``[min, max, 'log'|'linear', handle_style, roi_link_mode]``.
    Optional 5th entry ``roi_link_mode`` selects ROI coupling for any parameter name:
    ``None``, ``'link_size'``, ``'link_ratio'``, ``'link_angle'``, ``'link_box'``.
    """
    slider_init_dict = {
        'amp': [1e-3, 1000, 'log', 'triple_handle', None],
        'n_sersic': [1, 10, 'linear', 'triple_handle', None],
        'R_sersic': [1e-4, 10, 'log', 'triple_handle', 'link_size'],
        'q': [0.0, 1.0, 'linear', 'triple_handle', 'link_ratio'],
        'phi': [-np.pi, np.pi, 'linear', 'triple_handle', 'link_angle'],
        'sigma': [1e-4, 10, 'log', 'triple_handle', 'link_size'],
        'dpos': [1e-6, 0.2, 'log', 'single_handle', 'link_box'],
    }
    return slider_init_dict
    
def get_ROI_param_dict_light_model() :
    ROI_param_dict = {
        'CIRCLE1': [ 'GAUSSIAN', ['amp', 'sigma', 'dpos'] ],
        'CIRCLE2': [ 'SERSIC', ['amp', 'R_sersic', 'n_sersic', 'dpos'] ],
        'ELLIPSE': [ 'SERSIC_ELLIPSE_Q_PHI', ['amp', 'R_sersic', 'n_sersic', 'q', 'phi', 'dpos'] ],
        'POINT': [ 'SOURCE_POSITION', ['amp', 'dpos'] ]
    }
    return ROI_param_dict

def get_ROI_slider_init_dict_lens_model() :
    """Same shape as :func:`get_ROI_slider_init_dict_light_model` (optional 5th ``roi_link_mode``)."""
    slider_init_dict = {
        'sigma0': [1e-3, 1000, 'log', 'triple_handle', None],
        'Ra': [1e-3, 100, 'log', 'triple_handle', None],
        'Rs': [1e-3, 100, 'log', 'triple_handle', 'link_size'],
        'q': [-1.0, 1.0, 'linear', 'triple_handle', 'link_ratio'],
        'phi': [-1.0, 1.0, 'linear', 'triple_handle', 'link_angle'],
        'dpos': [1e-4, 10, 'log', 'single_handle', 'link_box'],
    }
    return slider_init_dict
    
def get_ROI_param_dict_lens_model() :
    ROI_param_dict = {
        'CIRCLE1': [ 'CIRCLE1', ['param'] ],
        'CIRCLE2': [ 'CIRCLE2', ['param'] ],
        'ELLIPSE': [ 'PJAFFE_ELLIPSE_POTENTIAL_Q_PHI', ['sigma0', 'Ra', 'Rs', 'q', 'phi', 'dpos'] ],
        'POINT': [ 'POINT', ['param'] ]
    }
    return ROI_param_dict



def make_lm_dict(imsim) :

    """ Get the lens model parameters from the rgb image panel's ROIs """
    # First, reset the lens model list and kwargs to just the base interpol model
    imsim.LensModel_list = [imsim.LensModel_list[0]]
    imsim.LensModel_kwargs = [imsim.LensModel_kwargs[0]]
    # Then, add the PJAFFE_ELLIPSE_POTENTIAL_Q_PHI model for each ellipse ROI on the rgb image plane
    for roi in imsim.image_plane_rgb.roi_list :
        if roi.type=='PJAFFE_ELLIPSE_POTENTIAL_Q_PHI' :
            imsim.LensModel_list.append('PJAFFE_ELLIPSE_POTENTIAL_Q_PHI')
            x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
            # This time, the intrinsic units of the rgb panel are pixels, so we need to convert the ROI's parameters to arcsec
            npix = imsim._crop_npix
            pix_scale = imsim.fits_image.pix_deg_scale * 3600.0
            center_x = (x_center - npix / 2.0) * pix_scale
            center_y = (npix / 2.0 - y_center) * pix_scale
            q = semi_minor / semi_major
            to_add = {
                'sigma0': roi.sliders['sigma0'].vmid,
                'Ra': roi.sliders['Ra'].vmid,
                'Rs': roi.sliders['Rs'].vmid,
                'q': q,
                'phi': angle,
                'center_x': center_x,
                'center_y': center_y,
            }
            imsim.LensModel_kwargs.append(to_add)

    imsim.LensModel = LensModel(lens_model_list=imsim.LensModel_list)
    # Might be useful at some point to shift the source plane widget by the new source center coordinates for the curent lens model
    imsim.source_center_coordinates_new = imsim.LensModel.ray_shooting(0.0, 0.0, imsim.LensModel_kwargs)

    """ Get the light model parameters from the source plane widget's ROIs """
    LightModel_source_list = []
    LightModel_source_kwargs = []
    PSModel_source_list = []
    PSModel_source_kwargs = []

    for roi in imsim.source_plane_widget.roi_list :
        if roi.type=='SERSIC_ELLIPSE_Q_PHI' :
            LightModel_source_list.append('SERSIC_ELLIPSE_Q_PHI')
            amp = roi.sliders['amp'].vmid
            n_sersic = roi.sliders['n_sersic'].vmid
            x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
            q = semi_minor / semi_major
            #e = (1-q) / (1+q)
            #e1 = e * np.cos(2*angle)
            #e2 = e * np.sin(2*angle)
            R_sersic = np.sqrt(semi_major * semi_minor)
            src_xr = x_center + imsim.source_center_coordinates[0]
            src_yr = y_center + imsim.source_center_coordinates[1]
            to_add = {'amp': amp, 'R_sersic': R_sersic, 'n_sersic': n_sersic, 'q': q, 'phi': angle, 'center_x': src_xr, 'center_y': src_yr}
            LightModel_source_kwargs.append(to_add)
        elif roi.type=='GAUSSIAN' :
            LightModel_source_list.append('GAUSSIAN')
            amp = roi.sliders['amp'].vmid
            x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
            src_xr = x_center + imsim.source_center_coordinates[0]
            src_yr = y_center + imsim.source_center_coordinates[1]
            sigma = abs(roi.size()[0])
            to_add = {'amp': amp, 'sigma': sigma, 'center_x': src_xr, 'center_y': src_yr}
            LightModel_source_kwargs.append(to_add)
        elif roi.type=='SERSIC' :
            LightModel_source_list.append('SERSIC')
            amp = roi.sliders['amp'].vmid
            n_sersic = roi.sliders['n_sersic'].vmid
            x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
            src_xr = x_center + imsim.source_center_coordinates[0]
            src_yr = y_center + imsim.source_center_coordinates[1]
            R_sersic = abs(roi.size()[0]) / 2.0
            to_add = {'amp': amp, 'R_sersic': R_sersic, 'n_sersic': n_sersic, 'center_x': src_xr, 'center_y': src_yr}
            LightModel_source_kwargs.append(to_add)
        elif roi.type=='SOURCE_POSITION' :
            PSModel_source_list.append('SOURCE_POSITION')
            amp = roi.sliders['amp'].vmid
            x_center, y_center, semi_major, semi_minor, angle = transform_ROI_params(roi)
            src_xr = x_center + imsim.source_center_coordinates[0]
            src_yr = y_center + imsim.source_center_coordinates[1]
            to_add = {'source_amp': amp, 'ra_source': src_xr, 'dec_source': src_yr} # point_amp?
            PSModel_source_kwargs.append(to_add)
    
    
    ######### Create the current lenstronomy model class #########
    lm_dict = {}
    lm_dict['models'] = {'lens_model_list': imsim.LensModel_list,
                         #'lens_light_model_list': imsim.LensModel_light_list,
                         'source_light_model_list': LightModel_source_list,
                         'point_source_model_list': PSModel_source_list}
    lm_dict['kwargs'] = {'kwargs_lens': imsim.LensModel_kwargs, 
                         'kwargs_source': LightModel_source_kwargs,
                         'kwargs_source_ps': PSModel_source_kwargs}
    return lm_dict


def make_lm_dict_opt(imsim, N_sigma=6.) :
    lm_dict_opt = make_lm_dict(imsim)
    
    """ Source Light models """
    roi_indices = np.where( [roi.type != 'SOURCE_POSITION' for roi in imsim.source_plane_widget.roi_list] )[0]
    
    kwargs_source_fixed = []
    kwargs_source_lower = []
    kwargs_source_upper = []
    kwargs_source_sigma = []
    kwargs_source_init = []
    for i, src in enumerate(lm_dict_opt['models']['source_light_model_list']) :
        kwargs = lm_dict_opt['kwargs']['kwargs_source'][i]
        kwargs_fixed = {}
        kwargs_lower = {}
        kwargs_upper = {}
        kwargs_sigma = {}
        kwargs_init = {}
        
        fixed_params = []#['center_x', 'center_y']                    
        opt_params = lm_dict_opt['kwargs']['kwargs_source'][i].keys()
        
        for param in fixed_params :
            opt_params.remove(param)
        for param in fixed_params :
            kwargs_fixed[param] = kwargs[param]
        for param in opt_params :
            if param in ['center_x', 'center_y'] :
                _, sbox = _roi_link_box_slider(imsim.source_plane_widget.roi_list[roi_indices[i]])
                dpos = sbox.vmid
                v = kwargs[param]
                kwargs_lower[param] = v - dpos
                kwargs_upper[param] = v + dpos
                kwargs_sigma[param] = 2*dpos/N_sigma
                kwargs_init[param] = v
            #elif param in ['e1', 'e2'] :
            #    kwargs_lower[param] = -1.
            #    kwargs_upper[param] = 1.
            #    kwargs_sigma[param] = 2/N_sigma
            elif param in ['q', 'phi'] :
                vmax = imsim.source_plane_widget.roi_list[ roi_indices[i] ].sliders[param].vmax
                vmin = imsim.source_plane_widget.roi_list[ roi_indices[i] ].sliders[param].vmin
                kwargs_lower[param] = vmin
                kwargs_upper[param] = vmax
                kwargs_sigma[param] = (vmax - vmin) / N_sigma
                kwargs_init[param] = kwargs[param]
            else :
                vmax = imsim.source_plane_widget.roi_list[ roi_indices[i] ].sliders[param].vmax
                vmin = imsim.source_plane_widget.roi_list[ roi_indices[i] ].sliders[param].vmin
                vmid = imsim.source_plane_widget.roi_list[ roi_indices[i] ].sliders[param].vmid
                kwargs_lower[param] = vmin
                kwargs_upper[param] = vmax
                kwargs_sigma[param] = (vmax + vmin)/N_sigma
                kwargs_init[param] = vmid
        
        
        kwargs_source_fixed.append(kwargs_fixed)
        kwargs_source_lower.append(kwargs_lower)
        kwargs_source_upper.append(kwargs_upper)
        kwargs_source_sigma.append(kwargs_sigma)
        kwargs_source_init.append(kwargs_init)


    """ Point source models """
    roi_indices = np.where( [roi.type == 'SOURCE_POSITION' for roi in imsim.source_plane_widget.roi_list] )[0]
    
    kwargs_ps_fixed = []
    kwargs_ps_lower = []
    kwargs_ps_upper = []
    kwargs_ps_sigma = []
    kwargs_ps_init = []
    for i, src in enumerate(lm_dict_opt['models']['point_source_model_list']) :
        kwargs = lm_dict_opt['kwargs']['kwargs_source_ps'][i]
        kwargs_fixed = {}
        kwargs_lower = {}
        kwargs_upper = {}
        kwargs_sigma = {}
        kwargs_init = {}
        
        fixed_params = []#['center_x', 'center_y']                    
        opt_params = lm_dict_opt['kwargs']['kwargs_source_ps'][i].keys()
        
        for param in fixed_params :
            opt_params.remove(param)
        for param in fixed_params :
            kwargs_fixed[param] = kwargs[param]
        for param in opt_params :
            if param in ['ra_source', 'dec_source'] :
                _, sbox = _roi_link_box_slider(imsim.source_plane_widget.roi_list[roi_indices[i]])
                dpos = sbox.vmid
                v = kwargs[param]
                kwargs_lower[param] = v - dpos
                kwargs_upper[param] = v + dpos
                kwargs_sigma[param] = 2*dpos/N_sigma
                kwargs_init[param] = v
            else : # Here only source_amp
                param_slider = 'amp' if param=='source_amp' else param
                vmax = imsim.source_plane_widget.roi_list[ roi_indices[i] ].sliders[param_slider].vmax
                vmin = imsim.source_plane_widget.roi_list[ roi_indices[i] ].sliders[param_slider].vmin
                vmid = imsim.source_plane_widget.roi_list[ roi_indices[i] ].sliders[param_slider].vmid
                kwargs_lower[param] = vmin
                kwargs_upper[param] = vmax
                kwargs_sigma[param] = (vmax + vmin)/N_sigma
                kwargs_init[param] = vmid
        
        
        kwargs_ps_fixed.append(kwargs_fixed)
        kwargs_ps_lower.append(kwargs_lower)
        kwargs_ps_upper.append(kwargs_upper)
        kwargs_ps_sigma.append(kwargs_sigma)
        kwargs_ps_init.append(kwargs_init)
        
    
    if False : # This produces an error if there is a SERSIC_ELLIPSE_Q_PHI light model
        param = Param(lm_dict_opt['models'],
                      kwargs_fixed_lens=imsim.LensModel_kwargs,
                      kwargs_fixed_source=kwargs_source_fixed,#self.imsim.fixed_source_kwargs,
                      #kwargs_fixed_lens_light=kwargs_fixed_lens_light,
                      kwargs_fixed_ps=kwargs_ps_fixed, 
                      kwargs_lower_lens=[{}],
                      kwargs_lower_source=kwargs_source_lower,
                      #kwargs_lower_lens_light=kwargs_lower_lens_light,
                      kwargs_lower_ps=kwargs_ps_lower,
                      kwargs_upper_lens=[{}],
                      kwargs_upper_source=kwargs_source_upper,
                      #kwargs_upper_lens_light=kwargs_upper_lens_light,
                      kwargs_upper_ps=kwargs_ps_upper,
                      #kwargs_lens_init=kwargs_lens,
                      #joint_lens_with_light: [[0, 0, ['center_x', 'center_y']]],
                      linear_solver=imsim.use_linear_solver)
        param.print_setting()
    
    lm_dict_opt['kwargs_fixed'] = {'kwargs_source': kwargs_source_fixed, 'kwargs_source_ps': kwargs_ps_fixed}
    lm_dict_opt['kwargs_lower'] = {'kwargs_source': kwargs_source_lower, 'kwargs_source_ps': kwargs_ps_lower}
    lm_dict_opt['kwargs_upper'] = {'kwargs_source': kwargs_source_upper, 'kwargs_source_ps': kwargs_ps_upper}
    lm_dict_opt['kwargs_sigma'] = {'kwargs_source': kwargs_source_sigma, 'kwargs_source_ps': kwargs_ps_sigma}
    lm_dict_opt['kwargs_init'] = {'kwargs_source': kwargs_source_init, 'kwargs_source_ps': kwargs_ps_init}
    
    return lm_dict_opt


def make_samples_dict(lm_dict_opt, samples, linear_solver=False) :
    param_list_full = samples[2]
    values = np.array( samples[1] ).T
    
    samples_dict = {'kwargs_source': [{} for model in lm_dict_opt['models']['source_light_model_list']],
                    'kwargs_source_ps': [{} for model in lm_dict_opt['models']['point_source_model_list']]}
    
    # This function just to match the parameter names used in the samples to those used in the lm_dict
    start = 0
    match_dict = {'amp_source_light': 'amp',
                  'R_sersic_source_light': 'R_sersic',
                  'n_sersic_source_light': 'n_sersic',
                  'sigma_source_light': 'sigma',
                  'center_x_source_light': 'center_x',
                  'center_y_source_light': 'center_y',
                  'e1_source_light': 'e1',
                  'e2_source_light': 'e2',
                  'q_source_light': 'q',
                  'phi_source_light': 'phi',
                  'source_amp': 'amp',
                  'ra_source': 'center_x',
                  'dec_source': 'center_y'}
    for i, model in enumerate( lm_dict_opt['models']['source_light_model_list'] ) :
        if model == 'GAUSSIAN' :
            n=3 if linear_solver else 4
        elif model == 'SERSIC' :
            n=4 if linear_solver else 5
        elif model in ('SERSIC_ELLIPSE', 'SERSIC_ELLIPSE_Q_PHI') :
            n=6 if linear_solver else 7
        param_list = param_list_full[ start:start+n ]
        for k, param in enumerate(param_list) :
            m = re.search(r'(\d+)$', param)
            if m :
                j = m.group(1)
                if int(j) != int(i) :
                    print('Debug --> attention, i != j: ', i, j)
                param_new = match_dict[param[:-len(j)]]
            else :
                param_new = match_dict[param]
            samples_dict['kwargs_source'][i][param_new] = np.array( values[start+k] )
        start = start+n
    
    for i, model in enumerate( lm_dict_opt['models']['point_source_model_list'] ) :
        n=2 if linear_solver else 3
        param_list = param_list_full[ start:start+n ]
        for k, param in enumerate(param_list) :
            param_new = match_dict[param]
            samples_dict['kwargs_source_ps'][i][param_new] = np.array( values[start+k] )
        start = start+n
        
    return samples_dict


