import numpy as np



def get_ROI_slider_init_dict_light_model() :
    """
    Slider specs: ``[min, max, 'log'|'linear', handle_style, roi_link_mode]``.
    Optional 5th entry ``roi_link_mode`` selects ROI coupling for any parameter name:
    ``None``, ``'link_size'``, ``'link_ratio'``, ``'link_angle'``, ``'link_box'``.
    """
    slider_init_dict = {
        'amp': [1e-3, 1000, 'log', 'triple_handle', None],
        'n_sersic': [0, 10, 'linear', 'triple_handle', None],
        'R_sersic': [1e-4, 10, 'log', 'triple_handle', 'link_size'],
        'q': [0.0, 1.0, 'linear', 'triple_handle', 'link_ratio'],
        'phi': [0., np.pi, 'linear', 'triple_handle', 'link_angle'],
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
        #'sigma0': [1e-3, 1000, 'log', 'triple_handle', None],
        'sigma0': [0.03, 100, 'log', 'triple_handle', None],
        #'Ra': [1e-4, 100, 'log', 'triple_handle', None],
        'Ra': [3e-6, 0.03, 'log', 'triple_handle', None],
        #'Rs': [1e-3, 100, 'log', 'triple_handle', 'link_size'],
        'Rs': [0.03, 100, 'log', 'triple_handle', 'link_size'],
        'q': [0.0, 1.0, 'linear', 'triple_handle', 'link_ratio'],
        'phi': [-np.pi/2, np.pi/2, 'linear', 'triple_handle', 'link_angle'],
        'dpos': [1e-5, 1, 'log', 'single_handle', 'link_box'],
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


