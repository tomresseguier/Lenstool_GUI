import copy
import numpy as np
import re




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



def make_samples_dict(lm_dict_opt, samples, linear_solver=False) :
    param_list_full = samples[2]
    values = np.array( samples[1] ).T
    
    samples_dict = {'kwargs_lens': [{} for model in lm_dict_opt['models']['lens_model_list']],
                    'kwargs_source': [{} for model in lm_dict_opt['models']['source_light_model_list']],
                    'kwargs_source_ps': [{} for model in lm_dict_opt['models']['point_source_model_list']]}
    
    # This function just to match the parameter names used in the samples to those used in the lm_dict
    start = 0
    match_dict = {'sigma0_lens': 'sigma0',
                  'Ra_lens': 'Ra',
                  'Rs_lens': 'Rs',
                  'q_lens': 'q',
                  'phi_lens': 'phi',
                  'center_x_lens': 'center_x',
                  'center_y_lens': 'center_y',
                  'amp_source_light': 'amp',
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
    
    for i, model in enumerate( lm_dict_opt['models']['lens_model_list'] ) :
        if model == 'INTERPOL' :
            continue
        if model == 'PJAFFE_ELLIPSE_POTENTIAL_Q_PHI' :
            n=7
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
            samples_dict['kwargs_lens'][i][param_new] = np.array( values[start+k] )
        start = start+n
        
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


