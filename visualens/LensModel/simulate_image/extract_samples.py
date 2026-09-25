import numpy as np
import re




def make_samples_dict(lm_dict_opt, samples, linear_solver=False) :
    if samples[0] in ('emcee', 'MCMC') :
        values = np.array( samples[1] ).T
    elif samples[0]=='PSO' :
        values = np.array( samples[1][1] ).T
    else :
        raise ValueError(f"Unknown sampling method: {samples[0]}")
    param_list_full = samples[2]

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
                  'source_amp': 'source_amp', # amp
                  'ra_source': 'ra_source', # center_x
                  'dec_source': 'dec_source'} # center_y
    
    for i, model in enumerate( lm_dict_opt['models']['lens_model_list'] ) :
        if model == 'INTERPOL' :
            continue
        if model == 'PJAFFE_ELLIPSE_POTENTIAL_Q_PHI' :
            n_params_fixed = len(lm_dict_opt['kwargs_fixed']['kwargs_lens'][i])
            n_params_sampled = 7 - n_params_fixed
        param_list = param_list_full[ start:start+n_params_sampled ]
        for k, param in enumerate(param_list) :
            m = re.search(r'(\d+)$', param) # Finds if the parameter name contains an index
            if m :
                j = m.group(1) # Finds the index from the parameter name
                if int(j) != int(i) :
                    print('Debug --> attention, i != j: ', i, j)
                param_new = match_dict[param[:-len(j)]]
            else :
                param_new = match_dict[param]
            samples_dict['kwargs_lens'][i][param_new] = np.array( values[start+k] )
        start = start+n_params_sampled
        
    for i, model in enumerate( lm_dict_opt['models']['source_light_model_list'] ) :
        if model == 'GAUSSIAN' :
            n_params = 3 if linear_solver else 4
            n_params_fixed = len(lm_dict_opt['kwargs_fixed']['kwargs_source'][i])
            n_params_sampled = n_params - n_params_fixed
        elif model == 'SERSIC' :
            n_params = 4 if linear_solver else 5
            n_params_fixed = len(lm_dict_opt['kwargs_fixed']['kwargs_source'][i])
            n_params_sampled = n_params - n_params_fixed
        elif model in ('SERSIC_ELLIPSE', 'SERSIC_ELLIPSE_Q_PHI') :
            n_params = 6 if linear_solver else 7
            n_params_fixed = len(lm_dict_opt['kwargs_fixed']['kwargs_source'][i])
            n_params_sampled = n_params - n_params_fixed
        param_list = param_list_full[ start:start+n_params_sampled ]
        for k, param in enumerate(param_list) :
            m = re.search(r'(\d+)$', param) # Finds if the parameter name contains an index
            if m :
                j = m.group(1) # Finds the index from the parameter name
                if int(j) != int(i) :
                    print('Debug --> attention, i != j: ', i, j)
                param_new = match_dict[param[:-len(j)]]
            else :
                param_new = match_dict[param]
            samples_dict['kwargs_source'][i][param_new] = np.array( values[start+k] )
        start = start+n_params_sampled
    
    for i, model in enumerate( lm_dict_opt['models']['point_source_model_list'] ) :
        n_params = 2 if linear_solver else 3
        n_params_fixed = len(lm_dict_opt['kwargs_fixed']['kwargs_source_ps'][i])
        n_params_sampled = n_params - n_params_fixed
        param_list = param_list_full[ start:start+n_params_sampled ]
        for k, param in enumerate(param_list) :
            param_new = match_dict[param]
            samples_dict['kwargs_source_ps'][i][param_new] = np.array( values[start+k] )
        start = start+n_params_sampled
        
    return samples_dict


