import copy




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
    
    


