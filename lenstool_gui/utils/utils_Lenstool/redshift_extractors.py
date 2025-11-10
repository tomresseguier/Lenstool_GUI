import os
import numpy as np
import pylenstool
import pandas as pd



def find_param_file(run_dir) :
    param_file_name = None
    for file_name in os.listdir(run_dir) :
        if file_name.endswith(".par") and not file_name.startswith("best") :
            print(file_name + ' found')
            param_file_name = file_name
    return param_file_name

def make_fixed_z_dict(run_dir) :
    param_file_name = find_param_file(run_dir)
    param_file_path = os.path.join(run_dir, param_file_name)
    
    fixed_z_dict = {}
    lenstool_file = pylenstool.lenstool_param_file(param_file_path)
    param_list = lenstool_file.get_parameter('image', 'z_m_limit')
    if type(param_list)==list :
        param_array = np.array(param_list)
        key_indices = np.where(param_array=='z_m_limit')[0]
        opt_indices = key_indices + 1
        image_name_indices = key_indices + 2
        fixed_z_image_name_indices = image_name_indices[np.where(param_array[opt_indices]=='0')[0]]
        fixed_z_image_names = param_array[fixed_z_image_name_indices]
        redshifts = param_array[fixed_z_image_name_indices + 2]
        for i, fixed_z_image_name in enumerate(fixed_z_image_names) :
            fixed_z_dict[fixed_z_image_name] = redshifts[i]
    return fixed_z_dict

def make_opt_z_dict(run_dir, best_file_name='best.par') :
    with open( os.path.join(run_dir, best_file_name) ) as best_par_file :
        lines = best_par_file.readlines()
        
    opt_redshift_key = '\tz_m_limit 1 '
    opt_redshifts_indexes = np.where( [line.startswith(opt_redshift_key) for line in lines] )[0]
    opt_redshifts_dict = {}
    #start = len('\tz_m_limit 1 A1 0 ')
    #end = len('\tz_m_limit 1 A1 0 1.478015')
    for idx in opt_redshifts_indexes :
        source_name = lines[idx].split()[2]
        print("source_name:", source_name)
        #source_name = lines[idx][len(opt_redshift_key):len(opt_redshift_key)+2]
        opt_redshift = float( lines[idx].split()[-3] )
        print("opt_redshift:", opt_redshift)
        #opt_redshift = float( lines[idx][start:end] )
        opt_redshifts_dict[source_name] = opt_redshift
    return opt_redshifts_dict

def make_source_z_dict(run_dir, use_family_name_only=False) :
    fixed_z_dict = make_fixed_z_dict(run_dir)
    opt_z_dict = make_opt_z_dict(run_dir)
    source_z_dict = {**fixed_z_dict, **opt_z_dict}
    if use_family_name_only :
        source_z_dict_bis = {}
        for key in source_z_dict.keys() :
            source_z_dict_bis[key[0]] = source_z_dict[key]
        source_z_dict = source_z_dict_bis
    return source_z_dict


















