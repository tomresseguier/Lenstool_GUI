import numpy as np
import subprocess
import shutil
import sys
import os


from ..utils_astro.get_cosmology import get_cosmo
from .redshift_extractors import make_source_z_dict
cosmo = get_cosmo()







def run_Lenstool_on_best(run_dir, best_file_name='best.par', env_name='lenstool_env') : #Planck_bullet_cluster
    cwd = os.getcwd()
    os.chdir(run_dir)
    #subprocess.run('conda init zsh', shell=True)
    #subprocess.run('conda activate ' + env_name, shell=True)
    #subprocess.run('lenstool ' + best_file_name + ' -n', shell=True)
    subprocess.run(['conda', 'run', '-n', env_name, 'lenstool', best_file_name, '-n'])
    os.chdir(cwd)


def make_magnifications_and_curves(run_dir, env_name='lenstool_env') :
    cwd = os.getcwd()
    os.chdir(run_dir)
    magnification_dir = './' + 'magnifications/'
    if not os.path.exists(magnification_dir) :
        os.mkdir(magnification_dir)
    curves_dir = './' + 'curves/'
    if not os.path.exists(curves_dir) :
        os.mkdir(curves_dir)
    
    if 'best_maps.par' not in os.listdir('./') or 'best_predictions.par' not in os.listdir('./') :
        best_files_maker('./')
    
    source_z_dict = make_source_z_dict('./')
    
    ############# make magnification maps #############
    with open('./best_maps.par') as best_par_file :
        lines = best_par_file.readlines()
    
    ampli_index = np.where( [line.startswith('\tampli') for line in lines] )[0][0]
    for source in source_z_dict.keys() :
        magnification_file_name = 'magnification_' + source + '.fits'
        
        redshift_cut_length = len(lines[ampli_index].split()[-2])
        ampli_str_cut_length = len(lines[ampli_index].split()[-1] + '\n')
        start = ampli_str_cut_length + redshift_cut_length + 1
        stop = ampli_str_cut_length + 1
        
        lines[ampli_index] = lines[ampli_index][:-start] + str(source_z_dict[source]) + lines[ampli_index][-stop:]
        lines[ampli_index] = lines[ampli_index][:-ampli_str_cut_length] + magnification_file_name + '\n'
                
        with open('./best_maps.par', 'w') as file:
            file.writelines(lines)
        
        subprocess.run(['conda', 'run', '-n', env_name, 'lenstool', 'best_maps.par', '-n'])
        #subprocess.run('conda activate Planck_bullet_cluster', shell=True)
        #subprocess.run('lenstool best_maps.par -n', shell=True)
        for file_name in os.listdir('./') :
            if file_name.startswith("magnification_") :
                shutil.move(os.path.join('./', file_name), os.path.join(magnification_dir, file_name))
            
    ############# make curves #############
    with open('./best_predictions.par') as best_par_file :
        lines = best_par_file.readlines()
    
    cline_index = np.where( [line.startswith('cline') for line in lines] )[0][0]
    curve_plan_index = np.where( [line.startswith('\tnplan') for line in lines] )[0][0]


    limitHigh_idx = np.where( ['limitHigh' in line.split() for line in lines] )[0]
    if len(limitHigh_idx)>0 :
        if float(lines[limitHigh_idx[0]].split()[1]) > 0.5 :
            lines[limitHigh_idx[0]] = '\tlimitHigh   0.5\n'
    else :
        lines.insert(curve_plan_index+1, '\tlimitHigh   0.5\n')
        
    limitLow_idx = np.where( ['limitLow' in line.split() for line in lines] )[0]
    if len(limitLow_idx)>0 :
        if float(lines[limitLow_idx[0]].split()[1]) > 0.1 :
            lines[limitLow_idx[0]] = '\tlimitLow   0.1\n'
    else :
        lines.insert(curve_plan_index+1, '\tlimitLow   0.1\n')
    
    step_idx = np.where( ['step' in line.split() for line in lines] )[0]
    if len(step_idx)>0 :
        print(lines[step_idx].split()[1])
        #if float(lines[step_idx].split()[1]) > 0.1 :
        #    lines[limitLow_idx] = '\tstep   0.1\n'
    else :
        lines.insert(curve_plan_index+1, '\tstep   0.1\n')
    
    
    for source in source_z_dict.keys() :
        lines[curve_plan_index] = lines[curve_plan_index][:len('\tnplan    1 ')] + str(source_z_dict[source]) + '\n'
        
        with open('./best_predictions.par', 'w') as file:
            file.writelines(lines)
        
        subprocess.run('lenstool best_predictions.par -n', shell=True)
        for file_name in os.listdir('./') :
            if file_name.startswith("ce.dat") :
                shutil.move(os.path.join('./', file_name), os.path.join(curves_dir, 'ce_' + source + '.dat'))
            if file_name.startswith("ci.dat") :
                shutil.move(os.path.join('./', file_name), os.path.join(curves_dir, 'ci_' + source + '.dat'))
    
    os.chdir(cwd)
    print('#################')
    

def best_files_maker(run_dir, run_Lenstool=True, best_file_name='best.par') :
    file_path = os.path.join(run_dir, best_file_name)
    with open(file_path, 'r') as file:
            lines = file.readlines()
    field_index = np.where( [ line.startswith('field') for line in lines ] )[0][0]
    maps_keywords_idx = np.where( [ line.startswith('\tmass') for line in lines ] )[0].tolist() + \
                        np.where( [ line.startswith('\tampli') for line in lines ] )[0].tolist() + \
                        np.where( [ line.startswith('\tshear') for line in lines ] )[0].tolist()
    predictions_keywords_idx = np.where( [ line.startswith('\timage') for line in lines ] )[0].tolist()
    
    #field_replacement_lines_predictions = ['\txmin -65.\n', '\txmax 40.\n', '\tymin -30.\n', '\tymax 30.\n']
    field_replacement_lines_predictions = ['\txmin -120.\n', '\txmax 80.\n', '\tymin -100.\n', '\tymax 100.\n']
    #field_replacement_lines_predictions = ['\txmin -100.\n', '\txmax 60.\n', '\tymin -50.\n', '\tymax 50.\n']
    field_replacement_lines_maps = ['\txmin -120.\n', '\txmax 80.\n', '\tymin -100.\n', '\tymax 100.\n']
    
    lines_predictions = lines.copy()
    for i in range(4) :
        lines_predictions[field_index+1+i] = field_replacement_lines_predictions[i]
    for i in maps_keywords_idx :
        lines_predictions.remove(lines[i])
    best_predictions_path = run_dir + 'best_predictions.par'
    with open(best_predictions_path, 'w') as file:
        file.writelines(lines_predictions)
    
    lines_maps = lines.copy()
    for i in range(4) :
        lines_maps[field_index+1+i] = field_replacement_lines_maps[i]
    for i in predictions_keywords_idx :
        lines_maps.remove(lines[i])    
    best_maps_path = run_dir + 'best_maps.par'
    with open(best_maps_path, 'w') as file:
        file.writelines(lines_maps)
    
    if run_Lenstool :
        run_Lenstool_on_best(run_dir, best_file_name='best_maps.par')
        run_Lenstool_on_best(run_dir, best_file_name='best_predictions.par')
    # TO DO: create 2 best.par files for maps (big field) and predictions (small)
    
    









def insert_lines(input_path, n):
    input_path = os.path.join(input_path, 'best.par')
    # Format n as string (with dot preserved)
    n_str = str(n)

    # New lines to insert
    new_lines = [
        f"\tampli\t   1 1000 {n_str} magnification_z{n_str}.fits\n",
        f"\tdpl\t   1 1000 {n_str} dx_z{n_str}.fits dy_z{n_str}.fits\n"
    ]

    # Read input file
    with open(input_path, 'r') as file:
        lines = file.readlines()

    # Find where to insert (after 'runmode')
    for i, line in enumerate(lines):
        if line.strip().startswith("runmode"):
            insert_index = i + 1
            break
    else:
        raise ValueError("No line starting with 'runmode' found.")

    # Insert new lines
    lines[insert_index:insert_index] = new_lines

    # Construct output path: insert _z{n} before extension
    base, ext = os.path.splitext(input_path)
    output_path = f"{base}_z{n_str}{ext}"

    # Write modified file
    with open(output_path, 'w') as file:
        file.writelines(lines)

    print(f"New file created: {output_path}")





