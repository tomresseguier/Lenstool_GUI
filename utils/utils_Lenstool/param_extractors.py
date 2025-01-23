import re
import pandas as pd
import numpy as np
import math
import os
import sys
import shutil
import astropy.units as u
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
import subprocess
import pylenstool
module_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(module_dir))
from utils_astro.set_cosmology import set_cosmo
from utils_astro.utils_general import world_to_relative
cosmo = set_cosmo()



def read_bayes_file(file_path, convert_to_kpc=True, z=None):
    columns = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("#"):
                col_name = line[1:].strip()
                columns.append(col_name)
            else:
                break
    
    df = pd.read_csv(file_path, delim_whitespace=True, comment='#', header=None)
    df.columns = columns
    
    if convert_to_kpc :
        if z is not None :
            ang_diam_dist_value = cosmo.angular_diameter_distance(z).to('kpc').value
            kpc_per_arcsec = ang_diam_dist_value * 1.*u.arcsec.to(u.rad)
            for col in df.columns :
                if 'rc (arcsec)' in col or 'rcut (arcsec)' in col :
                    df[col] = df[col]*kpc_per_arcsec
                    df = df.rename(columns={col: col[:-len('(arcsec)')] + '(kpc)'})
        else :
            raise ValueError("A redshift must be specified.")
    return df

def calculate_params_from_bayes(bayes_df) :
    param_dict = {}
    for col in bayes_df.select_dtypes(include='number').columns :
        lower, median, upper = np.percentile(bayes_df[col], [16, 50, 84])
        param_dict[col] = [lower, median, upper]
    return param_dict

def make_table_dict_from_bayes(bayes_df):
    param_dict = calculate_params_from_bayes(bayes_df)
    
    ### Extract the potential names ###
    pot_list = []
    for col in param_dict.keys() :
        pot_name = col.split()[0]
        if pot_name.startswith('O') or pot_name.startswith('Pot') :
            pot_list.append(pot_name)
    pots = np.unique(pot_list)
    
    correspondances = { 'x_centre': 'x',
                        'y_centre': 'y',
                        'ellipticity': 'emass',
                        'angle_pos': 'theta'}
    
    check_kpc = False
    for col in bayes_df.columns :
        if 'kpc' in col :
            check_kpc = True
    if check_kpc :
        correspondances['core_radius_kpc'] = 'rc'
        correspondances['cut_radius_kpc'] = 'rcut'
    else :
        correspondances['core_radius'] = 'rc'
        correspondances['cut_radius'] = 'rcut'
    correspondances['v_disp'] = 'sigma'
    
    final_pot_dict = {}
    for pot in pots :
        indiv_pot_dict = {}
        indiv_pot_dict['name'] = ['...' for i in range(len(correspondances.keys()))]
        indiv_pot_dict['value'] = ['...' for i in range(len(correspondances.keys()))]
        indiv_pot_dict['uncertainty1'] = ['...' for i in range(len(correspondances.keys()))]
        indiv_pot_dict['uncertainty2'] = ['...' for i in range(len(correspondances.keys()))]
        for col in param_dict.keys() :
            if col.split()[0]==pot :
                for i, name in enumerate(correspondances.keys()) :
                    indiv_pot_dict['name'][i] = name
                    if correspondances[name] in col.split() :
                        indiv_pot_dict['value'][i] = param_dict[col][1]
                        indiv_pot_dict['uncertainty1'][i] = param_dict[col][1] - param_dict[col][0]
                        indiv_pot_dict['uncertainty2'][i] = param_dict[col][2] - param_dict[col][1]
        final_pot_dict[pot] = pd.DataFrame(indiv_pot_dict)
    
    return final_pot_dict

def generate_latex_table(pot_dict, ref_coord=None) :
    if ref_coord is not None :
        ref_ra, ref_dec = ref_coord
    
    table_rows = []

    for key in pot_dict.keys() :
        df = pot_dict[key]
        table_row = key
        
        uncertainty_from_bayes = False
        for _, row in df.iterrows() :
            if type(row['uncertainty2']) is not str :
                if abs(row['uncertainty2'])>0. :
                    uncertainty_from_bayes = True
        
        
        for name in ['x_centre', 'y_centre', 'ellipticity', 'angle_pos', 'core_radius_kpc', 'cut_radius_kpc', 'v_disp'] :
            i_array = np.where(df['name']==name)[0]
            if len(i_array)!=0 :
                i = i_array[0]
                
                #if name=='x_centre' :
                #    value_ra = ref_ra - float(df['value'][i]) / 3600.0 / np.cos(np.deg2rad(ref_dec))
                #elif name=='y_centre' :
                #    value_dec = ref_dec + float(df['value'][i]) / 3600.0
                #else :
                #    value = float(df['value'][i])
                
                if uncertainty_from_bayes :
                    if type(df['value'][i]) is str :
                        table_row = table_row + " & ..."
                    else :
                        value = df['value'][i]
                        uncertainty1 = df['uncertainty1'][i]
                        uncertainty2 = df['uncertainty2'][i]
                        
                        
                        sig_fig1 = -math.floor(math.log10(abs(uncertainty1)))
                        if sig_fig1>=0 :
                            sig_fig1 = -math.floor(math.log10(abs( float(f"{uncertainty1:.{sig_fig1}f}") )))
                        if sig_fig1<0 :
                            sig_fig1 = 0
                        
                        sig_fig2 = -math.floor(math.log10(abs(uncertainty2)))
                        if sig_fig2>=0 :
                            sig_fig2 = -math.floor(math.log10(abs( float(f"{uncertainty2:.{sig_fig2}f}") )))
                        if sig_fig2<0 :
                            sig_fig2 = 0
                        
                        sig_fig = max(sig_fig1, sig_fig2)
                        
                        table_row = table_row + f" & ${value:.{sig_fig}f}" + "^{+" + f"{uncertainty2:.{sig_fig}f}" + "}_{-" + f"{uncertainty1:.{sig_fig}f}" + "}$"
                    
                else :
                    value = float(df['value'][i])
                    uncertainty = float(df['uncertainty1'][i])
                    
                    sig_fig = -math.floor(math.log10(abs(uncertainty)))
                    if sig_fig>=0 :
                        sig_fig = -math.floor(math.log10(abs( float(f"{uncertainty:.{sig_fig}f}") )))
                    if sig_fig<0 :
                        sig_fig = 0
                    
                    table_row = table_row + f" & {value:.{sig_fig}f}" + f" $\pm$ {uncertainty:.{sig_fig}f}"
                    
                    
            else :
                table_row = table_row + " & ..."
                
        table_row = table_row + " \\\\"
        table_rows.append(table_row)
        
    table_str = ( "\\begin{table*}[htb!]\n"
                  "\\begin{center}\n"
                  "    \\begin{tabular}{c c c c c c c c}\n"
                  "        \\hline\\hline\n"
                  "        Component & $\\Delta$ R.A. (\") & $\\Delta$ Dec (\") & ellipticity & $\\theta$ (deg) & $r_{\\rm core}$ (kpc) & $r_{\\rm cut}$ (kpc) & $\\sigma_0$ (km s$^{-1}$) \\\\\n"
                  "        \\hline\n" ) + \
                "        " + "\n        ".join(table_rows) + "\n" + \
                ( "        \\hline\n"
                  "    \\end{tabular}\n"
                  "    \\caption{Parameters of the lens model.}\n"
                  "    \\label{tab:SL_params}\n"
                  "\\end{center}\n"
                  "\\end{table*}" )
    
    return table_str

def make_param_latex_table(model_dir, convert_to_kpc=True, z=None) :
    bayes_file_path = os.path.join(model_dir, 'bayes.dat')
    bayes_df = read_bayes_file(bayes_file_path, convert_to_kpc=convert_to_kpc, z=z)
    final_pot_dict = make_table_dict_from_bayes(bayes_df)
    table_str = generate_latex_table(final_pot_dict)
    param_latex_table = os.path.join(model_dir, 'best_params.latex')
    open(param_latex_table, 'w').write(table_str)
    return table_str
    

    





def find_ref(bestopt_path) :
    with open(bestopt_path, 'r') as text_file :
        bestopt_full_string = text_file.read()
    pattern = re.compile(r"reference\s+\d+\s+([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)")
    match = pattern.search(bestopt_full_string)
    if match :
        ra = float(match.group(1))
        dec = float(match.group(2))
        return ra, dec
    else :
        raise ValueError("Reference RA and DEC not found in the file content.")

def bestopt_pot_param_extractor(bestopt_path) :
    with open(bestopt_path, 'r') as text_file :
        bestopt_full_string = text_file.read()
        
    pattern = re.compile( r"limit\s+(\w+)\s*"
                          r"((?:\s*\w+\s+\d+\s+[+-]?\d*\.?\d+\s+[+-]?\d*\.?\d+\s+[+-]?\d*\.?\d+\s*)+)"
                          r"end", re.DOTALL)
    
    pattern = re.compile(
    r"limit\s+(\w+)\s*"
    r"((?:\s*\w+\s+\d+\s+[+-]?\d*\.?\d+\s+[+-]?\d*\.?\d+\s+[+-]?\d*\.?\d+\s*)+)"
    r"end", re.DOTALL)
    
    sections = pattern.findall(bestopt_full_string)
    dataframes = {}
    for name, section_text in sections:
        param_pattern = re.compile(
            r"(\w+)\s+\d+\s+([+-]?\d*\.?\d+)\s+([+-]?\d*\.?\d+)\s+([+-]?\d*\.?\d+)"
        )
        matches = param_pattern.findall(section_text)
        
        dictionaries = []
        for match in matches:
            param_name = match[0].strip()
            value = match[1]
            uncertainty1 = match[2]
            uncertainty2 = match[3]
            dictionaries.append({
                'name': param_name,
                'value': value,
                'uncertainty1': uncertainty1,
                'uncertainty2': uncertainty2
            })
            
        df = pd.DataFrame(dictionaries)
        dataframes[name] = df
    
    for key in dataframes.keys() :
        print(f"Extracted values for limit {key}:")
        print(dataframes[key], "\n")
    
    ref_coord = find_ref(bestopt_path)
    
    return ref_coord, dataframes

def bestopt_latex_table(model_dir) :
    bestopt_path = os.path.join(model_dir, 'bestopt.par')
    ref_coord, bestopt_dataframes = bestopt_pot_param_extractor(bestopt_path)
    table_str = generate_latex_table(bestopt_dataframes, ref_coord=ref_coord)
    bestopt_latex_path = os.path.join(os.path.dirname(bestopt_path), 'bestopt.latex')
    open(bestopt_latex_path, 'w').write(table_str)
    return table_str



def make_param_str_from_dict(pot_dict, z_lens) :
    # Standard header
    param_str = """runmode\n\treference     3 179.488897 -10.766923\n\timage     1 PBC_arclets_to_predict.lenstool\n\tmass      4 1000 0.556000 mass.fits\n\tampli	   1 1000 3.000000 magnification.fits\n\tshear  1 1000 3.000000 shear.fits\n\tshearfield  1 25.000000 1.0 25\n\tend\ngrid\n\tnumber      30\n\tpolar     0\n\tnlens   156\n\tend\n"""
    # Map from our parameter names to Lenstool names
    param_map = {
        'x_centre': 'x_centre',
        'y_centre': 'y_centre', 
        'ellipticity': 'ellipticity',
        'angle_pos': 'angle_pos',
        'core_radius': 'core_radius',
        'core_radius_kpc': 'core_radius_kpc',
        'cut_radius': 'cut_radius',
        'cut_radius_kpc': 'cut_radius_kpc',
        'v_disp': 'v_disp'
    }
    
    # Add each potential
    for pot_name, pot_df in pot_dict.items():
        param_str += f"potential {pot_name}\n"
        param_str += "\tprofile       81\n"  # Assuming all are profile 81
        
        # Add each parameter
        for _, row in pot_df.iterrows():
            param_name = row['name']
            if param_name in param_map and row['value'] != '...':
                lenstool_name = param_map[param_name]
                value = float(row['value'])
                param_str += f"\t{lenstool_name}     {value:.6f}\n"
        
        # Add z_lens at the end of each potential
        param_str += f"\tz_lens     {z_lens:.4f}\n"
        param_str += "\tend\n"
    
    
    
    param_str += """cline\n\tnplan    1 2.000000\n\tdmax     0.000000\n\talgorithm   MARCHINGSQUARES\n\tlimitHigh   0.5\n\tlimitLow    0.100000\n\tend\ngrande\n\tiso         0 0 1.000000 0.000000 0.000000\n\tname        best\n\tprofile      0 0\n\tcontour     0 0\n\tlarge_dist  2.000000\n\tend\ncosmology\n\tmodel       1\n\tH0        70.000000\n\tomegaM    0.300000\n\tomegaX    0.700000\n\tomegaK    0.\n\twX        -1.000000\n\twa        0.000000\n\tend\nfield\n\txmin     -120.000000\n\txmax     80.000000\n\tymin     -100.000000\n\tymax     100.000000\n\tend\nfinish"""
    return param_str



def get_potfile_params(param_file_path) :
    param_file = pylenstool.lenstool_param_file(param_file_path)
    potfile_param_dict = {}
    potfile_name = param_file.get_parameter('potfile', 'filein')[-1]
    potfile_path = os.path.join(os.path.dirname(param_file_path), potfile_name)
    potfile_param_dict['potfile_path'] = potfile_path
    for name in ['type', 'zlens', 'mag0', 'corekpc'] :
        potfile_param_dict[name] = float(param_file.get_parameter('potfile', name)[1])
    for name in ['sigma', 'cutkpc', 'slope', 'vdslope'] :
        if int(param_file.get_parameter('potfile', name)[1])==0 :
            potfile_param_dict[name] = float(param_file.get_parameter('potfile', name)[2])
    return potfile_param_dict







def make_potfile_param_dict(param_file_path) :
    model_dir = os.path.dirname(param_file_path)
    bayes_df = read_bayes_file(os.path.join(model_dir, 'bayes.dat'), convert_to_kpc=False)
    final_pot_dict = make_table_dict_from_bayes(bayes_df)
    potfile_param_dict = get_potfile_params(param_file_path)
    Pot0_df = final_pot_dict['Pot0']
    for line in Pot0_df.iloc :
        if line['value']!='...' :
            potfile_param_dict[line['name']] = float(line['value'])
    return potfile_param_dict


def read_potfile(potfile_path) :
    with open(potfile_path, 'r') as file :
        lines = file.readlines()
    potfile_cat = Table(names=['id','ra','dec','a','b','theta','mag','lum'], dtype=['int', *['float',]*7])
    for line in lines :
        columns = line.split()
        if columns[0].isdigit() :
            potfile_cat.add_row( [col for col in columns] )
    return potfile_cat



def make_potfile_dict(param_file_path) :
    param_file = pylenstool.lenstool_param_file(param_file_path)
    ref_str = param_file.get_ref()
    ref = (float(ref_str[0]), float(ref_str[1]))
    
    potfile_param_dict = make_potfile_param_dict(param_file_path)
    potfile_cat = read_potfile(potfile_param_dict['potfile_path'])
        
    names = ['id','profile','x_centre','y_centre','ellipticity','angle_pos','core_radius','cut_radius','v_disp','mag','z_lens']
    potfile_cat_alt = Table(names=names, dtype=['int', *['float',]*(len(names)-1)])
    for i in range(len(potfile_cat)) :
        potfile_cat_alt.add_row([0] + [0.0] * (len(names)-1))
    
    potfile_cat_alt['id'] = potfile_cat['id']
    x, y = world_to_relative(potfile_cat['ra'], potfile_cat['dec'], reference=ref)
    potfile_cat_alt['x'] = x
    potfile_cat_alt['y'] = y
    a, b = potfile_cat['a'], potfile_cat['b']
    potfile_cat_alt['ellipticity'] = (a**2-b**2)/(a**2+b**2)
    potfile_cat_alt['angle_pos'] = potfile_cat['theta']
    
    
    
    return potfile_cat_alt



def get_lenstool_WCS(param_file_path) :
    param_file = pylenstool.lenstool_param_file(param_file_path)
    ref = param_file.get_ref()
    
    empty_param_str = "runmode\n\treference 3 " + ref[0] + " " + ref[1] + "\n\tsource 0\n\tmass 4 1000 0.556000 mass.fits\n\tend\npotentiel P01\n\tprofile 81\n\tx_centre 0.0\n\ty_centre 0.0\n\tz_lens 1.\n\tellipticite 0.0\n\tangle_pos 0.0\n\tcore_radius_kpc 10\n\tcut_radius_kpc 1500\n\tv_disp 600\n\tidentity P01\n\tend\nfini"
    
    temp_dir = os.path.join( os.path.dirname(param_file_path), 'TEMP/' )
    os.mkdir(temp_dir)
    empty_param_path = os.path.join(temp_dir, 'temp.par')
    with open(empty_param_path, 'w') as file :
        file.write(empty_param_str)
    
    cwd = os.getcwd()
    os.chdir(temp_dir)
    subprocess.run('conda activate Planck_bullet_cluster', shell=True)
    subprocess.run('lenstool temp.par -n', shell=True)
    with fits.open('./' + 'mass.fits') as hdulist :
        lenstool_WCS = WCS(hdulist[0].header)
    os.chdir(cwd)
    shutil.rmtree(temp_dir)
    
    return lenstool_WCS



param_file_path = '../../Lenstool_runs/PBC/RUN_059/PBC_strong.par'
make_potfile_dict(param_file_path)





