import re
import pandas as pd
from io import StringIO
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
from ..utils_astro.get_cosmology import get_cosmo
from ..utils_astro.utils_general import world_to_relative
cosmo = get_cosmo()



def extract_main_pot_names(param_file_path) :
    pot_names = []
    pattern = re.compile(r'^\s*potentiel\s+(\S+)', re.MULTILINE)
    with open(param_file_path, 'r') as file:
        content = file.read()
    pot_names = pattern.findall(content)
    return pot_names

###############################################################################
# Goes from bayes.dat file to dataFrame with all the parameters as columns and 
# their values in the different samples as lines
#
#                        Nsample  ln(Lhood)      ...    Pot0 sigma   Chi2
#                   0    1        -64393.427747  ...    195.359518   128880.763
#bayes.dat -->      1    1        -64400.548013  ...    195.365226   128895.003
#                   2    1        -64392.938091  ...    195.363315   128879.784
#                   3    1        -64392.622949  ...    195.359561   128879.153
###############################################################################
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

###############################################################################
#Goes from bayes dataFrame to dictionary with calculated statistics
#bayes_df --> param_dict
###############################################################################
def calculate_params_from_bayes(bayes_df) :
    param_dict = {}
    for col in bayes_df.select_dtypes(include='number').columns :
        lower, median, upper = np.percentile(bayes_df[col], [16, 50, 84])
        param_dict[col] = [lower, median, upper]
    return param_dict

###############################################################################
#Goes from bayes dataFrame to dictionary with calculated statistics sorted per 
#potential
#bayes_df --> pot_dict
###############################################################################
def make_table_dict_from_bayes(bayes_df) :
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
    
    pot_dict = {}
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
        pot_dict[pot] = pd.DataFrame(indiv_pot_dict)
        
    return pot_dict

###############################################################################
#Completes pot_dict with fixed values from param file
#param_file_path --> final_pot_dict
###############################################################################
def complete_pot_dict(param_file_path, convert_to_kpc=False, z=None) :
    model_dir = os.path.dirname(param_file_path)
    bayes_file_path = os.path.join(model_dir, 'bayes.dat')
    bayes_df = read_bayes_file(bayes_file_path, convert_to_kpc=convert_to_kpc, z=z)
    final_pot_dict = make_table_dict_from_bayes(bayes_df)
    param_file = pylenstool.lenstool_param_file(param_file_path)
    pot_names_alt = extract_main_pot_names(param_file_path)
    
    for i, pot in enumerate(final_pot_dict.keys()) :
        if pot!='Pot0' :
            for j, value in enumerate(final_pot_dict[pot]['value']) :
                if value=='...' :
                    param_line = param_file.get_parameter('potentiel ' + pot_names_alt[i], final_pot_dict[pot]['name'][j])
                    param_line_kpc = param_file.get_parameter('potentiel ' + pot_names_alt[i], final_pot_dict[pot]['name'][j] + '_kpc')
                    if param_line!=0 :
                        final_pot_dict[pot]['value'][j] = param_line[-1]
                    if param_line_kpc!=0 :
                        final_pot_dict[pot]['name'][j] = final_pot_dict[pot]['name'][j] + '_kpc'
                        final_pot_dict[pot]['value'][j] = param_line_kpc[-1]
    return final_pot_dict

###############################################################################
#Generates LaTeX string from dictionary with calculated statistics sorted per 
#potential
#final_pot_dict --> table_str
###############################################################################
def generate_latex_table(final_pot_dict, ref_coord=None) :
    if ref_coord is not None :
        ref_ra, ref_dec = ref_coord
    
    table_rows = []

    for key in final_pot_dict.keys() :
        df = final_pot_dict[key]
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

###############################################################################
#Saves LaTeX table as text file
#model_dir --> LaTeX table as text file
###############################################################################
def make_param_latex_table(param_file_path, convert_to_kpc=True, z=None) :
    model_dir = os.path.dirname(param_file_path)
    final_pot_dict = complete_pot_dict(param_file_path, convert_to_kpc=convert_to_kpc, z=z)
    table_str = generate_latex_table(final_pot_dict)
    param_latex_table = os.path.join(model_dir, 'best_params.latex')
    open(param_latex_table, 'w').write(table_str)
    return table_str

def find_ref(lenstool_file_path) :
    with open(lenstool_file_path, 'r') as text_file :
        bestopt_full_string = text_file.read()
    pattern = re.compile(r"reference\s+\d+\s+([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)")
    match = pattern.search(bestopt_full_string)
    if match :
        ra = float(match.group(1))
        dec = float(match.group(2))
        return ra, dec
    else :
        raise ValueError("Reference RA and DEC not found in the file content.")





###############################################################################
#Extracts the potfile and reference galaxy parameters from param file
#Extracts the optimized potfile and reference galaxy parameters from bayes.dat
#Joins both dictionaries
#param_file_path --> potfile_param_dict
###############################################################################
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
    
    final_pot_dict = complete_pot_dict(param_file_path)
    Pot0_df = final_pot_dict['Pot0']
    for line in Pot0_df.iloc :
        if line['value']!='...' :
            potfile_param_dict[line['name']] = float(line['value'])
    return potfile_param_dict

###############################################################################
#Reads potfile and puts galaxies in a catalog
#potfile_path --> potfile_cat
###############################################################################
def read_potfile(potfile_path) :
    with open(potfile_path, 'r') as file :
        lines = file.readlines()
    potfile_cat = Table(names=['id','ra','dec','a','b','theta','mag','lum'], dtype=['int', *['float',]*7])
    for line in lines :
        columns = line.split()
        if columns[0].isdigit() :
            potfile_cat.add_row( [col for col in columns] )
    return potfile_cat



def read_arclet(arclet_path) :
    with open(arclet_path, 'r') as file :
        lines = file.readlines()
    arclet_cat = Table(names=['id','ra','dec','a','b','theta','mag','lum'], dtype=['str', *['float',]*7])
    for line in lines :
        columns = line.split()
        arclet_cat.add_row( [col for col in columns] )
    return arclet_cat
###############################################################################
#Uses get_potfile_params and read_potfile to create a complete catalog of all 
#galaxies in the potfile with their optimized parameters according to the 
#scaling relation.
#param_file_path --> potfile_cat_opt
###############################################################################
def make_potfile_dict(param_file_path) :
    #param_file = pylenstool.lenstool_param_file(param_file_path)
    #ref_str = param_file.get_ref()
    #ref = (float(ref_str[0]), float(ref_str[1]))
    ref = find_ref(param_file_path)
    
    potfile_param_dict = get_potfile_params(param_file_path)
    profile_type = int(potfile_param_dict['type'])
    core_radius0 = potfile_param_dict['corekpc']
    cut_radius0 = potfile_param_dict['cut_radius']
    v_disp0 = potfile_param_dict['v_disp']
    mag0 = potfile_param_dict['mag0']
    z_lens = potfile_param_dict['zlens']
    
    slope = potfile_param_dict['slope']
    vdslope = potfile_param_dict['vdslope']
    
    
    potfile_cat = read_potfile(potfile_param_dict['potfile_path'])
        
    names = ['id','profile','x_centre','y_centre','ellipticity','angle_pos','core_radius','cut_radius','v_disp','mag','z_lens']
    potfile_cat_opt = Table(names=names, dtype=['int', *['float',]*(len(names)-1)])
    for i in range(len(potfile_cat)) :
        potfile_cat_opt.add_row([0] + [0.0] * (len(names)-1))
    
    potfile_cat_opt['id'] = potfile_cat['id']
    x, y = world_to_relative(potfile_cat['ra'], potfile_cat['dec'], ref)
    potfile_cat_opt['x_centre'] = x
    potfile_cat_opt['y_centre'] = y
    a, b = potfile_cat['a'], potfile_cat['b']
    potfile_cat_opt['ellipticity'] = (a**2-b**2)/(a**2+b**2)
    potfile_cat_opt['angle_pos'] = potfile_cat['theta']
    
    mag = potfile_cat['mag']
    potfile_cat_opt['mag'] = mag
    L_ratio = 10**((mag0-mag)/2.5)
    potfile_cat_opt['profile'] = profile_type
    potfile_cat_opt['core_radius'] = core_radius0 * L_ratio**0.5
    potfile_cat_opt['cut_radius'] = cut_radius0 * L_ratio**(2/slope)
    potfile_cat_opt['v_disp'] = v_disp0 * L_ratio**(1/vdslope)
    potfile_cat_opt['z_lens'] = z_lens
    return potfile_cat_opt

###############################################################################
#Formats potfile_cat_opt as a string according to the best file format.
#potfile_cat_opt --> opt_potfile_string
###############################################################################
def format_potfile_table(potfile_cat_opt) :
    output = []
    for row in potfile_cat_opt:        
        section = [
            f"potential       {row['id']}",
            f"\tprofile       {int(row['profile'])}",
            f"\tx_centre     {row['x_centre']:.6f}",
            f"\ty_centre     {row['y_centre']:.6f}",
            f"\tellipticity     {row['ellipticity']:.6f}",
            f"\tangle_pos       {row['angle_pos']:.6f}",
            f"\tcore_radius         {row['core_radius']:.6f}",
            f"\tcut_radius         {row['cut_radius']:.6f}",
            f"\tv_disp     {row['v_disp']:.6f}",
            f"\tmag\t  {row['mag']:.6f}",
            f"\tz_lens     {row['z_lens']:.4f}",
            "\tend"
        ]
        output.append("\n".join(section))
    
    opt_potfile_string = "\n".join(output)
    return opt_potfile_string

###############################################################################
#Extracts the main potentials' redshifts
#param_file_path --> main_pot_z_dict
###############################################################################
def get_main_pot_z(param_file_path) :
    param_file = pylenstool.lenstool_param_file(param_file_path)
    pot_names = extract_main_pot_names(param_file_path)
    main_pot_z_dict = {}
    for pot_name in pot_names :
        main_pot_z_dict[pot_name] = float(param_file.get_parameter('potentiel ' + pot_name, 'z_lens')[-1])
    return main_pot_z_dict


###############################################################################
#Creates the string containing the temporary best file
#param_file_path --> param_str
###############################################################################
def make_param_str(param_file_path) :
    model_dir = os.path.dirname(param_file_path)
    bayes_file_path = os.path.join(model_dir, 'bayes.dat')
    bayes_df = read_bayes_file(bayes_file_path, convert_to_kpc=False)
    final_pot_dict = complete_pot_dict(param_file_path)
    final_pot_dict.pop('Pot0')
    main_pot_z_dict = get_main_pot_z(param_file_path)
    main_pot_z_list = [main_pot_z_dict[key] for key in main_pot_z_dict.keys()]
    
    param_file = pylenstool.lenstool_param_file(param_file_path)
    ref_coord = param_file.get_ref()
    ref_coord_str = ref_coord[0] + ' ' + ref_coord[1]
    
    xmin = param_file.get_parameter('field', 'xmin')[1]
    xmax = param_file.get_parameter('field', 'xmax')[1]
    ymin = param_file.get_parameter('field', 'ymin')[1]
    ymax = param_file.get_parameter('field', 'ymax')[1]
    
    # Standard header
    param_str = """runmode\n\treference     3 """ + ref_coord_str + """\n\timage     1 arclets_to_predict.lenstool\n\tmass      4 1000 0.556000 mass.fits\n\tampli	   1 1000 3.000000 magnification.fits\n\tshear  1 1000 3.000000 shear.fits\n\tshearfield  1 25.000000 1.0 25\n\tend\ngrid\n\tnumber      30\n\tpolar     0\n\tnlens   156\n\tend\n"""
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
    for i, (pot_name, pot_df) in enumerate(final_pot_dict.items()) :
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
        param_str += f"\tz_lens     {main_pot_z_list[i]:.4f}\n"
        param_str += "\tend\n"
    
    
    potfile_cat_opt = make_potfile_dict(param_file_path)
    param_str += format_potfile_table(potfile_cat_opt)
    
    
    param_str += """cline\n\tnplan    1 2.000000\n\tdmax     0.000000\n\talgorithm   MARCHINGSQUARES\n\tlimitHigh   0.5\n\tlimitLow    0.100000\n\tend\ngrande\n\tiso         0 0 1.000000 0.000000 0.000000\n\tname        best\n\tprofile      0 0\n\tcontour     0 0\n\tlarge_dist  2.000000\n\tend\ncosmology\n\tmodel       1\n\tH0        70.000000\n\tomegaM    0.300000\n\tomegaX    0.700000\n\tomegaK    0.\n\twX        -1.000000\n\twa        0.000000\n\tend\nfield\n\txmin     """ + xmin + """\n\txmax     """ + xmax + """\n\tymin     """ + ymin + """\n\tymax     """ + ymax + """\n\tend\nfinish"""
    return param_str


###############################################################################
#Saves param_str as best.par file
#param_file_path --> best.par file
###############################################################################
def make_best_file_from_bayes(param_file_path) :
    model_dir = os.path.dirname(param_file_path)
    param_str = make_param_str(param_file_path)
    temp_best_file_path = os.path.join(model_dir, 'best_TEMP.par')
    open(temp_best_file_path, 'w').write(param_str)



def extract_magnification_from_dist(dist_path) :
    dist_text = open(dist_path, 'r').read()
    data_lines = [line for line in dist_text.splitlines() if line.strip() and not line.strip().startswith('#')]
    table_text = "\n".join(data_lines)
    columns = ['ID', 'X', 'Y', 'R', 'EPS', 'TAU', 'AMP', 'E_AMP', 'DMAG', 'TIME[days]', 'DTIME', 'PARITY']
    table_df = pd.read_csv(StringIO(table_text), delim_whitespace=True, names=columns)
    return table_df



































###############################################################################
# UNUSED. Not ideal as relies on running Lenstool
###############################################################################
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







"""
###############################################################################
# OBSOLETE
###############################################################################
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

###############################################################################
# OBSOLETE
###############################################################################
def bestopt_latex_table(model_dir) :
    bestopt_path = os.path.join(model_dir, 'bestopt.par')
    ref_coord, bestopt_dataframes = bestopt_pot_param_extractor(bestopt_path)
    table_str = generate_latex_table(bestopt_dataframes, ref_coord=ref_coord)
    bestopt_latex_path = os.path.join(os.path.dirname(bestopt_path), 'bestopt.latex')
    open(bestopt_latex_path, 'w').write(table_str)
    return table_str
"""

