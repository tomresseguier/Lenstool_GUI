import re
import pandas as pd
import numpy as np
import math
import os


def find_ref(bestopt_path) :
    with open(bestopt_path, 'r') as text_file :
        bestopt_full_string = text_file.read()
    pattern = re.compile(r"reference\s+\d+\s+([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)")
    match = pattern.search(bestopt_full_string)
    if match:
        ra = float(match.group(1))
        dec = float(match.group(2))
        return ra, dec
    else:
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

def generate_latex_table(ref_coord, bestopt_dataframes) :
    ref_ra, ref_dec = ref_coord
    
    table_rows = []

    for key in bestopt_dataframes.keys() :
        df = bestopt_dataframes[key]
        table_row = key
        
        for _, row in df.iterrows() :
            if float(row['uncertainty2'])!=0.0 :
                raise ValueError("Sorry, uniform prior uncertainties not yet implemented (to come).")
        
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

def bestopt_latex_table(bestopt_path) :
    ref_coord, bestopt_dataframes = bestopt_pot_param_extractor(bestopt_path)
    table_str = generate_latex_table(ref_coord, bestopt_dataframes)
    bestopt_latex_path = os.path.join(os.path.dirname(bestopt_path), 'bestopt.latex')
    open(bestopt_latex_path, 'w').write(table_str)
    return table_str








def make_dict_from_bayes(df) :
    param_dict = {}
    for col in df.select_dtypes(include='number').columns :
        lower, median, upper = np.percentile(df[col], [16, 50, 84])
        param_dict[col] = [lower, median, upper]
    return param_dict


def param_latex_table(df) :
    param_dict = make_dict_from_bayes(df)
    
    ### Extract the potential names ###
    pot_list = []
    for col in param_dict.keys() :
        pot_name = col.split()[0]
        if pot_name.startswith('O') :
            pot_list.append(pot_name)
    pots = np.unique(pot_list)
    
    correspondances = { 'x_centre': 'x', \
                        'y_centre': 'y', \
                        'ellipticity': 'emass', \
                        'angle_pos': 'theta', \
                        'core_radius_kpc': 'rc', \
                        'cut_radius_kpc': 'rcut', \
                        'v_disp': 'sigma'}
        
    pot_dict = {}
    pot_dict['name'] = pots
    for name in correspondances.keys() :
        pot_dict[name] = ['...', '...', '...']
    #pot_df = np.array(pot_dict)
    
    for i, pot in enumerate(pots) :
        for col in param_dict.keys() :
            if col.split()[0]==pot :
                for name in correspondances.keys() :
                    if correspondances[name] in col :
                        pot_dict[name][i] = param_dict[col]
        
        
        
    
    
    
    return table_str
















