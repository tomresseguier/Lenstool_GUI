import numpy as np
from scipy.spatial.distance import pdist, squareform



def read_text_file(file_path):
    lines = []
    with open(file_path, 'r') as file:
        for line in file:
            cleaned_line = line.strip()
            if not cleaned_line.startswith("#"):
                lines.append(cleaned_line)
    return lines


def find_close_coord(cat, d):
    ids = cat['id']
    x = cat['x']
    y = cat['y']
    
    points = np.vstack((x, y)).T
    distance_matrix = squareform(pdist(points))
    n = len(points)
    
    parent = list(range(n))
    
    def find(u):
        if parent[u] != u:
            parent[u] = find(parent[u])
        return parent[u]
    
    def union(u, v):
        root_u = find(u)
        root_v = find(v)
        if root_u != root_v:
            parent[root_u] = root_v

    for i in range(n):
        for j in range(i + 1, n):
            if distance_matrix[i, j] < d:
                union(i, j)
    
    groups = {}
    for i in range(n):
        root = find(i)
        if root in groups:
            groups[root].append(ids[i])
        else:
            groups[root] = [ids[i]]
    
    result = list(groups.values())
    
    return result


def remove_png_margins(im) :
    x_shape = im.shape[0]
    y_shape = im.shape[1]
    
    threshold = 3000.
    
    i=0
    check = False
    while not check :
        i+=1
        check = not np.sum(im[0:i,:,:])>i*y_shape*4.-threshold
        if check :
            print(np.sum(im[0:i,:,:]), i*y_shape*4.)
    left_x_margin = i-1
    
    i=0
    check = False
    while not check :
        i+=1
        check = not np.sum(im[-i-1:-1,:,:])>i*y_shape*4.-threshold
    right_x_margin = -i-1
    
    i=0
    check = False
    while not check :
        i+=1
        check = not np.sum(im[:,0:i,:])>i*x_shape*4.-threshold
    bottom_y_margin = i-1
    
    i=0
    check = False
    while not check :
        i+=1
        check = not np.sum(im[:,-i-1:-1,:])>i*x_shape*4.-threshold
    top_y_margin = -i-1
    
    im_cropped = im[left_x_margin:right_x_margin, bottom_y_margin:top_y_margin, :]
    
    return im_cropped


def make_colnames_dict(catalog, use_default_names=True):
    """
    Extracts column names for positions and shape parameters from an Astropy table.
    Parameters:
    - catalog: Astropy Table
        Input catalog containing astronomical data.
    Returns:
    - column_names: list
        List of column names for positions and shape parameters present in the catalog.
    """
    
    to_test_names_dict = {}
    to_test_names_dict['ra'] = ['ra', 'ALPHA_J2000', 'X_WORLD']
    to_test_names_dict['dec'] = ['dec', 'DELTA_J2000', 'Y_WORLD']
    #to_test_names_dict['x'] = ['X_IMAGE', 'x']
    #to_test_names_dict['y'] = ['Y_IMAGE', 'y']
    to_test_names_dict['a'] = ['a', 'A_IMAGE']
    to_test_names_dict['b'] = ['b', 'B_IMAGE']
    to_test_names_dict['theta'] = ['angle', 'theta', 'THETA_IMAGE']
    
    names_list = list(to_test_names_dict.keys())
    #names_list = ['ra', 'dec', 'x', 'y', 'a', 'b']
    names_dict = {}
    names_dict_default = {}
    for name in names_list :
        names_dict[name] = []
        names_dict_default[name] = None
    for name in names_list :
        for to_test_name in to_test_names_dict[name] :
            cat_colnames_lower = [col.lower() for col in catalog.colnames]
            
            #if 'colnames' in dir(catalog) :
            #    cat_colnames_lower = [col.lower() for col in catalog.colnames]
            #else :
            #    cat_colnames_lower = [col.lower() for col in catalog.columns.names]
            
            if to_test_name.lower() in cat_colnames_lower :
                col_idx = np.where( np.array(cat_colnames_lower) == to_test_name.lower() )[0][0]
                names_dict[name].append(catalog.colnames[col_idx])
                names_dict_default[name] = catalog.colnames[col_idx]
    
    print('Columns found in catalog: \n' + str(names_dict))
    yesno = 'y' if use_default_names else input('Columns to be used: \n' + str(names_dict_default) + \
                                                '\nKeep these names? (if no, user prompted to select other columns) [y][n]')
    if yesno == 'n' :
        for name in names_list :
            if len(names_dict[name]) > 1 :
                selected_name = input("Several columns found for name " + name + ": " + str(names_dict[name]) + ". Which one should be kept (if unit, should be image pixels)?")
                names_dict[name] = selected_name
            else :
                names_dict[name] = names_dict[name][0]
    else :
        names_dict = names_dict_default
    return names_dict

