from astropy.io import fits
from astropy.table import Table, vstack, hstack
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from tqdm import tqdm
from scipy import spatial

#from astropy.visualization import astropy_mpl_style
#plt.style.use(astropy_mpl_style)


def plot_galaxies(catalogs_paths) :
    fig, ax = plt.subplots()
    colors = ["blue", "teal", "turquoise", "cyan"]
    for i, catalog_path in enumerate(catalogs_paths) :
        catalog = Table.read(catalog_path, format='fits')
        ax.plot(catalog['ra'], catalog['dec'], 'x', color=colors[i])

def combine_catalogs_dumb(catalogs_paths) :
    catalogs = []
    for i, catalog_path in enumerate(catalogs_paths) :
        catalogs.append(Table.read(catalog_path, format='fits'))
        print("length of catalog ", i, ": ", catalogs[i]['ra'].shape)
    combined_catalog = vstack(catalogs)
    return combined_catalog

def match_cat2(cats, match_radius=0.5, fill_in_value = np.nan, return_match_mask=False, keep_all_col=False, column_to_transfer=None) :
    print('YIIHAA')
    catalogs = []
    for cat in cats :
        if isinstance(cat, str) :
            catalogs.append( Table.read(cat, format='fits') )
        else :
            catalogs.append(cat.copy())
    cat_receiver = catalogs[0]
    cat_giver = catalogs[1]
    
    ra_dec_names = [('RA', 'DEC'), ('ra', 'dec')]
    for ra_name, dec_name in ra_dec_names :
        if ra_name in cat_receiver.colnames :
            ra_dec_names_receiver = (ra_name, dec_name)
        if ra_name in cat_receiver.colnames :
            ra_dec_names_giver = (ra_name, dec_name)    
    
    coords_cat_receiver = SkyCoord(cat_receiver[ra_dec_names_receiver[0]], cat_receiver[ra_dec_names_receiver[1]], unit=u.deg)
    coords_cat_giver = SkyCoord(cat_giver[ra_dec_names_giver[0]], cat_giver[ra_dec_names_giver[1]], unit=u.deg)
    idx, d2d, _ = coords_cat_receiver.match_to_catalog_sky(coords_cat_giver, nthneighbor=1)
    match_mask = d2d < match_radius*u.arcsec
    matched_cat = cat_receiver.copy()
    print(( cat_giver.colnames if column_to_transfer is None else [column_to_transfer] ))
    for colname in ( cat_giver.colnames if column_to_transfer is None else [column_to_transfer] ) :
        if colname not in cat_receiver.colnames :
            print('Adding column: ' + colname)
            matched_cat.add_column( np.full(len(cat_receiver), fill_in_value), name=colname )
            for index in np.where(match_mask)[0] :
                matched_cat[index][colname] = cat_giver[idx[index]][colname]
        elif keep_all_col :
            extra_colname = colname + '_CAT2'
            print('Adding column: ' + extra_colname)
            if extra_colname not in matched_cat.colnames :
                matched_cat.add_column( np.full(len(cat_receiver), fill_in_value), name=extra_colname )
            else : matched_cat.replace_column( extra_colname, np.full(len(cat_receiver), fill_in_value) )
            for index in np.where(match_mask)[0] :
                matched_cat[index][extra_colname] = cat_giver[idx[index]][colname]
            
    if return_match_mask :
        return matched_cat, match_mask
    else :
        return matched_cat


def combine_catalogs_remove_double_detections(cats, match_radius=0.5, preference_criterion='area') :
    catalogs = []
    for cat in cats :
        if isinstance(cat, str) :
            catalogs.append( Table.read(cat, format='fits') )
        else :
            catalogs.append(cat.copy())
    
    def len_cat_list(cat_list) :
        count = 0
        for cat in cat_list :
            count += len(cat)
        return count
    initial_len = len_cat_list(catalogs)
    
    for i, catalog_new in enumerate(catalogs) :
        print('##############\n' + "step ", i+1)
        #print("length of new catalog at step ", i, " before removal: ", len(catalog_new))
        
        if i==0 :
            catalog_previous = catalog_new
        else :
            #print("length of previous catalog at step ", i, " before removal: ", len(catalog_previous))
            print("length of catalog", i, "before removal: ", len(catalog_previous))
            print("length of catalog", i+1, "before removal: ", len(catalog_new))
            coords_previous = SkyCoord(catalog_previous['RA'], catalog_previous['DEC'], unit=u.deg)
            coords_new = SkyCoord(catalog_new['RA'], catalog_new['DEC'], unit=u.deg)
            idx, d2d, _ = coords_previous.match_to_catalog_sky(coords_new, nthneighbor=1)
            separation_mask = d2d < match_radius*u.arcsec
            
            filled_columns_new = []
            empty_columns_new = []
            if catalog_new.mask is not None :
                for col_name in catalog_new.colnames :
                    if len(np.unique( catalog_new.mask[col_name] ))==1 and True in np.unique( catalog_new.mask[col_name] ) :
                        empty_columns_new.append(col_name)
                    else :
                        filled_columns_new.append(col_name)
            else :
                filled_columns_new = catalog_new.colnames
            
            filled_columns_previous = []
            empty_columns_previous = []
            if catalog_previous.mask is not None :
                for col_name in catalog_previous.colnames :
                    if len(np.unique( catalog_previous.mask[col_name] ))==1 and True in np.unique( catalog_previous.mask[col_name] ) :
                        empty_columns_previous.append(col_name)
                    else :
                        filled_columns_previous.append(col_name)
            else :
                filled_columns_previous = catalog_previous.colnames
            
            from_previous_to_new = []
            for col in filled_columns_previous :
                if col not in filled_columns_new :
                    from_previous_to_new.append(col)
            
            from_new_to_previous = []
            for col in filled_columns_new :
                if col not in filled_columns_previous :
                    from_new_to_previous.append(col)
                   
            #print('from_previous_to_new')
            #print(from_previous_to_new)
            #print('##########')
            
            #print('from_new_to_previous')
            #print(from_new_to_previous)
            #print('##########')
                    
            for col in from_new_to_previous :
                if col not in catalog_previous.colnames :
                    catalog_previous.add_column( np.full(len(catalog_previous), np.nan), name=col )
                    #print('adding column ' + col)
                for index in np.where(separation_mask)[0] :
                    catalog_previous[index][col] = catalog_new[idx[index]][col]
                #print('adding values from new to previous for column ' + col)
            
            for col in from_previous_to_new :
                if col not in catalog_new.colnames :
                    catalog_new.add_column( np.full(len(catalog_new), np.nan), name=col )
                    #print('adding column ' + col)
                for k, index in enumerate(idx[separation_mask]) :
                    catalog_new[index][col] = catalog_previous[separation_mask][k][col]
                #print('adding values from previous to new for column ' + col)
            
            
            #print('####### before vstack #######')
            #print(np.where(np.logical_not(np.isnan( catalog_new['magAB_F814W'])) & np.logical_not(np.isnan( catalog_new['magAB_F435W'])))[0].shape)
            #print('##############')
            
            ####################### Preference criterion ######################
            if preference_criterion=='area' :                
                areas_previous = catalog_previous["a"][separation_mask].data * catalog_previous["b"][separation_mask].data
                areas_new = catalog_new["a"][idx][separation_mask].data * catalog_new["b"][idx][separation_mask].data
                #to_keep = np.argmax( [catalog_previous["MAG_AUTO"][separation_mask].data, catalog_new["MAG_AUTO"][idx][separation_mask].data], axis=0 )
                to_keep = np.argmax( [areas_previous, areas_new], axis=0 )
            ###################################################################
            
            #print("len(np.where(d2d < 1.*u.arcsec)[0]) = ", len(np.where(d2d < 1.*u.arcsec)[0]))
            to_keep_previous = np.where(separation_mask)[0][ to_keep==0 ]
            to_remove_previous = np.where(separation_mask)[0][ to_keep==1 ]
            to_keep_new = idx[separation_mask][ to_keep==1 ]
            to_remove_new = idx[separation_mask][ to_keep==0 ]
            
            #print("to_remove_new:", to_remove_new)
            #print("to_keep_new:", to_keep_new)
            #print("to_remove_previous:", to_remove_previous)
            #print("to_keep_previous:", to_keep_previous)
            
            
            catalog_previous.remove_rows(to_remove_previous)
            catalog_new.remove_rows(to_remove_new)
            
            print("length of catalog", i, "after removal: ", len(catalog_previous))
            print("length of catalog", i+1, "after removal: ", len(catalog_new))
            
            catalog_previous = vstack( [catalog_previous, catalog_new], join_type='outer' )
            #print('######## after vstack ######')
            #print(np.where(np.logical_not(np.isnan( catalog_new['magAB_F814W'])) & np.logical_not(np.isnan( catalog_new['magAB_F435W'])))[0].shape)
            #print('##############')
            
        
        
    print('initial length: ' + str(initial_len))
    print('final length: ' + str(len(catalog_previous)))
    print( str(initial_len - len(catalog_previous)) + ' matched sources removed' )
    return catalog_previous


def combine_DIM_catalogs(cats) :
    catalogs = []
    for cat in cats :
        if isinstance(cat, str) :
            catalogs.append( Table.read(cat, format='fits') )
        else :
            catalogs.append(cat.copy())
    
    incomplete_col_dict = {}
    complete_col_dict = {}
    for i, cat in enumerate(catalogs) :
        for colname in cat.columns :
            if len( np.where(np.isnan(cat[colname]))[0] )==0 :
                if colname in complete_col_dict.keys() :
                    complete_col_dict[colname].append(i)
                else :
                    complete_col_dict[colname] = [i]
            else :
                if colname in incomplete_col_dict.keys() :
                    incomplete_col_dict[colname].append(i)
                else :
                    incomplete_col_dict[colname] = [i]
    
    combined_cat = Table()
    for colname in complete_col_dict.keys() :
        combined_cat[colname] = catalogs[complete_col_dict[colname][0]][colname]
    for colname in incomplete_col_dict.keys() :
        incomplete_cols = [ catalogs[i][colname] for i in incomplete_col_dict[colname] ]
        combined_cat[colname] = np.nanmean(incomplete_cols, axis=0)
    
    return combined_cat



def remove_object(cat, FWHM_to_radius=1) :
    '''
    From pyRRG
    '''
    if isinstance(cat, str) :
        cat = Table.read(cat, format='fits')
    
    print("num of objects in the catalogue:", len(cat))

    x = np.array(cat['X_IMAGE'])
    y = np.array(cat['Y_IMAGE'])
    FWHM = np.array(cat['FWHM_IMAGE'])*FWHM_to_radius
    mask = np.ones(len(x),dtype=bool)
    ori_order = np.arange(len(x))                     ##original index in rrg_catalogue

    temp=np.vstack((x,y,FWHM,ori_order))            ##temp: x, y, FWHM, ori_order
    temp=temp.T
    temp_sort=temp[temp[:,2].argsort()[::-1]]       ##sorting by FWHM, from big objects to small one
    pos=list(zip(temp_sort[:,0],temp_sort[:,1]))          ##x, y
    tree=spatial.KDTree(pos)

    for i in range(len(x)):
        if mask[i]==False:
            continue
        rm_obj=tree.query_ball_point(pos[i],2.0*temp_sort[i][2]) #find neighbour objects
        for index in rm_obj:
            if index==i:
                continue
            delta_x=temp_sort[index][0]-temp_sort[i][0]
            delta_y=temp_sort[index][1]-temp_sort[i][1]
            distance=np.sqrt((delta_x)**2+(delta_y)**2)
            if distance<(temp_sort[i][2]+temp_sort[index][2]): #remove overlapped FWHM ojects
                mask[index]=False

    temp_sort=temp_sort[mask]
    sort_data=temp_sort[temp_sort[:, 3].argsort()]  ##sorting by the original index
    cat=cat[sort_data[:, 3].astype(int)]

    print(("Num of objects after removing double-detection: %i" % len(cat)))
    return cat




"""
def combine_catalogs_remove_double_detections(catalogs_paths=None, catalog_list=None, match_radius=0.5) :
    catalogs = []
    if catalog_list is None :
        for i, catalog_path in enumerate(catalogs_paths) :
            catalogs.append( Table.read(catalog_path, format='fits') )
    else :
        for cat in catalog_list :
            catalogs.append(cat.copy())
    
    def len_cat_list(cat_list) :
        count = 0
        for cat in cat_list :
            count += len(cat)
        return count
    initial_len = len_cat_list(catalogs)
    
    catalogs_match = []
    for i, catalog in enumerate(catalogs) :
        print("length of catalogue at step ", i, " before removal: ", len(catalog))
        if i != 0 :
            for j in range(i) :
                coords_previous = SkyCoord(catalogs_match[i-1-j]['RA'], catalogs_match[i-1-j]['DEC'], unit=u.deg)
                coords_new = SkyCoord(catalog['RA'], catalog['DEC'], unit=u.deg)
                idx, d2d, _ = coords_new.match_to_catalog_sky(coords_previous, nthneighbor=1)
                separation_mask = d2d < match_radius*u.arcsec
                
                filled_columns_new = []
                empty_columns_new = []
                if catalog.mask is not None :
                    for col_name in catalog.colnames :
                        if catalog.mask[col_name][0] :
                            empty_columns_new.append(col_name)
                        else :
                            filled_columns_new.append(col_name)
                else :
                    filled_columns_new = catalog.colnames
                
                
                filled_columns_previous = []
                empty_columns_previous = []
                if catalogs_match[i-1-j].mask is not None :
                    for col_name in catalogs_match[i-1-j].colnames :
                        if catalogs_match[i-1-j].mask[col_name][0] :
                            empty_columns_previous.append(col_name)
                        else :
                            filled_columns_previous.append(col_name)
                else :
                    filled_columns_previous = catalogs_match[i-1-j].colnames
                
                
                from_previous_to_new = []
                for col in filled_columns_previous :
                    if col not in filled_columns_new :
                        from_previous_to_new.append(col)
                
                from_new_to_previous = []
                for col in filled_columns_new :
                    if col not in filled_columns_previous :
                        from_new_to_previous.append(col)
                       
                        
                print('from_previous_to_new')
                print(from_previous_to_new)
                print('##########')
                
                print('from_new_to_previous')
                print(from_new_to_previous)
                print('##########')
                        
                       
                for col in from_new_to_previous :
                    if col not in catalogs_match[i-1-j].colnames :
                        catalogs_match[i-1-j].add_column( np.full(len(catalogs_match[i-1-j]), np.nan), name=col )
                        print('adding column ' + col)
                    for k, index in enumerate(idx[separation_mask]) :
                        catalogs_match[i-1-j][index][col] = catalog[separation_mask][k][col]
                    print('adding values from new to previous for column ' + col)
                
                for col in from_previous_to_new :
                    if col not in catalog.colnames :
                        catalog.add_column( np.full(len(catalog), np.nan), name=col )
                        print('adding column ' + col)
                    for index in np.where(separation_mask)[0] :
                        catalog[index][col] = catalogs_match[i-1-j][idx[index]][col]
                    print('adding values from previous to new for column ' + col)
                
                
                #print('####### before vstack #######')
                #print(np.where(np.logical_not(np.isnan( catalog['magAB_F814W'])) & np.logical_not(np.isnan( catalog['magAB_F435W'])))[0].shape)
                #print('##############')
                
                areas_previous = catalogs_match[i-1-j]["a"][idx][separation_mask].data * catalogs_match[i-1-j]["b"][idx][separation_mask].data
                areas_new = catalog["a"][separation_mask].data * catalog["b"][separation_mask].data
                #to_keep = np.argmax( [catalogs_match[i-1-j]["MAG_AUTO"][idx][separation_mask].data, catalog["MAG_AUTO"][separation_mask].data], axis=0 )
                to_keep = np.argmax( [areas_previous, areas_new], axis=0 )
                #print("len(np.where(d2d < 1.*u.arcsec)[0]) = ", len(np.where(d2d < 1.*u.arcsec)[0]))
                to_keep_previous = idx[separation_mask][ to_keep==0 ]
                to_remove_previous = idx[separation_mask][ to_keep==1 ]
                to_keep_new = np.where(separation_mask)[0][ to_keep==1 ]
                to_remove_new = np.where(separation_mask)[0][ to_keep==0 ]
                
                #print("to_remove_new:", to_remove_new)
                #print("to_keep_new:", to_keep_new)
                #print("to_remove_previous:", to_remove_previous)
                #print("to_keep_previous:", to_keep_previous)
                catalog = vstack( [catalog[to_keep_new], catalogs_match[i-1-j][to_keep_previous]], join_type='outer' )
                #print('######## after vstack ######')
                #print(np.where(np.logical_not(np.isnan( catalog['magAB_F814W'])) & np.logical_not(np.isnan( catalog['magAB_F435W'])))[0].shape)
                #print('##############')
                catalogs_match[i-1-j].remove_rows(np.concatenate((to_keep_previous, to_remove_previous)))
                
        catalogs_match.append(catalog)
        print("length of catalogues at step ", i, " after removal: ", [len(cat['RA']) for cat in catalogs_match])
        
    print('initial length: ' + str(initial_len))
    combined_catalog = vstack(catalogs_match)
    print('final length: ' + str(len(combined_catalog)))
    print( str(initial_len - len(combined_catalog)) + ' matched sources removed' )
    return combined_catalog
"""






















