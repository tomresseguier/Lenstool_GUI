from astropy.table import Table, vstack
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy import spatial



def vstack_catalogs(catalogs_paths):
    """
    Vertically stacks multiple FITS catalogs into a single Astropy Table.

    Parameters
    ----------
    catalogs_paths : list of str
        A list of file paths to the FITS catalogs that need to be stacked.

    Returns
    -------
    astropy.table.Table
        A single table containing all the sources from the input catalogs.
    """
    catalogs = [Table.read(path, format='fits') for path in catalogs_paths]
    for i, catalog in enumerate(catalogs):
        print(f"length of catalog {i}: {len(catalog)}")
    return vstack(catalogs)


def match_cat2(cats, match_radius=0.5, fill_in_value=np.nan, return_match_idx=False, keep_all_col=False, column_to_transfer=None):
    """
    Matches two catalogs based on their celestial coordinates and transfers columns from one to the other.

    Parameters
    ----------
    cats : list of (str or astropy.table.Table)
        A list containing two catalogs. The first is the receiver catalog and the second is the giver.
    match_radius : float, optional
        The radius in arcseconds within which to consider sources as matched. Default is 0.5.
    fill_in_value : float, optional
        The value to use for non-matched sources in the transferred columns. Default is np.nan.
    return_match_idx : bool, optional
        If True, returns the indices of the matched sources in the giver catalog. Default is False.
    keep_all_col : bool, optional
        If True, appends '_CAT2' to the column name in case of a name conflict. Default is False.
    column_to_transfer : list of str, optional
        A list of column names to transfer. If None, all columns are transferred. Default is None.

    Returns
    -------
    astropy.table.Table or tuple
        The matched catalog. If `return_match_idx` is True, it returns a tuple of the catalog and the match indices.
    """
    catalogs = []
    for cat in cats:
        if isinstance(cat, str):
            catalogs.append(Table.read(cat, format='fits'))
        else:
            catalogs.append(cat.copy())
    cat_receiver, cat_giver = catalogs

    ra_dec_names = [('RA', 'DEC'), ('ra', 'dec')]
    ra_dec_names_receiver = next(((ra, dec) for ra, dec in ra_dec_names if ra in cat_receiver.colnames), (None, None))
    ra_dec_names_giver = next(((ra, dec) for ra, dec in ra_dec_names if ra in cat_giver.colnames), (None, None))

    if not ra_dec_names_receiver[0] or not ra_dec_names_giver[0]:
        raise ValueError("RA/DEC columns not found in one of the catalogs.")
    
    coords_cat_receiver = SkyCoord(cat_receiver[ra_dec_names_receiver[0]], cat_receiver[ra_dec_names_receiver[1]], unit=u.deg)
    coords_cat_giver = SkyCoord(cat_giver[ra_dec_names_giver[0]], cat_giver[ra_dec_names_giver[1]], unit=u.deg)
    
    idx, d2d, _ = coords_cat_receiver.match_to_catalog_sky(coords_cat_giver)
    match_mask = d2d < match_radius * u.arcsec
    
    matched_cat = cat_receiver.copy()
    
    columns_to_process = cat_giver.colnames if column_to_transfer is None else column_to_transfer

    for colname in columns_to_process:
        col_type = type(cat_giver[colname][0])
        if not np.isnan(fill_in_value) :
            fill_in_value_type = col_type if col_type!=np.ndarray else type(cat_giver[colname][0][0])
            fill_in_value_matched_type = fill_in_value_type(fill_in_value)
        else :
            fill_in_value_matched_type = fill_in_value
        if colname not in matched_cat.colnames:
            print(f'Adding column: {colname}')
            if col_type == np.ndarray :
                matched_cat[colname] = [np.full(len(cat_giver[colname][0]), fill_in_value_matched_type) for i in range(len(matched_cat))] #np.full(len(matched_cat), np.full(len(cat_giver[colname][0]), fill_in_value))
            else :
                matched_cat[colname] = np.full(len(matched_cat), fill_in_value_matched_type)
            matched_cat[colname][match_mask] = cat_giver[colname][idx[match_mask]]
        elif keep_all_col:
            new_colname = f"{colname}_CAT2"
            print(f'Adding column: {new_colname}')
            if col_type == np.ndarray :
                matched_cat[new_colname] = [np.full(len(cat_giver[colname][0]), fill_in_value_matched_type) for i in range(len(matched_cat))] #np.full(len(matched_cat), np.full(len(cat_giver[colname][0]), fill_in_value))
            else :
                matched_cat[new_colname] = np.full(len(matched_cat), fill_in_value_matched_type)
            matched_cat[new_colname][match_mask] = cat_giver[colname][idx[match_mask]]

    if return_match_idx:
        return matched_cat, idx
    return matched_cat


def _get_filled_and_empty_columns(catalog):
    if catalog.mask is None:
        return catalog.colnames, []
    
    empty_columns = [name for name in catalog.colnames if np.all(catalog.mask[name])]
    filled_columns = [name for name in catalog.colnames if name not in empty_columns]
    
    return filled_columns, empty_columns


def combine_catalogs_remove_double_detections(cats, match_radius=0.5, preference_criterion='area'):
    """
    Combines a list of catalogs, removing double detections based on a preference criterion.

    Parameters
    ----------
    cats : list of (str or astropy.table.Table)
        A list of catalogs to combine.
    match_radius : float, optional
        The matching radius in arcseconds. Default is 0.5.
    preference_criterion : str, optional
        The criterion to decide which detection to keep ('area' or other). Default is 'area'.

    Returns
    -------
    astropy.table.Table
        The combined catalog with double detections removed.
    """
    catalogs = [Table.read(cat, format='fits') if isinstance(cat, str) else cat.copy() for cat in cats]
    
    initial_len = sum(len(cat) for cat in catalogs)
    
    catalog_previous = catalogs[0]

    for i in range(1, len(catalogs)):
        catalog_new = catalogs[i]
        print(f'##############\nStep {i}')
        print(f"length of catalog {i-1} before removal: {len(catalog_previous)}")
        print(f"length of catalog {i} before removal: {len(catalog_new)}")

        coords_previous = SkyCoord(catalog_previous['RA'], catalog_previous['DEC'], unit=u.deg)
        coords_new = SkyCoord(catalog_new['RA'], catalog_new['DEC'], unit=u.deg)
        idx, d2d, _ = coords_previous.match_to_catalog_sky(coords_new)
        separation_mask = d2d < match_radius * u.arcsec

        filled_columns_previous, _ = _get_filled_and_empty_columns(catalog_previous)
        filled_columns_new, _ = _get_filled_and_empty_columns(catalog_new)

        from_previous_to_new = set(filled_columns_previous) - set(filled_columns_new)
        from_new_to_previous = set(filled_columns_new) - set(filled_columns_previous)

        for col in from_new_to_previous:
            if col not in catalog_previous.colnames:
                catalog_previous[col] = np.full(len(catalog_previous), np.nan)
            catalog_previous[col][separation_mask] = catalog_new[col][idx[separation_mask]]

        for col in from_previous_to_new:
            if col not in catalog_new.colnames:
                catalog_new[col] = np.full(len(catalog_new), np.nan)
            catalog_new[col][idx[separation_mask]] = catalog_previous[col][separation_mask]

        if preference_criterion == 'area':
            areas_previous = catalog_previous["a"][separation_mask] * catalog_previous["b"][separation_mask]
            areas_new = catalog_new["a"][idx[separation_mask]] * catalog_new["b"][idx[separation_mask]]
            to_keep = np.argmax([areas_previous, areas_new], axis=0)

        where_separation = np.where(separation_mask)[0]
        to_remove_previous = where_separation[to_keep == 1]
        to_remove_new = idx[separation_mask][to_keep == 0]

        catalog_previous.remove_rows(to_remove_previous)
        catalog_new.remove_rows(to_remove_new)

        print(f"length of catalog {i-1} after removal: {len(catalog_previous)}")
        print(f"length of catalog {i} after removal: {len(catalog_new)}")

        catalog_previous = vstack([catalog_previous, catalog_new], join_type='outer')

    final_len = len(catalog_previous)
    print(f'initial length: {initial_len}')
    print(f'final length: {final_len}')
    print(f'{initial_len - final_len} matched sources removed')
    return catalog_previous


def combine_DIM_catalogs(cats):
    """
    Combines Double Image Mode catalogs from SourceExtractor by taking the mean of columns with incomplete data.

    Parameters
    ----------
    cats : list of (str or astropy.table.Table)
        A list of catalogs to be combined.

    Returns
    -------
    astropy.table.Table
        The combined catalog.
    """
    catalogs = [Table.read(cat, format='fits') if isinstance(cat, str) else cat.copy() for cat in cats]
    
    col_info = {}
    for i, cat in enumerate(catalogs):
        for colname in cat.colnames:
            if colname not in col_info:
                col_info[colname] = {'complete': [], 'incomplete': []}
            
            if np.isnan(cat[colname]).any():
                col_info[colname]['incomplete'].append(i)
            else:
                col_info[colname]['complete'].append(i)

    combined_cat = Table()
    for colname, info in col_info.items():
        if info['complete']:
            combined_cat[colname] = catalogs[info['complete'][0]][colname]
        else:
            incomplete_cols = np.array([catalogs[i][colname] for i in info['incomplete']])
            combined_cat[colname] = np.nanmean(incomplete_cols, axis=0)
            
    return combined_cat


def remove_object(cat, FWHM_to_radius=1):
    """
    Removes objects that are too close to each other, keeping the one with the larger FWHM.

    Parameters
    ----------
    cat : str or astropy.table.Table
        The input catalog.
    FWHM_to_radius : float, optional
        A factor to convert FWHM to a radius for overlap checking. Default is 1.

    Returns
    -------
    astropy.table.Table
        The catalog with overlapping objects removed.
    """
    if isinstance(cat, str):
        cat = Table.read(cat, format='fits')

    print(f"num of objects in the catalogue: {len(cat)}")

    x = np.array(cat['X_IMAGE'])
    y = np.array(cat['Y_IMAGE'])
    fwhm = np.array(cat['FWHM_IMAGE']) * FWHM_to_radius
    
    ori_order = np.arange(len(x))
    
    # Sort by FWHM in descending order
    sort_indices = np.argsort(fwhm)[::-1]
    x_sorted, y_sorted, fwhm_sorted, ori_order_sorted = x[sort_indices], y[sort_indices], fwhm[sort_indices], ori_order[sort_indices]

    pos = np.column_stack((x_sorted, y_sorted))
    tree = spatial.KDTree(pos)
    
    keep_mask = np.ones(len(x), dtype=bool)
    
    for i in range(len(x)):
        if not keep_mask[i]:
            continue
            
        # Find neighbors within 2 * FWHM
        neighbors_idx = tree.query_ball_point(pos[i], 2.0 * fwhm_sorted[i])
        
        for index in neighbors_idx:
            if i == index:
                continue
            
            distance = np.hypot(pos[index, 0] - pos[i, 0], pos[index, 1] - pos[i, 1])
            if distance < (fwhm_sorted[i] + fwhm_sorted[index]):
                keep_mask[index] = False

    final_indices = ori_order_sorted[keep_mask]
    
    cat_filtered = cat[final_indices]
    
    print(f"Num of objects after removing double-detection: {len(cat_filtered)}")
    return cat_filtered




















