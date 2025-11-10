import os
import sys
import numpy as np
from lenstool_gui.fits_image import fits_image


DATA_dir = 'path/to/data/'



# Red sequence
im = fits_image(DATA_dir + "macs0308_rgb.fits")
im.boost()
phot_cat_path = DATA_dir + "macs0308_phot-eazy.cat"
im.import_catalog(phot_cat_path)


im.imported_cat.plot()
im.imported_cat.plot_column('z_phot')

im.extract_sources()
im.sources.cat
im.sources.plot()

im.imported_cat.transfer_col('a', which_cat='sources')
im.imported_cat.transfer_col('b', which_cat='sources')
im.imported_cat.transfer_col('theta', which_cat='sources')

im.imported_cat.cat['a'] = im.imported_cat.cat['a_CAT2']
im.imported_cat.cat['b'] = im.imported_cat.cat['b_CAT2']
im.imported_cat.cat['theta'] = im.imported_cat.cat['theta_CAT2']


im.imported_cat.plot()



def add_magnitude_column(catalog, flux_col):
    flux = catalog[flux_col]
    valid_flux = flux > 0
    magnitude = np.full(len(flux), np.nan)  # Initialize with NaNs
    magnitude[valid_flux] = -2.5 * np.log10(flux[valid_flux])
    catalog[flux_col[:-len('flux')]+'mag'] = magnitude
    return catalog

add_magnitude_column(im.imported_cat.cat, 'f200w_flux')
add_magnitude_column(im.imported_cat.cat, 'f105w_flux')

im.imported_cat.mag_colnames = im.imported_cat.cat.colnames[-2:]
im.imported_cat.plot_selection_panel()
im.imported_cat.plot()



im.imported_cat.export_to_potfile()




selection_mask = im.imported_cat.seselection_mask
RS_catalog = im.imported_cat.cat[selection_mask]






# Multiple images
im.plot_image()
im.import_catalog(phot_cat_path)
im.imported_cat.plot()








