import os
import sys
import numpy as np
from visualens import Visualens


DATA_dir = 'path/to/data/'



# Red sequence
vl = Visualens()
vl.import_image(DATA_dir + "macs0308_rgb.fits")
vl.image.boost()
phot_cat_path = DATA_dir + "macs0308_phot-eazy.cat"
vl.import_catalog(phot_cat_path)


vl.catalog.plot()
vl.catalog.plot_column('z_phot')

vl.image.extract_sources()
vl.image.sources.cat
vl.image.sources.plot()

vl.catalog.transfer_col('a', which_cat='sources')
vl.catalog.transfer_col('b', which_cat='sources')
vl.catalog.transfer_col('theta', which_cat='sources')

vl.catalog.cat['a'] = vl.catalog.cat['a_CAT2']
vl.catalog.cat['b'] = vl.catalog.cat['b_CAT2']
vl.catalog.cat['theta'] = vl.catalog.cat['theta_CAT2']


vl.catalog.plot()



def add_magnitude_column(catalog, flux_col):
    flux = catalog[flux_col]
    valid_flux = flux > 0
    magnitude = np.full(len(flux), np.nan)  # Initialize with NaNs
    magnitude[valid_flux] = -2.5 * np.log10(flux[valid_flux])
    catalog[flux_col[:-len('flux')]+'mag'] = magnitude
    return catalog

add_magnitude_column(vl.catalog.cat, 'f200w_flux')
add_magnitude_column(vl.catalog.cat, 'f105w_flux')

vl.catalog.mag_colnames = vl.catalog.cat.colnames[-2:]
vl.catalog.plot_selection_panel()
vl.catalog.plot()



vl.catalog.export_to_potfile()




selection_mask = vl.catalog.seselection_mask
RS_catalog = vl.catalog.cat[selection_mask]




# Multiple images
vl.image.plot_image()
vl.import_catalog(phot_cat_path)
vl.catalog.plot()
