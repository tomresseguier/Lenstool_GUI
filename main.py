import argparse
import os
from fits_image import fits_image
from photometry_catalog import photometry_catalog
#from astropy import units as u

#from Lenstool_GUI.pyextract import pysex



parser = argparse.ArgumentParser(prog='Lenstool_GUI',
                                 description='A GUI for Lenstool',
                                 epilog='#######################')


parser.add_argument('FITS_image')
parser.add_argument('photometry_catalog')
#parser.add_argument('-w', '--weight_file')      # option that takes a value
#weight_file = args.weight_file


args = parser.parse_args()


### Check that the files exist ###
if os.path.exists( os.path.join( os.path.dirname(os.path.abspath(__file__)), args.FITS_image ) ):
    image_path = os.path.join( os.path.dirname(os.path.abspath(__file__)), args.FITS_image )
elif os.path.exists( args.FITS_image ) :
    image_path = args.FITS_image
else :
    print('FITS image not found')

if os.path.exists( os.path.join( os.path.dirname(os.path.abspath(__file__)), args.photometry_catalog ) ):
    photometry_catalog_path = os.path.join( os.path.dirname(os.path.abspath(__file__)), args.photometry_catalog )
elif os.path.exists( args.photometry_catalog ) :
    photometry_catalog_path = args.photometry_catalog
else :
    print('Photometry catalog not found')
##################################






RGB_image_path = '../../spt0615/DATA/v7-wisps/color/' + 'spt0615_color_rgb.fits'
F444W_image_path = '/Users/Margaux/Desktop/RESEARCH/spt0615/DATA/background_subtracted/j061548m5745-grizli-v7.0-f444w-clear_drc_sci_bkgsub.fits'

photometry_catalog_path = '../../spt0615/DATA/v7-deblend/catalogs/' + 'spt0615_phot-eazy.cat'
mult_file_path = '../../spt0615/spt0615_process/SL_models/' + 'spt0615_mult_images.lenstool'

potfile_path = '../../spt0615/spt0615_process/SL_models/galcat.cat'


image = fits_image(RGB_image_path)
#image.plot_image()
image.import_catalog(photometry_catalog_path)
image.imported_cat.plot()
image.imported_cat.plot_selection_panel()
image.imported_cat.plot()

image.extract_sources(F444W_image_path)
image.import_multiple_images(mult_file_path)
image.multiple_images.plot()
image.load_potfile(potfile_path)
image.imported_cat.plot(color=[0,1,0])
image.potfile.plot(color=[0,0,0])






RGB_image_path = '/Users/Margaux/Desktop/RESEARCH/Planck_bullet_cluster/DATA/Treated_data/images/total/' + 'hlsp_16429_hst_acs-wfc3-60mas_plckg282+49_total_drz.fits'
mult_file_path = '/Users/Margaux/Desktop/RESEARCH/Planck_bullet_cluster/PBC_process/Lenstool/SL_only/RUN_045_forme4/' + 'PBC_mult_images.lenstool'

image2 = fits_image(RGB_image_path)
image2.import_multiple_images(mult_file_path)
image2.multiple_images.plot()



#photometry_catalog = photometry_catalog(photometry_catalog_path)


### Extract galaxy shapes ###
#pysex.run()

