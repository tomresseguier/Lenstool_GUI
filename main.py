import argparse
import os
from fits_image import fits_image
from photometry_catalog import photometry_catalog

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






image_path = '../../spt0615/DATA/v7-wisps/color/' + 'spt0615_color_rgb.fits'
photometry_catalog_path = '../../spt0615/DATA/v7-deblend/catalogs/' + 'spt0615_phot-eazy.cat'



image = fits_image(image_path)
photometry_catalog = photometry_catalog(photometry_catalog_path)


### Extract galaxy shapes ###
#pysex.run()

