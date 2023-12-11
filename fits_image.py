import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Polygon, Circle
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import *
from astropy.coordinates import SkyCoord
from astropy.visualization import ZScaleInterval, ImageNormalize

#import sys
#module_dir = os.path.dirname(os.path.abspath(__file__))
#sys.path.append(module_dir)

from source_extraction.source_extract import source_extract



class fits_image() :
    def __init__(self, image_path) :
        self.image_path = image_path
        def open_image(image_path) :
            ### OPENING THE FITS FILE ###
            with fits.open(image_path) as hdus :
                ### Finding the header ###
                for hdu in hdus :
                    if isinstance(hdu.header, fits.header.Header) :
                        header = hdu.header
                        break
                wcs = WCS(hdus[0].header)
                ### Finding the pixel scale ###
                #if 'CD1_1' in hdus[0].header.keys() :
                #    CD1_1 = hdus[0].header['CD1_1']
                #    CD1_2 = hdus[0].header['CD1_2']
                #    pix_deg_scale = np.sqrt(CD1_1**2+CD1_2**2)
                #elif 'CDELT1' in hdus[0].header.keys() :
                #    pix_deg_scale = abs(hdus[0].header['CDELT1'])
                #else :
                #    pix_deg_scale = input('Pixel scale not found in header. Please provide manually in degrees:')
                pix_deg_scale = np.sqrt(wcs.pixel_scale_matrix[0, 0]**2+wcs.pixel_scale_matrix[0, 1]**2)
                ### Finding the data ###
                data = []
                for hdu in hdus :
                    if isinstance(hdu.data, np.ndarray) :
                        data.append(hdu.data)
                if len(data)<=2 :
                    image = data[0]
                if len(data)>=3 :
                    data_red = data[0]
                    data_green = data[1]
                    data_blue = data[2]
                    image = np.dstack((data_red, data_green, data_blue))
            return pix_deg_scale, wcs, image
        self.pix_deg_scale, self.wcs, self.image = open_image(image_path)
        self.sources = None
        self.fig = None
        self.ax = None
        self.multiple_images = None
        self.galaxy_selection = None
    
    def plot_image(self, fig=None, wcs_projection=True, units='pixel', pos=111, make_axes_labels=True) :
        if fig is None :
            fig = plt.figure()
        if wcs_projection :
            ax = fig.add_subplot(pos, projection=self.wcs)
            ax.coords.grid(True, color='black', ls='dotted')
        else :
            ax = fig.add_subplot(pos)
            if units=='pixel' or units=='pixels' or units=='image' :
                scaling = 1
            if units=='arcsec' :
                scaling = pix_deg_scale*60*60
            if units=='arcmin' :
                scaling = pix_deg_scale*60
        if make_axes_labels and wcs_projection :
            ax.coords[0].set_axislabel('Right ascension')
            ax.coords[1].set_axislabel('Declination')
        elif make_axes_labels and not wcs_projection :
            ax.set_xlabel('x (' + units + ')')
            ax.set_ylabel('y (' + units + ')')
        elif not make_axes_labels and wcs_projection :
            ax.coords[0].set_axislabel(' ')
            ax.coords[1].set_axislabel(' ')
        else :
            ax.set_xlabel(' ')
            ax.set_ylabel(' ')
        if wcs_projection :
            ax.imshow(self.image, origin="lower")
        if not wcs_projection :
            ax.imshow(self.image, origin='lower', extent=[0, self.image.shape[1]*scaling, 0, self.image.shape[0]*scaling])
        #ax.figure.tight_layout()
        return fig, ax        
    
    def extract_sources(self) :
        self.sources = source_extract(self.image_path, weight_path=None, zero_point=None,
                       out_dir='PWD', outfile_name='SExtractor_cat', return_sources=True)
        return self.sources
    
    def select_multiple_images




    
    
    