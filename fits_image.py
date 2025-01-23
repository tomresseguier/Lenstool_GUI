import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Polygon, Circle, Rectangle
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import astropy.constants as c
from reproject import reproject_interp
#from astropy.visualization.wcsaxes import *
from astropy.coordinates import SkyCoord
from astropy.visualization import ZScaleInterval, ImageNormalize
from matplotlib.patches import Ellipse, Polygon, Circle
from tqdm import tqdm
import requests
import types

import PyQt5
from PyQt5.QtWidgets import QMainWindow, QWidget, QHBoxLayout
import pyqtgraph as pg
from PyQt5.QtWidgets import QGraphicsEllipseItem, QGraphicsScene, QGraphicsView, QGraphicsSceneMouseEvent, QApplication
#from PyQt5.QtCore import Qt
from pyqtgraph.Qt import QtCore
from astropy.table import Table
import pandas as pd

###############################################################################
import sys
#module_dir = os.path.dirname(os.path.abspath(__file__))
#sys.path.append(module_dir)

from .source_extraction.source_extract import source_extract, source_extract_DIM
from .source_extraction.match_cat import run_match
#sys.path.append(module_dir + "/utils")
from .utils.utils_plots.plot_utils_general import *
from .utils.utils_Qt.selectable_classes import *
from .utils.utils_Qt.utils_general import *
from .utils.utils_general.utils_general import find_close_coord
from .utils.utils_plots.plt_framework import plt_framework
from .utils.utils_Lenstool.redshift_extractors import make_source_z_dict, find_param_file
from .utils.utils_Lenstool.param_extractors import read_potfile, read_bayes_file, make_param_latex_table
#from utils_general.utils import flux_muJy_to_magAB
from .utils.utils_astro.set_cosmology import set_cosmo
cosmo = set_cosmo()
###############################################################################

"""
from source_extraction.source_extract import source_extract
from utils.utils_plots.plot_utils_general import *
from utils.utils_Qt.selectable_classes import *
from utils.utils_Qt.utils_general import *
from utils.utils_general.utils import flux_muJy_to_magAB
"""


pg.setConfigOption('imageAxisOrder', 'row-major')




class fits_image :
    def __init__(self, image_path) :
        self.image_path = image_path
        if os.path.isfile(self.image_path[:-8] + 'wht.fits') :
            print("Weight file found: " + self.image_path[:-8] + 'wht.fits')
            self.weight_path = self.image_path[:-8] + 'wht.fits'
        else :
            self.weight_path = None
        self.image_data, self.pix_deg_scale, self.orientation, self.wcs, self.header = self.open_image()
        ### Following lines useless, just to see the possible attributes of the class ###
        self.sources = None
        self.fig = None
        self.ax = None
        self.multiple_images = None
        self.galaxy_selection = None
        self.potfile = None
        self.imported_cat = None
        self.qt_plot = None
        #self.qt_plot = self.plot_image()
        self.qtItems_dict = {'sources': None,
                             'potfile_cat': None,
                             'imported_cat': None,
                             'multiple_images': None}
        self.ax = None
        self.redshift = None
    
    def open_image(self) :
        with fits.open(self.image_path) as hdus :
            if isinstance(hdus, fits.hdu.hdulist.HDUList) :
                image_hdus = []
                for hdu in hdus :
                    if hdu.is_image and isinstance(hdu.data, np.ndarray) :
                        image_hdus.append(hdu)
                print(f"{len(image_hdus)} image hdus found.")
                
                if len(image_hdus)==0 :
                    print('No image found in FITS file.')
                
                true_image_hdus = []
                for i, image_hdu in enumerate(image_hdus):
                    if image_hdu.data is not None and len(image_hdu.data.shape) >= 2:
                        print(f"HDU {i} contains image data with shape {hdu.data.shape}")
                        true_image_hdus.append(image_hdu)
                print(f"{len(true_image_hdus)} non empty images found.")
                
                
                keyword = None
                if True in np.unique(['EXTNAME' in true_image_hdu.header for true_image_hdu in true_image_hdus]) :
                    keyword = 'EXTNAME'
                elif True in np.unique(['FILETYPE' in true_image_hdu.header for true_image_hdu in true_image_hdus]) :
                    keyword = 'FILETYPE'
                if keyword is None :
                    print("No 'SCI' extname found.")
                    selected_hdus = true_image_hdus
                else :
                    #sci_hdus = []
                    selected_hdus = []
                    wht_hdus = []
                    for true_image_hdu in true_image_hdus :
                        if keyword in true_image_hdu.header :
                            if true_image_hdu.header[keyword]=='SCI' :
                                #sci_hdus.append(true_image_hdu)
                                selected_hdus.append(true_image_hdu)
                                print('Science image found.')
                            if true_image_hdu.header[keyword]=='WHT' :
                                wht_hdus.append(true_image_hdu)
                                print('Weight image found.')
                
                
                
                if len(selected_hdus)==3 :
                    x_sizes = [selected_hdu.data.shape[0] for selected_hdu in selected_hdus]
                    y_sizes = [selected_hdu.data.shape[1] for selected_hdu in selected_hdus]
                    if len(np.unique(x_sizes))==1 and len(np.unique(y_sizes))==1 :
                        print("Assuming RGB data.")
                    data_red = selected_hdus[0].data
                    data_green = selected_hdus[1].data
                    data_blue = selected_hdus[2].data
                    image = np.dstack((data_red, data_green, data_blue))
                else :
                    print("Using first hdu.")
                    image = selected_hdus[0].data
                
                
                wcs = WCS(selected_hdus[0].header)
                header = selected_hdus[0].header
                
                if 'ORIENTAT' in header :
                    orientation = header['ORIENTAT']
                elif 'CD1_1' in header and 'CD1_2' in header and 'CD2_1' in header and 'CD2_2' in header :
                    cd = np.array([[header['CD1_1'], header['CD1_2']], [header['CD2_1'], header['CD2_2']]])
                    #det = np.linalg.det(cd)
                    #sign = np.sign(det)
                    orientation = np.arctan2(cd[1,0], cd[1,1])
                elif 'PC1_1' in header and 'PC1_2' in header and 'PC2_1' in header and 'PC2_2' in header :
                    cd = np.array([[header['PC1_1'], header['PC1_2']], [header['PC2_1'], header['PC2_2']]])
                    #det = np.linalg.det(cd)
                    #sign = np.sign(det)
                    orientation = np.arctan2(cd[1,0], cd[1,1])
                else :
                    orientation = None
                
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
                
                return image, pix_deg_scale, np.rad2deg(orientation), wcs, header
            
            else :
                print('Unable to extract image data from FITS file')
                return None
                
        
    
    def plot_image(self) :
        #if self.qt_plot is None :
        #to_plot = self.image_data
        #to_plot = np.transpose(self.image_data, axes=[1,0,2])
        to_plot = np.flip(self.image_data, axis=0)
        
        #self.qt_plot = pg.image(to_plot)
        
        self.qt_plot = pg.ImageView()
        self.qt_plot.setImage(to_plot)
        self.qt_plot.autoLevels()
        
        self.image_widget_layout = QHBoxLayout()
        self.image_widget_layout.addWidget(self.qt_plot)
        
        #self.image_widget = QWidget()
        self.image_widget = DragWidget(self.qt_plot)
        self.image_widget.setLayout(self.image_widget_layout)
        
        self.window = QMainWindow()
        self.window.setWindowTitle(os.path.basename(self.image_path))
        self.window.setCentralWidget(self.image_widget)
        self.window.show()
        
        #win = pg.GraphicsLayoutWidget()
        #win.setWindowTitle(os.path.basename(image.image_path))
        #win.show()
        
        #win_element = win.addPlot()
        #win_element.addItem(image.qt_plot, row=0, col=0)
        ########################################
        
        """
        self.qt_plot = pg.ImageItem(to_plot)
        self.window = pg.GraphicsLayoutWidget()
        plot_element = self.window.addPlot(title=os.path.basename(self.image_path))
        plot_element.addItem(self.qt_plot)
        self.window.show()
        """
        
        return self.qt_plot
    
    def plot_image_mpl(self, wcs_projection=True, units='pixel', pos=111, make_axes_labels=True, make_grid=True, crop=None, replace_image=True, extra_pad=None) :
        fig, ax = plot_image_mpl(self.image_data, wcs=self.wcs, wcs_projection=wcs_projection, units=units, pos=pos, \
                                 make_axes_labels=make_axes_labels, make_grid=make_grid, crop=crop, extra_pad=extra_pad)
        plot_NE_arrows(ax, self.wcs)
        if replace_image :
            self.fig, self.ax = fig, ax
        return fig, ax
    
    def set_weight(self, weight_path) :
        self.weight_path = weight_path
    
    def extract_sources(self, image_path=None, weight_path=None, DIM_ref_path=None, rerun=False, reproject=True) :
        if image_path is None :
            image_path = self.image_path
            weight_path = self.weight_path
            
        out_dir = os.path.join( os.path.dirname(image_path), 'source_extraction/' )
        if not os.path.exists(out_dir) :
            os.mkdir(out_dir)
        
        outfile_name = 'SExtractor_cat.fits' if DIM_ref_path is None else 'SExtractor_cat_DIM.fits'
        out_path = os.path.join( out_dir, outfile_name )
        if os.path.isfile(out_path) and not rerun :
            print('Previous SExtractor catalog found.')
            with fits.open(out_path) as hdu :
                self.sources_all = Table(hdu[1].data)
        elif DIM_ref_path is None :
            self.sources_all = source_extract(image_path, weight_path=weight_path, pixel_scale=self.pix_deg_scale*3600, zero_point=None, out_dir=out_dir,
                                              outfile_name=outfile_name, return_sources=True)
        else :
            if type(DIM_ref_path) is list :
                if reproject :
                    reprojected_image_path = self.reproject(DIM_ref_path[0], image_path)
                    reprojected_weight_path = self.reproject(DIM_ref_path[0], weight_path)
                    reprojected_image_path = [reprojected_image_path, reprojected_weight_path]
                else :
                    reprojected_image_path = [image_path, weight_path]
            else :
                if reproject :
                    reprojected_image_path = self.reproject(DIM_ref_path, image_path)
                else :
                    reprojected_image_path = image_path
            
            self.sources_all = source_extract_DIM(DIM_ref_path, reprojected_image_path, pixel_scale=self.pix_deg_scale*3600, zero_point=None, out_dir=out_dir,
                                                  outfile_name=outfile_name, return_sources=True)
        
        # This next part doesn't make sense here as the purpose of extract_sources() is to extract sources from the imported image,
        # but the code should be added to import_cat()/make_catalog()
        # THE WAY TO CONVERT ANGLES IS MORE COMPLICATED!!!
        x, y = self.world_to_image(self.sources_all['RA'], self.sources_all['DEC'], unit='deg')
        if x[0] != self.sources_all['X_IMAGE'][0] :
            #print('Catalog sextracted from different image: replacing X_IMAGE, Y_IMAGE and THETA_IMAGE columns with current image coordinates.')
            #self.sources_all['X_IMAGE'], self.sources_all['Y_IMAGE'] = x, y
            print('Catalog SExtracted from different image: keeping X_IMAGE, Y_IMAGE and THETA_IMAGE and adding x, y, theta columns from current image coordinates.')
            self.sources_all.add_column(x, name='x')
            self.sources_all.add_column(y, name='y')
            ref_image_angle = ( np.arctan2(self.wcs.wcs.get_pc()[1, 0], self.wcs.wcs.get_pc()[0, 0]) %np.pi ) * 360/np.pi
            print('Overall angle of imported image: ' + str(ref_image_angle))
            #self.sources_all.add_column(self.sources_all['THETA_WORLD'] + ref_image_angle, name='prout')
        
        #mask_mag = self.sources_all['MAG_AUTO']<-10.
        #mask = self.sources['KRON_RADIUS']
        #mask_galstar = self.sources_all['CLASS_STAR']<0.4
        #mask_size = self.sources_all['A_IMAGE']*self.sources_all['B_IMAGE']*np.pi>1000.
        #mask = mask_mag & mask_galstar & mask_size
        
        #self.sources = self.sources_all #[mask]
        self.make_photometry(self.sources_all)
        self.sources = self.make_catalog(self.sources_all)
        return str(len(self.sources.cat)) + ' sources found.'
    
    def reproject(self, ref_image_path, image_path) :
        reprojected_image_path = image_path[:-len('.fits')] + '_reprojected.fits'
        if os.path.isfile(reprojected_image_path) :
            print('Previous reprojected image found.')
        else :
            with fits.open(ref_image_path) as hdu :
                reference_header = hdu[0].header
            with fits.open(image_path) as hdu :
                print('Reprojecting image ' + image_path + ' onto reference ' + ref_image_path)
                reprojected_data, footprint = reproject_interp(hdu[0], reference_header)
                fits.writeto(reprojected_image_path, reprojected_data, reference_header)
        return reprojected_image_path
    
    def make_photometry(self, cat) :
        print("################")
        print("Figuring out the photometry:")
        print("Assuming instrument HST/ACS")
        print("PHOTFLAM = " + str(self.header['PHOTFLAM']))
        print("Pivot wavelength = " + str(self.header['PHOTPLAM']))
        print("################")
        flux_lambda = self.header['PHOTFLAM'] * cat['FLUX_ISO'] * u.erg/u.cm**2/u.s/u.AA #FLUX_AUTO, FLUX_ISO, FLUX_APER
        magST = -2.5*np.log10(flux_lambda.value) - 21.1
        pivot_wavelength = self.header['PHOTPLAM'] * u.AA
        flux_nu = flux_lambda.to(u.erg/u.cm**2/u.s/u.Hz, u.spectral_density(pivot_wavelength))
        magAB = -2.5*np.log10(flux_nu.value) - 48.6
        if 'FILTER2' in self.header.keys() :
            print("filter " + self.header['FILTER2'] + " found")
            cat.add_column(magAB, name='magAB_' + self.header['FILTER2'])
            cat.add_column(magST, name='magST_' + self.header['FILTER2'])
        else :
            cat.add_column(magAB, name='magAB')
            cat.add_column(magST, name='magST')
        return 'Magnitudes calculated'
    
    def select_multiple_images(self) :
        return 'in progress'
    
    def import_multiple_images(self, mult_file_path, unit_is_pixel=False):
        
        def make_masks_and_colors(cat, which='all', filled_markers=False) :
            if which=='all':
                families = np.unique([cat[i]['id'][0] for i in range(len(cat))])
                which = families
            to_plot_mask = []
            for symbol in which:
                to_plot_mask.append([cat[i]['id'].startswith(symbol) for i in range(len(cat))])
            if filled_markers:
                colors = make_palette(len(which), 1, alpha=0.5)
            else:
                colors = make_palette(len(which), 1, alpha=0)
            colors_dict = {}
            for i, mask in enumerate(to_plot_mask):
                for multiple_image in cat[mask]:
                    colors_dict[multiple_image['id']] = colors[i]
            return to_plot_mask, colors_dict, colors
        
        
        multiple_images = Table(names=['id','ra','dec','a','b','theta','z','mag'], dtype=['str',*['float',]*7])
        with open(mult_file_path, 'r') as mult_file:
            for line in mult_file:
                cleaned_line = line.strip()
                if not cleaned_line.startswith("#") and len(cleaned_line)>0 :
                    split_line = cleaned_line.split()
                    row = [split_line[0]]
                    for element in split_line[1:8] :
                        row.append(float(element))
                    multiple_images.add_row(row)
        
        multiple_images['theta'] = multiple_images['theta'] - self.orientation
        
        def plot_multiple_images(self, which='all', size=40, marker='o', filled_markers=False, colors=None, mpl=False, fontsize=9,
                               make_thumbnails=False, square_size=150, margin=50, distance=200, savefig=False, square_thumbnails=True,
                               boost=[2,1.5,1], linewidth=1.7, text_color='white', text_alpha=0.5):
            
            # Get colors for each family
            to_plot_mask, colors_dict, default_colors = make_masks_and_colors(self.cat, which, filled_markers)
            if colors is not None:
                # Override default colors if custom colors provided
                colors_dict = {}
                for i, mask in enumerate(to_plot_mask):
                    for multiple_image in self.cat[mask]:
                        colors_dict[multiple_image['id']] = colors[i]
            else :
                colors = default_colors
            
            cat_contains_ellipse_params = len(np.unique(self.cat['a']))!=1
            count = 0
            for i, mask in enumerate(to_plot_mask) :
                for multiple_image in self.cat[mask] :
                    # Remove the *1000
                    if not cat_contains_ellipse_params :
                        a, b = size, size
                    else :
                        a, b = multiple_image['a'], multiple_image['b']
                    ellipse = self.plot_one_object(multiple_image['x'], multiple_image['y'], a, b,
                                                   multiple_image['theta'], count, color=colors[i])
                    self.qtItems[count] = ellipse
                    count += 1
                    
                    if mpl :
                        font = {'size':fontsize, 'family':'DejaVu Sans'}
                        plt.rc('font', **font)
                        self.plot_one_galaxy_mpl(multiple_image['x'], multiple_image['y'], a, b, multiple_image['theta'], color=colors[i][:3],
                                                 text=multiple_image['id'], linewidth=linewidth, text_color=text_color, text_alpha=text_alpha)
                        #self.plot_one_galaxy_mpl(multiple_image['x'], multiple_image['y'], a, b, multiple_image['theta'], color=colors[i][:3], text=multiple_image['id'])
            
            if make_thumbnails :
                if boost is not None :
                    adjusted_image = adjust_contrast(self.image_data, boost[0], pivot=boost[1])
                    adjusted_image = adjust_luminosity(adjusted_image, boost[2])
                else :
                    adjusted_image = self.image_data
                
                group_list = find_close_coord(self.cat, distance)
                for group in group_list :
                    
                    x_array = [self.cat[np.where(self.cat['id']==name)[0][0]]['x'] for name in group]
                    y_array = [self.cat[np.where(self.cat['id']==name)[0][0]]['y'] for name in group]
                    
                    x_pix = (np.max(x_array) + np.min(x_array)) / 2
                    y_pix = (np.max(y_array) + np.min(y_array)) / 2
                    
                    half_side = square_size // 2
                    
                    x_min = round( max( min( np.min(x_array) - margin, x_pix - half_side ), 0) )
                    x_max = round( min( max( np.max(x_array) + margin, x_pix + half_side ), self.image_data.shape[1]) )
                    y_min = round( max( min( np.min(y_array) - margin, y_pix - half_side ), 0) )
                    y_max = round( min( max( np.max(y_array) + margin, y_pix + half_side ), self.image_data.shape[0]) )
                    
                    if square_thumbnails :
                        x_side_size = x_max - x_min
                        y_side_size = y_max - y_min
                        if x_side_size!=y_side_size :
                            demi_taille_unique = round( max(x_side_size, y_side_size)/2 )
                            x_pix = round( (x_max + x_min)/2 )
                            y_pix = round( (y_max + y_min)/2 )
                            x_min = x_pix - demi_taille_unique
                            x_max = x_pix + demi_taille_unique
                            y_min = y_pix - demi_taille_unique
                            y_max = y_pix + demi_taille_unique
                    
                    plt_framework(image=True, figsize=3, drawscaler=1.2)
                    font = {'size':9, 'family':'DejaVu Sans'}
                    plt.rc('font', **font)
                    
                    
                    #cropped_image = self.image_data[y_min:y_max, x_min:x_max, :]
                    fig, ax = plot_image_mpl(adjusted_image, wcs=None, wcs_projection=False, units='pixel',
                                             pos=111, make_axes_labels=False, make_grid=False, crop=[x_min, x_max, y_min, y_max])
                    
                    for multiple_image_id in group :
                        multiple_image = self.cat[np.where(self.cat['id']==multiple_image_id)[0][0]]
                        color = colors_dict[multiple_image_id]
                        if not cat_contains_ellipse_params :
                            a, b = 75, 75
                        else :
                            a, b = multiple_image['a'], multiple_image['b']
                        self.plot_one_galaxy_mpl(multiple_image['x']-x_min, multiple_image['y']-y_min, a, b, multiple_image['theta'],
                                                 color=color[:3], text=multiple_image['id'], ax=ax, linewidth=linewidth, text_color=text_color, text_alpha=text_alpha)
                    
                    ax.axis('off')
                    #plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
                    
                    fig.show()
                        
                    plot_scale_bar(ax, deg_per_pix=self.pix_deg_scale, unit='arcsec',
                                   length=1 , color='white', linewidth=2, text_offset=0.01)
                    if savefig :
                        fig.savefig(os.path.join(os.path.dirname(self.image_path), 'mult_' + group[0]), bbox_inches='tight', pad_inches=0)
                        
                    plt_framework(full_tick_framework=True, ticks='out', image=True, width='full', drawscaler=0.8, tickscaler=0.5, minor_ticks=False)
        
        def plot_multiple_images_column(self, text_column, which='all'):
            if text_column not in self.cat.colnames:
                print(f"Column '{text_column}' not found in catalog")
                return
            if not hasattr(self, 'text_items'):
                self.text_items = []
            for text_item in self.text_items:
                self.qt_plot.removeItem(text_item)
            self.text_items.clear()
            
            to_plot_mask, colors_dict, default_colors = make_masks_and_colors(self.cat, which='all')
            colors = []
            for default_color in default_colors :
                colors.append(list(np.array(default_color)*255)[:3])
            print(colors)
            
            for i, mask in enumerate(to_plot_mask) :
                for multiple_image in self.cat[mask] :
                    
                    text = str(multiple_image[text_column])
                    text_item = pg.TextItem(text, color=colors[i])
                    
                    x = multiple_image['x']
                    y = self.image_data.shape[0] - multiple_image['y']  # Flip y to match PyQtGraph convention
                    semi_major = multiple_image['a']
                    semi_minor = multiple_image['b']
                    offset = max(semi_major, semi_minor)
                    text_item.setPos(x + offset/2, y - offset/2)
                    
                    font = PyQt5.QtGui.QFont()
                    font.setPointSize(15)
                    text_item.setFont(font)
                    
                    self.qt_plot.addItem(text_item)
                    self.text_items.append(text_item)
        
        
        self.multiple_images = self.make_catalog(multiple_images, unit_is_pixel=unit_is_pixel)
        self.multiple_images.plot = types.MethodType(plot_multiple_images, self.multiple_images)
        self.multiple_images.plot_column = types.MethodType(plot_multiple_images_column, self.multiple_images)
        return self.multiple_images.cat
    
    def load_potfile(self, potfile_path) :
        potfile_cat = read_potfile(potfile_path)
        self.potfile = self.make_catalog(cat=potfile_cat)
        return self.potfile.cat
    
    def open_cat(self, cat_path) :
        if '.fits' in cat_path :
            cat = Table.read(cat_path, format='fits')
        else :
            with open(cat_path, 'r') as raw_cat :
                first_line = raw_cat.readlines()[0]
            if len(first_line.split()) > len(first_line.split(',')) :
                cat_df = pd.read_csv(cat_path, delim_whitespace=True)[1:].apply(pd.to_numeric, errors='coerce')
            else :
                cat_df = pd.read_csv(cat_path)[1:].apply(pd.to_numeric, errors='coerce')
            cat = Table.from_pandas(cat_df)
        return cat
    
    def import_catalog(self, cat, color=[1., 1., 0], mag_colnames=['magAB_F814W', 'magAB_F435W']) :
        ref_path = None
        if isinstance(cat, str) :
            ref_path = cat
            cat = self.open_cat(cat)
        elif isinstance(cat, list) :
            ref_path = os.path.dirname(cat[0])
            run_match(cat[0], cat[1])
            for i in range(len(cat)-3) :
                run_match('matched_A_B.fits', cat[i+2])
            matched_cat = run_match('matched_A_B.fits', cat[-1])
            cat = Table(matched_cat[1].data)
            os.remove('matched_A_B.fits')
            
        self.imported_cat = self.make_catalog(cat=cat, make_selection_panel=True, color=color, mag_colnames=mag_colnames, ref_path=ref_path)
        return self.imported_cat.cat
    
    ################## Transform catalog into catalog class ###################
    
    def make_colnames_dict(self, catalog, use_default_names=True):
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
    
    def make_uniform_names_cat(self, cat, use_default_names=True, unit_is_pixel=False) :
        uniform_names_cat = cat.copy()
        colnames_dict = self.make_colnames_dict(cat, use_default_names=use_default_names)
        
        print('Column names to be used:')
        print(colnames_dict)
        
        for colname in colnames_dict.keys() :
            if colnames_dict[colname] is not None :
                
                if colname in uniform_names_cat.columns :
                    uniform_names_cat[colname] = uniform_names_cat[colnames_dict[colname]]
                else :
                    #uniform_names_cat.rename_column(colnames_dict[colname], colname)
                    uniform_names_cat.add_column( uniform_names_cat[colnames_dict[colname]], name=colname )
                    
                
        if colnames_dict['a'] != 'A_IMAGE' and colnames_dict['b'] != 'B_IMAGE' :
            yesno = 'y' if unit_is_pixel else input("ellipticity parameters " + str(colnames_dict['a']) \
                                                        + " and " + str(colnames_dict['b']) + " in pixels? [y][arcsec][deg]")
            if yesno == 'deg' :
                uniform_names_cat.replace_column( 'a', uniform_names_cat['a']/(self.pix_deg_scale) )
                uniform_names_cat.replace_column( 'b', uniform_names_cat['b']/(self.pix_deg_scale) )
            if yesno == 'arcsec' :
                uniform_names_cat.replace_column( 'a', uniform_names_cat['a']/(self.pix_deg_scale*3600) )
                uniform_names_cat.replace_column( 'b', uniform_names_cat['b']/(self.pix_deg_scale*3600) )
        
        #if colnames_dict['x']==None :
        x, y = self.world_to_image(uniform_names_cat['ra'], uniform_names_cat['dec'], unit='deg')
        uniform_names_cat['x'] = x
        uniform_names_cat['y'] = y
        
        yesno = 'y'
        if colnames_dict['a'] is not None and not use_default_names :
            yesno = input("'a', 'b' and 'theta' columns found in catalog. Use them as ellipticity parameters (if not, sources will be shown as circles)? [y] or [n]")
        if colnames_dict['a'] is None or yesno != 'y' :
            size = 40.
            uniform_names_cat['a'] = np.full(len(uniform_names_cat), size)
            uniform_names_cat['b'] = np.full(len(uniform_names_cat), size)
            uniform_names_cat['theta'] = np.full(len(uniform_names_cat), 0.)
        return uniform_names_cat
    
    def make_catalog(self, cat=None, cat_path=None, make_selection_panel=False, color=[1., 1., 0.], \
                     mag_colnames=['magAB_F814W', 'magAB_F435W'], ref_path=None, unit_is_pixel=False) :
        if cat_path is not None :
            cat = self.open_cat(cat_path)
        uniform_names_cat = self.make_uniform_names_cat(cat, unit_is_pixel=unit_is_pixel)
        if self.qt_plot is None :
            self.plot_image()
        return self.catalog(uniform_names_cat, self.image_data, self.qt_plot, window=self.window, make_selection_panel=make_selection_panel, \
                            image_path=self.image_path, image_widget = self.image_widget, image_widget_layout=self.image_widget_layout, \
                            color=color, mag_colnames=mag_colnames, ref_path=ref_path, mpl_fig=self.fig, mpl_ax=self.ax, pix_deg_scale=self.pix_deg_scale)
        
    ###########################################################################
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    class catalog :
        def __init__(self, cat, image_data, qt_plot, window=None, make_selection_panel=False, image_path=None, image_widget=None, image_widget_layout=None, \
                     color=[1., 1., 0.], mag_colnames=['magAB_F814W', 'magAB_F435W'], ref_path=None, mpl_fig=None, mpl_ax=None, pix_deg_scale=None) :
            self.cat = cat
            self.image_data = image_data
            self.qt_plot = qt_plot
            self.window = window
            self.qtItems = np.empty(len(cat), dtype=PyQt5.QtWidgets.QGraphicsEllipseItem)
            #self.qtItems = np.empty(len(cat), dtype=utils.utils_classes.selectable_ellipse.SelectableEllipse)
            self.color = color
            self.selection_mask = np.full(len(cat), False)
            self.selection_regions = []
            self.make_selection_panel = make_selection_panel
            self.image_path = image_path
            self.image_widget = image_widget
            self.image_widget_layout = image_widget_layout
            self.RS_widget = None
            self.x_axis_cleaned = np.full(len(cat), None)
            self.y_axis_cleaned = np.full(len(cat), None)
            self.mag_colnames = mag_colnames
            self.ref_path = ref_path
            self.mpl_ax = mpl_ax
            self.mpl_fig = mpl_fig
            self.pix_deg_scale = pix_deg_scale
        
        def make_mask_naninf(self, xy_axes=None) :
            #mag_F444W = flux_muJy_to_magAB(self.cat['f444w_tot_0'])
            #mag_F090W = flux_muJy_to_magAB(self.cat['f090w_tot_0'])
            #x_axis = mag_F444W
            #y_axis = mag_F090W - mag_F444W
            
            if xy_axes is None :
                mag_F814W = self.cat[self.mag_colnames[0]]
                mag_F435W = self.cat[self.mag_colnames[1]]
                x_axis = mag_F814W
                y_axis = mag_F435W - mag_F814W
            else :
                x_axis = self.cat[xy_axes[0]]
                y_axis = self.cat[xy_axes[1]]
            
            nan_mask = np.logical_not(np.isnan(x_axis)) & np.logical_not(np.isnan(y_axis))
            inf_mask = (x_axis!=np.inf) & (y_axis!=np.inf)
            #extremes_mask = (x_axis>0) & (x_axis<50)
            self.mask_naninf = nan_mask & inf_mask
            self.x_axis_cleaned = x_axis[self.mask_naninf]
            self.y_axis_cleaned = y_axis[self.mask_naninf]
            
            self.clear()
            self.cat = self.cat[self.mask_naninf]
            self.qtItems = np.empty(len(self.cat), dtype=PyQt5.QtWidgets.QGraphicsEllipseItem)
            self.selection_mask = np.full(len(self.cat), False)
        
        def plot(self, scale=1., color=None, text_column=None) :
            x = self.cat['x']
            y = self.cat['y']
            semi_major = self.cat['a'] * scale
            semi_minor = self.cat['b'] * scale
            angle = self.cat['theta']
            for i in tqdm(range(len(semi_major))) :
                ellipse = self.plot_one_object(x[i], y[i], semi_major[i], semi_minor[i], angle[i], i, color=color)
                self.qtItems[i] = ellipse
            
            # Add text labels if requested
            if text_column is not None:
                self.plot_column(text_column, color=color)
        
        def clear(self) :
            # Clear ellipses
            for i in tqdm( range(len(self.qtItems)) ) :
                self.qt_plot.removeItem(self.qtItems[i])
            
            # Clear text items if they exist
            if hasattr(self, 'text_items'):
                self.clear_column()
        
        def clear_column(self) :
            for text_item in self.text_items:
                self.qt_plot.removeItem(text_item)
            self.text_items.clear()
        
        def clear_selection(self) :
            self.selection_mask[np.full(len(self.cat), True)] = False
            self.selection_regions.clear()
            self.clear()
            self.plot()
        
        def plot_one_object(self, x, y, semi_major, semi_minor, angle, idx, color=None) :
            if color is None :
                color = list(np.array(self.color)*255)
            else :
                color = list(np.array(color)*255)
            #make the flip to accomosate pyqtgraph's strange plotting conventions
            y = self.image_data.shape[0] - y
            angle = -angle
            #####################################################################
            ellipse = SelectableEllipse(x-semi_major/2, y-semi_minor/2, semi_major, semi_minor, idx, self.selection_mask, \
                                        self.qtItems, color, scatter_pos=(self.x_axis_cleaned[idx], self.y_axis_cleaned[idx]), RS_widget=self.RS_widget)
            ellipse.setTransformOriginPoint( PyQt5.QtCore.QPointF(x, y) )
            ellipse.setRotation(angle)
            self.qt_plot.addItem(ellipse)
            return ellipse
        
        def plot_selection_panel(self, xy_axes=None) :
            if self.make_selection_panel :
                self.make_mask_naninf(xy_axes=xy_axes)
                
                self.RS_widget, self.selection_ROI = plot_panel(self.x_axis_cleaned, self.y_axis_cleaned, self.image_widget_layout, self.qt_plot)
                
                data=(self.x_axis_cleaned, self.y_axis_cleaned)
                #self.selection_mask = np.full(len(self.mag_F444W_cleaned), False)
                self.selectable_scatter = SelectableScatter(self.RS_widget, self.selection_ROI, data, self.selection_mask, \
                                                            qtItems=self.qtItems, color=list(np.array(self.color)*255))
                
        def make_image_ROI(self) :
            center_y = self.image_data.shape[0]/2
            center_x = self.image_data.shape[1]/2
            self.image_ROI = ellipse_maker_ROI([center_x-200, center_y-100], [400, 200], self.qt_plot, self.window, self.cat)
            make_handles(self.image_ROI)
            self.qt_plot.addItem(self.image_ROI)
            
        def make_cleaner_ROI(self) :
            self.image_widget.cat = self
            self.select_sources = SelectSources(self.cat, self.qt_plot, self.image_widget.current_ROI, self.selection_mask, self.selection_regions, \
                                                window=self.window, qtItems=self.qtItems, color=list(np.array(self.color)*255))
            
        def save_selection_mask(self, path=None) :
            self.selection_mask_path = self.make_path(path, self.ref_path, 'selection_mask.npy')
            np.save(self.selection_mask_path, self.selection_mask)
            
        def load_selection_mask(self, path=None) :
            self.selection_mask_path = self.make_path(path, self.ref_path, 'selection_mask.npy')
            self.selection_mask = np.load(self.selection_mask_path)
            
        def save_selection_regions(self, path=None) :
            self.selection_regions_path = self.make_path(path, self.image_path, 'selection_regions.npy')
            np.save(self.selection_regions_path, self.selection_regions)
            
        def load_selection_regions(self, path=None, name='selection_regions.npy') :
            self.selection_regions_path = self.make_path(path, self.image_path, name)
            self.selection_regions = np.load(self.selection_regions_path).tolist()
            
            size_y = self.qt_plot.image.shape[0]
            for rect_params in self.selection_regions :
                indiv_mask = InRectangle(self.cat['x'], size_y - self.cat['y'], rect_params)
                self.selection_mask[indiv_mask] = True
            
            fig, ax = plt.subplots()
            ax.axis('equal')
            size = max(self.qt_plot.image.shape[0], self.qt_plot.image.shape[1])
            #ax.invert_yaxis()
            ax.set_ylim([size+2000, -2000])
            ax.set_xlim([-4000, size+4000])
            
            for i, rect_params in enumerate(self.selection_regions) :
                x0, y0, a, b, angle = rect_params
                #x1, y1 = x0, y0
                #x2, y2 = x0 + a*np.cos(angle), y0 + a*np.sin(angle)
                #x3, y3 = x0 + a*np.cos(angle) - b*np.sin(angle),  y0 + a*np.sin(angle) + b*np.cos(angle)
                #x4, y4 = x0 - b*np.sin(angle), y0 + b*np.cos(angle)
                #ax.plot([x1, x2, x3, x4, x1], size_y-np.array([y1, y2, y3, y4, y1]), c='b')
                #ax.axis('equal')
                ax.add_patch( Rectangle((x0, y0), a, b, angle=angle*180/np.pi, alpha=0.4 ))
                ax.text(x0, y0, str(i))
                fig.show()
                plt.pause(0.05)
                
        def make_path(self, path, ref_path, name) :
            if path is None :#and self.ref_path is not None :
                #to_return = os.path.join(os.path.dirname(ref_path), name)
                to_return = os.path.join(os.path.dirname(ref_path), os.path.basename(ref_path).split('.')[0] + '_' + name)
            elif os.path.isdir(path) :
                to_return = os.path.join(path, name)
            elif os.path.isdir(os.path.dirname(path)) :
                to_return = path
            return to_return
            
            
        def plot_one_galaxy_mpl(self, x, y, a, b, theta, color=[1,1,1], text=None, ax=None, linewidth=1., text_color='white', text_alpha=0.5) :
            edgecolor = list(color).copy()
            edgecolor.append(1)
            facecolor = edgecolor.copy()
            facecolor[-1] = 0
            ellipse = Ellipse( (x, y), a, b, angle=theta, facecolor=facecolor, edgecolor=edgecolor, lw=linewidth )
            if ax is None :
                self.mpl_ax.add_artist(ellipse)
                if text is not None :
                    #self.mpl_ax.text(x-1.5*b*np.abs(np.sin(theta)), y-1.5*b*np.abs(np.cos(theta)), text, color=edgecolor[:3], \
                    #                 ha='right', va='top')
                    offset = 0.85
                    theta_modulo = theta%180 * np.pi/180
                    if theta_modulo<np.pi/2 :
                        x_text, y_text = x+offset*b*np.abs(np.sin(theta_modulo)), y-offset*b*np.abs(np.cos(theta_modulo))
                        horizontalalignment, verticalalignment = 'left', 'top'
                    else :
                        x_text, y_text = x-offset*b*np.abs(np.sin(theta_modulo)), y-offset*b*np.abs(np.cos(theta_modulo))
                        horizontalalignment, verticalalignment = 'right', 'top'
                    self.mpl_ax.text( x_text, y_text, text, c=text_color, alpha=1, fontsize=15, \
                                      ha=horizontalalignment, va=verticalalignment, \
                                      bbox=dict(facecolor=edgecolor[:3], alpha=text_alpha, edgecolor='none') )
            else :
                ax.add_artist(ellipse)
                if text is not None :
                    offset = 0.85
                    theta_modulo = theta%180 * np.pi/180
                    if theta_modulo<np.pi/2 :
                        x_text, y_text = x+offset*b*np.abs(np.sin(theta_modulo)), y-offset*b*np.abs(np.cos(theta_modulo))
                        horizontalalignment, verticalalignment = 'left', 'top'
                    else :
                        x_text, y_text = x-offset*b*np.abs(np.sin(theta_modulo)), y-offset*b*np.abs(np.cos(theta_modulo))
                        horizontalalignment, verticalalignment = 'right', 'top'
                    #ax.text(x_text, y_text, text, color=edgecolor[:3], \
                    #        ha=horizontalalignment, va=verticalalignment)
                    ax.text( x_text, y_text, text, c=text_color, alpha=1, fontsize=15, \
                             ha=horizontalalignment, va=verticalalignment, \
                             bbox=dict(facecolor=edgecolor[:3], alpha=text_alpha, edgecolor='none') )
        
        def export_to_mult_file(self, file_path=None) :
            if file_path is None :
                file_path = os.path.join(os.path.dirname(self.image_path), 'mult.lenstool')
            
            sub_cat = self.cat[self.selection_mask]
            
            header = "#REFERENCE 0\n## id   RA      Dec        a         b         theta     z         mag\n"
            with open(file_path, 'w') as f :
                f.write(header)
                for index, row in enumerate(sub_cat) :
                    if 'THETA_WORLD' in sub_cat.colnames :
                        line = (f"{row['id']:<3}  {row['ra']:10.6f}  {row['dec']:10.6f}  "
                                f"{row['a']:8.6f}  {row['b']:8.6f}  {row['THETA_WORLD']:8.6f}  "
                                f"{row['zb']:8.6f}  {row['f814w_mag']:8.6f}\n")
                    else :
                        line = (f"{row['id']:<3}  {row['ra']:10.6f}  {row['dec']:10.6f}  "
                                "0.0  0.0  0.0  0.0  0.0\n")
                    f.write(line)
        
        def plot_column(self, text_column, color=None):
            """
            Add text labels from a specified column to existing plotted ellipses.
            
            Parameters:
            -----------
            text_column : str
                Name of the column in the catalog to use for labels
            color : list or None
                RGB color for the text. If None, uses the same color as the ellipses
            """
            if text_column not in self.cat.colnames:
                print(f"Column '{text_column}' not found in catalog")
                return
        
            if color is None:
                color = list(np.array(self.color[:3])*255)
            else:
                color = list(np.array(color[:3])*255)
        
            # Store text items to prevent garbage collection
            if not hasattr(self, 'text_items'):
                self.text_items = []
        
            # Clear existing text items if any
            for text_item in self.text_items:
                self.qt_plot.removeItem(text_item)
            self.text_items.clear()
        
            # Add new text labels
            for i in range(len(self.cat)):
                text = str(self.cat[text_column][i])
                text_item = pg.TextItem(text, color=color)
        
                # Get ellipse position and size for offset calculation
                x = self.cat['x'][i]
                y = self.image_data.shape[0] - self.cat['y'][i]  # Flip y to match PyQtGraph convention
                semi_major = self.cat['a'][i]
                semi_minor = self.cat['b'][i]
        
                # Position text slightly offset from the ellipse
                offset = max(semi_major, semi_minor)
                text_item.setPos(x + offset/2, y - offset/2)
        
                # Set font
                font = PyQt5.QtGui.QFont()
                font.setPointSize(15)
                text_item.setFont(font)

                self.qt_plot.addItem(text_item)
                self.text_items.append(text_item)

    
    
    
    
    def world_to_image(self, ra, dec, unit='deg') :
        coord = SkyCoord(ra, dec, unit=unit)
        image_coord = WCS.world_to_pixel(self.wcs, coord)
        if len(image_coord[0].shape)==0 :
            image_coord = (image_coord[0]*1., image_coord[1]*1.)
        return image_coord
    
    def clear_Items(self) :
        for key in self.qtItems_dict.keys() :
            if self.qtItems_dict[key] is not None :
                for i in tqdm( range(len(self.qtItems_dict[key])) ) :
                    self.qt_plot.removeItem(self.qtItems_dict[key][i])
    
    
    
    
    
    def import_lenstool(self, model_dir) :
        self.lt_dir = model_dir
        
        ### Add correct optimized/fixed redshifts ###
        if self.multiple_images is not None :
            source_z_dict = make_source_z_dict(model_dir, use_family_name_only=True)
            for i in range(len(self.multiple_images.cat)) :
                self.multiple_images.cat['z'][i] = source_z_dict[ self.multiple_images.cat['id'][i][0] ]
            
            magnification_dir = os.path.join(model_dir, 'magnifications')
            if os.path.isdir(magnification_dir) :
                if 'magnification' not in self.multiple_images.cat.colnames :
                    
                    magnification_files = os.listdir(magnification_dir)
                    magnification_dict = {}
                    for name in source_z_dict.keys() :
                        idx = np.where( np.array([ magnification_file[:-len('.fits')-1][-1] \
                                                   for magnification_file in magnification_files ]) == name)[0][0]
                        
                        magnification_path = os.path.join( magnification_dir, magnification_files[idx])
                        with fits.open(magnification_path) as magnification_fits :
                            magnification_dict[name] = magnification_fits[0].data, WCS(magnification_fits[0].header)
                    
                    magnification_column = np.zeros(len(self.multiple_images.cat))
                    for i in range(len(self.multiple_images.cat)) :
                        family_name = self.multiple_images.cat['id'][i][0]
                        world_coord = SkyCoord( self.multiple_images.cat[i]['ra'], self.multiple_images.cat[i]['dec'], unit='deg' )
                        pix_coord_array = WCS.world_to_pixel(magnification_dict[family_name][1], world_coord)
                        pix_coord = ( round(pix_coord_array[0]*1.), round(pix_coord_array[1]*1.) )
                        magnification = magnification_dict[family_name][0][pix_coord[1], pix_coord[0]]
                        magnification_column[i] = magnification
                        
                    self.multiple_images.cat.add_column(magnification_column, name='magnification')
                
        if self.redshift is None :
            self.redshift = float(input("Please specify the redshift of the lens: "))
        self.lt_bayes = read_bayes_file( os.path.join(model_dir, 'bayes.dat'), convert_to_kpc=True, z=self.redshift )
        self.lt_latex_table_str = make_param_latex_table(model_dir, convert_to_kpc=True, z=self.redshift)
        def print_lt_latex_table() :
            print(self.lt_latex_table_str)
        self.lt_latex_table = print_lt_latex_table
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    def make_hand_made_catalogue(self) :
        empty_cat_dict = {'x': [], 'y': [], 'a': [], 'b': [], 'theta': []}
        empty_cat = Table(empty_cat_dict)
        self.hand_made_catalogue = self.catalog(empty_cat, self.image_data, self.qt_plot, window=self.window, make_selection_panel=False, image_widget_layout=self.image_widget_layout)
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    """
    def plot_sources(self, color=[1., 1., 0.], alpha=0.5, scale=1.) :
        if self.sources is None :
            self.extract_sources()
        if self.qt_plot is None :
            self.plot_image()
        x = self.sources['X_IMAGE']
        y = self.sources['Y_IMAGE']
        semi_major = self.sources['A_IMAGE'] * scale
        semi_minor = self.sources['B_IMAGE'] * scale
        angle = self.sources['THETA_IMAGE']
        
        qtItems_local = np.empty(len(x), dtype=PyQt5.QtWidgets.QGraphicsEllipseItem)
        for i in tqdm(range(len(semi_major))) :
            ellipse = self.plot_one_source(x[i], y[i], semi_major[i], semi_minor[i], angle[i], color=color)
            qtItems_local[i] = ellipse
        self.qtItems_dict['sources'] = qtItems_local
    
    def plot_one_source(self, x, y, semi_major, semi_minor, angle, color=[1., 1., 0]) :
        color = list(np.array(color)*255)
        #make the flip to accomosate pyqtgraph's strange plotting conventions
        y = self.image_data.shape[0] - y
        angle = -angle
        #####################################################################
        ellipse = PyQt5.QtWidgets.QGraphicsEllipseItem(x-semi_major/2, y-semi_minor/2, semi_major, semi_minor)
        ellipse.setTransformOriginPoint( PyQt5.QtCore.QPointF(x, y) ) 
        #ellipse.setTransform( PyQt5.QtGui.QTransform().rotate(angle[i]) )
        ellipse.setRotation(angle)
        ellipse.setPen( pg.mkPen(color + [255]) )
        ellipse.setBrush( pg.mkBrush(color + [127]) )
        self.qt_plot.addItem(ellipse)
        return ellipse
    
    def plot_catalog(self, cat_path=None, cat=None, color=[1., 0., 1.], scale=1.) :
        if self.imported_cat is None :
            self.import_catalog(self, cat_path=cat_path, cat=cat)
        uniform_names_cat = self.make_uniform_names_cat(self.imported_cat.cat)
        if self.qt_plot is None :
            self.plot_image()
        x = uniform_names_cat['x']
        y = uniform_names_cat['y']
        semi_major = uniform_names_cat['a'] * scale
        semi_minor = uniform_names_cat['b'] * scale
        angle = uniform_names_cat['theta']
        
        qtItems_local = np.empty(len(x), dtype=PyQt5.QtWidgets.QGraphicsEllipseItem)
        for i in tqdm(range(len(x))) :
            ellipse = self.plot_one_source(x[i], y[i], semi_major[i], semi_minor[i], angle[i], color=color)
            qtItems_local[i] = ellipse
        self.qtItems_dict['imported_cat'] = qtItems_local
    
    def plot_multiple_images(self, which='all', size=80, alpha=0.7, marker='o', filled_markers=False, plotter='qt') :
        if plotter=='qt' and self.qt_plot is None :
            self.plot_image()
        if plotter=='mpl' and self.ax is None :
            self.plot_image_mpl()
        if self.multiple_images is None :
            print('No multiple images to be plotted.')
        else :
            if which=='all' :
                families = np.unique( [self.multiple_images[i]['n'][0] for i in range(len(self.multiple_images))] )
                which = families
            to_plot_mask = []
            for symbol in which :
                to_plot_mask.append( [self.multiple_images[i]['n'].startswith(symbol) for i in range(len(self.multiple_images))] )
            
            colors = make_palette(len(which), 1, alpha=1)
            if filled_markers :
                facecolors = make_palette(len(which), 1, alpha=0.3)
            else :
                facecolor = [[0,0,0,0] for i in range(len(which))]
            
            if plotter=='qt' :
                qtItems_local = []
                for i, mask in enumerate(to_plot_mask) :
                    for multiple_image in self.multiple_images[mask] :
                        x, y = self.world_to_image( multiple_image['x'], multiple_image['y'] )
                        # Remove the *1000
                        ellipse = self.plot_one_source(x, y, multiple_image['a']*1000, multiple_image['b']*1000, multiple_image['theta'], color=colors[i][:3])
                        qtItems_local.append(ellipse)
                qtItems_dict['multiple_images'] = np.array(qtItems_local)
            if plotter=='mpl' :
                for i, mask in enumerate(to_plot_mask) :
                    x, y = self.world_to_image( self.multiple_images['x'][mask], self.multiple_images['y'][mask] )
                    self.ax.scatter( x, y, sizes=[size], marker=marker, color=colors[i], facecolor=facecolor[i], label=which[i] )
                    self.ax.legend()
    """
    
    
    
    """
    def plot_sources_mpl(self, color='blue', alpha=0.5, scale=1., facecolor=None) :
        if self.sources is None :
            self.extract_sources()
        if self.ax is None :
            self.plot_image_mpl()
        semi_major = self.sources['A_IMAGE'] * scale
        semi_minor = self.sources['B_IMAGE'] * scale
        angle = self.sources['THETA_IMAGE']
        
        for i in tqdm(range(len(semi_major))):
            if facecolor is None :
                ellipse = Ellipse( (self.sources['X_IMAGE'][i], self.sources['Y_IMAGE'][i]), 2*semi_major[i], 2*semi_minor[i], angle=angle[i], \
                                   alpha=alpha, color=color )
            else :
                edgecolor = facecolor.copy()
                edgecolor[-1] = min(1, facecolor[-1]+0.1)
                ellipse = Ellipse( (self.sources['X_IMAGE'][i], self.sources['Y_IMAGE'][i]), 2*semi_major[i], 2*semi_minor[i], angle=angle[i], \
                                   facecolor=facecolor, edgecolor=edgecolor )
            self.ax.add_artist(ellipse)
        plt.show()
        return self.fig, self.ax
    """
    
    
    
    def plot_sub_region(self, ra, dec, size=3):
        """
        Plots a square region around given RA and Dec coordinates.

        Parameters:
        ra (float): Right Ascension of the center in degrees.
        dec (float): Declination of the center in degrees.
        size (float): Size of the square region in arcseconds (default is 10).
        """
        
        x_center, y_center = self.world_to_image(ra, dec, unit='deg')
                
        size_pix = int( size / (self.pix_deg_scale*3600) / 2 )
        
        fig, axs = plt.subplots(1,3)
        for i in range(len(x_center)) :
            x_min = int(x_center[i]) - size_pix
            x_max = int(x_center[i]) + size_pix
            y_min = int(y_center[i]) - size_pix
            y_max = int(y_center[i]) + size_pix
            region = self.image_data[y_min:y_max, x_min:x_max, :]
            axs[i].imshow(region, origin='lower')
            axs[i].axis('off')
            
        return fig, axs
        
        
        
        
        




