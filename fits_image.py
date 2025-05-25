import os
import glob
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

from .catalog import catalog
from .lenstool_model import lenstool_model
from .source_extraction.source_extract import source_extract, source_extract_DIM
from .source_extraction.match_cat import run_match
#sys.path.append(module_dir + "/utils")
from .utils.utils_plots.plot_utils_general import *
from .utils.utils_Qt.selectable_classes import *
from .utils.utils_Qt.utils_general import *
from .utils.utils_general.utils_general import find_close_coord, make_colnames_dict
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
        self.image_data, self.pix_deg_scale, self.orientation, self.wcs, self.header = self.open_image(self.image_path)
        ### Following lines useless, just to see the possible attributes of the class ###
        self.sources = None
        self.fig = None
        self.ax = None
        self.multiple_images = None
        self.galaxy_selection = None
        self.imported_cat = None
        self.qt_plot = None
        #self.qt_plot = self.plot_image()
        self.qtItems_dict = {'sources': None,
                             'potfile_cat': None,
                             'imported_cat': None,
                             'multiple_images': None}
        self.ax = None
        self.redshift = None
        self.boosted_image_path = self.image_path[:-5] + '_boosted.fits'
        if os.path.isfile(self.boosted_image_path) :
            print("Boosted image found: " + self.boosted_image_path)
            self.boosted_image, _, _, _, _ = self.open_image(self.boosted_image_path)
        else :
            self.boosted_image = None
        self.boosted = False
        self.extra_qt_plots = []
        self.extra_windows = []
        
        self.plot_image()
    
    def open_image(self, image_path) :
        with fits.open(image_path) as hdus :
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
                print("keyword is", keyword)
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
                    if len(selected_hdus)==0 and len(wht_hdus)==0 :
                        selected_hdus = true_image_hdus
                
                
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
                elif 'CD1_1' in header and 'CD2_2' in header :
                    if 'CD1_2' in header and 'CD2_1' in header :
                        cd = np.array([[header['CD1_1'], header['CD1_2']], [header['CD2_1'], header['CD2_2']]])
                        #det = np.linalg.det(cd)
                        #sign = np.sign(det)
                        orientation = np.arctan2(cd[1,0], cd[1,1])
                    else :
                        orientation = 0.0
                elif 'PC1_1' in header and 'PC2_2' in header :
                    if 'PC1_2' in header and 'PC2_1' in header :
                        cd = np.array([[header['PC1_1'], header['PC1_2']], [header['PC2_1'], header['PC2_2']]])
                        #det = np.linalg.det(cd)
                        #sign = np.sign(det)
                        orientation = np.arctan2(cd[1,0], cd[1,1])
                    else :
                        orientation = 0.0                    
                else :
                    orientation = None
                orientation = np.rad2deg(orientation) if orientation is not None else None
                
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
                
                return image, pix_deg_scale, orientation, wcs, header
            
            else :
                print('Unable to extract image data from FITS file')
                return None
    
    def create_qt_plot(self) :
        to_plot = np.flip(self.image_data, axis=0) if not self.boosted else np.flip(self.boosted_image, axis=0)
        
        #qt_plot = pg.image(to_plot)
        qt_plot = pg.ImageView()
        qt_plot.setImage(to_plot)
        #qt_plot.autoLevels()
        
        image_widget_layout = QHBoxLayout()
        image_widget_layout.addWidget(qt_plot)
        
        #self.image_widget = QWidget()
        image_widget = DragWidget(qt_plot)
        image_widget.setLayout(image_widget_layout)
        
        window = QMainWindow()
        window.setWindowTitle(os.path.basename(self.image_path))
        window.setCentralWidget(image_widget)
        window.show()
        return qt_plot, image_widget_layout, image_widget, window
    
    def plot_image(self) :
        #to_plot = self.image_data
        #to_plot = np.transpose(self.image_data, axes=[1,0,2])
        if self.qt_plot is None or not self.qt_plot.isVisible() :
            print('Creating main window...')
            self.qt_plot, self.image_widget_layout, self.image_widget, self.window = self.create_qt_plot()
            print('Done')
            to_return = self.qt_plot
        else :
            print('Creating secondary window...')
            extra_qt_plot, _, _, extra_window = self.create_qt_plot()
            extra_window.setWindowTitle(os.path.basename(self.image_path) + ' (' + str(len(self.extra_qt_plots)+2) + ')')
            extra_window.show()
            self.extra_qt_plots.append(extra_qt_plot)
            self.extra_windows.append(extra_window)
            print('Done')
            to_return = extra_qt_plot
        return to_return
    
    def boost(self, boost=[2,1.5,1]) :
        if self.boosted_image is None :
            print('Adjusting contrast...')
            adjusted_image = adjust_contrast(self.image_data, boost[0], pivot=boost[1])
            print('Adjusting luminosity...')
            self.boosted_image = adjust_luminosity(adjusted_image, boost[2])
            print('Writing to memory...')
            hdul = fits.HDUList([fits.PrimaryHDU()] + [ fits.ImageHDU(data=self.boosted_image[:,:,i], header=self.wcs.to_header(), name=name) for i, name in enumerate(['RED','GREEN','BLUE']) ])
            hdul.writeto(self.boosted_image_path, overwrite=True)
        if not self.boosted :
            print('Plotting...')
            self.qt_plot.setImage(np.flip(self.boosted_image, axis=0))
            #self.qt_plot.autoLevels()
            self.boosted = True
            for extra_qt_plot in self.extra_qt_plots :
                extra_qt_plot.setImage(np.flip(self.boosted_image, axis=0))
            print('Done')
    
    def unboost(self) :
        if self.boosted :
            print('Plotting...')
            self.qt_plot.setImage(np.flip(self.image_data, axis=0))
            #self.qt_plot.autoLevels()
            self.boosted = False
            for extra_qt_plot in self.extra_qt_plots :
                extra_qt_plot.setImage(np.flip(self.image_data, axis=0))
            print('Done')
    
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
    
    
    def make_catalog(self, cat, color=[1., 1., 0], mag_colnames=['magAB_F814W', 'magAB_F435W'], units=None) :
        if self.qt_plot is None :
            self.plot_image()
        #to_return = catalog(cat, self.image_data, self.wcs, self.qt_plot, window=self.window, image_path=self.image_path, 
        #                            image_widget = self.image_widget, image_widget_layout=self.image_widget_layout, color=color, 
        #                            mag_colnames=mag_colnames, mpl_fig=self.fig, mpl_ax=self.ax, 
        #                            pix_deg_scale=self.pix_deg_scale, units=units)
        to_return = catalog(cat, self, color=color, mag_colnames=mag_colnames, units=units)
        return to_return
    
    def import_catalog(self, cat, color=[1., 1., 0], mag_colnames=['magAB_F814W', 'magAB_F435W'], units='pixel') :
        self.imported_cat = self.make_catalog(cat, color=color, mag_colnames=mag_colnames, units=units)
        #print(self.imported_cat)
        
    ###########################################################################
    
    
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
        #self.lt_dir = model_dir
        self.lt = lenstool_model(model_dir, self)
        """
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
    """
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
        
        
        
        
        




