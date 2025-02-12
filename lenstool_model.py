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
from .utils.utils_astro.utils_general import relative_to_world
#from utils_general.utils import flux_muJy_to_magAB
from .utils.utils_astro.set_cosmology import set_cosmo
cosmo = set_cosmo()









def import_multiple_images(self, mult_file_path, fits_image, units=None, AttrName='mult', filled_markers=False):
    
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
    
    multiple_images['theta'] = multiple_images['theta'] - fits_image.orientation
    
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
    
    
    setattr(self, AttrName, fits_image.make_catalog(multiple_images, units=units))
    
    to_plot_mask, colors_dict, default_colors = make_masks_and_colors(getattr(self, AttrName).cat, which='all', filled_markers=filled_markers)
    colors_dict_families = {}
    for name in colors_dict.keys() :
        family = name[:-1]
        colors_dict_families[family] = colors_dict[name]
    
    getattr(self, AttrName).plot = types.MethodType(plot_multiple_images, getattr(self, AttrName))
    getattr(self, AttrName).plot_column = types.MethodType(plot_multiple_images_column, getattr(self, AttrName))
    getattr(self, AttrName).color = colors_dict_families
    #return getattr(self, AttrName).cat



def plot_curves(self, fits_image, which_critcaus='critical', join=False, size=2) :
    curves_dir = os.path.join(self.model_dir, 'curves')
    curve_paths = glob.glob(os.path.join(curves_dir, "*.dat"))
    curve_dict = {}
    for name in self.families :
        i = np.where( [name in os.path.basename(path) for path in curve_paths] )[0][0]
        curve_path = curve_paths[i]
        curve_dict[name] = curve_path
        
        file = open(curve_path, 'r')
        all_lines = file.readlines()
        lines = all_lines[1:]
        
        ra_ref = float( all_lines[0].split()[-2] )
        dec_ref = float( all_lines[0].split()[-1] )
        
        if which_critcaus=='critical' :
            delta_ra = np.array( [ float( lines[i].split()[1] ) for i in range(len(lines)) ] )
            delta_dec = np.array( [ float( lines[i].split()[2] ) for i in range(len(lines)) ] )
            ra, dec = relative_to_world(delta_ra, delta_dec, reference=(ra_ref, dec_ref))
            x, y = fits_image.world_to_image(ra, dec)
            
            y = fits_image.image_data.shape[0] - y
            
            shorten_indices = np.linspace(0, len(x) - 1, 10000, dtype=int)
            x = x[shorten_indices]
            y = y[shorten_indices]
            
            #if join :
            #    x, y = rearrange_points(x, y)
            color = np.round(self.mult.color[name]*255).astype(int)
            color[3] = 255
            scatter = pg.ScatterPlotItem(x, y, pen=None, brush=pg.mkBrush(color), size=size)
            fits_image.qt_plot.addItem(scatter)
                 
        if which_critcaus=='caustic' :
            delta_ra = np.array( [ float( lines[i].split()[3] ) for i in range(len(lines)) ] )
            delta_dec = np.array( [ float( lines[i].split()[4] ) for i in range(len(lines)) ] )
            ra, dec = relative_to_world(delta_ra, delta_dec, reference=(ra_ref, dec_ref))
            x, y = fits_image.world_to_image(ra, dec)
            
            #if join :
            #    x, y = rearrange_points(x, y)
            scatter = pg.ScatterPlotItem(x, y, pen=None, brush=pg.mkBrush(255, 0, 0, 255), size=5)
            fits_image.qt_plot.addItem(scatter)
        
    return curve_dict, scatter
    
    
    
    
    
    
    




class lenstool_model :
    def __init__(self, model_path, fits_image) :
        self.fits_image = fits_image
        self.model_dir = model_path if os.path.isdir(model_path) else os.path.dirname(model_path)
        all_par_file_paths = glob.glob(os.path.join(self.model_dir, "*.par"))
        all_cat_file_paths = glob.glob(os.path.join(self.model_dir, "*.lenstool"))
        
        if os.path.isfile(model_path) :
            self.param_file_path = model_path
        else :
            for file_path in all_par_file_paths :
                if not os.path.basename(file_path).startswith('best') :
                    with open(file_path, 'r') as file :
                        for line in file :
                            stripped_line = line.strip()
                            if stripped_line and not stripped_line.startswith('#') :  # Skip empty and comment lines
                                if stripped_line.startswith('runmode') :
                                    self.param_file_path = file_path
        
        all_par_file_names = [ os.path.basename(file_path) for file_path in all_par_file_paths ]
        
        self.has_run = 'best.par' in all_par_file_names
        self.best_file_path = os.path.join(self.model_dir, 'best.par') if 'best.par' in all_par_file_names else None
        self.bestopt_file_path = os.path.join(self.model_dir, 'bestopt.par') if 'bestopt.par' in all_par_file_names else None
        potfile_paths_list = glob.glob(os.path.join(self.model_dir, "*potfile*.lenstool"))
        self.potfile_path = potfile_paths_list[0] if len(potfile_paths_list)>=1 else None
        
        potfile_Table = read_potfile(self.potfile_path)
        self.potfile = fits_image.make_catalog(potfile_Table, color=[1.,0.,0.], units='arcsec')
        
        mult_idx = np.where(['mult' in name for name in all_cat_file_paths])[0]
        if len(mult_idx)==1 :
            mult_file_path = all_cat_file_paths[mult_idx[0]]
            import_multiple_images(self, mult_file_path, fits_image, units='pixel')
            self.families = np.unique( [s[:-1] for s in self.mult.cat['id']] )
        else :
            self.mult = None
            self.families = None
        
        #dot_all_paths = glob.glob(os.path.join(self.model_dir, "*.all"))
        arclets_path = os.path.join(self.model_dir, 'image.dat')
        if os.path.isfile(arclets_path) :
            import_multiple_images(self, arclets_path, fits_image, AttrName='image', units='pixel')
        
        
            
        
        
    def select_multiple_images(self) :
        return 'in progress'
    
    def plot_curves(self) :
        plot_curves(self, self.fits_image)
    
    
    #def relative_to_mosaic_pixel
        
        
        
        
        
        
        