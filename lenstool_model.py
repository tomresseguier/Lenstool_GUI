import os
import glob
import shutil
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Polygon, Circle, Rectangle
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import astropy.constants as c
from reproject import reproject_interp
from collections import defaultdict
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
import pyqtgraph.exporters
from PyQt5.QtWidgets import QGraphicsEllipseItem, QGraphicsScene, QGraphicsView, QGraphicsSceneMouseEvent, QApplication
from PyQt5.QtCore import QRectF
#from PyQt5.QtCore import Qt
from pyqtgraph.Qt import QtCore
from astropy.table import Table
import pandas as pd
import lenstool

###############################################################################
import sys
#module_dir = os.path.dirname(os.path.abspath(__file__))
#sys.path.append(module_dir)

from .source_extraction.source_extract import source_extract, source_extract_DIM
from .source_extraction.match_cat import run_match
#sys.path.append(module_dir + "/utils")
from .utils.utils_astro.cat_manip import match_cat2
from .utils.utils_plots.plot_utils_general import *
from .utils.utils_Qt.selectable_classes import *
from .utils.utils_Qt.utils_general import *
from .utils.utils_general.utils_general import find_close_coord, make_colnames_dict
from .utils.utils_plots.plt_framework import plt_framework
from .utils.utils_Lenstool.redshift_extractors import make_source_z_dict, find_param_file
from .utils.utils_Lenstool.param_extractors import read_potfile, read_bayes_file, make_param_latex_table

from .utils.utils_Lenstool.file_makers import best_files_maker, make_magnifications_and_curves                  # This import is problematic. The two functions run Lenstool
                                                                                                                # and are therefore dependent on my own install.
from .utils.utils_Lenstool.operations import make_image_to_source
from .utils.utils_astro.utils_general import relative_to_world
#from utils_general.utils import flux_muJy_to_magAB
from .utils.utils_astro.set_cosmology import set_cosmo
cosmo = set_cosmo()






def make_which_colors(self, filled_markers=False, saturation=None) :
    which = self.families if self.which=='all' else self.which
    
    if filled_markers:
        colors = make_palette(len(which), 1, alpha=0.5, sat_fixed=saturation)
    else:
        colors = make_palette(len(which), 1, alpha=0, sat_fixed=saturation)
    
    which_colors_dict = {}
    for i, name in enumerate(which) :
        which_colors_dict[name] = colors[i]
    #for i, mask in enumerate(to_plot_mask):
    #    for multiple_image in cat[mask]:
    #        which_colors_dict[multiple_image['id']] = colors[i]
    return which_colors_dict


def make_full_color_function(families) :
    n_families = len(families)
    def make_full_color_dict(filled_markers=False, saturation=None) :
        if filled_markers:
            colors = make_palette(n_families, 1, alpha=0.5, sat_fixed=saturation)
        else:
            colors = make_palette(n_families, 1, alpha=0, sat_fixed=saturation)
        full_colors_dict = {}
        for i, family in enumerate(families) :
            full_colors_dict[family] = colors[i]
        return full_colors_dict
    return make_full_color_dict


def import_multiple_images(self, mult_file_path, fits_image, units=None, AttrName='mult', filled_markers=False, saturation=None) :
    multiple_images = Table(names=['id','family','ra','dec','a','b','theta','z','mag'], dtype=['str','str',*['float',]*7])
    with open(mult_file_path, 'r') as mult_file:
        for line in mult_file:
            cleaned_line = line.strip()
            if not cleaned_line.startswith("#") and len(cleaned_line)>0 :
                split_line = cleaned_line.split()
                row = [split_line[0], '---'] #split_line[0][:-1]
                for element in split_line[1:8] :
                    row.append(float(element))
                multiple_images.add_row(row)
    multiple_images['theta'] = multiple_images['theta'] - fits_image.orientation
    multiple_images['family'] = find_families(multiple_images['id'])
    
    setattr(self, AttrName, fits_image.make_catalog(multiple_images, units=units))
    
    local_families = np.unique( find_families(getattr(self, AttrName).cat['id']) ).tolist()
    #local_families = np.unique([ im_id[:-1] for im_id in getattr(self, AttrName).cat['id'] ]).tolist()
    self.families = np.unique(self.families + local_families).tolist()
    self.which = self.families.copy()
    
    self.mult_colors = make_full_color_function(self.families)
    
    #if AttrName=='mult' :
        #self.families = np.unique( find_families(getattr(self, AttrName).cat['id']) )
        #self.mult_colors = make_full_color_function(self.families)
        #self.which = self.families.tolist()
    
    def make_to_plot_masks() :
        to_plot_masks = {}
        #for i, name in enumerate(self.which) :
        #    other_names_mask = np.full(len(self.which), True)
        #    other_names_mask[i] = False
        #    other_names = np.array(self.which)[other_names_mask]
        #    ambiguous_names = []
        #    for other_name in other_names :
        #        if other_name.startswith(name) :
        #           ambiguous_names.append(other_name)
        #    to_plot_mask = np.full(len(getattr(self, AttrName).cat), False)
        #    for j, im_id in enumerate(getattr(self, AttrName).cat['id']) :
        #        if im_id.startswith(name) and True not in [ im_id.startswith(ambiguous_name) for ambiguous_name in ambiguous_names ] :
        #            to_plot_mask[j] = True
        #    to_plot_masks[name] = to_plot_mask
        for family in self.which :
            to_plot_masks[family] = getattr(self, AttrName).cat['family'] == family
        return to_plot_masks
    def make_overall_mask() :
        overall_mask = np.full(len(getattr(self, AttrName).cat), False)
        for mask in make_to_plot_masks().values() :
            overall_mask = np.logical_or(overall_mask, mask)
        return overall_mask
    getattr(self, AttrName).masks = make_to_plot_masks
    getattr(self, AttrName).mask = make_overall_mask
    
    lenstool_model = self
    
    def plot_multiple_images(self, size=1, marker='o', filled_markers=filled_markers, colors=None, mpl=False, fontsize=9,
                             make_thumbnails=False, square_size=150, margin=50, distance=200, savefig=False, square_thumbnails=True,
                             boost=[2,1.5,1], linewidth=1.7, text_color='white', text_alpha=0.5, saturation=saturation) :
        self.clear()
        saturation = fits_image.lt.saturation if saturation is None else saturation
        self.saturation = saturation
        
        if colors is not None :
            colors_dict = {}
            for i, family in enumerate(lenstool_model.which) :
                colors_dict[family] = colors[i]
        else :
            colors_dict = lenstool_model.mult_colors(filled_markers=filled_markers, saturation=saturation)
        
        
        cat_contains_ellipse_params = len(np.unique(self.cat['a']))!=1
        count = 0
        for name, mask in self.masks().items() :
            for multiple_image in self.cat[mask] :
                # Remove the *1000
                if not cat_contains_ellipse_params :
                    a, b = size*40, size*40
                else :
                    a, b = multiple_image['a']*size, multiple_image['b']*size
                ellipse = self.plot_one_object(multiple_image['x'], multiple_image['y'], a, b,
                                               multiple_image['theta'], count, color=colors_dict[name], linewidth=linewidth)
                self.qtItems[count] = ellipse
                count += 1
                
                if mpl :
                    font = {'size':fontsize, 'family':'DejaVu Sans'}
                    plt.rc('font', **font)
                    self.plot_one_galaxy_mpl(multiple_image['x'], multiple_image['y'], a, b, multiple_image['theta'], color=colors_dict[name][:3],
                                             text=multiple_image['id'], linewidth=linewidth, text_color=text_color, text_alpha=text_alpha)
                    #self.plot_one_galaxy_mpl(multiple_image['x'], multiple_image['y'], a, b, multiple_image['theta'], color=colors_dict[name][:3], text=multiple_image['id'])
        
        
        if make_thumbnails :
            if boost is not None :
                adjusted_image = adjust_contrast(self.fits_image.image_data, boost[0], pivot=boost[1])
                adjusted_image = adjust_luminosity(adjusted_image, boost[2])
            else :
                adjusted_image = self.fits_image.image_data
            
            if group_images :
                group_list = find_close_coord(self.cat[self.mask()], distance)
            else :
                group_list = [[name] for name in self.cat[self.mask()]['id']]
            
            for group in group_list :
                
                x_array = [self.cat[np.where(self.cat['id']==name)[0][0]]['x'] for name in group]
                y_array = [self.cat[np.where(self.cat['id']==name)[0][0]]['y'] for name in group]
                
                x_pix = (np.max(x_array) + np.min(x_array)) / 2
                y_pix = (np.max(y_array) + np.min(y_array)) / 2
                
                half_side = square_size // 2
                
                x_min = round( max( min( np.min(x_array) - margin, x_pix - half_side ), 0) )
                x_max = round( min( max( np.max(x_array) + margin, x_pix + half_side ), self.fits_image.image_data.shape[1]) )
                y_min = round( max( min( np.min(y_array) - margin, y_pix - half_side ), 0) )
                y_max = round( min( max( np.max(y_array) + margin, y_pix + half_side ), self.fits_image.image_data.shape[0]) )
                
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
                
                
                #cropped_image = self.fits_image.image_data[y_min:y_max, x_min:x_max, :]
                fig, ax = plot_image_mpl(adjusted_image, wcs=None, wcs_projection=False, units='pixel',
                                         pos=111, make_axes_labels=False, make_grid=False, crop=[x_min, x_max, y_min, y_max])
                
                for multiple_image_id in group :
                    multiple_image = self.cat[np.where(self.cat['id']==multiple_image_id)[0][0]]
                    color = colors_dict[multiple_image_id[:-1]]
                    if not cat_contains_ellipse_params :
                        a, b = 75, 75
                    else :
                        a, b = multiple_image['a'], multiple_image['b']
                    self.plot_one_galaxy_mpl(multiple_image['x']-x_min, multiple_image['y']-y_min, a, b, multiple_image['theta'],
                                             color=color[:3], text=multiple_image['id'], ax=ax, linewidth=linewidth, text_color=text_color, text_alpha=text_alpha)
                
                ax.axis('off')
                #plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
                
                fig.show()
                    
                plot_scale_bar(ax, deg_per_pix=self.fits_image.pix_deg_scale, unit='arcsec',
                               length=1 , color='white', linewidth=2, text_offset=0.01)
                if savefig :
                    fig.savefig(os.path.join(os.path.dirname(self.fits_image.image_path), 'mult_' + group[0]), bbox_inches='tight', pad_inches=0)
                    
                plt_framework(full_tick_framework=True, ticks='out', image=True, width='full', drawscaler=0.8, tickscaler=0.5, minor_ticks=False)
    
    
    def plot_multiple_images_column(self, text_column, which='all'):
        if text_column not in self.cat.colnames:
            print(f"Column '{text_column}' not found in catalog")
            return
        if not hasattr(self, 'text_items'):
            self.text_items = []
        for text_item in self.text_items:
            self.fits_image.qt_plot.removeItem(text_item)
        self.text_items.clear()
        
        colors_dict = lenstool_model.mult_colors(filled_markers=False, saturation=self.saturation)
        
        for name, mask in self.masks().items() :
            for multiple_image in self.cat[mask] :
                
                text = str(multiple_image[text_column])
                text_item = pg.TextItem( text, color=list( np.array(colors_dict[name])*255 )[:3] )
                
                x = multiple_image['x']
                y = self.fits_image.image_data.shape[0] - multiple_image['y']  # Flip y to match PyQtGraph convention
                semi_major = multiple_image['a']
                semi_minor = multiple_image['b']
                offset = max(semi_major, semi_minor)
                text_item.setPos(x + offset/2, y - offset/2)
                
                font = PyQt5.QtGui.QFont()
                font.setPointSize(15)
                text_item.setFont(font)
                
                self.fits_image.qt_plot.addItem(text_item)
                self.text_items.append(text_item)
    
    
    getattr(self, AttrName).plot = types.MethodType(plot_multiple_images, getattr(self, AttrName))
    getattr(self, AttrName).plot_column = types.MethodType(plot_multiple_images_column, getattr(self, AttrName))
    
    def transfer_ids(self, id_name='id') :
        if fits_image.imported_cat is not None :
            if id_name in fits_image.imported_cat.cat.colnames :
                temp_cat = match_cat2([self.cat, fits_image.imported_cat.cat], keep_all_col=True, fill_in_value=-1, column_to_transfer=id_name)
                if id_name in self.cat.colnames :
                    id_name = id_name + '_CAT2'
                self.cat[id_name] = temp_cat[id_name]
                print('###############\nColumn ' + id_name + ' added.\n###############')
            else :
                print(id_name + ' not found in imported_cat')
        else :
            print('No imported_cat')
    
    getattr(self, AttrName).transfer_ids = types.MethodType(transfer_ids, getattr(self, AttrName))
    
    


def export_thumbnails(self, group_images=True, square_thumbnails=True, square_size=150, margin=50, distance=200, export_dir=None, boost=True, make_broad_view=True, broad_view_params=None) :
    export_dir = os.path.join(os.path.dirname(self.fits_image.image_path), 'mult_thumbnails') if export_dir is None else os.path.abspath(os.path.join(export_dir, 'mult_thumbnails'))
    if os.path.isdir( os.path.dirname( os.path.dirname(export_dir) ) ) and not os.path.isdir( os.path.dirname(export_dir) ) :
        os.mkdir(os.path.dirname(export_dir))
    if not os.path.isdir(export_dir) :
        os.mkdir(export_dir)
    
    if not self.fits_image.boosted and boost :
        self.fits_image.boost()
    
    if group_images :
        group_list = find_close_coord(self.cat[self.mask()], distance)
    else :
        group_list = [[name] for name in self.cat[self.mask()]['id']]
    
    for group in group_list :
        
        x_array = [self.cat[np.where(self.cat['id']==name)[0][0]]['x'] for name in group]
        y_array = [self.cat[np.where(self.cat['id']==name)[0][0]]['y'] for name in group]
        
        x_pix = (np.max(x_array) + np.min(x_array)) / 2
        y_pix = (np.max(y_array) + np.min(y_array)) / 2
        
        half_side = square_size // 2
        
        x_min = round( max( min( np.min(x_array) - margin, x_pix - half_side ), 0) )
        x_max = round( min( max( np.max(x_array) + margin, x_pix + half_side ), self.fits_image.image_data.shape[1]) )
        y_min = round( max( min( np.min(y_array) - margin, y_pix - half_side ), 0) )
        y_max = round( min( max( np.max(y_array) + margin, y_pix + half_side ), self.fits_image.image_data.shape[0]) )
        
        #if square_thumbnails :
        x_side_size = x_max - x_min
        y_side_size = y_max - y_min
        demi_taille_unique = round( max(x_side_size, y_side_size)/2 )
        if x_side_size!=y_side_size :
            x_pix = round( (x_max + x_min)/2 )
            y_pix = round( (y_max + y_min)/2 )
            x_min = x_pix - demi_taille_unique
            x_max = x_pix + demi_taille_unique
            y_min = y_pix - demi_taille_unique
            y_max = y_pix + demi_taille_unique
        
        zoom_rect = QRectF(x_min, self.fits_image.image_data.shape[0] - y_max, demi_taille_unique*2, demi_taille_unique*2)
        self.fits_image.qt_plot.getView().setRange(zoom_rect)        
        
        
        thumbnail_path = os.path.join( export_dir, 'mult_' + group[0] + '.png' )
        print('Creating ' + thumbnail_path)
        exporter = pg.exporters.ImageExporter(self.fits_image.qt_plot.view)
        exporter.export(thumbnail_path)
        print('Done')
        
    
    ##### Adding broad view #####
    if make_broad_view :
        # broad_view_params = [ [x_min, x_max], [y_min, y_max] ]
        if broad_view_params is not None :
            x = broad_view_params[0][0]
            y = self.fits_image.image_data.shape[0] - broad_view_params[1][1]
            x_width = broad_view_params[0][1]-broad_view_params[0][0]
            y_width = broad_view_params[1][1]-broad_view_params[1][0]
            zoom_rect = QRectF(x, y, x_width, y_width)
        else :
            zoom_rect = QRectF(0, 0, self.fits_image.image_data.shape[1], self.fits_image.image_data.shape[0])
        self.fits_image.qt_plot.getView().setRange(zoom_rect)        
        
        broadview_filename = 'broadview'
        for name in self.fits_image.lt.which :
            broadview_filename += '_' + name
        broadview_path = os.path.join( export_dir, broadview_filename + '.png' )
        print('Creating ' + broadview_path)
        exporter = pg.exporters.ImageExporter(self.fits_image.qt_plot.view)
        exporter.export(broadview_path)
        print('Done')
    ##############################

        



class curves :
    def __init__(self, curves_dir, lenstool_model, fits_image, which_critcaus='critical', join=False, size=2) :
        self.dir = curves_dir
        self.paths = glob.glob(os.path.join(curves_dir, "*.dat"))
        self.lenstool_model = lenstool_model
        self.fits_image = fits_image
        self.size = size
        
        self.qtItems = {}
        for name in lenstool_model.families :
            self.qtItems[name] = None
        
        self.coords = {}
        for name in lenstool_model.families :
            curve_mask = np.array([name in os.path.basename(path) for path in self.paths])
            if True in curve_mask :
                lines = []
                for curve_path in np.array(self.paths)[curve_mask] :
                    file = open(curve_path, 'r')
                    all_lines = file.readlines()
                    lines += all_lines[1:]
                
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
                         
                if which_critcaus=='caustic' :
                    delta_ra = np.array( [ float( lines[i].split()[3] ) for i in range(len(lines)) ] )
                    delta_dec = np.array( [ float( lines[i].split()[4] ) for i in range(len(lines)) ] )
                    ra, dec = relative_to_world(delta_ra, delta_dec, reference=(ra_ref, dec_ref))
                    x, y = fits_image.world_to_image(ra, dec)
                    
                    y = fits_image.image_data.shape[0] - y
                    
                    shorten_indices = np.linspace(0, len(x) - 1, 10000, dtype=int)
                    x = x[shorten_indices]
                    y = y[shorten_indices]
                    
                    #if join :
                    #    x, y = rearrange_points(x, y)
                
                self.coords[name] = (x, y)
        
    def plot(self) :
        self.clear()
        for name in self.lenstool_model.which :
            color = np.round(self.lenstool_model.mult_colors(saturation=self.lenstool_model.saturation)[name]*255).astype(int)
            color[3] = 255
            
            x, y = self.coords[name]
            
            scatter = pg.ScatterPlotItem(x, y, pen=None, brush=pg.mkBrush(color), size=self.size)
            self.qtItems[name] = scatter
            self.fits_image.qt_plot.addItem(scatter)
            
    def clear(self) :
        for name, qtItem in self.qtItems.items() :
            if qtItem is not None :
                self.fits_image.qt_plot.removeItem(qtItem)
                self.qtItems[name] = None
                
    
    
    


class lenstool_model :
    def __init__(self, model_path, fits_image) :
        self.safe_mode = False
        self.fits_image = fits_image
        self.saturation = 1
        self.model_dir = model_path if os.path.isdir(model_path) else os.path.dirname(model_path)
        all_par_file_paths = glob.glob(os.path.join(self.model_dir, "*.par"))
        #all_cat_file_paths = glob.glob(os.path.join(self.model_dir, "*.lenstool"))
        
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
        
        if self.potfile_path is not None :
            potfile_Table = read_potfile(self.potfile_path)
            self.potfile = fits_image.make_catalog(potfile_Table, color=[1.,0.,0.], units='arcsec')
        else :
            self.potfile = None
        
        self.families = []
        self.which = []
        
        mult_path_list = glob.glob( os.path.join(self.model_dir, "*mult*.lenstool") )
        if len(mult_path_list)==1 :
            mult_file_path = mult_path_list[0]
            import_multiple_images(self, mult_file_path, fits_image, units='pixel', filled_markers=True)
        else :
            self.mult = None
        
        arclets_path_list = glob.glob( os.path.join(self.model_dir, "*arclet*.lenstool") )
        if len(arclets_path_list)==1 :
            arclets_path = arclets_path_list[0]
            import_multiple_images(self, arclets_path, fits_image, AttrName='arclets', units='pixel', filled_markers=False)
        else :
            self.arclets = None
        
        if len(arclets_path_list)==1 and len(mult_path_list)==1 :
            import_multiple_images(self, mult_file_path, fits_image, units='pixel', filled_markers=True)
            import_multiple_images(self, arclets_path, fits_image, AttrName='arclets', units='pixel', filled_markers=False)
        
        #dot_all_paths = glob.glob(os.path.join(self.model_dir, "*.all"))
        predicted_images_path = os.path.join(self.model_dir, 'image.dat')
        if os.path.isfile(predicted_images_path) :
            import_multiple_images(self, predicted_images_path, fits_image, AttrName='image', units='pixel', filled_markers=False)
        else :
            self.image = None
        
        curves_dir = os.path.join(self.model_dir, 'curves')
        if os.path.isdir(curves_dir) :
            self.curves = curves(curves_dir, self, fits_image, which_critcaus='critical', join=False, size=2)
        else :
            self.curves = None
        
        #self.reference = [179.4888967, -10.7669233]
        self.reference = [47.2332780, 26.7604953]
    
    
    def SafeMode(self) :
        if self.safe_mode :
            print("Already in safe directory")
        elif self.model_dir.endswith('_safe/') :
            print("Safe directory already selected, moving to it")
            os.chdir(self.model_dir)
            self.safe_mode = True
        else :
            safe_dir = self.model_dir.rstrip('/') + '_safe/'
            if os.path.exists(safe_dir) :
                print("Safe directory already exists, moving to it")
                self.__init__(safe_dir, self.fits_image)
                os.chdir(safe_dir)
                self.safe_mode = True
            else :
                os.makedirs(safe_dir, exist_ok=False)
                for item in os.listdir(self.model_dir):
                    s = os.path.join(self.model_dir, item)
                    d = os.path.join(safe_dir, item)
                    if os.path.isdir(s):
                        shutil.copytree(s, d, dirs_exist_ok=True)
                    else:
                        shutil.copy2(s, d)
                self.__init__(safe_dir, self.fits_image)
                os.chdir(safe_dir)
                self.safe_mode = True
        print("Now in " + os.getcwd())
    
    
    def select_multiple_images(self) :
        return 'in progress'
    
    def plot(self, which=None) :
        if which is not None :
            self.set_which(which)
        if self.mult is not None :
            self.mult.plot()
            self.mult.plot_column('id')
        if self.image is not None :
            self.image.plot()
            self.image.plot_column('id')
        if self.curves is not None :
            self.curves.plot()
    
    def clear(self) :
        if self.mult is not None :
            self.mult.clear()
        if self.arclets is not None :
            self.arclets.clear()
        if self.image is not None :
            self.image.clear()
        if self.curves is not None :
            self.curves.clear()
            
    def set_which(self, *names) :
        if names[0]=='all' :
            self.which = self.families.tolist()
        elif isinstance(names[0], list) :
            self.which = names[0]
        else :
            self.which = list(names)
        print("Images to plot are now ", self.which)
        #self.clear()
        #self.plot()
        
    def make_files(self) :
        best_files_maker(self.model_dir)
        make_magnifications_and_curves(self.model_dir)
    #def relative_to_mosaic_pixel
    
    
    def export_thumbnails(self, group_images=True, square_thumbnails=True, square_size=150, margin=50, distance=200, export_dir=None, boost=True, make_broad_view=True, broad_view_params=None) :
        export_thumbnails(self.mult, group_images=group_images, square_thumbnails=square_thumbnails, square_size=square_size, margin=margin, \
                          distance=distance, export_dir=export_dir, boost=boost, make_broad_view=make_broad_view, broad_view_params=broad_view_params)
        
    
    def world_to_relative(self, ra, dec) :
        ref = SkyCoord(self.reference[0], self.reference[1], unit='deg')
        world_radec = SkyCoord(ra, dec, unit='deg')
        relative_coord = ( (world_radec.ra - ref.ra)*np.cos(ref.dec.rad), world_radec.dec - ref.dec )
        return -relative_coord[0].arcsec, relative_coord[1].arcsec
    
    def relative_to_world(self, xr, yr) :
        ref = SkyCoord(self.reference[0], self.reference[1], unit='deg')
        dec = ref.dec.deg + yr*u.arcsec.to('deg')
        ra = ref.ra.deg - xr*u.arcsec.to('deg') / np.cos(dec*u.deg.to('rad'))
        return ra, dec
    
    def make_webpage(self) :
        print('in progress')
        
    
    
    def set_lt_z(self, z) :
        if not self.safe_mode :
            print("Moving to safe directory")
            self.copysafe()
        
        self.lt_z = z
        print(self.best_file_path)
        print(os.getcwd())
        self.lt = lenstool.Lenstool( os.path.basename(self.best_file_path) )
        
        self.dx_map = None
    
    
    def start_im2source(self) :
        if self.dx_map is None :
            self.dx_map, self.dy_map, self.dmap_wcs = self.lt.g_dpl(1000, self.lt_z)
        
        transform_coords_radec = make_image_to_source(self.dx_map, self.dy_map, self.dmap_wcs)
        
        def transform_coords(x, y) :
            ra, dec = self.fits_image.image_to_world(x, y)
            
            ra_source, dec_source = transform_coords_radec(ra, dec)
            x_source, y_source = self.fits_image.world_to_image(ra_source, dec_source)
            return x_source, y_source, ra_source, dec_source
            
        self.transform_coords = transform_coords
        
        # Set up label if needed (optional)
        if not hasattr(self, 'transform_label'):
            self.transform_label = pg.TextItem(anchor=(0, 1), color='w')
            self.fits_image.qt_plot.addItem(self.transform_label)
        
        # Create scatter point for transformed location
        self.transformed_point = pg.ScatterPlotItem(size=10, brush='r')
        self.fits_image.qt_plot.addItem(self.transformed_point)
        
        self.images_scatter = pg.ScatterPlotItem(size=10, brush='g')
        self.fits_image.qt_plot.addItem(self.images_scatter)
    
        # Ensure view does not auto-range when updating
        self.fits_image.qt_plot.getView().enableAutoRange(pg.ViewBox.XAxis, False)
        self.fits_image.qt_plot.getView().enableAutoRange(pg.ViewBox.YAxis, False)
    
        # Mouse tracking setup
        def mouse_moved(evt):
            pos = evt[0]
            if self.fits_image.qt_plot.getView().sceneBoundingRect().contains(pos):
                mouse_point = self.fits_image.qt_plot.getView().mapSceneToView(pos)
                x, y_flipped = mouse_point.x(), mouse_point.y()
                x, y = x, self.fits_image.image_data.shape[0] - y_flipped
                try:
                    x_source, y_source, ra_source, dec_source = self.transform_coords(x, y)
                    
                    #xr_source, yr_source = self.world_to_relative(ra_source, dec_source)
                    source = Table( rows=[('test', ra_source, dec_source, 1, 1, 0, self.lt_z, 25)],names=['n','x','y','a','b','theta','z','mag'], dtype=['str',*['float',]*7] )
                    self.lt.set_sources(source)
                    self.lt.e_lensing()
                    image_cat = self.lt.get_images()
                    x_images = []
                    y_images = []
                    for image in image_cat :
                        ra_image, dec_image = self.relative_to_world(image['x'], image['y'])
                        x_image, y_image = self.fits_image.world_to_image(ra_image, dec_image)
                        x_images.append(x_image)
                        y_images.append(y_image)
                    self.images_scatter.setData( x_images, self.fits_image.image_data.shape[0] - np.array(y_images) )
                    
                    self.transformed_point.setData([x_source], [self.fits_image.image_data.shape[0] - y_source])
                    self.transform_label.setText(f"({x:.2f}, {y:.2f}) → ({x_source:.2f}, {y_source:.2f})")
                    self.transform_label.setPos(x, y_flipped)
                except Exception as e:
                    self.transform_label.setText(f"Error: {e}")
                    self.transform_label.setPos(x, y_flipped)
    
        self._transform_proxy = pg.SignalProxy(self.fits_image.qt_plot.getView().scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)
    
        # Key press handling
        def keyPressEvent(event):
            if event.key() == Qt.Key_Escape:
                self.fits_image.qt_plot.removeItem(self.transformed_point)
                self.fits_image.qt_plot.removeItem(self.transform_label)
                self._transform_proxy.disconnect()
                del self._transform_proxy
                del self.transformed_point
                del self.transform_label
                del self.images_scatter
                self.fits_image.window.keyPressEvent = self._original_keyPressEvent
    
        # Save existing keyPressEvent so we can restore it
        self._original_keyPressEvent = self.fits_image.window.keyPressEvent
        self.fits_image.window.keyPressEvent = keyPressEvent
    
    
    
    
    




def find_families(image_ids):
    # Step 1: Initial guess by chopping last character
    id_to_family = {img_id: img_id[:-1] for img_id in image_ids}
    
    # Step 2: Group by these tentative families
    family_groups = defaultdict(list)
    for img_id, fam in id_to_family.items():
        family_groups[fam].append(img_id)

    # Step 3: Merge singleton families if their name starts with another family name
    updated = True
    while updated:
        updated = False
        singletons = {fam for fam, ids in family_groups.items() if len(ids) == 1}
        for fam in list(singletons):
            for target in family_groups:
                if fam != target and fam.startswith(target):
                    family_groups[target].extend(family_groups[fam])
                    del family_groups[fam]
                    updated = True
                    break
            if updated:
                break

    # Step 4: Merge families with 'alt' in original IDs if the ID starts with another family name
    for fam in list(family_groups):
        for img_id in family_groups[fam]:
            if 'alt' in img_id:
                for target in family_groups:
                    if fam != target and img_id.startswith(target):
                        family_groups[target].extend(family_groups[fam])
                        del family_groups[fam]
                        break
                break  # Only need to check one 'alt' image to trigger a merge

    # Step 5: Build final output mapping
    final_map = {}
    for fam, ids in family_groups.items():
        for img_id in ids:
            final_map[img_id] = fam

    return [final_map[img_id] for img_id in image_ids]

    








