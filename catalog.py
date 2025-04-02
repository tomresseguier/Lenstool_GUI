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
#from utils_general.utils import flux_muJy_to_magAB
from .utils.utils_astro.set_cosmology import set_cosmo
cosmo = set_cosmo()






def open_cat(cat_path) :
    if '.fits' in cat_path :
        cat = Table.read(cat_path, format='fits')
    else :
        with open(cat_path, 'r') as raw_cat :
            first_line, second_line = raw_cat.readlines()[0:2]
            start_line = 1 if second_line.startswith('--') else 0
        if len(first_line.split()) > len(first_line.split(',')) :
            cat_df = pd.read_csv(cat_path, delim_whitespace=True)[start_line:].apply(pd.to_numeric, errors='coerce')
        else :
            cat_df = pd.read_csv(cat_path)[start_line:].apply(pd.to_numeric, errors='coerce')
        cat = Table.from_pandas(cat_df)
    return cat

    
def make_uniform_names_cat(cat, self) :
    uniform_names_cat = cat.copy()
    colnames_dict = make_colnames_dict(cat, use_default_names=self.use_default_names)
    
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
        if self.units is None :
            self.units = input("ellipticity parameters " + str(colnames_dict['a']) \
                               + " and " + str(colnames_dict['b']) + " in pixels? [y][arcsec][deg]")
        if self.units == 'deg' :
            uniform_names_cat.replace_column( 'a', uniform_names_cat['a']/(self.fits_image.pix_deg_scale) )
            uniform_names_cat.replace_column( 'b', uniform_names_cat['b']/(self.fits_image.pix_deg_scale) )
        if self.units == 'arcsec' :
            uniform_names_cat.replace_column( 'a', uniform_names_cat['a']/(self.fits_image.pix_deg_scale*3600) )
            uniform_names_cat.replace_column( 'b', uniform_names_cat['b']/(self.fits_image.pix_deg_scale*3600) )
    
    #if colnames_dict['x']==None :
    x, y = self.fits_image.world_to_image(uniform_names_cat['ra'], uniform_names_cat['dec'], unit='deg')
    uniform_names_cat['x'] = x
    uniform_names_cat['y'] = y
    
    yesno = 'y'
    if colnames_dict['a'] is not None and not self.use_default_names :
        yesno = input("'a', 'b' and 'theta' columns found in catalog. Use them as ellipticity parameters (if not, sources will be shown as circles)? [y] or [n]")
    if colnames_dict['a'] is None or yesno != 'y' :
        size = 40.
        uniform_names_cat['a'] = np.full(len(uniform_names_cat), size)
        uniform_names_cat['b'] = np.full(len(uniform_names_cat), size)
        uniform_names_cat['theta'] = np.full(len(uniform_names_cat), 0.)
    return uniform_names_cat


def initialize_catalog(cat, self) :
    ##### Combine or just open the catalog #####
    ref_path = None
    if isinstance(cat, str) :
        ref_path = cat
        cat = open_cat(cat)
    elif isinstance(cat, list) :
        ref_path = os.path.dirname(cat[0])
        run_match(cat[0], cat[1])
        for i in range(len(cat)-3) :
            run_match('matched_A_B.fits', cat[i+2])
        matched_cat = run_match('matched_A_B.fits', cat[-1])
        cat = Table(matched_cat[1].data)
        os.remove('matched_A_B.fits')
    
    ##### Standardize the column names #####
    uniform_names_cat = make_uniform_names_cat(cat, self)
    return uniform_names_cat, ref_path
    







class catalog :
    def __init__(self, cat, fits_image, color=[1., 1., 0.], mag_colnames=['magAB_F814W', 'magAB_F435W'], use_default_names=True, units=None) :
        self.fits_image = fits_image
        
        self.mag_colnames = mag_colnames
        self.use_default_names = use_default_names
        self.units = units
        
        self.cat, self.ref_path = initialize_catalog(cat, self)
        self.qtItems = np.empty(len(self.cat), dtype=PyQt5.QtWidgets.QGraphicsEllipseItem)
        #self.qtItems = np.empty(len(self.cat), dtype=utils.utils_classes.selectable_ellipse.SelectableEllipse)
        self.color = color
        self.selection_mask = np.full(len(self.cat), False)
        self.selection_regions = []
        self.RS_widget = None
        self.x_axis_cleaned = np.full(len(self.cat), None)
        self.y_axis_cleaned = np.full(len(self.cat), None)
    
    
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
        self.clear()
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
            if self.qtItems[i] is not None :
                self.fits_image.qt_plot.removeItem(self.qtItems[i])
                self.qtItems[i] = None
        
        # Clear text items if they exist
        if hasattr(self, 'text_items'):
            self.clear_column()
    
    def clear_column(self) :
        for text_item in self.text_items:
            self.fits_image.qt_plot.removeItem(text_item)
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
        y = self.fits_image.image_data.shape[0] - y
        angle = -angle
        #####################################################################
        ellipse = SelectableEllipse(x-semi_major/2, y-semi_minor/2, semi_major, semi_minor, idx, self.selection_mask, \
                                    self.qtItems, color, scatter_pos=(self.x_axis_cleaned[idx], self.y_axis_cleaned[idx]), RS_widget=self.RS_widget)
        ellipse.setTransformOriginPoint( PyQt5.QtCore.QPointF(x, y) )
        ellipse.setRotation(angle)
        self.fits_image.qt_plot.addItem(ellipse)
        return ellipse
    
    def plot_selection_panel(self, xy_axes=None) :
        self.make_mask_naninf(xy_axes=xy_axes)
        
        self.RS_widget, self.selection_ROI = plot_panel(self.x_axis_cleaned, self.y_axis_cleaned, self.fits_image.image_widget_layout, self.fits_image.qt_plot)
        
        data=(self.x_axis_cleaned, self.y_axis_cleaned)
        #self.selection_mask = np.full(len(self.mag_F444W_cleaned), False)
        self.selectable_scatter = SelectableScatter(self.RS_widget, self.selection_ROI, data, self.selection_mask, \
                                                    qtItems=self.qtItems, color=list(np.array(self.color)*255))
            
    def make_image_ROI(self) :
        center_y = self.fits_image.image_data.shape[0]/2
        center_x = self.fits_image.image_data.shape[1]/2
        self.image_ROI = ellipse_maker_ROI([center_x-200, center_y-100], [400, 200], self.fits_image.qt_plot, self.fits_image.window, self.cat)
        make_handles(self.image_ROI)
        self.fits_image.qt_plot.addItem(self.image_ROI)
        
    def make_cleaner_ROI(self) :
        self.fits_image.image_widget.cat = self
        self.select_sources = SelectSources(self.cat, self.fits_image.qt_plot, self.fits_image.image_widget.current_ROI, self.selection_mask, self.selection_regions, \
                                            window=self.fits_image.window, qtItems=self.qtItems, color=list(np.array(self.color)*255))
        
    def save_selection_mask(self, path=None) :
        self.selection_mask_path = self.make_path(path, self.ref_path, 'selection_mask.npy')
        np.save(self.selection_mask_path, self.selection_mask)
        print("Selection mask saved at " + self.selection_mask_path)
        
    def load_selection_mask(self, path=None) :
        self.selection_mask_path = self.make_path(path, self.ref_path, 'selection_mask.npy')
        self.selection_mask = np.load(self.selection_mask_path)
        
    def save_selection_regions(self, path=None) :
        self.selection_regions_path = self.make_path(path, self.fits_image.image_path, 'selection_regions.npy')
        np.save(self.selection_regions_path, self.selection_regions)
        
    def load_selection_regions(self, path=None, name='selection_regions.npy') :
        self.selection_regions_path = self.make_path(path, self.fits_image.image_path, name)
        self.selection_regions = np.load(self.selection_regions_path).tolist()
        
        size_y = self.fits_image.qt_plot.image.shape[0]
        for rect_params in self.selection_regions :
            indiv_mask = InRectangle(self.cat['x'], size_y - self.cat['y'], rect_params)
            self.selection_mask[indiv_mask] = True
        
        fig, ax = plt.subplots()
        ax.axis('equal')
        size = max(self.fits_image.qt_plot.image.shape[0], self.fits_image.qt_plot.image.shape[1])
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
            self.fits_image.mpl_ax.add_artist(ellipse)
            if text is not None :
                #self.fits_image.mpl_ax.text(x-1.5*b*np.abs(np.sin(theta)), y-1.5*b*np.abs(np.cos(theta)), text, color=edgecolor[:3], \
                #                 ha='right', va='top')
                offset = 0.85
                theta_modulo = theta%180 * np.pi/180
                if theta_modulo<np.pi/2 :
                    x_text, y_text = x+offset*b*np.abs(np.sin(theta_modulo)), y-offset*b*np.abs(np.cos(theta_modulo))
                    horizontalalignment, verticalalignment = 'left', 'top'
                else :
                    x_text, y_text = x-offset*b*np.abs(np.sin(theta_modulo)), y-offset*b*np.abs(np.cos(theta_modulo))
                    horizontalalignment, verticalalignment = 'right', 'top'
                self.fits_image.mpl_ax.text( x_text, y_text, text, c=text_color, alpha=1, fontsize=15, \
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
            file_path = os.path.join(os.path.dirname(self.fits_image.image_path), 'mult.lenstool')
        
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
        print('Selected sources exported at ' + file_path)
                
    
    def export_to_potfile(self, file_path=None) :
        cat = self.cat[self.selection_mask] if True in self.selection_mask else self.cat
        ref_col = self.mag_colnames[0]
        sort_array = np.argsort(cat[ref_col])
        sorted_cat = cat[sort_array]
        
        lines = []
        lines.append('#REFERENCE 0\n')
        lines.append('## id   RA   Dec        a        b        theta     mag       lum\n')
        
        for i, galaxy in enumerate(sorted_cat) :
            lines.append( '%d %f %f %f %f %f %f 0.\n' % (i+1, galaxy['ra'], galaxy['dec'], galaxy['a'], galaxy['b'], galaxy['theta']-self.fits_image.orientation, galaxy[ref_col]) )
        
        if file_path is None :
            if self.ref_path is not None :
                file_path = os.path.join( os.path.dirname(self.ref_path), 'exported_potfile.lenstool')
            else :
                file_path = os.path.join( os.path.dirname(self.fits_image.image_path), 'exported_potfile.lenstool')
        print('Exporting selected sources to ' + file_path)
        with open(file_path, 'w') as file:
            file.writelines(lines)
        
    
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
            self.fits_image.qt_plot.removeItem(text_item)
        self.text_items.clear()
    
        # Add new text labels
        for i in range(len(self.cat)):
            text = str(self.cat[text_column][i])
            text_item = pg.TextItem(text, color=color)
    
            # Get ellipse position and size for offset calculation
            x = self.cat['x'][i]
            y = self.fits_image.image_data.shape[0] - self.cat['y'][i]  # Flip y to match PyQtGraph convention
            semi_major = self.cat['a'][i]
            semi_minor = self.cat['b'][i]
    
            # Position text slightly offset from the ellipse
            offset = max(semi_major, semi_minor)
            text_item.setPos(x + offset/2, y - offset/2)
    
            # Set font
            font = PyQt5.QtGui.QFont()
            font.setPointSize(15)
            text_item.setFont(font)

            self.fits_image.qt_plot.addItem(text_item)
            self.text_items.append(text_item)
            
    
    
    def export_thumbnails(self, mask=None, group_images=True) :
        if mask is None :
            mask = self.selection_mask
        ### TO DO ###







