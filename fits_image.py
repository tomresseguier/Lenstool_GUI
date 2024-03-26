import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse, Polygon, Circle
from astropy.io import fits
from astropy.wcs import WCS
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
#from PyQt5.QtWidgets import QGraphicsEllipseItem, QGraphicsScene, QGraphicsView, QGraphicsSceneMouseEvent, QApplication
#from PyQt5.QtCore import Qt
from pyqtgraph.Qt import QtCore
from astropy.table import Table
import lenstool
import pandas as pd

#import sys
#module_dir = os.path.dirname(os.path.abspath(__file__))
#sys.path.append(module_dir)

from source_extraction.source_extract import source_extract
from utils.utils_plots.plot_utils_general import *
from utils.utils_Qt.selectable_classes import *
from utils.utils_Qt.utils_general import *
from utils.utils_general.utils import flux_muJy_to_magAB


pg.setConfigOption('imageAxisOrder', 'row-major')




class fits_image :
    def __init__(self, image_path) :
        self.image_path = image_path
        self.pix_deg_scale, self.wcs, self.image_data = self.open_image()
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
    
    def open_image(self) :
        with fits.open(self.image_path) as hdus :
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
    
    def plot_image(self) :
        #if self.qt_plot is None :
        #to_plot = self.image_data
        #to_plot = np.transpose(self.image_data, axes=[1,0,2])
        to_plot = np.flip(self.image_data, axis=0)
        
        self.qt_plot = pg.image(to_plot)
        
        ############### NEW CODE ###############
        self.window = QMainWindow()
        self.window.setWindowTitle(os.path.basename(self.image_path))
        self.image_widget_layout = QHBoxLayout()
        
        self.image_widget_layout.addWidget(self.qt_plot)
        #self.qt_plot.setSizePolicy(pg.QtWidgets.QSizePolicy.Fixed, pg.QtWidgets.QSizePolicy.Expanding)
        
        self.image_widget = QWidget()
        self.image_widget.setLayout(self.image_widget_layout)
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
    
    def extract_sources(self, image_path=None, rerun=False) :
        if image_path is None :
            self.sources_all = source_extract(self.image_path, weight_path=None, zero_point=None, out_dir='PWD',
                                              outfile_name='SExtractor_cat.fits', return_sources=True, rerun=rerun)
        else :
            self.sources_all = source_extract(image_path, weight_path=None, zero_point=None, out_dir='PWD',
                                              outfile_name='SExtractor_cat.fits', return_sources=True, rerun=rerun)
        
        x, y = self.world_to_image(self.sources_all['RA'], self.sources_all['DEC'], unit='deg')
        if x[0] != self.sources_all['X_IMAGE'][0] :
            print('Catalog sextracted from different image: replacing X_IMAGE, Y_IMAGE and THETA_IMAGE columns with current image coordinates.')
            self.sources_all['X_IMAGE'], self.sources_all['Y_IMAGE'] = x, y
            #self.wcs.get_pc()[0, 1]
            ref_image_angle = ( np.arctan2(self.wcs.wcs.get_pc()[1, 0], self.wcs.wcs.get_pc()[0, 0]) %np.pi ) * 360/np.pi
            self.sources_all['THETA_IMAGE'] = self.sources_all['THETA_WORLD'] + ref_image_angle
        
        #mask_mag = self.sources_all['MAG_AUTO']<-10.
        #mask = self.sources['KRON_RADIUS']
        #mask_galstar = self.sources_all['CLASS_STAR']<0.4
        #mask_size = self.sources_all['A_IMAGE']*self.sources_all['B_IMAGE']*np.pi>1000.
        #mask = mask_mag & mask_galstar & mask_size
        
        #self.sources = self.sources_all #[mask]
        self.sources = self.make_catalog(self.sources_all)
        return str(len(self.sources.cat)) + ' sources found.'
    
    def select_multiple_images(self) :
        return 'in progress'
    
    def import_multiple_images(self, mult_file_path) :
        #lt = lenstool.Lenstool()
        multiple_images = Table(names=['id','ra','dec','a','b','theta','z','mag'], dtype=['str',*['float',]*7])
        with open(mult_file_path, 'r') as mult_file:
            for line in mult_file:
                cleaned_line = line.strip()
                if not cleaned_line.startswith("#") and len(cleaned_line)>0 :
                    split_line = cleaned_line.split()
                    row = [split_line[0]]
                    for element in split_line[1:8] :
                        row.append(float(element))
        #            if row[-2]==0. :
        #                row[-2] = -0.01 #set a non zero redshift to keep the Lenstool wrapper happy
                    multiple_images.add_row(row)
        #lt.set_sources(multiple_images)
        
        self.multiple_images = self.make_catalog(multiple_images)
        
        def plot_multiple_images(self, which='all', size=80, alpha=0.7, marker='o', filled_markers=False) :
            if which=='all' :
                families = np.unique( [self.cat[i]['id'][0] for i in range(len(self.cat))] )
                which = families
            to_plot_mask = []
            for symbol in which :
                to_plot_mask.append( [self.cat[i]['id'].startswith(symbol) for i in range(len(self.cat))] )
            
            colors = make_palette(len(which), 1, alpha=1)
            if filled_markers :
                facecolors = make_palette(len(which), 1, alpha=0.3)
            else :
                facecolor = [[0,0,0,0] for i in range(len(which))]
            
            count = 0
            for i, mask in enumerate(to_plot_mask) :
                for multiple_image in self.cat[mask] :
                    # Remove the *1000
                    if len(np.unique(multiple_image['a']))==1 :
                        a, b = 75, 75
                    else :
                        a, b = multiple_image['a'], multiple_image['b']
                    ellipse = self.plot_one_object(multiple_image['x'], multiple_image['y'], \
                                                   a, b, \
                                                   multiple_image['theta'], count, color=colors[i][:3])
                    self.qtItems[count] = ellipse
                    count += 1
        
        self.multiple_images.plot = types.MethodType(plot_multiple_images, self.multiple_images)
        
        #self.multiple_images = multiple_images
        return self.multiple_images.cat
    
    def load_potfile(self, potfile_path) :
        with open(potfile_path, 'r') as file :
            lines = file.readlines()
        potfile_cat = Table(names=['id','ra','dec','a','b','theta','mag','lum'], dtype=['int', *['float',]*7])
        for line in lines :
            columns = line.split()
            if columns[0].isdigit() :
                potfile_cat.add_row( [col for col in columns] )
        
        self.potfile = self.make_catalog(cat=potfile_cat)
        return self.potfile.cat
    
    def open_cat(self, cat_path) :
        with open(cat_path, 'r') as raw_cat :
            first_line = raw_cat.readlines()[0]
        if len(first_line.split()) > len(first_line.split(',')) :
            cat_df = pd.read_csv(cat_path, delim_whitespace=True)[1:].apply(pd.to_numeric, errors='coerce')
        else :
            cat_df = pd.read_csv(cat_path)[1:].apply(pd.to_numeric, errors='coerce')
        cat = Table.from_pandas(cat_df)
        return cat
    
    def import_catalog(self, cat_path=None, cat=None) :
        if cat==None :
            cat = self.open_cat(cat_path)
        self.imported_cat = self.make_catalog(cat=cat, make_selection_panel=True)
        return self.imported_cat.cat
    
    ################## Transform catalog into catalog class ###################
    
    def make_colnames_dict(self, catalog):
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
        to_test_names_dict['ra'] = ['ra', 'ALPHA_J2000']
        to_test_names_dict['dec'] = ['dec', 'DELTA_J2000']
        to_test_names_dict['x'] = ['x', 'X_IMAGE']
        to_test_names_dict['y'] = ['y', 'Y_IMAGE']
        to_test_names_dict['a'] = ['a', 'A_IMAGE']
        to_test_names_dict['b'] = ['b', 'B_IMAGE']
        to_test_names_dict['theta'] = ['theta', 'THETA_IMAGE', 'angle']
        
        names_list = list(to_test_names_dict.keys())
        #names_list = ['ra', 'dec', 'x', 'y', 'a', 'b']
        names_dict = {}
        for name in names_list :
            names_dict[name] = None
        for name in names_list :
            for to_test_name in to_test_names_dict[name] :
                cat_colnames_lower = [col.lower() for col in catalog.colnames]
                
                #if 'colnames' in dir(catalog) :
                #    cat_colnames_lower = [col.lower() for col in catalog.colnames]
                #else :
                #    cat_colnames_lower = [col.lower() for col in catalog.columns.names]
                
                if to_test_name.lower() in cat_colnames_lower :
                    col_idx = np.where( np.array(cat_colnames_lower) == to_test_name.lower() )[0][0]
                    #names_dict[name] = to_test_name
                    names_dict[name] = catalog.colnames[col_idx]
        return names_dict
    
    def make_uniform_names_cat(self, cat) :
        uniform_names_cat = cat.copy()
        colnames_dict = self.make_colnames_dict(cat)
        
        print('Column names to be used:')
        print(colnames_dict)
        
        for colname in colnames_dict.keys() :
            if colnames_dict[colname] is not None :
                uniform_names_cat.rename_column(colnames_dict[colname], colname)
                if colnames_dict[colname] == 'a' or colnames_dict[colname] == 'b' :
                    uniform_names_cat.replace_column( colname, uniform_names_cat[colname]/(self.pix_deg_scale*3600)*10 )
        if colnames_dict['x']==None :
            x, y = self.world_to_image(uniform_names_cat['ra'], uniform_names_cat['dec'], unit='deg')
            uniform_names_cat['x'] = x
            uniform_names_cat['y'] = y
        if colnames_dict['a']==None :
            size = 40.
            uniform_names_cat['a'] = np.full(len(uniform_names_cat), size)
            uniform_names_cat['b'] = np.full(len(uniform_names_cat), size)
            uniform_names_cat['theta'] = np.full(len(uniform_names_cat), 0.)
        return uniform_names_cat
    
    def make_catalog(self, cat=None, cat_path=None, make_selection_panel=False) :
        if cat_path is not None :
            cat = self.open_cat(cat_path)
        uniform_names_cat = self.make_uniform_names_cat(cat)
        if self.qt_plot is None :
            self.plot_image()
        return self.catalog(uniform_names_cat, self.image_data, self.qt_plot, make_selection_panel=make_selection_panel, image_widget_layout=self.image_widget_layout)
        
    ###########################################################################
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    class catalog :
        def __init__(self, cat, image_data, qt_plot, make_selection_panel=False, image_widget_layout=None) :
            self.cat = cat
            self.image_data = image_data
            self.qt_plot = qt_plot
            self.qtItems = np.empty(len(cat), dtype=PyQt5.QtWidgets.QGraphicsEllipseItem)
            #self.qtItems = np.empty(len(cat), dtype=utils.utils_classes.selectable_ellipse.SelectableEllipse)
            self.color = [1., 1., 0]
            self.selection_mask = np.full(len(cat), False)
            self.make_selection_panel = make_selection_panel
            self.image_widget_layout = image_widget_layout
            self.RS_widget = None
            self.x_axis_cleaned = np.full(len(cat), None)
            self.y_axis_cleaned = np.full(len(cat), None)
        
        def make_mask_naninf(self) :
            mag_F444W = flux_muJy_to_magAB(self.cat['f444w_tot_0'])
            mag_F090W = flux_muJy_to_magAB(self.cat['f090w_tot_0'])
            x_axis = mag_F444W
            y_axis = mag_F090W - mag_F444W
            
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
        
        def plot(self, scale=1.) :
            x = self.cat['x']
            y = self.cat['y']
            semi_major = self.cat['a'] * scale
            semi_minor = self.cat['b'] * scale
            angle = self.cat['theta']
            for i in tqdm(range(len(semi_major))) :
                ellipse = self.plot_one_object(x[i], y[i], semi_major[i], semi_minor[i], angle[i], i)
                self.qtItems[i] = ellipse
        
        def clear(self) :
            for i in tqdm( range(len(self.qtItems)) ) :
                self.qt_plot.removeItem(self.qtItems[i])
        
        def plot_one_object(self, x, y, semi_major, semi_minor, angle, idx) :
            color = list(np.array(self.color)*255)
            #make the flip to accomosate pyqtgraph's strange plotting conventions
            y = self.image_data.shape[0] - y
            angle = -angle
            #####################################################################
            ellipse = SelectableEllipse(x-semi_major/2, y-semi_minor/2, semi_major, semi_minor, idx, self.selection_mask, \
                                        self.qtItems, color, scatter_pos=(self.x_axis_cleaned[idx], self.y_axis_cleaned[idx]), RS_widget=self.RS_widget)
            #ellipse = PyQt5.QtWidgets.QGraphicsEllipseItem(x-semi_major/2, y-semi_minor/2, semi_major, semi_minor)
            ellipse.setTransformOriginPoint( PyQt5.QtCore.QPointF(x, y) ) 
            #ellipse.setTransform( PyQt5.QtGui.QTransform().rotate(angle[i]) )
            ellipse.setRotation(angle)
            self.qt_plot.addItem(ellipse)
            return ellipse
        
        def plot_selection_panel(self) :
            if self.make_selection_panel :
                self.make_mask_naninf()
                
                self.RS_widget, self.selection_ROI = plot_panel(self.x_axis_cleaned, self.y_axis_cleaned, self.image_widget_layout, self.qt_plot)
                
                data=(self.x_axis_cleaned, self.y_axis_cleaned)
                #self.selection_mask = np.full(len(self.mag_F444W_cleaned), False)
                self.selectable_scatter = SelectableScatter(self.RS_widget, self.selection_ROI, data, self.selection_mask, qtItems=self.qtItems, color=list(np.array(self.color)*255))
                
                
            
        #def select(self) :
        #    class MouseClickHandler(QtCore.QObject) :
        #        clicked = QtCore.Signal(float, float)
        #        def __init__(self, plot_item):
        #            super(MouseClickHandler, self).__init__()
        #            self.plot_item = plot_item
        #            self.plot_item.scene().sigMouseClicked.connect(self.mouse_click_event)
        #        def mouse_click_event(self, event):
        #            if event.double():
        #                # Ignore double-click events
        #                return
        #            # Map the mouse click position to the plot coordinates
        #            pos = self.plot_item.vb.mapSceneToView(event.scenePos())
        #            x, y = pos.x(), pos.y()
        #            # Emit the clicked signal with the coordinates
        #            self.clicked.emit(x, y)            
        #    # Create a handler for mouse clicks
        #    click_handler = MouseClickHandler(self.qt_plot)
        #    
        #    def on_mouse_click(x, y):
        #        print(f"Mouse clicked at ({x}, {y})")
        #    
        #    # Connect the handler's signal to a custom slot
        #    click_handler.clicked.connect(on_mouse_click)
    
    
    
    
    
    
    
    
    
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    def plot_image_mpl(self, wcs_projection=True, units='pixel', pos=111, make_axes_labels=True) :
        if self.ax is None :
            self.fig = plt.figure()
        if wcs_projection :
            self.ax = self.fig.add_subplot(pos, projection=self.wcs)
            self.ax.coords.grid(True, color='black', ls='dotted')
        else :
            self.ax = self.fig.add_subplot(pos)
            if units=='pixel' or units=='pixels' or units=='image' :
                scaling = 1
            if units=='arcsec' :
                scaling = pix_deg_scale*60*60
            if units=='arcmin' :
                scaling = pix_deg_scale*60
        if make_axes_labels and wcs_projection :
            self.ax.coords[0].set_axislabel('Right ascension')
            self.ax.coords[1].set_axislabel('Declination')
        elif make_axes_labels and not wcs_projection :
            self.ax.set_xlabel('x (' + units + ')')
            self.ax.set_ylabel('y (' + units + ')')
        elif not make_axes_labels and wcs_projection :
            self.ax.coords[0].set_axislabel(' ')
            self.ax.coords[1].set_axislabel(' ')
        else :
            self.ax.set_xlabel(' ')
            self.ax.set_ylabel(' ')
        if wcs_projection :
            self.ax.imshow(self.image_data, origin="lower")
        if not wcs_projection :
            self.ax.imshow(self.image_data, origin='lower', extent=[0, self.image_data.shape[1]*scaling, 0, self.image_data.shape[0]*scaling])
        #ax.figure.tight_layout()
        return self.fig, self.ax
    
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
    
        
        




