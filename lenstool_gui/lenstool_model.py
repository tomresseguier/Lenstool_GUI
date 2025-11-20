import os
import glob
import shutil
import numpy as np
from matplotlib import pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from tqdm import tqdm
import copy

import PyQt5
import pyqtgraph as pg
import pyqtgraph.exporters
from PyQt5.QtCore import Qt
from astropy.table import Table
import pickle
import lenstool
import pylenstool

###############################################################################
from .utils.utils_astro.utils_general import world_to_relative, relative_to_world
from .utils.utils_lenstool_model import import_multiple_images, import_sources, export_thumbnails, curves
from .utils.utils_simulate_image import save_lenstronomy_model, load_lenstronomy_model, make_LENSTRONOMY_plot
from .utils.utils_plots.plot_utils_general import plot_corner
from .utils.utils_general.utils_general import extract_line
from .utils.utils_general.sort_points import break_curves
from .utils.utils_Lenstool.param_extractors import read_potfile, make_best_file_from_bayes, make_param_latex_table, read_bayes_file
from .simulate_image import lenstronomy_model

from .utils.utils_Lenstool.file_makers import best_files_maker, make_magnifications_and_curves                  # This import is problematic. The two functions run Lenstool
                                                                                                                # and are therefore dependent on my own install.
from .utils.utils_Lenstool.operations import make_image_to_source, MakeFunctionFromMap





class lenstool_model :
    def __init__(self, model_path, fits_image) :
        self.safe_mode = False
        self.fits_image = fits_image
        self.saturation = 1
        self.model_dir = model_path if os.path.isdir(model_path) else os.path.dirname(model_path)
        all_par_file_paths = glob.glob(os.path.join(self.model_dir, "*.par"))
        #all_cat_file_paths = glob.glob(os.path.join(self.model_dir, "*.lenstool"))
        self.param_file_path = None
        
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
        
        if self.param_file_path is not None :
            param_file = pylenstool.lenstool_param_file(self.param_file_path)
            ref_coord = param_file.get_ref()
            self.reference = ( float(ref_coord[0]), float(ref_coord[1]) )
        
        all_par_file_names = [ os.path.basename(file_path) for file_path in all_par_file_paths ]
        
        self.has_run = 'best.par' in all_par_file_names
        self.best_file_path = os.path.join(self.model_dir, 'best.par') if 'best.par' in all_par_file_names else None
        self.bayes_file_path = os.path.join(self.model_dir, 'bayes.dat') if 'bayes.dat' in os.listdir(self.model_dir) else None
        if self.best_file_path is None and 'bayes.dat' in os.listdir(self.model_dir) :
            yesno = input('bayes.dat file found, create best_TEMP.par from bayes file? ([Y]/n)')
            if yesno.lower() in ['y', ''] :
                make_best_file_from_bayes(self.param_file_path)
                self.best_file_path = os.path.join(self.model_dir, 'best_TEMP.par')
        
        self.bestopt_file_path = os.path.join(self.model_dir, 'bestopt.par') if 'bestopt.par' in all_par_file_names else None
        potfile_paths_list = glob.glob(os.path.join(self.model_dir, "*potfile*.lenstool"))
        self.potfile_path = potfile_paths_list[0] if len(potfile_paths_list)>=1 else None
        
        if self.potfile_path is not None :
            potfile_Table = read_potfile(self.potfile_path)
            self.potfile = fits_image.make_catalog(potfile_Table, color=[1.,0.,0.], units='arcsec')
        else :
            self.potfile = None
        
        self.families = []
        self.broad_families = []
        self.which = []
        
        mult_path_list = glob.glob( os.path.join(self.model_dir, "*mult*.lenstool") )
        if len(mult_path_list)==1 :
            print("Multiple images file found: '" + mult_path_list[0] + "'")
            self.mult_file_path = mult_path_list[0]
            import_multiple_images(self, self.mult_file_path, fits_image, units='pixel', filled_markers=True)
        else :
            self.mult = None
        
        arclets_path_list = glob.glob( os.path.join(self.model_dir, "*arclet*.lenstool") )
        if len(arclets_path_list)==1 :
            arclets_path = arclets_path_list[0]
            import_multiple_images(self, arclets_path, fits_image, AttrName='arclets', units='pixel', filled_markers=False)
        else :
            self.arclets = None
        
        if len(arclets_path_list)==1 and len(mult_path_list)==1 :
            import_multiple_images(self, self.mult_file_path, fits_image, units='pixel', filled_markers=True)
            import_multiple_images(self, arclets_path, fits_image, AttrName='arclets', units='pixel', filled_markers=False)
        
        #dot_all_paths = glob.glob(os.path.join(self.model_dir, "*.all"))
        predicted_images_path = os.path.join(self.model_dir, 'image.dat')
        if os.path.isfile(predicted_images_path) :
            import_multiple_images(self, predicted_images_path, fits_image, AttrName='image', units='pixel', filled_markers=False)
            import_multiple_images(self, predicted_images_path, fits_image, AttrName='image_filtered', units='pixel', filled_markers=False)
            self.filter_image()
        else :
            self.image = None
            self.image_filtered = None
        
        curves_dir = os.path.join(self.model_dir, 'curves')
        if os.path.isdir(curves_dir) :
            self.curves = curves(curves_dir, self, fits_image, which_critcaus='critical', join=False, size=2)
        else :
            self.curves = None
        
        predicted_sources_path = os.path.join(self.model_dir, 'source.dat')
        if os.path.isfile(predicted_sources_path) :
            import_sources(self, predicted_sources_path, fits_image, AttrName='source', units='pixel', filled_markers=False)
        
        self.curve_plot = None
        
        self.lt = None
        #######################################################################
        self.dpl_maps_path = os.path.join(self.model_dir, 'dpl_maps.pkl')
        if os.path.exists(self.dpl_maps_path) :
            print('dpl_maps.pkl found')
            with open(self.dpl_maps_path, 'rb') as f :
                self.dpl_maps = pickle.load(f)
        else :
            self.dpl_maps = {}
        #######################################################################
        self.lt_curves_path = os.path.join(self.model_dir, 'lt_curves.pkl')
        if os.path.exists(self.lt_curves_path) :
            print('lt_curves.pkl found')
            with open(self.lt_curves_path, 'rb') as f :
                self.lt_curves = pickle.load(f)
        else :
            self.lt_curves = {}
        #######################################################################
        self.lt_magnification_maps_path = os.path.join(self.model_dir, 'lt_magnification_maps.pkl')
        if os.path.exists(self.lt_magnification_maps_path) :
            print('lt_magnification_maps.pkl found')
            with open(self.lt_magnification_maps_path, 'rb') as f :
                self.lt_magnification_maps = pickle.load(f)
        else :
            self.lt_magnification_maps = {}
        #######################################################################
        self.lt_caustics_path = os.path.join(self.model_dir, 'lt_caustics.pkl')
        if os.path.exists(self.lt_caustics_path) :
            print('lt_caustics.pkl found')
            with open(self.lt_caustics_path, 'rb') as f :
                self.lt_caustics = pickle.load(f)
        else :
            self.lt_caustics = {}
        
        self.z_lens = None
        self.magnification_res = 1000
        self.magnification_line_ax = None
        self.previous_state_current_ROI = None
        self.LENSTRONOMY_fixed_source_kwargs = []
        
    
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
            self.mult.plot(marker='o', filled_markers=False, scale=1.25)#size=1.5, linewidth=2, filled_markers=False)
            self.mult.plot_column('id')
        if self.image is not None :
            #self.image.plot(marker='x', filled_markers=True, scale=1)
            self.image_filtered.plot(marker='s', filled_markers=True, scale=0.5)
            self.image.saturation = 1.
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
        if self.image_filtered is not None :
            self.image_filtered.clear()
        if self.curves is not None :
            self.curves.clear()
        if self.curve_plot is not None :
            self.fits_image.qt_image.removeItem(self.curve_plot)
            
    def set_which(self, *names) :
        if names[0]=='all' :
            self.which = self.broad_families.tolist()
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
    
    
    def export_thumbnails(self, group_images=True, square_thumbnails=True, square_size=150, margin=50, distance=200, export_dir=None, boost=True, make_broad_view=True, broad_view_params=None) :
        export_thumbnails(self.mult, group_images=group_images, square_thumbnails=square_thumbnails, square_size=square_size, margin=margin, \
                          distance=distance, export_dir=export_dir, boost=boost, make_broad_view=make_broad_view, broad_view_params=broad_view_params)
        
    
    def world_to_relative(self, ra, dec) :
        return world_to_relative(ra, dec, self.reference)
    
    def relative_to_world(self, xr, yr) :
        return relative_to_world(xr, yr, self.reference)
    
    def make_webpage(self) :
        print('in progress')
        
    def make_latex(self) :
        latex_str = make_param_latex_table(self.param_file_path, convert_to_kpc=True, z=self.z_lens)
        return latex_str
    
    
    def set_lt_z(self, z, color=[255,100,255], recompute=False) :
        if not self.safe_mode :
            self.SafeMode()
        #if self.curve_plot is not None :
        #    self.fits_image.qt_image.removeItem(self.curve_plot)
        
        self.lt_z = z
        print(self.best_file_path)
        print(os.getcwd())
        if self.lt==None :
            self.lt = lenstool.Lenstool( os.path.basename(self.best_file_path) )
        
        #self.lt.set_grid(50, 0)
        
        
        
        ######## Curves ########
        if z not in self.lt_curves.keys() or recompute :
            self.compute_lt_curve(z)
        self.lt_curve_coords_image = self.lt_curves[z]
        self.lt_caustic_coords_image = self.lt_caustics[z]
        self.plot_lt_curve(color=color)
        
        ######## Magnification ########
        if z not in self.lt_magnification_maps.keys() or recompute :
            print('Computing magnification map (can take a little while)...')
            self.lt_magnification_maps[z] = self.lt.g_ampli(1, self.magnification_res, self.lt_z)
            print('done')
            with open(self.lt_magnification_maps_path, 'wb') as f:
                pickle.dump(self.lt_magnification_maps, f)
        self.magnification_map, self.magnification_wcs = self.lt_magnification_maps[z]            
        
        ######## Displacement maps ########
        if self.lt_z not in self.dpl_maps.keys() or recompute :
            self.compute_lt_dpl()
        self.dx_map, self.dy_map, self.dmap_wcs = self.dpl_maps[self.lt_z]
        
        mmap, wcs = self.lt_magnification_maps[z]
        self.get_magnification = MakeFunctionFromMap(mmap, wcs)
        
    def compute_lt_dpl(self, z=None, npix=2000) :
        if z==None :
            z = self.lt_z
        print('Computing displacement maps (can take a little while)...')
        self.dpl_maps[self.lt_z] = self.lt.g_dpl(npix, z)
        print('done')
        with open(self.dpl_maps_path, 'wb') as f:
            pickle.dump(self.dpl_maps, f)
    
    
    def start_im2source(self) :
        self.lt.set_grid(50, 0)
        
        self.transform_coords_radec = make_image_to_source(self.dx_map, self.dy_map, self.dmap_wcs)
        
        def transform_coords(x, y) :
            ra, dec = self.fits_image.image_to_world(x, y)
            
            ra_source, dec_source = self.transform_coords_radec(ra, dec)
            x_source, y_source = self.fits_image.world_to_image(ra_source, dec_source)
            return x_source, y_source, ra_source, dec_source
            
        self.transform_coords = transform_coords
        
        # Set up label if needed (optional)
        if not hasattr(self, 'transform_label'):
            self.transform_label = pg.TextItem(anchor=(0, 1), color='w')
            self.fits_image.qt_image.addItem(self.transform_label)
        
        # Create scatter point for transformed location
        self.transformed_point = pg.ScatterPlotItem(size=10, brush='r')
        self.fits_image.qt_image.addItem(self.transformed_point)
        
        self.images_scatter = pg.ScatterPlotItem(size=10, symbol='+', brush='g')
        self.fits_image.qt_image.addItem(self.images_scatter)
    
        # Ensure view does not auto-range when updating
        self.fits_image.qt_image.getView().enableAutoRange(pg.ViewBox.XAxis, False)
        self.fits_image.qt_image.getView().enableAutoRange(pg.ViewBox.YAxis, False)
    
        # Mouse tracking setup
        def mouse_moved(evt):
            pos = evt[0]
            if self.fits_image.qt_image.getView().sceneBoundingRect().contains(pos):
                mouse_point = self.fits_image.qt_image.getView().mapSceneToView(pos)
                x, y_flipped = mouse_point.x(), mouse_point.y()
                x, y = x, self.fits_image.image_data.shape[0] - y_flipped
                try:
                    x_source, y_source, ra_source, dec_source = self.transform_coords(x, y)
                    xr_source, yr_source = self.world_to_relative(ra_source, dec_source)
                    
                    #if int( importlib.metadata.version('lenstool').split('.')[1] )>=6 :
                        ### For Lenstool version 8.6.3 ###
                    #print(str(ra_source)[3:], str(dec_source)[3:])
                    source = Table( rows=[('test', ra_source, dec_source, 1, 1, 0, self.lt_z, 25)],names=['n','x','y','a','b','theta','z','mag'], dtype=['str',*['float',]*7] )
                    #else :
                        ### For Lenstool version ?? ###
                        #source = Table( rows=[('test', xr_source, yr_source, 1, 1, 0, self.lt_z, 25)],names=['n','x','y','a','b','theta','z','mag'], dtype=['str',*['float',]*7] )
                    
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
                        
                        #ellipse = QGraphicsEllipseItem(x_image, self.fits_image.image_data.shape[0] - y_image, image['a'], image['b'])
                        #ellipse.setTransformOriginPoint( PyQt5.QtCore.QPointF(x_image, self.fits_image.image_data.shape[0] - y_image) )
                        #ellipse.setRotation(-image['theta'])
                    self.images_scatter.setData( x_images, self.fits_image.image_data.shape[0] - np.array(y_images) )
                    
                    self.transformed_point.setData([x_source], [self.fits_image.image_data.shape[0] - y_source])
                    self.transform_label.setText(f"({x:.2f}, {y:.2f}) → ({x_source:.2f}, {y_source:.2f})")
                    self.transform_label.setPos(x, y_flipped)
                    
                    self._last_transform_coords = {'x': x, 'y': y, 'x_source': x_source, 'y_source': y_source}
                except Exception as e:
                    self.transform_label.setText(f"Error: {e}")
                    self.transform_label.setPos(x, y_flipped)
        
        self._transform_proxy = pg.SignalProxy(self.fits_image.qt_image.getView().scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)
        
        
        
        
        self.doubleclick_image_marker = pg.ScatterPlotItem(size=12, symbol='s', brush='g', pen='g')
        self.doubleclick_source_marker = pg.ScatterPlotItem(size=12, symbol='+', brush='r', pen='r')
        self.fits_image.qt_image.addItem(self.doubleclick_image_marker)
        self.fits_image.qt_image.addItem(self.doubleclick_source_marker)
        
        self.source_markers_x = []
        self.source_markers_y = []
        
        def mouse_clicked(evt):
            if evt.double():
                if hasattr(self, '_last_transform_coords'):
                    coords = self._last_transform_coords
                    x, y = coords['x'], coords['y']
                    x_source, y_source = coords['x_source'], coords['y_source']
                    
                    self.source_markers_x.append(x_source)
                    self.source_markers_y.append(self.fits_image.image_data.shape[0] - y_source)
                    
                    self.doubleclick_image_marker.setData([x], [self.fits_image.image_data.shape[0] - y])
                    #self.doubleclick_source_marker.setData([x_source], [self.fits_image.image_data.shape[0] - y_source])
                    self.doubleclick_source_marker.setData(self.source_markers_x, self.source_markers_y)
        self._doubleclick_connection = self.fits_image.qt_image.scene.sigMouseClicked.connect(mouse_clicked)
        
        
        def keyPressEvent(event):
            if event.key() == Qt.Key_Escape or event.key() == Qt.Key_Space :
                self.stop_im2source()
        
        self._original_keyPressEvent = self.fits_image.window.keyPressEvent
        self.fits_image.window.keyPressEvent = keyPressEvent
        
    def stop_im2source(self) :
        self.fits_image.qt_image.removeItem(self.images_scatter)
        self.fits_image.qt_image.removeItem(self.transformed_point)
        self.fits_image.qt_image.removeItem(self.transform_label)
        self._transform_proxy.disconnect()
        del self._transform_proxy
        del self.images_scatter
        del self.transformed_point
        del self.transform_label
        self.fits_image.window.keyPressEvent = self._original_keyPressEvent
        if hasattr(self, 'doubleclick_image_marker'):
            self.fits_image.qt_image.removeItem(self.doubleclick_image_marker)
            del self.doubleclick_image_marker
        if hasattr(self, 'doubleclick_source_marker'):
            self.fits_image.qt_image.removeItem(self.doubleclick_source_marker)
            del self.doubleclick_source_marker
        if hasattr(self, '_doubleclick_connection'):
            self.fits_image.qt_image.scene.sigMouseClicked.disconnect(self._doubleclick_connection)
            del self._doubleclick_connection
        print('Interactive source & images viewer closed.')
    
    
    def start_magnification(self) :
        return None
    
    def add_lensing_columns(self, cat=None, which_cat='imported_cat', index=None) :
        if cat is None :
            if index is not None :
                cat = self.fits_image.imported_cat_list[index].cat
            else :
                cat = getattr(self.fits_image, which_cat, None).cat
        
        check = False
        for name in ['magnification', 'convergence', 'shear', 'tangential_magnification']:
            if name in cat.colnames :
                check = True
        if check :
            yesno = input('Lensing column already exists. Overwrite? [Y]/n')
        
        if yesno.lower() in ['y', ''] :
            mu_col = np.full(len(cat), np.nan)
            gamma_col = np.full(len(cat), np.nan)
            kappa_col = np.full(len(cat), np.nan)
            tmu_col = np.full(len(cat), np.nan)
            rmu_col = np.full(len(cat), np.nan)
            
            initial_field = self.lt.get_field([])
            
            z_colname = None
            for name in ['z_spec', 'zspec', 'z_phot', 'zphot'] :
                if name in cat.colnames :
                    z_colname = name
                    break
            if z_colname is None :
                print('Redshift column not found in catalog. Using current source redshift = ' + str(self.lt_z))
                
            for i in tqdm(range(len(cat))) :
                if z_colname is not None :
                    z = cat[z_colname][i]
                    if z>self.z_lens :
                        #print(str(cat['id'][i]) + ': computing lensing maps at redshift ' + str(z))
                        xr, yr = self.world_to_relative(cat['ra'][i], cat['dec'][i])
                        #delta = 1.
                        delta = self.fits_image.pix_deg_scale * 3600 / 2
                        self.lt.set_field([xr-delta, xr+delta, yr-delta, yr+delta])
                        
                        #npix = 11
                        npix = 2
                        mmap, wcs = self.lt.g_ampli(1, npix, z)
                        #get_magnification = MakeFunctionFromMap(mmap, wcs)
                        #mu_col[i] = get_magnification(cat['ra'][i], cat['dec'][i])
                        mu_col[i] = np.mean(mmap)
                        
                        kappa, wcs = self.lt.g_mass(1, npix, self.z_lens, z)
                        #get_convergence = MakeFunctionFromMap(kappa, wcs)
                        #kappa_col[i] = get_convergence(cat['ra'][i], cat['dec'][i]) 
                        kappa_col[i] = np.mean(kappa)
                        
                        #shear, wcs = self.lt.g_ampli(6, npix, z)
                        #get_shear = MakeFunctionFromMap(shear, wcs)
                        #gamma_col[i] = get_shear(cat['ra'][i], cat['dec'][i])
                        #gamma_col[i] = np.mean(shear)
                        gamma_col[i] = ( (1-kappa_col[i])**2 - 1/mu_col[i] )**0.5
                    #else :
                        #print(str(cat['id'][i]) + ': redshift ' + str(z) + ' lower than lens redshift --> NaN')
                else :
                    mu_col[i] = self.get_magnification(cat['ra'][i], cat['dec'][i])
                    kappa_col[i] = self.get_convergence(cat['ra'][i], cat['dec'][i])     
                    gamma_col[i] = self.get_shear(cat['ra'][i], cat['dec'][i])
                tmu_col[i] = 1 / (1 - kappa_col[i] - gamma_col[i] )
                rmu_col[i] = 1 / (1 - kappa_col[i] + gamma_col[i] )
            
            columns_to_add = [mu_col, kappa_col, gamma_col, tmu_col, rmu_col]
            names = ['magnification', 'convergence', 'shear', 'tangential_magnification', 'radial_magnification']
            for i, name in enumerate(names) :
                if name in cat.colnames :
                    cat.replace_column(name, columns_to_add[i])
                else :
                    cat.add_column(columns_to_add[i], name=name)
            self.lt.set_field(initial_field)
    
    def compute_lt_curve(self, z=None) :
        if z==None :
            z = self.lt_z
        print('Computing critical curve (can take a little while)...')
        self.lt_curve = self.lt.criticnew(zs=z, limitHigh=1, limitLow=0.05) #limitHigh=0.5, limitLow=0.1
        
        ni = len(self.lt_curve[0])
        ne = len(self.lt_curve[1])
        lt_curve_xr = np.zeros(ni + ne)
        lt_curve_yr = np.zeros(ni + ne)
        for i in range(ni) :
            lt_curve_xr[i] = self.lt_curve[0][i].I.x
            lt_curve_yr[i] = self.lt_curve[0][i].I.y
        for i in range(ne) :
            lt_curve_xr[ni+i] = self.lt_curve[1][i].I.x
            lt_curve_yr[ni+i] = self.lt_curve[1][i].I.y
            
        lt_curve_ra, lt_curve_dec = self.relative_to_world(lt_curve_xr, lt_curve_yr)
        lt_curve_x, lt_curve_y = self.fits_image.world_to_image(lt_curve_ra, lt_curve_dec)
        self.lt_curve_coords_image = [lt_curve_x, self.fits_image.image_data.shape[0] - lt_curve_y]
        self.lt_curve_coords_world = [lt_curve_ra, lt_curve_dec]
        self.lt_curve_coords_relative = [lt_curve_xr, lt_curve_yr]
        print('done')
        
        self.lt_curves[z] = self.lt_curve_coords_image
        with open(self.lt_curves_path, 'wb') as f:
            pickle.dump(self.lt_curves, f)
        
        
        ###### Caustics ######
        ni = len(self.lt_curve[0])
        ne = len(self.lt_curve[1])
        lt_caustic_xr = np.zeros(ni + ne)
        lt_caustic_yr = np.zeros(ni + ne)
        for i in range(ni) :
            lt_caustic_xr[i] = self.lt_curve[0][i].S.x
            lt_caustic_yr[i] = self.lt_curve[0][i].S.y
        for i in range(ne) :
            lt_caustic_xr[ni+i] = self.lt_curve[1][i].S.x
            lt_caustic_yr[ni+i] = self.lt_curve[1][i].S.y
            
        lt_caustic_ra, lt_caustic_dec = self.relative_to_world(lt_caustic_xr, lt_caustic_yr)
        lt_caustic_x, lt_caustic_y = self.fits_image.world_to_image(lt_caustic_ra, lt_caustic_dec)
        self.lt_caustic_coords_image = [lt_caustic_x, self.fits_image.image_data.shape[0] - lt_caustic_y]
        self.lt_caustic_coords_world = [lt_caustic_ra, lt_caustic_dec]
        self.lt_caustic_coords_relative = [lt_caustic_xr, lt_caustic_yr]
        print('done')
        
        self.lt_caustics[z] = self.lt_caustic_coords_image
        with open(self.lt_caustics_path, 'wb') as f:
            pickle.dump(self.lt_caustics, f)
        
        self.plot_lt_curve()
    
    def plot_lt_curve(self, color=[255, 0, 255], which='critical') :
        if self.curve_plot is not None :
            self.fits_image.qt_image.removeItem(self.curve_plot)
            #del self.curve_plot
        
        if which=='critical' :
            coords = self.lt_curve_coords_image
        elif which=='caustic' :
            coords = self.lt_caustic_coords_image
        
        self.lt_curve_coords_image_sorted = break_curves(coords)
        #self.lt_curve_coords_image_sorted = sort_points(coords, distance_threshold=1.0/(self.fits_image.pix_deg_scale*3600), angle_threshold=np.pi)
        
        self.curve_plot = pg.PlotDataItem()
        #self.curve_plot = pg.ScatterPlotItem()
        self.curve_plot.setPen( color=color+[255], width=4.0001 )
        #self.curve_plot.setBrush( color=color+[255], width=3 )
        self.curve_plot.setData(self.lt_curve_coords_image_sorted[0], self.lt_curve_coords_image_sorted[1])
        #self.curve_scatter = pg.ScatterPlotItem(size=1, brush='g')
        #self.curve_scatter.setData(coords[0], coords[1])
        self.fits_image.qt_image.addItem(self.curve_plot)
        
    def plot_bayes(self) :
        if self.z_lens is None :
            self.z_lens = float(input('redshift of lens?'))
        self.df = read_bayes_file(self.bayes_file_path, z=self.z_lens)
        
        # Extract the numeric columns (skip non-numeric, zero-range etc.)
        self.df_param_only = self.df.select_dtypes(include='number')
        for i, col in enumerate(self.df_param_only.columns) :
            col_min = np.min(self.df_param_only[col])
            col_max = np.max(self.df_param_only[col])
            if col_min==col_max or col=='Chi2' or col=='Nsample' or col=='ln(Lhood)':
                del self.df_param_only[col]
        
        plot_corner(self.df_param_only)
        corr_matrix = self.df_param_only.corr()
        self.fig_cov, self.ax_cov = plt.subplots()
        cax = self.ax_cov.imshow(corr_matrix, cmap='PuOr')
        cbar = self.fig_cov.colorbar(cax, ax=self.ax_cov)
        self.ax_cov.set_xticks(np.arange(len(corr_matrix.columns)))
        self.ax_cov.set_yticks(np.arange(len(corr_matrix.index)))
        self.ax_cov.set_xticklabels(corr_matrix.columns, rotation=45, ha='right')
        self.ax_cov.set_yticklabels(corr_matrix.index)
    
    def filter_image(self, threshold_arcsec=1) :        
        to_remove = []
        for i, image in enumerate(self.image_filtered.cat) :
            ref_mask = self.mult.cat['id']==image['id']
            if len(np.unique(ref_mask))==1 and not np.unique(ref_mask)[0] :
                d = 0
            else :
                ref = self.mult.cat[ np.where(ref_mask)[0][0] ]
                d = ( (ref['x'] - image['x'])**2 + (ref['y'] - image['y'])**2 )**0.5
            if d==0 :
                to_remove.append(i)
                
        self.image_filtered.cat.remove_rows(to_remove)
        
        
        if False :
            threshold_pix = threshold_arcsec / 3600 / self.fits_image.pix_deg_scale
            
            N = len(self.image_filtered.cat)
            distance_matrix = np.zeros((N, N))
            for i in range(N) :
                for j in range(N) :
                    im_i = self.image_filtered.cat[i]
                    im_j = self.image_filtered.cat[j]
                    distance_matrix[i, j] = ( (im_i['x'] - im_j['x'])**2 + (im_i['y'] - im_j['y'])**2 )**0.5
                    #if i==j :
                    #    distance_matrix[i, j] = np.nan
            to_group_matrix = np.zeros((N, N))
            #to_groug_matrix[ np.logical_and(distance_matrix<threshold_pix, distance_matrix!=0.) ] = 1
            to_group_matrix[ distance_matrix<threshold_pix ] = 1.
            
            def find_related_groups(matrix):
                N = len(matrix)
                visited = [False] * N
                groups = []
                def dfs(node, group):
                    visited[node] = True
                    group.append(node)
                    for neighbor in range(N):
                        if matrix[node][neighbor] == 1 and not visited[neighbor]:
                            dfs(neighbor, group)
                for i in range(N):
                    if not visited[i]:
                        group = []
                        dfs(i, group)
                        groups.append(group)
                return groups
            groups = find_related_groups(to_group_matrix)
            
            to_remove = []
            for i, group in enumerate(groups) :
                x_mean = np.mean(self.image_filtered.cat['x'][group])
                y_mean = np.mean(self.image_filtered.cat['y'][group])
                self.image_filtered.cat[group[0]]['x'] = x_mean
                self.image_filtered.cat[group[0]]['y'] = y_mean
                to_remove += list(np.array(group)[1:])
            self.image_filtered.cat.remove_rows(to_remove)
    
    
    def start_extract_magnification_line(self) :
        self.doubleclick_magnification_marker = pg.ScatterPlotItem(size=12, symbol='x', brush='b', pen='b')
        self.source_magnification_marker = pg.ScatterPlotItem(size=8, symbol='o', brush='y', pen='y')
        self.fits_image.qt_image.addItem(self.doubleclick_magnification_marker)
        self.fits_image.qt_image.addItem(self.source_magnification_marker)
        self.magnification_markers_x = []
        self.magnification_markers_y = []
        self.magnification_source_markers_x = []
        self.magnification_source_markers_y = []
        self.magnification_temp_SkyCoords = []
        
        def mouse_clicked(evt):
            if evt.double():
                pos = evt.scenePos()
                if self.fits_image.qt_image.getView().sceneBoundingRect().contains(pos):
                    if len(self.magnification_temp_SkyCoords)==2 :
                        self.magnification_markers_x = []
                        self.magnification_markers_y = []
                        self.magnification_source_markers_x = []
                        self.magnification_source_markers_y = []
                        self.magnification_temp_SkyCoords = []
                        self.doubleclick_magnification_marker.setData([], [])
                        self.source_magnification_marker.setData([], [])
                    
                    mouse_point = self.fits_image.qt_image.getView().mapSceneToView(pos)
                    x, y_flipped = mouse_point.x(), mouse_point.y()
                    x, y = x, self.fits_image.image_data.shape[0] - y_flipped
                    ra, dec = self.fits_image.image_to_world(x, y)
                    
                    self.magnification_markers_x.append(x)
                    self.magnification_markers_y.append(self.fits_image.image_data.shape[0] - y)
                    self.magnification_temp_SkyCoords.append(SkyCoord(ra, dec, unit='deg'))
                    self.doubleclick_magnification_marker.setData(self.magnification_markers_x, self.magnification_markers_y)
                    
                    start = WCS.world_to_pixel(self.magnification_wcs, self.magnification_temp_SkyCoords[0])
                    self.magnification_line_start = (start[0]*1., start[1]*1.)
                    if len(self.magnification_temp_SkyCoords)==2 :
                        end = WCS.world_to_pixel(self.magnification_wcs, self.magnification_temp_SkyCoords[1])
                        self.magnification_line_end = (end[0]*1., end[1]*1.)
                        self.magnification_line = extract_line( self.magnification_line_start, self.magnification_line_end, self.magnification_map )
                        magnification_wcs = self.lt_magnification_maps[self.lt_z][1]
                        cd = magnification_wcs.wcs.cdelt[np.newaxis, :] * magnification_wcs.wcs.pc
                        deg_per_pix = np.sqrt((cd**2).sum(axis=0))[0]
                        self.magnification_line[0] = np.array(self.magnification_line[0]) * deg_per_pix * 3600 #x axis in arcsec
                        if self.magnification_line_ax==None :
                            print('Creating new magnification plot')
                            self.magnification_line_fig, self.magnification_line_ax = plt.subplots()
                            self.magnification_line_ax.set_yscale('log')
                        self.magnification_line_ax.clear()
                        self.magnification_line_ax.plot(self.magnification_line[0], np.abs(self.magnification_line[1]))
                        #self.magnification_line_fig.show()
            if evt.button()==PyQt5.QtCore.Qt.MiddleButton :
                pos = evt.scenePos()
                print(pos)
                mouse_point = self.fits_image.qt_image.getView().mapSceneToView(pos)
                x, y_flipped = mouse_point.x(), mouse_point.y()
                x, y = x, self.fits_image.image_data.shape[0] - y_flipped
                self.magnification_source_markers_x.append(x)
                self.magnification_source_markers_y.append(self.fits_image.image_data.shape[0] - y)
                self.source_magnification_marker.setData(self.magnification_source_markers_x, self.magnification_source_markers_y)
                
                distance = ( (self.magnification_markers_x[0] - x)**2 + (self.magnification_markers_y[0] - y_flipped)**2 )**0.5 * self.fits_image.pix_deg_scale*3600 #in arcsec
                print("self.magnification_line_ax", self.magnification_line_ax)
                self.magnification_line_ax.plot(np.full(10, distance), np.linspace(0, np.max(self.magnification_line[1]), 10), ls='--', c='grey')
                
                
        
        self._doubleclick_connection = self.fits_image.qt_image.scene.sigMouseClicked.connect(mouse_clicked)
        
        def keyPressEvent(event):
            #print('Hand selection stopped.')
            if event.key() == Qt.Key_Escape :
                if hasattr(self, 'doubleclick_magnification_marker'):
                    self.fits_image.qt_image.removeItem(self.doubleclick_magnification_marker)
                    del self.doubleclick_magnification_marker
                if hasattr(self, 'source_magnification_marker'):
                    self.fits_image.qt_image.removeItem(self.source_magnification_marker)
                    del self.source_magnification_marker
                if hasattr(self, '_doubleclick_connection'):
                    self.fits_image.qt_image.scene.sigMouseClicked.disconnect(self._doubleclick_connection)
                    del self._doubleclick_connection
                self.magnification_temp_SkyCoords = []
                self.fits_image.window.keyPressEvent = self._original_keyPressEvent
                print('Magnification line extraction stopped.')
        
        self._original_keyPressEvent = self.fits_image.window.keyPressEvent
        self.fits_image.window.keyPressEvent = keyPressEvent
    
    
    def send_to_source_plane(self) :
        for row in self.fits_image.imported_cat.cat :
            row['ra'], row['dec'] = self.transform_coords_radec(row['ra'], row['dec'])
            row['x'], row['y'] = self.fits_image.world_to_image(row['ra'], row['dec'])
    
    
    
    
    
    
    def start_simulate_image(self) :
        self.lm = lenstronomy_model(self.fits_image)
    
        
    
    def save_lenstronomy_model(self) :
        self.LENSTRONOMY_source_model = save_lenstronomy_model(os.path.join(self.model_dir, "lenstronomy_model.pkl"),
                                                                self.LENSTRONOMY_models,
                                                                self.LENSTRONOMY_result_kwargs,
                                                                self.LENSTRONOMY_center_world)
        
    def load_lenstronomy_model(self) :
        self.LENSTRONOMY_source_model_full = load_lenstronomy_model(os.path.join(self.model_dir, "lenstronomy_model.pkl"), self.LENSTRONOMY_center_world)
        
        self.models = copy.deepcopy(self.LENSTRONOMY_models)
        self.models['source_light_model_list'] = self.LENSTRONOMY_source_model_full['models']['source_light_model_list']
        
        self.kwargs = {'kwargs_lens': self.LENSTRONOMY_LensModel_kwargs}
        self.kwargs['kwargs_source'] = self.LENSTRONOMY_source_model_full['results']['kwargs_source']
        
        make_LENSTRONOMY_plot(self, self.models, self.kwargs)









