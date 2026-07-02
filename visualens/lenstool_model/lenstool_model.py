import os
import glob
import shutil
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from tqdm import tqdm
import warnings
import gc
import PyQt5
import pyqtgraph as pg
import pyqtgraph.exporters
from PyQt5.QtCore import Qt
import pickle
try:
    import lenstool
    _HAS_LENSTOOL = True
except ImportError :
    _HAS_LENSTOOL = False
    warnings.warn(
        "It looks like Lenstool is not installed. It is required to enable lens modeling functions.\n"
        "You can install it with conda:\n\n"
        "    conda install conda-forge::lenstool\n",
        ImportWarning,
        stacklevel=2)

###############################################################################
from ..utils.utils_astro.utils_general import world_to_relative, relative_to_world
from ..utils.utils_plots.plot_utils_general import plot_corner
from ..utils.utils_general.utils_general import extract_line
from ..utils.utils_general.sort_points import break_curves
from .simulate_image.simulate_image import image_simulator
from .im2source import start_im2source, stop_im2source
from .utils.operations import MakeFunctionFromMap
from .utils.utils_general import import_multiple_images, export_thumbnails, get_lenstool_file_path, import_lenstool_files
from .utils.param_extractors import read_potfile, make_best_file_from_bayes, make_param_latex_table, read_bayes_file, parse_lenstool_parameter_file, write_single_sample_best_file
from ..utils.utils_Qt.utils_general import transform_rectangle

from .utils.file_makers import best_files_maker, make_magnifications_and_curves                  # This import is problematic. The two functions run Lenstool
                                                                                                 # and are therefore dependent on my own install.




class lenstool_model :
    def __init__(self, model_path, fits_image, use_wrapper=None, compute_predictions=True, verbose=True) :
        self.fits_image = fits_image
        self.reference = None
        self.saturation = 1.
        self.lt = None
        self.z_lens = None
        self.mult = None
        self.source = None
        self.image = None
        self.image_filtered = None
        self.arclets = None
        self.potfile = None
        self.curve_plot = None
        self.curves = None
        self._safe_mode = False
        self._compute_predictions = compute_predictions
        self.verbose = verbose
        self.has_run = False
        
        # Get model directory
        self.model_dir = model_path if os.path.isdir(model_path) else os.path.dirname(model_path)
        
        # Get parameter file
        self.param_file_path = None
        self.best_file_path = None
        if os.path.isfile(model_path) :
            name = os.path.basename(model_path)
            if name.startswith('best') and name.endswith('.par') :
                self.best_file_path = model_path
                self.has_run = True
            else :
                self.param_file_path = model_path
        
        all_par_file_paths = glob.glob(os.path.join(self.model_dir, "*.par"))
        all_par_file_names = [ os.path.basename(file_path) for file_path in all_par_file_paths ]
        if self.param_file_path is None :
            for file_path in all_par_file_paths :
                if not os.path.basename(file_path).startswith('best') :
                    with open(file_path, 'r') as file :
                        for line in file :
                            stripped_line = line.strip()
                            if stripped_line and not stripped_line.startswith('#') :  # Skip empty and comment lines
                                if stripped_line.startswith('runmode') :
                                    self.param_file_path = file_path
                                    
        self.param = parse_lenstool_parameter_file(self.param_file_path) if self.param_file_path is not None else None
        
        # Get best and bayes file paths if they exist
        if self.best_file_path is None :
            if 'best.par' in all_par_file_names :
                self.best_file_path = os.path.join(self.model_dir, 'best.par')
                self.has_run = True
            else :
                l = [ file_name.startswith('best') and file_name.endswith('.par') for file_name in all_par_file_names ]
                if True in l :
                    i = np.where(l)[0][0]
                    self.has_run = True
                    self.best_file_path = all_par_file_paths[i]
        
        self.param_best = parse_lenstool_parameter_file(self.best_file_path) if self.best_file_path is not None else None
        self.bayes_file_path = os.path.join(self.model_dir, 'bayes.dat') if 'bayes.dat' in os.listdir(self.model_dir) else None
                
        if self.param is not None :
            # Get the multiple images and potfile galaxies from file names in parameter file 
            self.mult_path = os.path.join(self.model_dir, self.param['image']['multfile'][1]) if 'image' in self.param else None
            self.potfile_path = os.path.join(self.model_dir, self.param['potfile']['filein'][1]) if 'potfile' in self.param else None
            self.load_potfile(self.potfile_path)
        else :
            # Look for multiple image file and potfile without getting their file names from parameter file
            self.mult_path = get_lenstool_file_path(self.model_dir, 'mult')
            self.potfile_path = get_lenstool_file_path(self.model_dir, 'potfile')
            self.load_potfile(self.potfile_path)
            if self.potfile is not None :
                self._vprint("Potfile was found. Current instance now has potfile plotting capabilities.")
            if self.mult is None and self.potfile is None :
                self._vprint("No valid multiple image or potfile files found")
        
        if self.mult_path is not None :
            self._vprint("Multiple image file found.")
            self.families = []
            self.broad_families = []
            self.which = []
            import_multiple_images(self, self.mult_path, self.fits_image, units='pixel', filled_markers=True)
        
        # Get the lens redshift if unique
        if self.param is not None  or self.param_best is not None :
            param_dict = self.param if self.param is not None else self.param_best
            redshifts = []
            for name in param_dict :
                if 'potential' in name or 'potentiel' in name :
                    redshifts.append(param_dict[name]['z_lens'])
            if len(np.unique(redshifts))==1 :
                self.z_lens = redshifts[0]
                self._vprint(f"Lens redshift fixed at z={self.z_lens}")
            else :
                self._vprint("Several different redshift values found: " + str(np.unique(redshifts)))
                self.z_lens = np.max(redshifts)
                self._vprint("Setting lenstool_model.z_lens to the furthest lens' redshift: " + str(self.z_lens))
        
        # Checks if Lenstool files were found and if so ask user if they want to use Lenstool's wrapper and its capabilities
        # If so, moves to the model's directory (required by the Lenstool wrapper)
        # Check which file to use: best.par, best_TEMP.par from bayes.dat, or parameter file  
        self._FileToUse = None
        if self.has_run :
            # Use best file
            self._FileToUse = os.path.basename(self.best_file_path)
            
            if self.bayes_file_path is not None :
                self.samples_df_full = read_bayes_file(self.bayes_file_path, z=self.z_lens) if self.bayes_file_path is not None else None
                # Extract the numeric columns (skip non-numeric, zero-range etc.)
                self.samples_df = self.samples_df_full.copy()
                for col in self.samples_df.columns :
                    if col=='Chi2' or col=='Nsample' or col=='ln(Lhood)':
                        del self.samples_df[col]
                    
        elif self.param_file_path is not None :
            if 'bayes.dat' in os.listdir(self.model_dir) :
                # Option to create temporary best file from bayes if optimization hasn't finished
                yesno = input('bayes.dat file found, create best_TEMP.par from bayes file to be used in Lenstool wrapper? ([Y]/n)')
                if yesno.lower() in ['y', ''] :
                    make_best_file_from_bayes(self.param_file_path)
                    self.best_file_path = os.path.join(self.model_dir, 'best_TEMP.par')
                    self._FileToUse = os.path.basename(self.best_file_path)
                else :
                    self._vprint("Loading parameter file only. Limited capabilities available.")
                    self._FileToUse = os.path.basename(self.param_file_path)
            else :
                self._vprint("Only parameter file was found (looks like lenstool optimization hasn't run yet). Limited capabilities available.")
                self._FileToUse = os.path.basename(self.param_file_path)
        else :
            self._vprint("No best file or parameter file found."
                         "Make sure best file has name 'best.par' and parameter file has extension '.par'")
            import_lenstool_files(self)
        
        
        # Check if user wants to use Lenstool wrapper or not
        if self._FileToUse is not None :
            if not _HAS_LENSTOOL :
                self._vprint("\n------\n"
                             "It looks like Lenstool is not installed. Only limited lens modeling functions will be available.\n"
                             "You can install Lenstool with conda:\n\n"
                             "    conda install conda-forge::lenstool\n"
                             "------\n")
                yesno = 'n'
            elif not self.model_dir.endswith('_safe/') :
                self._vprint(f"\nA valid lenstool file was found: {self._FileToUse}")
                if use_wrapper is None :
                    self._vprint("Lens modeling tools rely on Lenstool and may modify some of the output files present in your model's directory (image.all, source.dat for example).\n")
                    yesno = input(f"Continue? Or make a copy of your model's current directory? ({os.path.basename(self.model_dir[:-1])+'_safe'})\n\nYes/copy/no [Y]/c/n\n-- to skip this message use .import_lenstool(dir, use_wrapper=True) --\n")
                elif use_wrapper :
                    yesno = 'y'
                else :
                    yesno = 'n'
            else :
                yesno = 'c'
            
            # Load files according to files found and user input
            if yesno.lower()=='n' or yesno.lower()=='no' :
                self._vprint("Lens model loaded with limited capabilities.")
                self.lt = None
                #Functionalities without Lenstool's wrapper
                if param_dict is not None :
                    self.reference = tuple(param_dict['runmode']['reference'][1:]) if 'runmode' in param_dict else None
                import_lenstool_files(self)
            else :
                if yesno.lower()=='c' or yesno.lower()=='copy' :
                    self.SafeMode()
                    self._vprint(f"Loading {self._FileToUse}")
                    self.lt = lenstool.Lenstool( self._FileToUse )
                else :
                    self._vprint('--------------------')
                    self._vprint("Moving to " + self.model_dir)
                    self._vprint('--------------------')
                    os.chdir(self.model_dir)
                    self._vprint(f"Loading {self._FileToUse}")
                    self.lt = lenstool.Lenstool( self._FileToUse )
                self.reference = (self.lt.M.ref_ra, self.lt.M.ref_dec)
        
        
        # some useful initializations
        self.curve_plot = None
        self.magnification_res = 1000
        self.magnification_line_ax = None
        self.previous_state_current_ROI = None
        self.LENSTRONOMY_fixed_source_kwargs = []
        
        # load maps if some have been saved to save compute time
        self.load_saved_maps()
        
        # compute sources and images from multiple images
        if self.mult is not None and self.lt is not None and compute_predictions :
            self.add_lensing_columns(cat=self.mult.cat)
            self.compute_sources_and_images()
            
            
    def compute_sources_and_images(self, ngrid=256, restrict_to_mult_field=False) :
        if self.source is not None :
            self.source.clear()
        if self.image is not None :
            self.image.clear()
        
        self._vprint("Computing sources")
        source = self.mult.cat.copy()
        source.rename_column('ra', 'ra_image')
        source.rename_column('dec', 'dec_image')
        source.rename_column('ra_source', 'ra')
        source.rename_column('dec_source', 'dec')
        
        import_multiple_images(self, source, self.fits_image, AttrName='source', units='pixel', filled_markers=True)
        
        # Format source catalog in the Lenstool format
        source.rename_column('id', 'n')
        source['x'], source['y'] = source['ra'], source['dec'] #self.world_to_relative(source['ra'], source['dec'])
        
        for colname in source.colnames :
            if colname not in ['n','x','y','a','b','z','theta','mag'] :
                source.remove_column(colname)
        
        to_remove = []
        for i, s in enumerate(source) :
            if np.isnan(s['z']) or s['z'] <= self.z_lens :
                to_remove.append(i)

        source.remove_rows(to_remove)
        self._vprint('done')
        
        self._vprint("Computing images")
        _initial_ngrid_value = self.lt.G.ngrid
        self.lt.set_grid(ngrid, 0)
        
        if restrict_to_mult_field :
            xr, yr = self.world_to_relative(self.mult.cat['ra'], self.mult.cat['dec'])
            _initial_field = self.lt.get_field([])
            self.lt.set_field([np.min(xr)-1., np.max(xr)+1., np.min(yr)-1., np.max(yr)+1.])
        
        self.lt.set_sources(source)
        self.lt.e_lensing()
        image = self.lt.get_images()
        
        image['ra'], image['dec'] = self.relative_to_world(image['x'], image['y'])
        image['x'], image['y'] = self.fits_image.world_to_image(image['ra'], image['dec'])
        
        image.rename_column('n', 'id')
        
        cols_to_add = [[], [], []]
        
        for row in image :
            i = np.where(self.mult.cat['id']==row['id'])[0][0]
            for j, colname in enumerate(['family','broad_family', 'confidence']) :
                cols_to_add[j].append(self.mult.cat[colname][i])
        for j, colname in enumerate(['family','broad_family', 'confidence']) :
            image.add_column(cols_to_add[j], name=colname)
        
        import_multiple_images(self, image, self.fits_image, AttrName='image', units='pixel', filled_markers=True)
        import_multiple_images(self, image, self.fits_image, AttrName='image_filtered', units='pixel', filled_markers=True)
        self.filter_image()
        
        self.lt.set_grid(_initial_ngrid_value, 0)
        if restrict_to_mult_field :
            self.lt.set_field(_initial_field)
        self._vprint('done')
    
    def _vprint(self, *args, **kwargs) :
        if self.verbose :
            print(*args, **kwargs)

    def world_to_relative(self, ra, dec) :
        return world_to_relative(ra, dec, self.reference)
    
    def relative_to_world(self, xr, yr) :
        return relative_to_world(xr, yr, self.reference)
    
    def SafeMode(self) :
        if self._safe_mode :
            self._vprint("Already in safe directory")
        elif self.model_dir.endswith('_safe/') :
            self._vprint("Safe directory already selected, moving to it")
            os.chdir(self.model_dir)
            self._safe_mode = True
        else :
            safe_dir = self.model_dir.rstrip('/') + '_safe/'
            if os.path.exists(safe_dir) :
                self._vprint("Safe directory already exists, moving to it")
                self.__init__(safe_dir, self.fits_image)
                os.chdir(safe_dir)
                self._safe_mode = True
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
                self._safe_mode = True
        self._vprint("Now in " + os.getcwd())
    
    def load_potfile(self, path) :
        if self.potfile is not None :
            self.potfile.clear()
        if path is not None :        
            potfile_Table = read_potfile(path)
            self.potfile = self.fits_image.make_catalog(potfile_Table, color=[1.,0.,0.], units='arcsec', verbose=self.verbose)
    
    def load_saved_maps(self) :
        self.convergence_maps_path = os.path.join(self.model_dir, 'convergence_maps.pkl')
        if os.path.exists(self.convergence_maps_path) :
            self._vprint('convergence_maps.pkl found')
            with open(self.convergence_maps_path, 'rb') as f :
                self.convergence_maps = pickle.load(f)
        else :
            self.convergence_maps = {}
        #######################################################################
        self.dpl_maps_path = os.path.join(self.model_dir, 'dpl_maps.pkl')
        if os.path.exists(self.dpl_maps_path) :
            self._vprint('dpl_maps.pkl found')
            with open(self.dpl_maps_path, 'rb') as f :
                self.dpl_maps = pickle.load(f)
        else :
            self.dpl_maps = {}
        #######################################################################
        self.lt_curves_path = os.path.join(self.model_dir, 'lt_curves.pkl')
        if os.path.exists(self.lt_curves_path) :
            self._vprint('lt_curves.pkl found')
            with open(self.lt_curves_path, 'rb') as f :
                self.lt_curves = pickle.load(f)
        else :
            self.lt_curves = {}
        #######################################################################
        self.lt_magnification_maps_path = os.path.join(self.model_dir, 'lt_magnification_maps.pkl')
        if os.path.exists(self.lt_magnification_maps_path) :
            self._vprint('lt_magnification_maps.pkl found')
            with open(self.lt_magnification_maps_path, 'rb') as f :
                self.lt_magnification_maps = pickle.load(f)
        else :
            self.lt_magnification_maps = {}
        #######################################################################
        self.lt_caustics_path = os.path.join(self.model_dir, 'lt_caustics.pkl')
        if os.path.exists(self.lt_caustics_path) :
            self._vprint('lt_caustics.pkl found')
            with open(self.lt_caustics_path, 'rb') as f :
                self.lt_caustics = pickle.load(f)
        else :
            self.lt_caustics = {}
    
    #def select_multiple_images(self) :
    #    print('function to be implemented in the future')
    
    def plot(self, which=None) :
        if which is not None :
            self.set_which(which)
        if self.mult is not None :
            self.mult.plot(marker='o', filled_markers=False, scale=1.25)#size=1.5, linewidth=2, filled_markers=False)
            self.mult.plot_column('id')
        if self.image is not None :
            #self.image.plot(marker='x', filled_markers=True, scale=1)
            self.image.saturation = 1.
            #self.image.plot_column('id')
        if self.image_filtered is not None :
            self.image_filtered.plot(marker='x', filled_markers=True, scale=0.5)
            #self.image_filtered.saturation = 1.
            self.image.plot_column('id')
        if self.curves is not None :
            self.curves.plot()
    
    def clear(self) :
        if self.mult is not None :
            self.mult.clear()
        if self.image is not None :
            self.image.clear()
        if self.image_filtered is not None :
            self.image_filtered.clear()
        if self.arclets is not None :
            self.arclets.clear()
        if self.curves is not None :
            self.curves.clear()
        if self.curve_plot is not None :
            self.fits_image.qt_image.removeItem(self.curve_plot)
            
    def set_which(self, *names) :
        if names[0]=='all' :
            self.which = self.broad_families if type(self.broad_families) is list else self.broad_families.tolist()
        elif isinstance(names[0], list) :
            self.which = names[0]
        else :
            self.which = list(names)
        self._vprint("Images to plot are now ", self.which)
        #self.clear()
        #self.plot()
        
    def __make_files(self) : #You can probably remove this function now that we are using Lenstool's wrapper
        best_files_maker(self.model_dir)
        make_magnifications_and_curves(self.model_dir)
    
    
    def export_thumbnails(self, group_images=True, square_thumbnails=True, square_size=150, margin=50, distance=200, export_dir=None, boost=True, make_broad_view=True, broad_view_params=None) :
        export_thumbnails(self.mult, group_images=group_images, square_thumbnails=square_thumbnails, square_size=square_size, margin=margin, \
                          distance=distance, export_dir=export_dir, boost=boost, make_broad_view=make_broad_view, broad_view_params=broad_view_params)
    
    #def make_webpage(self) :
    #    print('function to be implemented in the future')
        
    def make_latex(self) :
        latex_str = make_param_latex_table(self.param_file_path, convert_to_kpc=True, z=self.z_lens)
        return latex_str
    
    
    def set_lt_z(self, z, color=[255,100,255], recompute=False) :
        #if self.curve_plot is not None :
        #    self.fits_image.qt_image.removeItem(self.curve_plot)
        self.lt_z = z
        self._vprint(self.best_file_path)
        self._vprint(os.getcwd())
        #self.lt.set_grid(50, 0)

        ######## Curves ########
        if z not in self.lt_curves.keys() or recompute :
            self.compute_lt_curve(z)
        self.lt_curve_coords_relative = self.lt_curves[z]
        self.lt_caustic_coords_relative = self.lt_caustics[z]
        self._curves_add_all_coords()
        self.plot_lt_curve(color=color)
        
        ######## Magnification ########
        if z not in self.lt_magnification_maps.keys() or recompute :
            self._vprint('Computing magnification map (can take a little while)...')
            self.lt_magnification_maps[z] = self.lt.g_ampli(1, self.magnification_res, self.lt_z)
            self._vprint('done')
            with open(self.lt_magnification_maps_path, 'wb') as f:
                pickle.dump(self.lt_magnification_maps, f)
        self.magnification_map, self.magnification_wcs = self.lt_magnification_maps[z]            
        
        ######## Convergence ########
        if self.lt_z not in self.convergence_maps.keys() or recompute :
            self.compute_lt_convergence()
        else :
            self.convergence_map, self.convergence_map_wcs = self.convergence_maps[self.lt_z]
        
        ######## Displacement maps ########
        if self.lt_z not in self.dpl_maps.keys() or recompute :
            self.compute_lt_dpl()
        else :
            self.dx_map, self.dy_map, self.dmap_wcs = self.dpl_maps[self.lt_z]
        
        mmap, wcs = self.lt_magnification_maps[z]
        self.get_magnification = MakeFunctionFromMap(mmap, wcs)
        
    def compute_lt_convergence(self, z=None, npix=2000) :
        z = self.lt_z if z is None else z
        self._vprint('Computing convergence map (can take a little while)...')
        self.convergence_maps[z] = self.lt.g_mass(1, npix, self.z_lens, z)
        self.convergence_map, self.convergence_map_wcs = self.convergence_maps[z]
        self._vprint('done')
        with open(self.convergence_maps_path, 'wb') as f:
            pickle.dump(self.convergence_maps, f)
        
    def compute_lt_dpl(self, z=None, npix=2000) :
        z = self.lt_z if z is None else z
        self._vprint('Computing displacement maps (can take a little while)...')
        self.dpl_maps[z] = self.lt.g_dpl(npix, z)
        self.dx_map, self.dy_map, self.dmap_wcs = self.dpl_maps[z]
        self._vprint('done')
        with open(self.dpl_maps_path, 'wb') as f:
            pickle.dump(self.dpl_maps, f)
    
    
    def start_im2source(self) :
        start_im2source(self)
        
    def stop_im2source(self) :
        stop_im2source(self)
    
    def start_magnification(self) :
        return None
    
    def add_lensing_columns(self, cat=None, which_cat='imported_cat', index=None, z_source=None, overwrite=None) :
        if cat is None :
            if index is not None :
                cat = self.fits_image.imported_cat_list[index].cat
            else :
                cat = getattr(self.fits_image, which_cat, None).cat
        
        lensing_columns = ['magnification', 'convergence', 'shear', 'gamma1', 'gamma2', 'time', 'tangential_magnification', 'radial_magnification', 'ra_source', 'dec_source']
        check = False
        existing_cols = []
        for name in lensing_columns :
            if name in cat.colnames :
                check = True
                existing_cols.append(name)
        yesno = 'y' if overwrite is None else 'y' if overwrite else 'n'
        if check and overwrite is None :
            yesno = input(f"Following lensing columns already exist: {existing_cols}. Overwrite? [Y]/n")
        
        if yesno.lower() in ['y', ''] :
            mu_col = np.full(len(cat), np.nan)
            gamma_col = np.full(len(cat), np.nan)
            kappa_col = np.full(len(cat), np.nan)
            tmu_col = np.full(len(cat), np.nan)
            rmu_col = np.full(len(cat), np.nan)
            time_col = np.full(len(cat), np.nan)
            gamma1_col = np.full(len(cat), np.nan)
            gamma2_col = np.full(len(cat), np.nan)
            
            ra_source_col = np.full(len(cat), np.nan)
            dec_source_col = np.full(len(cat), np.nan)
            
            initial_field = self.lt.get_field([])
            
            z_colname = None
            for name in ['z', 'z_spec', 'zspec', 'z_phot', 'zphot'] :
                if name in cat.colnames :
                    z_colname = name
                    break
            if z_colname is None and z_source is None :
                self._vprint('Redshift column not found in catalog. Using current source redshift = ' + str(self.lt_z))
                z_source = self.lt_z
                
            for i in tqdm(range(len(cat))) :
                z = cat[z_colname][i] if z_source is None else z_source
                if z>self.z_lens : # Also returns False if z is np.nan
                    #print(str(cat['id'][i]) + ': computing lensing maps at redshift ' + str(z))
                    ra, dec = cat['ra'][i], cat['dec'][i]
                    world_coord = SkyCoord(ra, dec, unit='deg')
                    xr, yr = self.world_to_relative(ra, dec)
                    #delta = 1.
                    delta = self.fits_image.pix_deg_scale * 3600 / 2
                    self.lt.set_field([xr-delta, xr+delta, yr-delta, yr+delta])
                    
                    #npix = 11
                    npix = 2
                    ampli, wcs = self.lt.g_ampli(1, npix, z)
                    mu_col[i] = np.mean(ampli)
                    
                    #if False : # this crashes Lenstool if 
                    try :
                        kappa, wcs = self.lt.g_mass(1, npix, self.z_lens, z)
                    except ZeroDivisionError :
                        self.lt.set_field([xr-1., xr+1., yr-1., yr+1.])
                        kappa, wcs = self.lt.g_mass(1, 3, self.z_lens, z)
                        self.lt.set_field([xr-delta, xr+delta, yr-delta, yr+delta])
                    
                    kappa_col[i] = np.mean(kappa)
                    #gamma_col[i] = ( (1-kappa_col[i])**2 - 1/mu_col[i] )**0.5
                    
                    dx, dy, wcs = self.lt.g_dpl(npix, z)
                    dx = np.mean(dx)
                    dy = np.mean(dy)                    
                    ra_source_col[i] = ra + dx /3600 /np.cos( dec * np.pi/180 )
                    dec_source_col[i] = dec - dy /3600

                    time, wcs = self.lt.g_time(1, npix, z)
                    time_col[i] = np.mean(time)

                    shear, wcs = self.lt.g_shear(1, npix, z)
                    gamma_col[i] = np.mean(shear)

                    shear1, wcs = self.lt.g_shear(3, npix, z)
                    gamma1_col[i] = np.mean(shear1)

                    shear2, wcs = self.lt.g_shear(4, npix, z)
                    gamma2_col[i] = np.mean(shear2)

                #else :
                    #print(str(cat['id'][i]) + ': redshift ' + str(z) + ' lower than lens redshift --> NaN')
                
                tmu_col[i] = 1 / (1 - kappa_col[i] - gamma_col[i] )
                rmu_col[i] = 1 / (1 - kappa_col[i] + gamma_col[i] )
            
            columns_to_add = [mu_col, kappa_col, gamma_col, gamma1_col, gamma2_col, time_col, tmu_col, rmu_col, ra_source_col, dec_source_col]
            for i, name in enumerate(lensing_columns) :
                if name in cat.colnames :
                    cat.replace_column(name, columns_to_add[i])
                else :
                    cat.add_column(columns_to_add[i], name=name)
            self.lt.set_field(initial_field)
    
    def compute_lt_curve(self, z=None, limitHigh=0.5, limitLow=0.1) :
        if z==None :
            z = self.lt_z
        self._vprint('Computing critical curve (can take a little while)...')
        self.lt_curve = self.lt.criticnew(zs=z, limitHigh=limitHigh, limitLow=limitLow) #limitHigh=1., limitLow=0.05
        
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
        self._vprint('done')
        
        self.lt_curves[z] = self.lt_curve_coords_relative
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
        self._vprint('done')
        
        self.lt_caustics[z] = self.lt_caustic_coords_relative
        with open(self.lt_caustics_path, 'wb') as f:
            pickle.dump(self.lt_caustics, f)
        
        self.plot_lt_curve()
    
    def _curves_add_all_coords(self) :
        lt_curve_xr, lt_curve_yr = self.lt_curve_coords_relative
        lt_curve_ra, lt_curve_dec = self.relative_to_world(lt_curve_xr, lt_curve_yr)
        lt_curve_x, lt_curve_y = self.fits_image.world_to_image(lt_curve_ra, lt_curve_dec)
        
        self.lt_curve_coords_world = [lt_curve_ra, lt_curve_dec]
        self.lt_curve_coords_image = [lt_curve_x, self.fits_image.image_data.shape[0] - lt_curve_y]
        
        lt_caustic_xr, lt_caustic_yr = self.lt_caustic_coords_relative
        lt_caustic_ra, lt_caustic_dec = self.relative_to_world(lt_caustic_xr, lt_caustic_yr)
        lt_caustic_x, lt_caustic_y = self.fits_image.world_to_image(lt_caustic_ra, lt_caustic_dec)
        
        self.lt_caustic_coords_world = [lt_caustic_ra, lt_caustic_dec]
        self.lt_caustic_coords_image = [lt_caustic_x, self.fits_image.image_data.shape[0] - lt_caustic_y]
    
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
        plot_corner(self.samples_df_param_only)
        corr_matrix = self.samples_df_param_only.corr()
        self.fig_cov, self.ax_cov = plt.subplots()
        cax = self.ax_cov.imshow(corr_matrix, cmap='PuOr')
        cbar = self.fig_cov.colorbar(cax, ax=self.ax_cov)
        self.ax_cov.set_xticks(np.arange(len(corr_matrix.columns)))
        self.ax_cov.set_yticks(np.arange(len(corr_matrix.index)))
        self.ax_cov.set_xticklabels(corr_matrix.columns, rotation=45, ha='right')
        self.ax_cov.set_yticklabels(corr_matrix.index)
    
    def filter_image(self, threshold_arcsec=0.1) :
        if self.mult is not None :
            threshold_pix = threshold_arcsec / 3600 / self.fits_image.pix_deg_scale
            to_remove = []
            for i, image in enumerate(self.image_filtered.cat) :
                ref_mask = self.mult.cat['id']==image['id']
                #if not ref_mask.any() :
                #    d = 0
                #else :
                ref = self.mult.cat[ np.where(ref_mask)[0][0] ]
                d = ( (ref['x'] - image['x'])**2 + (ref['y'] - image['y'])**2 )**0.5
                if d<threshold_pix :
                    to_remove.append(i)
            self.image_filtered.cat.remove_rows(to_remove)
        
        ### Grouping of images with similar positions, not used anymore ###
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
                        magnification_wcs = self.magnification_wcs
                        cd = magnification_wcs.wcs.cdelt[np.newaxis, :] * magnification_wcs.wcs.pc
                        deg_per_pix = np.sqrt((cd**2).sum(axis=0))[0]
                        self.magnification_line[0] = np.array(self.magnification_line[0]) * deg_per_pix * 3600 #x axis in arcsec
                        if True : #self.magnification_line_ax==None :
                            self.magnification_line_distances = []
                            
                            plt.close()
                            self._vprint('Creating new magnification plot')
                            self.magnification_line_fig, self.magnification_line_ax = plt.subplots()
                            self.magnification_line_ax.set_yscale('log')
                            
                            #self.magnification_line_ax_xlim = self.magnification_line_ax.get_xlim()
                            #self.magnification_line_ax_ylim = self.magnification_line_ax.get_ylim()
                            
                        self.magnification_line_ax.clear()
                        self.magnification_line_ax.plot(self.magnification_line[0], np.abs(self.magnification_line[1]))
                        #self.magnification_line_fig.show()
            elif evt.button()==PyQt5.QtCore.Qt.MiddleButton :
                pos = evt.scenePos()
                self._vprint(pos)
                mouse_point = self.fits_image.qt_image.getView().mapSceneToView(pos)
                x, y_flipped = mouse_point.x(), mouse_point.y()
                x, y = x, self.fits_image.image_data.shape[0] - y_flipped
                self.magnification_source_markers_x.append(x)
                self.magnification_source_markers_y.append(self.fits_image.image_data.shape[0] - y)
                self.source_magnification_marker.setData(self.magnification_source_markers_x, self.magnification_source_markers_y)
                
                distance = ( (self.magnification_markers_x[0] - x)**2 + (self.magnification_markers_y[0] - y_flipped)**2 )**0.5 * self.fits_image.pix_deg_scale*3600 #in arcsec
                self.magnification_line_distances.append(distance)
                
                if True : # remove this when plot update available
                    plt.close()
                    self._vprint('Creating new magnification plot')
                    self.magnification_line_fig, self.magnification_line_ax = plt.subplots()
                    self.magnification_line_ax.set_yscale('log')
                    self.magnification_line_ax.plot(self.magnification_line[0], np.abs(self.magnification_line[1]))
                    
                    #self.magnification_line_ax_xlim = self.magnification_line_ax.get_xlim()
                    #self.magnification_line_ax_ylim = self.magnification_line_ax.get_ylim()
                    self.magnification_line_ax.set_xlim(self.magnification_line_ax.get_xlim())
                    self.magnification_line_ax.set_ylim(self.magnification_line_ax.get_ylim())
                
                for distance in self.magnification_line_distances :
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
                self.fits_image.qt_image.keyPressEvent = self._original_keyPressEvent
                self._vprint('Magnification line extraction stopped.')
        
        self._original_keyPressEvent = self.fits_image.qt_image.keyPressEvent
        self.fits_image.qt_image.keyPressEvent = keyPressEvent
    
    
    def send_to_source_plane(self) :
        for row in self.fits_image.imported_cat.cat :
            row['ra'], row['dec'] = self.transform_coords_radec(row['ra'], row['dec'])
            row['x'], row['y'] = self.fits_image.world_to_image(row['ra'], row['dec'])
    
    
    def start_simulate_image(self, which_filter=None, throttle_mode=0) :
        self.imsim = image_simulator(self.fits_image, which_filter=which_filter, throttle_mode=throttle_mode)
        
    
    def compute_mass_map(self, z=None, npix=1000) :
        z = self.lt_z if z is None else z
        self.mass_map, self.mass_map_wcs = self.lt.g_mass(1, npix, self.z_lens, z)
        fig, ax = plt.subplots()
        ax.imshow(np.arctan(self.mass_map), origin='lower')
    
    def compute_magnification_map(self, z=None, npix=1000) :
        z = self.lt_z if z is None else z
        self.magnification_map, self.magnification_wcs = self.lt.g_ampli(1, npix, z)
        fig, ax = plt.subplots()
        ax.imshow(np.arctan(np.abs(self.magnification_map)), origin='lower')
    
    def compute_time_map(self, z=None, npix=1000) :
        z = self.lt_z if z is None else z
        self.time_map, self.time_wcs = self.lt.g_time(1, npix, z)
        fig, ax = plt.subplots()
        ax.imshow(self.time_map, origin='lower')
        
    def compute_shear_map(self, z=None, npix=1000, which='gamma') :
        z = self.lt_z if z is None else z
        ishear = 1 if which=='gamma' else 3 if which=='gamma1' else 4 if which=='gamma2' else None
        if ishear is None :
            raise ValueError('Invalid shear type: ' + which)
        self.shear_map, self.shear_wcs = self.lt.g_shear(ishear, npix, z)
        fig, ax = plt.subplots()
        ax.imshow(self.shear_map, origin='lower')
        
    def compute_mass(self, ra, dec, z=None) :
        initial_field = self.lt.get_field([])
        z = self.lt_z if z is None else z
        xr, yr = self.world_to_relative(ra, dec)
        delta = 0.0001 # 0.1 milli arcsec
        self.lt.set_field([xr-delta, xr+delta, yr-delta, yr+delta])
        npix = 2
        kappa, _ = self.lt.g_mass(1, npix, self.z_lens, z)
        self.lt.set_field(initial_field)
        return np.mean(kappa)
    
    def compute_magnification(self, ra, dec, z=None) :
        initial_field = self.lt.get_field([])
        z = self.lt_z if z is None else z
        xr, yr = self.world_to_relative(ra, dec)
        delta = 0.0001 # 0.1 milli arcsec
        self.lt.set_field([xr-delta, xr+delta, yr-delta, yr+delta])
        npix = 2
        ampli, _ = self.lt.g_ampli(1, npix, z)
        self.lt.set_field(initial_field)
        return np.mean(ampli)
    
    def compute_time(self, ra, dec, z=None) :
        initial_field = self.lt.get_field([])
        z = self.lt_z if z is None else z
        xr, yr = self.world_to_relative(ra, dec)
        delta = 0.0001 # 0.1 milli arcsec
        self.lt.set_field([xr-delta, xr+delta, yr-delta, yr+delta])
        npix = 2
        time, _ = self.lt.g_time(1, npix, z)
        self.lt.set_field(initial_field)
        return np.mean(time)
    
    def compute_shear(self, ra, dec, z=None, which='gamma') :
        initial_field = self.lt.get_field([])
        z = self.lt_z if z is None else z
        ishear = 1 if which=='gamma' else 3 if which=='gamma1' else 4 if which=='gamma2' else None
        xr, yr = self.world_to_relative(ra, dec)
        delta = 0.0001 # 0.1 milli arcsec
        self.lt.set_field([xr-delta, xr+delta, yr-delta, yr+delta])
        npix = 2
        shear, _ = self.lt.g_shear(ishear, npix, z)
        self.lt.set_field(initial_field)
        return np.mean(shear)
        
    def compute_displacement(self, ra, dec, z=None) :
        initial_field = self.lt.get_field([])
        z = self.lt_z if z is None else z
        xr, yr = self.world_to_relative(ra, dec)
        delta = 0.0001 # 0.1 milli arcsec
        self.lt.set_field([xr-delta, xr+delta, yr-delta, yr+delta])
        npix = 2
        dx, dy, _ = self.lt.g_dpl(npix, z)
        self.lt.set_field(initial_field)
        return np.mean(dx), np.mean(dy)
    
    
    def set_field(self) :
        # Careful, this function only works when the ROI is flat and image in world frame
        
        self.ROI = self.fits_image.image_widget.current_ROI
        x0 = self.ROI.getState()['pos'][0]
        y0 = self.ROI.getState()['pos'][1]
        a = self.ROI.getState()['size'][0]
        b = self.ROI.getState()['size'][1]
        angle = self.ROI.getState()['angle'] *np.pi/180
        
        x0, y0, a, b, angle = transform_rectangle(x0, y0, a, b, angle) #x0, y0 at the top left
        size_y = self.fits_image.image_data.shape[0]
        y0 = size_y-y0 #Counting pixels from bottom instead of top
        
        ra_left, dec_top = self.fits_image.image_to_world(x0, y0)
        ra_right, dec_bottom = self.fits_image.image_to_world(x0 + a, y0 - b)
        
        xr_left, yr_top = self.world_to_relative(ra_left, dec_top)
        xr_right, yr_bottom = self.world_to_relative(ra_right, dec_bottom)
        
        self._previous_field = self.lt.get_field([])
        self.lt.set_field([xr_left, xr_right, yr_bottom, yr_top])
        
    
    def compute_uncertainties(self, nsamples=None, recompute_samples=False) :
        self.samples_dir = os.path.join(self.model_dir, 'samples')
        if not os.path.exists(self.samples_dir) :
            self._vprint('Creating samples directory')
            os.makedirs(self.samples_dir)
        
        lensing_columns = [ 'magnification', 
                            'convergence', 
                            'shear', 
                            'gamma1', 
                            'gamma2', 
                            'time', 
                            'tangential_magnification', 
                            'radial_magnification', 
                            'ra_source', 
                            'dec_source',
                            'z_opt' ]
        
        samples_dict_path = os.path.join(self.samples_dir, 'samples_dict.pkl')
        if os.path.exists(samples_dict_path) and not recompute_samples :
            self._vprint('Loading samples dictionary from ' + samples_dict_path)
            with open(samples_dict_path, 'rb') as f :
                self.samples_dict = pickle.load(f)

        else :
            if os.path.exists(samples_dict_path) :
                self._vprint('Renaming previous samples dictionary to ' + samples_dict_path.replace('.pkl', '_previous.pkl'))
                shutil.move(samples_dict_path, samples_dict_path.replace('.pkl', '_previous.pkl'))

            if os.path.exists(os.path.join(self.samples_dir, os.path.basename(self.mult_path))) :
                self._vprint('Overwriting mult file in samples directory')
            else :
                self._vprint('Copying mult file to samples directory')
            shutil.copy2(self.mult_path, self.samples_dir)
            
            nsamples = len(self.samples_df.index) if nsamples is None else nsamples
            
            self.samples_dict = {}
            for imID in self.mult.cat['id'] :
                self.samples_dict[imID] = {col: np.full(nsamples, np.nan) for col in lensing_columns}
                
            for i in tqdm(range(nsamples)) :
                sample_file_path = write_single_sample_best_file(self, i)
                sample_lt = import_lenstool(sample_file_path, self.fits_image, use_wrapper=True, compute_predictions=False, verbose=False)

                """ Sample the redshift as well for those that were optimized """
                for im in sample_lt.mult.cat :
                    #if not np.isnan(im['z_opt']) :
                    broad_family_members_mask = np.array([ member['broad_family']==im['broad_family'] for member in sample_lt.mult.cat ])
                    for member in sample_lt.mult.cat[broad_family_members_mask] :
                        for col in self.samples_df.columns :
                            if member['id'] == col[len('Redshift of '):] :
                                self._vprint('Sampling redshift of ' + im['id'] + ': ' + str(im['z']) + ' --> ' + str(self.samples_df[col][i]))
                                im['z'] = self.samples_df[col][i]
                                self.samples_dict[im['id']]['z_opt'][i] = self.samples_df[col][i]
                                
                """ Compute magnification etc. for the sample """
                sample_lt.add_lensing_columns(cat=sample_lt.mult.cat)

                for im in sample_lt.mult.cat :
                    for col in lensing_columns :
                        if col != 'z_opt' :
                            self.samples_dict[im['id']][col][i] = im[col]
                    
                del sample_lt
                gc.collect()
                os.remove(sample_file_path)

                if i % 10 == 0 :
                    pickle.dump(self.samples_dict, open(samples_dict_path, 'wb'))
            
            pickle.dump(self.samples_dict, open(samples_dict_path, 'wb'))
        
        for col in lensing_columns :
            col_16_percentile = np.full(len(self.mult.cat), np.nan)
            col_84_percentile = np.full(len(self.mult.cat), np.nan)
            col_50_percentile = np.full(len(self.mult.cat), np.nan)

            name = f'{col}_16_percentile'
            if name in self.mult.cat.colnames :
                self.mult.cat.replace_column(name, col_16_percentile)
            else :
                self.mult.cat.add_column(col_16_percentile, name=name)

            name = f'{col}_84_percentile'
            if name in self.mult.cat.colnames :
                self.mult.cat.replace_column(name, col_84_percentile)
            else :
                self.mult.cat.add_column(col_84_percentile, name=name)

            name = f'{col}_50_percentile'
            if name in self.mult.cat.colnames :
                self.mult.cat.replace_column(name, col_50_percentile)
            else :
                self.mult.cat.add_column(col_50_percentile, name=name)

            for im in self.mult.cat :
                im[f'{col}_16_percentile'] = np.percentile(self.samples_dict[im['id']][col], 16)
                im[f'{col}_84_percentile'] = np.percentile(self.samples_dict[im['id']][col], 84)
                im[f'{col}_50_percentile'] = np.percentile(self.samples_dict[im['id']][col], 50)
        

    def compute_map_as(self, fits_file_path, which='dpl') :
        with fits.open(fits_file_path) as hdul :
            data_hdu = None
            for hdu in hdul :
                if getattr(hdu, 'data', None) is not None :
                    data_hdu = hdu
                    break
            if data_hdu is None :
                raise ValueError(f"No image data found in FITS file: {fits_file_path}")
            
            target_header = data_hdu.header
            target_data = data_hdu.data
            target_wcs = WCS(target_header)
        
        ndim = target_data.ndim
        if ndim > 2 :
            ny, nx = target_data.shape[-2], target_data.shape[-1]
        else :
            ny, nx = target_data.shape
        
        which_norm = which.lower().strip()
        if which_norm in ['displacement', 'disp'] :
            which_norm = 'dpl'
        elif which_norm in ['magnification', 'mu', 'ampli'] :
            which_norm = 'magnification'
        elif which_norm in ['convergence', 'kappa', 'mass'] :
            which_norm = 'convergence'
        elif which_norm in ['shear', 'gamma'] :
            which_norm = 'shear'
        
        valid_which = ['dpl', 'magnification', 'shear', 'convergence']
        if which_norm not in valid_which :
            raise ValueError(
                f"Invalid 'which' value: {which}. "
                f"Allowed values are: {valid_which}"
            )
        
        if which_norm == 'dpl' :
            map_data = np.full((2, ny, nx), np.nan, dtype=float)
        else :
            map_data = np.full((ny, nx), np.nan, dtype=float)
        self._vprint(f"Computing '{which_norm}' map on input pixel grid ({ny} x {nx})...")
        for iy in tqdm(range(ny), disable=not self.verbose) :
            for ix in range(nx) :
                if ndim == 3 :
                    ra, dec, _ = target_wcs.pixel_to_world_values(ix, iy, 0.)
                else :
                    ra, dec = target_wcs.pixel_to_world_values(ix, iy)
                    
                if which_norm == 'dpl' :
                    dx, dy = self.compute_displacement(ra, dec)
                    map_data[0, iy, ix] = dx
                    map_data[1, iy, ix] = dy
                elif which_norm == 'magnification' :
                    map_data[iy, ix] = self.compute_magnification(ra, dec)
                elif which_norm == 'shear' :
                    map_data[iy, ix] = self.compute_shear(ra, dec, which='gamma')
                elif which_norm == 'convergence' :
                    map_data[iy, ix] = self.compute_mass(ra, dec)
        
        out_header = target_wcs.to_header()
        out_header['BUNIT'] = 'arcsec' if which_norm == 'dpl' else ''
        out_header['MAPTYPE'] = which_norm
                
        input_dir = os.path.dirname(fits_file_path)
        input_name = os.path.basename(fits_file_path)
        root_name, _ = os.path.splitext(input_name)
        out_path = os.path.join(input_dir, f"{root_name}_{which_norm}_map.fits")
        
        if which_norm == 'dpl' :
            fits.PrimaryHDU(data=map_data[0], header=out_header).writeto(out_path.replace('.fits', '_dx.fits'), overwrite=True)
            fits.PrimaryHDU(data=map_data[1], header=out_header).writeto(out_path.replace('.fits', '_dy.fits'), overwrite=True)
        else :
            fits.PrimaryHDU(data=map_data, header=out_header).writeto(out_path, overwrite=True)
        self._vprint(f"Saved map to {out_path}")
        
        arcsec_per_pixel = target_wcs.wcs.cdelt[0] * 3600
        def plot_source_plane_grid(dpl_maps) :
            # Careful: this function only works when x, y match ra, dec axes
            fig, ax = plt.subplots()
            x, y = np.meshgrid(np.arange(0, dpl_maps[0].shape[1], 1), np.arange(0, dpl_maps[0].shape[0], 1))
            x_im = x * arcsec_per_pixel
            y_im = y * arcsec_per_pixel
            x_src = x_im - dpl_maps[0]
            y_src = y_im - dpl_maps[1]
    
            ax.scatter(x_src, y_src, marker = '.', sizes=[0.5], alpha=0.5)
            ax.set_xlabel('x (arcsec)')
            ax.set_ylabel('y (arcsec)')
            ax.set_title('Data grid in the source plane')
            plt.show()
        
        if which_norm == 'dpl' :
            plot_source_plane_grid(map_data)
        
        return map_data, target_wcs


def import_lenstool(model_dir, fits_image, use_wrapper=None, compute_predictions=False, verbose=True) :
    return lenstool_model(model_dir, fits_image, use_wrapper=use_wrapper, compute_predictions=compute_predictions, verbose=verbose)











