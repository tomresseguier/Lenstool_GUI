import os
import re
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
from tqdm import tqdm

import PyQt5
from PyQt5.QtWidgets import QMainWindow, QSplitter
import pyqtgraph as pg
from PyQt5.QtCore import Qt
from astropy.table import Table

from .catalog import catalog
from .lenstool_model.lenstool_model import lenstool_model
from .source_extraction.source_extract import source_extract, source_extract_DIM
from .utils.utils_fits.utils_fits_image import open_image
from .utils.utils_plots.plot_utils_general import *
from .utils.utils_Qt.selectable_classes import *
from .utils.utils_Qt.drag_widgets import DragWidget
from .utils.utils_Qt.utils_general import *
from .utils.utils_astro.get_cosmology import get_cosmo




pg.setConfigOption('imageAxisOrder', 'row-major')



class fits_image :
    def __init__(self, image_path, main_window=None, plot_image=True) :
        self.image_path = image_path
        # If an external QMainWindow is provided (e.g. from the GUI), use it
        # instead of spawning a new independent window.
        self.main_window = main_window  # type: ignore[assignment]
        if os.path.isfile(self.image_path[:-8] + 'wht.fits') :
            print("Weight file found: " + self.image_path[:-8] + 'wht.fits')
            self.weight_path = self.image_path[:-8] + 'wht.fits'
        else :
            self.weight_path = None
        self.image_data, self.pix_deg_scale, self.orientation, self.wcs, self.header = open_image(self.image_path)
        self.sources = None
        self.fig = None
        self.ax = None
        self.multiple_images = None
        self.galaxy_selection = None
        self.imported_cat = None
        self.imported_cat_list = []
        self.qt_image = None
        self.qtItems_dict = {'sources': None,
                             'potfile_cat': None,
                             'imported_cat': None,
                             'multiple_images': None}
        self.ax = None
        self.redshift = None
        self.boosted_image_path = self.image_path[:-5] + '_boosted.fits'
        if os.path.isfile(self.boosted_image_path) :
            print("Boosted image found: " + self.boosted_image_path)
            self.boosted_image, _, _, _, _ = open_image(self.boosted_image_path)
        else :
            self.boosted_image = None
        self.boosted = False
        self.qt_image_list = []
        self.qt_window_list = []
        if self.orientation==None :
            self.orientation = 0.
        self.cosmo = get_cosmo()

        if plot_image :
            self.plot_image()
        self.filters = None
    
    
    def set_cosmo(self, cosmo_name) :
        self.cosmo = get_cosmo(cosmo_name)
    
    def create_qt_image(self) :
        to_plot = np.flip(self.image_data, axis=0) if not self.boosted else np.flip(self.boosted_image, axis=0)
        
        #qt_image = pg.image(to_plot)
        qt_image = pg.ImageView()
        qt_image.setImage(to_plot)
        #qt_image.autoLevels()
        
        image_widget = DragWidget(qt_image)
        
        qt_layout = QSplitter(Qt.Horizontal)
        qt_layout.addWidget(image_widget)
        
        if self.main_window is None:
            window = QMainWindow()
            window.setWindowTitle(os.path.basename(self.image_path))
            window.setCentralWidget(qt_layout)
            window.show()
        else:
            # Reuse the provided main window.
            window = self.main_window
            # Replace any existing central widget.
            window.setCentralWidget(qt_layout)
            window.setWindowTitle(os.path.basename(self.image_path))
            # Ensure the window is visible (may already be shown).
            window.show()

        return qt_image, qt_layout, image_widget, window
    
    def _flush_closed_windows(self) :
        if self.qt_image is not None :
            if not self.qt_image.isVisible() :
                self.window.close()
                self.qt_image, self.qt_layout, self.image_widget, self.window = None, None, None, None
        for window in self.qt_window_list :
            if not window.isVisible() :
                window.close()
        self.qt_image_list = [ image for window, image in zip(self.qt_window_list, self.qt_image_list) if window.isVisible() ]
        self.qt_window_list = [ window for window in self.qt_window_list if window.isVisible() ]
        return

    def link_zoom(self) :
        self._view_range_sync_state()['zoom'] = True
        self._connect_view_range_sync()

    def unlink_zoom(self) :
        self._view_range_sync_state()['zoom'] = False
        self._disconnect_view_range_sync_if_idle()

    def link_pan(self) :
        self._view_range_sync_state()['pan'] = True
        self._connect_view_range_sync()

    def unlink_pan(self) :
        self._view_range_sync_state()['pan'] = False
        self._disconnect_view_range_sync_if_idle()

    def link_zoom_and_pan(self) :
        sync = self._view_range_sync_state()
        sync['zoom'] = True
        sync['pan'] = True
        self._connect_view_range_sync()

    def unlink_zoom_and_pan(self) :
        sync = self._view_range_sync_state()
        sync['zoom'] = False
        sync['pan'] = False
        self._disconnect_view_range_sync_if_idle()

    def plot_image(self) :
        #to_plot = self.image_data
        #to_plot = np.transpose(self.image_data, axes=[1,0,2])
        self._flush_closed_windows()
        if self.qt_image is None :
            print('Creating main window...')
            self.qt_image, self.qt_layout, self.image_widget, self.window = self.create_qt_image()
            self.qt_window_list = [self.window] + self.qt_window_list
            self.qt_image_list = [self.qt_image] + self.qt_image_list
            print('Done')
        else :
            print('Creating secondary window...')
            extra_qt_image, _, _, extra_window = self.create_qt_image()
            extra_window.setWindowTitle(os.path.basename(self.image_path) + ' (' + str(len(self.qt_image_list)+1) + ')')
            extra_window.show()
            self.qt_image_list.append(extra_qt_image)
            self.qt_window_list.append(extra_window)
            print('Done')
        if hasattr(self, '_view_range_sync') and (self._view_range_sync['zoom'] or self._view_range_sync['pan']) :
            self._connect_view_range_sync()
        return
    
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
            self.qt_image.setImage(np.flip(self.boosted_image, axis=0))
            #self.qt_image.autoLevels()
            self.boosted = True
            for extra_qt_image in self.qt_image_list :
                extra_qt_image.setImage(np.flip(self.boosted_image, axis=0))
            print('Done')
    
    def unboost(self) :
        if self.boosted :
            print('Plotting...')
            self.qt_image.setImage(np.flip(self.image_data, axis=0))
            #self.qt_image.autoLevels()
            self.boosted = False
            for extra_qt_image in self.qt_image_list :
                extra_qt_image.setImage(np.flip(self.image_data, axis=0))
            print('Done')
    
    def set_weight(self, weight_path) :
        self.weight_path = weight_path
    
    def extract_sources(self, image_path=None, weight_path=None, DIM_ref_path=None, rerun=False, reproject=True) :
        print("Source extraction not yet available. Work in progress.")
        if False :
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
    
    
    def make_catalog(self, cat, color=[1., 1., 0], units=None, verbose=True) :
        if self.qt_image is None :
            self.plot_image()
        #to_return = catalog(cat, self.image_data, self.wcs, self.qt_image, window=self.window, image_path=self.image_path, 
        #                            image_widget = self.image_widget, qt_layout=self.qt_layout, color=color, 
        #                            mag_colnames=mag_colnames, mpl_fig=self.fig, mpl_ax=self.ax, 
        #                            pix_deg_scale=self.pix_deg_scale, units=units)
        to_return = catalog(cat, self, color=color, units=units, verbose=verbose)
        return to_return
    
    def import_catalog(self, cat, color=None, units='pixel') :
        self.imported_cat = self.make_catalog(cat, color=color, units=units)
        self.imported_cat_list.append(self.imported_cat)
        
    ###########################################################################
    
    
    def world_to_image(self, ra, dec, unit='deg') :
        coord = SkyCoord(ra, dec, unit=unit)
        image_coord = WCS.world_to_pixel(self.wcs, coord)
        if len(image_coord[0].shape)==0 :
            image_coord = (image_coord[0]*1., image_coord[1]*1.)
        return image_coord
    
    def image_to_world(self, x, y, unit='deg') :
        world_coord = WCS.pixel_to_world(self.wcs, x, y)
        return world_coord.ra.deg, world_coord.dec.deg
    
    def clear_Items(self) :
        for key in self.qtItems_dict.keys() :
            if self.qtItems_dict[key] is not None :
                for i in tqdm( range(len(self.qtItems_dict[key])) ) :
                    self.qt_image.removeItem(self.qtItems_dict[key][i])
    
    
        
    def import_lenstool(self, model_dir, use_wrapper=None, compute_predictions=True, verbose=True) :
        #self.lt_dir = model_dir
        #if hasattr(self, 'lt'):
        #    del self.lt
        self.lt = lenstool_model(model_dir, self, use_wrapper=use_wrapper, compute_predictions=compute_predictions, verbose=verbose)
    
    
    
    def start_hand_select(self) :
        cat_dict = {'id': [], 'ra': [], 'dec': [], 'x': [], 'y': [], 'a': [], 'b': [], 'theta': []}
        cat = Table(cat_dict)
        self.hand_made_cat = catalog(cat, self, units='pixel')
        
        def mouse_clicked(evt):
            if evt.double():
                pos = evt.scenePos()
                if self.qt_image.getView().sceneBoundingRect().contains(pos):
                    mouse_point = self.qt_image.getView().mapSceneToView(pos)
                    x, y_flipped = mouse_point.x(), mouse_point.y()
                    x, y = x, self.image_data.shape[0] - y_flipped
                    ra, dec = self.image_to_world(x, y)
                    self.hand_made_cat.cat.add_row( {'id': [len(self.hand_made_cat.cat) + 1], 'ra': [ra], 'dec': [dec], 'x': [x], 'y': [y], 'a': [10], 'b': [10], 'theta': [0]} )
                    self.hand_made_cat.qtItems.append(PyQt5.QtWidgets.QGraphicsEllipseItem())
                    self.hand_made_cat.plot(color=[1,1,1,0])
                    
        self._doubleclick_connection = self.qt_image.scene.sigMouseClicked.connect(mouse_clicked)
        
        def keyPressEvent(event):
            #print('Hand selection stopped.')
            if event.key() == Qt.Key_Escape or event.key() == Qt.Key_Space :
                if hasattr(self, '_doubleclick_connection'):
                    self.qt_image.scene.sigMouseClicked.disconnect(self._doubleclick_connection)
                    del self._doubleclick_connection
                self.hand_made_cat.clear()
                self.window.keyPressEvent = self._original_keyPressEvent
                print('Hand selection stopped.')
        
        self._original_keyPressEvent = self.window.keyPressEvent
        self.window.keyPressEvent = keyPressEvent
    
    
    
    def plot_image_mpl(self, wcs_projection=True, units='pixel', pos=111, make_axes_labels=True, make_grid=True, crop=None, replace_image=True, extra_pad=None) :
        fig, ax = plot_image_mpl(self.image_data, wcs=self.wcs, wcs_projection=wcs_projection, units=units, pos=pos, \
                                 make_axes_labels=make_axes_labels, make_grid=make_grid, crop=crop, extra_pad=extra_pad)
        plot_NE_arrows(ax, self.wcs)
        if replace_image :
            self.fig, self.ax = fig, ax
        return fig, ax
    
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
    
    def load_filters(self, filter_dir=None):
        if filter_dir==None :
            listdir = os.listdir( os.path.dirname(self.image_path) )
            listdir_lower = [ name.lower() for name in os.listdir( os.path.dirname(self.image_path) ) ]
            indices = np.where([ 'filter' in name for name in listdir_lower ])[0]
            filters = []
            if len(indices)>0 :
                for i in indices :
                    filter_dir = os.path.join(os.path.dirname(self.image_path), listdir[i])
                    print('Looking for filters in directory: ' + filter_dir)
                    f = load_filters(filter_dir)
                    filters.append(f)
                l = np.array( [ len(f) for f in filters ] )
                idx = np.argmax(l)
                self.filters = filters[idx]
                print( str(len(self.filters)) + ' filters loaded from directory: ' + listdir[indices[idx]] )
            else :
                raise ValueError('No filter directory found in image directory.')
        else :
            self.filters = load_filters(filter_dir)
    
    def add_scale_bar(self, z=None, unit='arcmin', length=1, color='white', linewidth=4, text_offset=0.01, position=['bottom', 'right']):
        """Add a floating scale bar that stays anchored in the chosen viewport corner.

        The bar is a UIGraphicsItem (pg.ScaleBar) so it lives in screen/viewport
        coordinates: it never pans with the image and its on-screen pixel width
        automatically updates as the user zooms in or out, always representing the
        same angular or physical scale.

        Parameters
        ----------
        z : float, optional
            Source redshift, required when *unit* is a physical length ('pc', 'kpc', 'Mpc').
        unit : str
            One of 'arcsec', 'arcmin', 'pc', 'kpc', 'Mpc'.
        length : float
            Desired scale-bar length in *unit*.
        color : str or tuple
            Bar and label colour accepted by pyqtgraph (e.g. 'white', (255,255,255)).
        linewidth : int
            Thickness of the bar rectangle in screen pixels.
        text_offset : float
            Unused legacy parameter (kept for API compatibility).
        position : list of str
            Two-element list with vertical ('top'/'bottom') and horizontal
            ('left'/'right') placement of the bar, e.g. ['bottom', 'left'].
        """
        print(self.cosmo)

        if unit == 'arcsec':
            length_deg = length / 3600.
            label = f'{length}"'
        elif unit == 'arcmin':
            length_deg = length / 60.
            label = f"{length}'"
        elif unit in ['pc', 'kpc', 'Mpc']:
            if z is None:
                raise ValueError("Redshift z must be provided for physical units.")
            ang_diam_dist = self.cosmo.angular_diameter_distance(z).to(unit).value
            length_rad = length / ang_diam_dist
            length_deg = np.rad2deg(length_rad)
            # Angular equivalent shown as a subtitle under the physical label.
            # Pick arcsec when the value is below 60, arcmin otherwise.
            length_arcsec = length_deg * 3600.
            if length_arcsec < 60.:
                ang_fmt = f'{length_arcsec:.1f}"' if length_arcsec < 10. else f'{length_arcsec:.0f}"'
            else:
                length_arcmin = length_deg * 60.
                ang_fmt = f"{length_arcmin:.1f}'" if length_arcmin < 10. else f"{length_arcmin:.0f}'"
            label = f'{length} {unit}\n{ang_fmt}'
        else:
            raise ValueError(f"Unknown unit: {unit}")

        # Bar length expressed in image pixels (the data coordinate unit of the ViewBox)
        bar_length_pix = length_deg / self.pix_deg_scale

        # Compute anchor parameters before any item manipulation
        # UIGraphicsItem.anchor(itemPos, parentPos, offset):
        #   itemPos  – which point of the bar item to pin (0=top/left, 1=bottom/right)
        #   parentPos – which point of the ViewBox to pin to (same convention)
        #   offset    – additional screen-pixel nudge (positive x → right, positive y → down)
        margin = 20  # screen pixels from the viewport edge
        vertical   = 'bottom' if 'bottom' in position else 'top'
        horizontal = 'left'   if 'left'   in position else 'right'
        vy = 1 if vertical   == 'bottom' else 0
        vx = 0 if horizontal == 'left'   else 1
        ox = margin  if horizontal == 'left'   else -margin
        oy = -margin if vertical   == 'bottom' else  margin

        # If a bar already exists, update it in-place rather than destroying it.
        # Destroying via setParentItem(None) leaves the bar's viewRangeChanged slot
        # still connected to the ViewBox; when the view later updates it fires on the
        # now-parentless item and raises "Cannot anchor; parent is not set."
        if hasattr(self, '_scale_bar') and self._scale_bar is not None:
            self._scale_bar.size   = bar_length_pix
            self._scale_bar._width = linewidth
            self._scale_bar.brush  = pg.mkBrush(color)
            self._scale_bar.pen    = pg.mkPen(None)
            self._scale_bar.text.setText(label)
            self._scale_bar.text.setColor(color)
            self._scale_bar.anchor((vx, vy), (vx, vy), offset=(ox, oy))
            self._scale_bar.update()
            return self._scale_bar

        view = self.qt_image.getView()

        # pg.ScaleBar is a UIGraphicsItem: it renders at a fixed screen position and
        # redraws its pixel width automatically whenever the ViewBox zoom changes so
        # that it always represents bar_length_pix data-space units.
        bar = pg.ScaleBar(
            size=bar_length_pix,
            width=linewidth,
            brush=pg.mkBrush(color),
            pen=pg.mkPen(None),   # no outline around the filled rectangle
            suffix='',
        )
        bar.text.setText(label)
        bar.text.setColor(color)
        bar.setParentItem(view)
        bar.anchor((vx, vy), (vx, vy), offset=(ox, oy))

        self._scale_bar = bar
        return bar
    




    ########## View range synchronization functions##########
    def _all_viewboxes(self) :
        self._flush_closed_windows()
        return [iv.getView() for iv in self.qt_image_list]

    def _view_range_sync_state(self) :
        if not hasattr(self, '_view_range_sync') :
            self._view_range_sync = {'zoom': False, 'pan': False, 'connections': {}, 'guard': False}
        return self._view_range_sync

    def _connect_view_range_sync(self) :
        sync = self._view_range_sync_state()
        vbs = self._all_viewboxes()
        for vb in list(sync['connections'].keys()) :
            if vb not in vbs :
                try :
                    vb.sigRangeChanged.disconnect(sync['connections'][vb])
                except (TypeError, RuntimeError) :
                    pass
                del sync['connections'][vb]
        for vb in vbs :
            if vb not in sync['connections'] :
                sync['connections'][vb] = vb.sigRangeChanged.connect(
                    lambda _vb=vb: self._sync_linked_view_range(_vb)
                )

    def _sync_linked_view_range(self, source_vb) :
        sync = self._view_range_sync_state()
        if sync['guard'] or (not sync['zoom'] and not sync['pan']) :
            return
        vbs = self._all_viewboxes()
        if len(vbs) < 2 :
            return
        sync['guard'] = True
        try :
            xrange, yrange = source_vb.viewRange()
            cx = 0.5 * (xrange[0] + xrange[1])
            cy = 0.5 * (yrange[0] + yrange[1])
            wx = xrange[1] - xrange[0]
            wy = yrange[1] - yrange[0]
            for vb in vbs :
                if vb is source_vb :
                    continue
                ox0, ox1 = vb.viewRange()[0]
                oy0, oy1 = vb.viewRange()[1]
                if sync['zoom'] and sync['pan'] :
                    new_x, new_y = xrange, yrange
                elif sync['zoom'] :
                    occx = 0.5 * (ox0 + ox1)
                    occy = 0.5 * (oy0 + oy1)
                    new_x = (occx - 0.5 * wx, occx + 0.5 * wx)
                    new_y = (occy - 0.5 * wy, occy + 0.5 * wy)
                else :
                    owx = ox1 - ox0
                    owy = oy1 - oy0
                    new_x = (cx - 0.5 * owx, cx + 0.5 * owx)
                    new_y = (cy - 0.5 * owy, cy + 0.5 * owy)
                vb.setRange(xRange=new_x, yRange=new_y, padding=0)
        finally :
            sync['guard'] = False

    def _disconnect_view_range_sync_if_idle(self) :
        sync = self._view_range_sync_state()
        if sync['zoom'] or sync['pan'] :
            return
        for vb, conn in list(sync['connections'].items()) :
            try :
                vb.sigRangeChanged.disconnect(conn)
            except (TypeError, RuntimeError) :
                pass
        sync['connections'].clear()






class FilterImage(fits_image) :
    def __init__(self, image_path, main_window=None, plot_image=False) :
        super().__init__(image_path, main_window=main_window, plot_image=plot_image)
        self.filter = self._get_filter_name()
        self.wavelength = self._get_pivot_wavelength()
        psf_path = self._get_psf_path()
        if psf_path is not None :
            self.psf = PSF(psf_path)
        else :
            self.psf = None
            
    def _get_filter_name(self):
        """Extract filter name from header or filename."""
        # Try header first
        for key in ['FILTER', 'FILTER1', 'FILTER2']:
            if key in self.header and isinstance(self.header[key], str):
                val = self.header[key].strip().upper()
                if re.match(r'F\d{3,4}[WMN]*', val):
                    return val

        # Try filename
        m = re.search(r'F\d{3,4}[WMN]*', os.path.basename(self.image_path).upper())
        if m:
            return m.group(0)

        print(f"Could not determine filter name for {os.path.basename(self.image_path)}.")
        return "UNKNOWN"
    
    def _get_pivot_wavelength(self):
        lam = None
        f_lower = self.filter.lower()
        
        for key in ['PHOTPLAM', 'PIVOTWL', 'WAVELENGTH']:
            if key in self.header:
                lam = self.header[key]
                break

        if lam is None:
            lam = np.nan
            print(f"Could not determine wavelength for {self.filter} ({self.image_path}).")

        return lam
    
    def _get_psf_path(self):
        # Search first in parent of image dir, then in image dir
        search_dirs = [
            os.path.dirname(os.path.dirname(self.image_path)),
            os.path.dirname(self.image_path)
        ]

        for base in search_dirs:
            if not os.path.isdir(base):
                continue
            try:
                listdir = os.listdir(base)
            except Exception:
                continue

            listdir_lower = [name.lower() for name in listdir]
            indices = [i for i, name in enumerate(listdir_lower) if 'psf' in name]
            if len(indices) == 0:
                continue

            i = indices[0]
            psf_dir = os.path.join(base, listdir[i])
            print('Looking for psf in directory: ' + psf_dir)
            possible_psf_files = []
            for f in glob.glob(os.path.join(psf_dir, '*.fits')):
                if self.filter.lower() in os.path.basename(f).lower():
                    possible_psf_files.append(f)
            if len(possible_psf_files) > 0:
                print(f"PSF file found for filter {self.filter}: {os.path.basename(possible_psf_files[0])}")
                return possible_psf_files[0]
            else:
                print(f"No PSF file found for filter {self.filter} in {psf_dir}.")

        return None


class PSF :
    def __init__(self, psf_path) :
        with fits.open(psf_path) as hdu :
            self.data = hdu[0].data
            self.wcs = WCS(hdu[0].header)
    
    def plot(self) :
        fig, ax = plt.subplots()
        ax.imshow( np.log(self.data) , origin='lower')
        

def load_filters(path):
    """
    Load all FITS files in `path` and return a dict {filter_name: FilterImage},
    ordered by increasing pivot wavelength.
    If multiple files correspond to the same filter, selects the one most likely
    to be the science image (e.g., containing 'SCI' in header or filename).
    """
    fits_files = [f for f in os.listdir(path) if f.lower().endswith('.fits')]
    filter_groups = {}

    # --- Group files by detected filter name ---
    for fname in fits_files:
        full_path = os.path.join(path, fname)
        try:
            # Just read header to find filter (faster)
            with fits.open(full_path) as hdul:
                hdr = hdul[0].header
                filt = None
                for key in ['FILTER', 'FILTER1', 'FILTER2']:
                    if key in hdr and isinstance(hdr[key], str):
                        val = hdr[key].strip().upper()
                        if re.match(r'F\d{3,4}[WMN]*', val):
                            filt = val
                            break
                if filt is None:
                    m = re.search(r'F\d{3,4}[WMN]*', fname.upper())
                    if m:
                        filt = m.group(0)
        except Exception:
            print(f"Could not read header for {fname}. Skipping.")
            continue

        if filt is None:
            continue

        filter_groups.setdefault(filt, []).append(full_path)

    # --- For each filter, pick best file if multiple ---
    selected_files = {}
    for filt, files in filter_groups.items():
        if len(files) == 1:
            selected_files[filt] = files[0]
        else:
            sci_candidates = []
            for f in files:
                try:
                    with fits.open(f) as hdul:
                        hdr = hdul[0].header
                        for key in ['EXTNAME', 'FILETYPE', 'IMAGETYP']:
                            if key in hdr and 'SCI' in str(hdr[key]).upper():
                                sci_candidates.append(f)
                                break
                    if 'SCI' in os.path.basename(f).upper():
                        sci_candidates.append(f)
                except Exception:
                    continue

            if sci_candidates:
                selected_files[filt] = sci_candidates[0]
                print(f"Selected {os.path.basename(sci_candidates[0])} for {filt}")
            else:
                print(f"Multiple files for {filt}, picking first arbitrarily.")
                selected_files[filt] = files[0]

    # --- Create FilterImage instances ---
    filters_list = [FilterImage(fpath, plot_image=False) for fpath in selected_files.values()]
    filters_list.sort(key=lambda f: np.inf if np.isnan(f.wavelength) else f.wavelength)

    return {f.filter: f for f in filters_list}



