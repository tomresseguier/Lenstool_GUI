import numpy as np
import os
import pyqtgraph as pg
from PyQt5.QtWidgets import QMainWindow, QSplitter, QWidget, QVBoxLayout, QHBoxLayout
from PyQt5.QtCore import Qt
from pyqtgraph.Qt import QtWidgets

from lenstronomy.Data.pixel_grid import PixelGrid     
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.Data.psf import PSF
#from lenstronomy.Plots import lens_plot
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.Util import util

from ...utils.utils_astro.utils_general import world_to_relative
from ...utils.utils_general.sort_points import break_curves
from ...utils.utils_Qt.drag_widgets import DragPlotWidget_special
from ...utils.utils_Qt.utils_general import transform_rectangle
from .event_filters import SourceFilter, ImageFilter
from .lenstronomy_model import lenstronomy_model
from .utils import make_lm_dict_opt
from .simulate import simulate
from .optimize import optimize



class image_simulator :
    def __init__(self, fits_image, which_filter=None, throttle_mode=0, use_linear_solver=False, dpl_resolution=200) :
        self.fits_image = fits_image
        self.z_source = fits_image.lt.lt_z
        self.ROI = fits_image.image_widget.current_ROI
        #self.fits_image.qt_image.addItem(ROI)
        self.models=None
        self.use_linear_solver = use_linear_solver
        
        if not hasattr(self, 'previous_state_current_ROI') :
            self.previous_state_current_ROI = None
        
        ######### Get the ROI's coordinates #########
        x0 = self.ROI.getState()['pos'][0]
        y0 = self.ROI.getState()['pos'][1]
        a = self.ROI.getState()['size'][0]
        b = self.ROI.getState()['size'][1]
        angle = self.ROI.getState()['angle'] *np.pi/180
        
        x0, y0, a, b, angle = transform_rectangle(x0, y0, a, b, angle) #x0, y0 at the top left, a and b > 0
        y0 = fits_image.image_data.shape[0] - y0 # Counting pixels from bottom instead of top

        # Here we transform the ROI's parameters into a square area with corners corresponding to integer pixel coordinates.
        x_center, y_center = x0 + a/2, y0 - b/2
        square_size = max(a, b)
        square_size_floor = int(square_size)
        _x, _y = abs(x_center - round(x_center)), abs(y_center - round(y_center))
        if _x + _y > 0.5 : # The center will be the center of a pixel --> odd number of pixels
            rule = 'odd'
            x_center_pix, y_center_pix = int(x_center) + 0.5, int(y_center) + 0.5
            square_size_pix = square_size_floor + 1 if square_size_floor%2 == 0 else square_size_floor # odd number
        else : # The center will be the corner of a pixel --> even number of pixels
            rule = 'even'
            x_center_pix, y_center_pix = round(x_center), round(y_center)
            square_size_pix = square_size_floor if square_size_floor%2 == 0 else square_size_floor + 1 # even number
        
        # Store the exact pixel-crop origin used for image_data so interactive clicks
        # can be mapped back to full-image coordinates without sub-pixel offsets.
        self._crop_npix = square_size_pix
        self._crop_x0 = int( x_center_pix - square_size_pix/2 ) # x_center_pix - square_size_pix/2 is float with 0 decimal
        self._crop_y0 = int( y_center_pix + square_size_pix/2 )
        self.anchor = (self._crop_x0, self._crop_y0)
        self.center_world = fits_image.image_to_world(x_center_pix, y_center_pix)

        xr_center, yr_center = fits_image.lt.world_to_relative(self.center_world[0], self.center_world[1])
        square_size_arcsec = square_size_pix * fits_image.pix_deg_scale*3600
        self._SquareOfInterest_side_arcsec = square_size_arcsec
        self._SquareOfInterest_xr_bottomleft = xr_center - square_size_arcsec / 2
        self._SquareOfInterest_yr_bottomleft = yr_center - square_size_arcsec / 2
        new_field = [self._SquareOfInterest_xr_bottomleft, self._SquareOfInterest_xr_bottomleft + square_size_arcsec, 
                     self._SquareOfInterest_yr_bottomleft, self._SquareOfInterest_yr_bottomleft + square_size_arcsec]
        
        ######### Set Lenstool field corresponding to ROI's position and calculate the local displacement map #########
        if fits_image.lt.lt is not None :
            self._initial_field = fits_image.lt.lt.get_field([])
            fits_image.lt.lt.set_field(new_field)
            if self.previous_state_current_ROI is None or self.previous_state_current_ROI != self.ROI.getState() or fits_image.lt.dpl_maps[fits_image.lt.lt_z][0].shape[0] != dpl_resolution :
                fits_image.lt.compute_lt_dpl(npix=dpl_resolution)
                fits_image.lt.compute_lt_convergence(npix=dpl_resolution)
        else :
            if fits_image.lt.lt_z not in fits_image.lt.dpl_maps :
                raise KeyError('Displacement maps for redshift {} have not been calculated yet.'.format(fits_image.lt.lt_z))
            else :
                print('Lenstool instance not initialized. Using saved displacement maps instead. \
                      Careful: saved maps only cover the area you selected if you previously computed them for this specific area!')
        
        #-------------- Lens model Lenstronomy definitions --------------#
        
        deltaPix = abs( fits_image.lt.dpl_maps[fits_image.lt.lt_z][2].wcs.cdelt[0]*3600 )
        x_grid_interp, y_grid_interp = util.make_grid(dpl_resolution, deltaPix)
        x_axes, y_axes = util.get_axes(x_grid_interp, y_grid_interp)
        
        self.LensModel_kwargs = [{'grid_interp_x': x_axes,
                                  'grid_interp_y': y_axes,
                                  'f_': fits_image.lt.convergence_maps[fits_image.lt.lt_z][0],
                                  'f_x': fits_image.lt.dpl_maps[fits_image.lt.lt_z][0],
                                  'f_y': fits_image.lt.dpl_maps[fits_image.lt.lt_z][1]}]
        
        self.LensModel_list = ['INTERPOL']
        self.LensModel = LensModel(lens_model_list=self.LensModel_list)
        
        #m=4
        #f, ax = plt.subplots(1, 1, figsize=(10, 10), sharex=False, sharey=False)
        #lens_plot.lens_model_plot(ax, lensModel=self.lens_model, kwargs_lens=kwargs_lens,
        #                          with_caustics=True, fast_caustic=True, coord_inverse=False, numPix=round(dpl_resolution/m), deltaPix=deltaPix*m)
        
        #-------------- Qt graphics definitions --------------#
        
        #This version of the code tends to freeze unfortunately
        self.source_plane_widget = DragPlotWidget_special(throttle_mode=throttle_mode)
        self.source_plane_widget.setTitle('Source plane')
        #source_plane_widget.setAspectLocked(lock=True, ratio=1)
        
        self._imsim_srcplane_layout = QSplitter(Qt.Vertical)
        self._imsim_implane_layout = QWidget()
        
        self._imsim_srcplane_layout.addWidget(self.source_plane_widget)
        
        self._imsim_layout = QSplitter(Qt.Horizontal)
        self._imsim_layout.addWidget(self._imsim_srcplane_layout)
        self._imsim_layout.addWidget(self._imsim_implane_layout)
        
        self.window = QMainWindow()
        self.window.setWindowTitle('Lenstronomy image simulator')
        self.window.setCentralWidget(self._imsim_layout)
        
        
        #-------------- Plot critical & caustic curves --------------#

        self.source_center_coordinates = self.LensModel.ray_shooting(0, 0, self.LensModel_kwargs)
        
        if fits_image.lt.lt is not None and (self.previous_state_current_ROI is None or self.previous_state_current_ROI != self.ROI.getState()) :
            fits_image.lt.compute_lt_curve()
        fits_image.lt.plot_lt_curve(color=[0, 255, 0], which='caustic') #to create self.lt_curve_coords_image_sorted
        
        caustic_plot = pg.PlotDataItem()
        caustic_plot.setPen( color=[0,255,0,255], width=4.0001 )
        coords = break_curves(fits_image.lt.lt_caustic_coords_relative, distance_threshold=8.*fits_image.pix_deg_scale*3600) #sort_points(self.lt_caustic_coords_relative, distance_threshold=1.0, angle_threshold=np.pi)
        x = coords[0] - xr_center - self.source_center_coordinates[0]
        y = coords[1] - yr_center - self.source_center_coordinates[1]
        caustic_plot.setData(x, y)
        self.source_plane_widget.addItem(caustic_plot)
        self.source_plane_widget.setXRange(-3, 3)
        self.source_plane_widget.setYRange(-3, 3)
        
        self._lt_curve_coords_relative_broken = break_curves(fits_image.lt.lt_curve_coords_relative, distance_threshold=8.*fits_image.pix_deg_scale*3600) #sort_points(self.lt_curve_coords_relative, distance_threshold=1.0, angle_threshold=np.pi)
        
        self.filter_source = SourceFilter(self)
        self.source_plane_widget.installEventFilter(self.filter_source)
        
        
        #-------------- Plot source positions from multiple image catalog --------------#
        SX = []
        SY = []
        XR = []
        YR = []
        for mult in fits_image.lt.mult.cat :
            xr, yr = world_to_relative( mult['ra'], mult['dec'], self.center_world )
            XR.append(xr)
            YR.append(yr)
            sx, sy = self.LensModel.ray_shooting(xr, yr, self.LensModel_kwargs)
            SX.append(sx - self.source_center_coordinates[0])
            SY.append(sy - self.source_center_coordinates[1])
        
        source_coord_plot = pg.ScatterPlotItem()
        source_coord_plot.setData(SX, SY)
        self.source_plane_widget.addItem(source_coord_plot)
        
        
        #-------------- Set the PixelGrid for lenstronomy --------------#
        transform_pix2angle = np.array([[1, 0], [0, 1]]) * fits_image.pix_deg_scale*3600
        
        npix = self._crop_npix
        self.PixelGrid_kwargs = {'nx': npix,
                                    'ny': npix,  # number of pixels per axis
                                    'ra_at_xy_0': - square_size_arcsec / 2,
                                    'dec_at_xy_0': - square_size_arcsec / 2,
                                    'transform_pix2angle': transform_pix2angle} 
        self.PixelGrid = PixelGrid(**self.PixelGrid_kwargs)
        
        
        #-------------- DATA --------------#
        if which_filter is not None :
            if which_filter in fits_image.filters.keys() :
                print('Using filter ' + which_filter)
                self.individual_filter = fits_image.filters[which_filter]
            elif type(which_filter) is int :
                if which_filter>2 or which_filter<0 :
                    raise('Filter index must belong to [0,1,2].')
                c = ['red', 'green', 'blue']
                which_filter_str = str(c[which_filter])
                print('Using ' + which_filter_str + ' image.')
                if fits_image.filters is None :
                    fits_image.filters = {}
                if which_filter_str not in fits_image.filters :
                    exptime = fits_image.header['EXPTIME'] if 'EXPTIME' in fits_image.header else 1000.
                    self.individual_filter = filter_lite(fits_image.image_data[:,:,which_filter], exptime=exptime)
                    fits_image.filters[which_filter_str] = self.individual_filter
                else :
                    self.individual_filter = fits_image.filters[which_filter_str]
            else :
                print('Filter ' + which_filter + " has not been imported. Import individual filters with fits_image.load_filters(). Filter fits files have to be in a folder named 'filters' placed in the same directory as your main RGB image.")
                raise ValueError('Filter ' + which_filter + " has not been imported. Import individual filters with fits_image.load_filters(). Filter fits files have to be in a folder named 'filters' placed in the same directory as your main RGB image.")
        elif fits_image.filters is not None and not all(i in [0, 'red', 'green', 'blue'] for i in fits_image.filters) :
            filter_list = list(fits_image.filters.keys())
            filter_list[0] = [filter_list[0]]
            filter_list_str = str(filter_list)
            filter_list_str = filter_list_str.replace(", ", "/")
            filter_list_str = filter_list_str.replace("'", "")
            which_filter = input('Select filter: ' + filter_list_str)
            which_filter = filter_list[0][0] if which_filter=='' else which_filter
            print('Using filter ' + which_filter)
            self.individual_filter = fits_image.filters[which_filter]
        else :
            print('\n------------------')
            print('\nUsing main image as individual filter. \nCareful: image simulator needs flux units when calculating the residuals. Make sure the main image has such units when running the optimization, or import individual filters separately with fits_image.load_filters()')
            which_filter = 'red' if len(fits_image.image_data.shape)==3 else 0
            if fits_image.filters is not None :
                if which_filter in fits_image.filters :
                    self.individual_filter = fits_image.filters[which_filter]
            else :
                exptime = fits_image.header['EXPTIME'] if 'EXPTIME' in fits_image.header else 1000.
                self.individual_filter = filter_lite(fits_image.image_data[:,:,0], exptime=exptime) if which_filter=='red'\
                                         else filter_lite(fits_image.image_data, exptime=exptime)
                fits_image.filters = {which_filter: self.individual_filter}
        
        #-------------- ImageData --------------#
        if not hasattr(self.individual_filter, 'rms') :
            print('Calculating RMS...')
            self.individual_filter.rms = np.std(self.individual_filter.image_data)
            print('done: ' + str(self.individual_filter.rms))
            if not hasattr(self.individual_filter, 'wcs') : # test if is instance of filter_lite instead of full filter class
                fits_image.rms = self.individual_filter.rms
        if 'EXPTIME' not in self.individual_filter.header :
            print('\nCareful: EXPTIME not found in header. Using arbitrary value EXPTIME=1ks.')
            exptime = 1000.
        else :
            exptime = self.individual_filter.header['EXPTIME']

        self.image_data = self.individual_filter.image_data[ self._crop_y0 - self._crop_npix : self._crop_y0, self._crop_x0 : self._crop_x0 + self._crop_npix ]
        self.ImageData_kwargs = {'image_data': self.image_data,
                                'background_rms': self.individual_filter.rms,
                                'exposure_time': exptime,
                                'transform_pix2angle': transform_pix2angle,
                                'ra_at_xy_0': - square_size_arcsec / 2,
                                'dec_at_xy_0': - square_size_arcsec / 2}
        self.ImageData = ImageData(**self.ImageData_kwargs)
        
        #-------------- PSF --------------#
        if self.individual_filter.psf is not None :
            self.PSF_kwargs = {'psf_type': 'PIXEL',
                                            'kernel_point_source': self.individual_filter.psf.data,
                                            #'truncation': 35,
                                            #'point_source_supersampling_factor': 1,
                                            'pixel_size': fits_image.pix_deg_scale*3600}
        else :
            self.PSF_kwargs = {'psf_type': 'GAUSSIAN',
                                            'fwhm': fits_image.pix_deg_scale*3600*2,
                                            'truncation': 5,
                                            'pixel_size': fits_image.pix_deg_scale*3600}
        self.PSF = PSF(**self.PSF_kwargs)
        
        self.kwargs_numerics = {'supersampling_factor': 8, #ideally, supersampling_factor=16 for light source model, but 8 is ok. Doesn't matter for point source model.
                                'supersampling_convolution': True}
        
        self.sigma = 0.00036

        #-------------- Image plane: observed / simulated / residual / RGB (linked zoom & pan) --------------#
        self.image_plane_observed = pg.ImageView()
        self.image_plane_simulated = pg.ImageView()
        self.image_plane_residual = pg.ImageView()
        self.image_plane_rgb = pg.ImageView()
        self.image_plane_plot = self.image_plane_simulated  # backward compatibility

        sim0 = np.zeros_like(self.image_data)
        rgb_data = fits_image.image_data[ self._crop_y0 - self._crop_npix : self._crop_y0, self._crop_x0 : self._crop_x0 + self._crop_npix ]
        self.simulated_image = sim0
        self.residual_image = self.image_data - sim0
        self.rgb_image_data = rgb_data
        self.image_plane_observed.setImage(self.image_data[::-1, :])
        self.image_plane_simulated.setImage(sim0[::-1, :])
        self.image_plane_residual.setImage(self.residual_image[::-1, :])
        self.image_plane_rgb.setImage(rgb_data[::-1, ...])

        x_cc = (self._lt_curve_coords_relative_broken[0] - self._SquareOfInterest_xr_bottomleft) / (fits_image.pix_deg_scale * 3600)
        y_cc = self.image_data.shape[0] - (self._lt_curve_coords_relative_broken[1] - self._SquareOfInterest_yr_bottomleft) / (fits_image.pix_deg_scale * 3600)
        pen_cc = pg.mkPen(color=[255, 0, 255, 255], width=4.0001)
        self.critical_curve_plots = []
        for iv in (self.image_plane_observed, self.image_plane_simulated, self.image_plane_residual, self.image_plane_rgb):
            cc = pg.PlotDataItem()
            cc.setPen(pen_cc)
            cc.setData(x_cc, y_cc)
            iv.addItem(cc)
            self.critical_curve_plots.append(cc)

            mult_plot = pg.ScatterPlotItem()
            mult_plot.setData(np.array(XR)/fits_image.pix_deg_scale/3600 + self._crop_npix/2, self._crop_npix/2 - np.array(YR)/fits_image.pix_deg_scale/3600)
            iv.addItem(mult_plot)
        self.critical_curve_plot = self.critical_curve_plots[1]

        vb0 = self.image_plane_observed.getView()
        vb1 = self.image_plane_simulated.getView()
        vb2 = self.image_plane_residual.getView()
        vb3 = self.image_plane_rgb.getView()
        vb1.setXLink(vb0)
        vb1.setYLink(vb0)
        vb2.setXLink(vb0)
        vb2.setYLink(vb0)
        vb3.setXLink(vb0)
        vb3.setYLink(vb0)

        # Fixed proportions: top half = Simulated | Observed (each half width); bottom half = Residual (full width).
        # Outer horizontal split (source | image plane) remains adjustable via QSplitter.
        implane_v = QVBoxLayout(self._imsim_implane_layout)
        implane_v.setContentsMargins(0, 0, 0, 0)
        implane_v.setSpacing(2)

        top_row = QWidget()
        top_h = QHBoxLayout(top_row)
        top_h.setContentsMargins(0, 0, 0, 0)
        top_h.setSpacing(2)

        sim_col = QWidget()
        sim_v = QVBoxLayout(sim_col)
        sim_v.setContentsMargins(0, 0, 0, 0)
        sim_v.setSpacing(2)
        sim_v.addWidget(QtWidgets.QLabel("Simulated"))
        sim_v.addWidget(self.image_plane_simulated, stretch=1)

        obs_col = QWidget()
        obs_v = QVBoxLayout(obs_col)
        obs_v.setContentsMargins(0, 0, 0, 0)
        obs_v.setSpacing(2)
        obs_v.addWidget(QtWidgets.QLabel("Observed"))
        obs_v.addWidget(self.image_plane_observed, stretch=1)

        top_h.addWidget(sim_col, stretch=1)
        top_h.addWidget(obs_col, stretch=1)

        bottom_row = QWidget()
        bottom_h = QHBoxLayout(bottom_row)
        bottom_h.setContentsMargins(0, 0, 0, 0)
        bottom_h.setSpacing(2)

        res_col = QWidget()
        res_v = QVBoxLayout(res_col)
        res_v.setContentsMargins(0, 0, 0, 0)
        res_v.setSpacing(2)
        res_v.addWidget(QtWidgets.QLabel("Residual (obs − sim)"))
        res_v.addWidget(self.image_plane_residual, stretch=1)

        rgb_col = QWidget()
        rgb_v = QVBoxLayout(rgb_col)
        rgb_v.setContentsMargins(0, 0, 0, 0)
        rgb_v.setSpacing(2)
        rgb_v.addWidget(QtWidgets.QLabel("RGB"))
        rgb_v.addWidget(self.image_plane_rgb, stretch=1)

        bottom_h.addWidget(res_col, stretch=1)
        bottom_h.addWidget(rgb_col, stretch=1)

        implane_v.addWidget(top_row, stretch=1)
        implane_v.addWidget(bottom_row, stretch=1)

        self.image_plane_views = {
            'observed': self.image_plane_observed,
            'simulated': self.image_plane_simulated,
            'residual': self.image_plane_residual,
            'rgb': self.image_plane_rgb,
        }

        self.filter_image = ImageFilter(self)
        for _iv in (self.image_plane_observed, self.image_plane_simulated, self.image_plane_residual, self.image_plane_rgb):
            _iv.installEventFilter(self.filter_image)

        self.window.show()
        
        # Initial width ratio source:image = 3:4
        x_left, x_right = self._imsim_layout.sizes()
        total_size = x_left + x_right
        ratio = 3/7
        x_left_new = round(total_size * ratio)
        x_right_new = total_size - x_left_new
        self._imsim_layout.setSizes([x_left_new, x_right_new])

        # Initialize image-plane y-range to full image height for all linked views.
        vb0.setYRange(0, self.image_data.shape[0], padding=0)
        
        #fits_image.lt.lt.set_field(self._initial_field)
        self.previous_state_current_ROI = self.ROI.getState()
    
    
    def simulate(self, N_sigma=6.) :
        self.lm_dict_opt = make_lm_dict_opt(self, N_sigma=N_sigma)
        self.lm_current = lenstronomy_model(self.lm_dict_opt, self)
        simulate(self)
    
    def optimize(self, N_sigma=6., fitting_kwargs_list=None) :
        self.simulate(N_sigma=6.)
        optimize(self, fitting_kwargs_list=fitting_kwargs_list)
    
    #def save(self) :
    #    self.model_local = format_lm_local(self.models, self.result_kwargs)
    #    self.source_model.save( path=os.path.join(self.fits_image.lt.model_dir, "lenstronomy_model.pkl") )
        
    def load(self, path=None) :
        if path is None :
            self.lm_imported = lenstronomy_model( os.path.join(self.fits_image.lt.model_dir, "lenstronomy_model.pkl"), self )
        else :
            self.lm_imported = lenstronomy_model( path, self )
        
        #self.imported_models = copy.deepcopy(self.models)
        #self.imported_models = {'lens_model_list': self.LensModel_list}
        #self.imported_models['source_light_model_list'] = self.imported_source_model['models']['source_light_model_list']
        
        #self.imported_kwargs = {'kwargs_lens': self.LensModel_kwargs}
        #self.imported_kwargs['kwargs_source'] = self.imported_source_model['results']['kwargs_source']
        




class filter_lite() :
    def __init__(self, image_data, exptime=1000.) :
        self.image_data = image_data
        self.header = {'EXPTIME': exptime}
        self.psf = None


#class SpecialWindow(QMainWindow):
#    def __init__(self, lt, parent=None):
#        super().__init__(parent)
#        #self.lt = lt
#        
#    def closeEvent(self, event):
#        print("lenstool instance not re-initialized")
#        #self.lt.lt = lenstool.Lenstool( self.lt._FileToUse )
#        #self.lt.lt.set_field(self.lt.imsim._initial_field)
#        #event.accept()



