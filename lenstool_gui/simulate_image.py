import numpy as np
import pyqtgraph as pg
from PyQt5.QtWidgets import QMainWindow, QSplitter
from PyQt5.QtCore import Qt
from pyqtgraph.Qt import QtWidgets

from lenstronomy.Data.pixel_grid import PixelGrid     
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.Data.psf import PSF
#from lenstronomy.Plots import lens_plot
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.Util import util

from .utils.utils_astro.utils_general import world_to_relative
from .utils.utils_general.sort_points import break_curves
from .utils.utils_simulate_image import SourceFilter, ImageFilter
from .utils.utils_Qt.drag_widgets import DragPlotWidget_special
from .utils.utils_Qt.utils_general import transform_rectangle



class lenstronomy_model :
    def __init__(self, fits_image) :
        self.fits_image = fits_image
        self.ROI = fits_image.image_widget.current_ROI
        #self.fits_image.qt_image.addItem(ROI)
        
        if not hasattr(self, 'previous_state_current_ROI') :
            self.previous_state_current_ROI = None
        
        x0 = self.ROI.getState()['pos'][0]
        y0 = self.ROI.getState()['pos'][1]
        a = self.ROI.getState()['size'][0]
        b = self.ROI.getState()['size'][1]
        angle = self.ROI.getState()['angle'] *np.pi/180
        
        x0, y0, a, b, angle = transform_rectangle(x0, y0, a, b, angle) #x0, y0 at the top left
        size_y = fits_image.image_data.shape[0]
        y0 = size_y-y0 #Counting pixels from bottom instead of top
        self.anchor = (x0, y0)
        
        ra_topleft, dec_topleft = fits_image.image_to_world(x0, y0)
        ra_center, dec_center = fits_image.image_to_world(x0 + a/2, y0 - b/2)
        self.center_world = ra_center, dec_center
        
        xr_topleft, yr_topleft = fits_image.lt.world_to_relative(ra_topleft, dec_topleft)
        xr_center, yr_center = fits_image.lt.world_to_relative(ra_center, dec_center)
        
        offset_arcsec_x = xr_topleft - xr_center
        offset_arcsec_y = - (yr_topleft - yr_center) #because ROI's anchor corner is top left instead of bottom left
        
        offset_arcsec = np.max(np.abs( [offset_arcsec_x, offset_arcsec_y] ))
        self._SquareOfInterest_side_arcsec = offset_arcsec * 2
        self._SquareOfInterest_xr_bottomleft = xr_center - offset_arcsec
        self._SquareOfInterest_yr_bottomleft = yr_center - offset_arcsec
        new_field = [xr_center-offset_arcsec, xr_center+offset_arcsec, yr_center-offset_arcsec, yr_center+offset_arcsec]
        
        numPix = 100
        initial_field = fits_image.lt.lt.get_field([])
        fits_image.lt.lt.set_field(new_field)
        if self.previous_state_current_ROI is None or self.previous_state_current_ROI != self.ROI.getState() or fits_image.lt.dpl_maps[6.2][0].shape[0] != numPix :
            fits_image.lt.compute_lt_dpl(npix=numPix)
        
        #-------------- Lens model Lenstronomy definitions --------------#
        
        deltaPix = abs( fits_image.lt.dpl_maps[6.2][2].wcs.cdelt[0]*3600 )
        x_grid_interp, y_grid_interp = util.make_grid(numPix, deltaPix)
        x_axes, y_axes = util.get_axes(x_grid_interp, y_grid_interp)
        
        self.LensModel_kwargs = [{'grid_interp_x': x_axes,
                                            'grid_interp_y': y_axes,
                                            'f_x': fits_image.lt.dpl_maps[6.2][0],
                                            'f_y': fits_image.lt.dpl_maps[6.2][1]}]
        self.LensModel_list = ['INTERPOL']
        self.LensModel = LensModel(lens_model_list=self.LensModel_list)
        
        #m=4
        #f, ax = plt.subplots(1, 1, figsize=(10, 10), sharex=False, sharey=False)
        #lens_plot.lens_model_plot(ax, lensModel=self.lens_model, kwargs_lens=kwargs_lens,
        #                          with_caustics=True, fast_caustic=True, coord_inverse=False, numPix=round(numPix/m), deltaPix=deltaPix*m)
        
        #-------------- Qt graphics definitions --------------#
        
        self.source_plane_widget = DragPlotWidget_special()
        self.source_plane_widget.setTitle('Source plane')
        #source_plane_widget.setAspectLocked(lock=True, ratio=1)
        
        self.to_plot = fits_image.image_data[round(y0-b):round(y0),round(x0):round(x0+a),:][::-1,:,:]
        self.image_plane_plot = pg.ImageView()
        self.image_plane_plot.setImage(self.to_plot)        
        
        self._imsim_srcplane_layout = QSplitter(Qt.Vertical)
        self._imsim_implane_layout = QSplitter(Qt.Vertical)
        
        self._imsim_srcplane_layout.addWidget(self.source_plane_widget)
        label = QtWidgets.QLabel("Image plane")
        self._imsim_implane_layout.addWidget(label)
        self._imsim_implane_layout.addWidget(self.image_plane_plot)
        
        self._imsim_layout = QSplitter(Qt.Horizontal)
        self._imsim_layout.addWidget(self._imsim_srcplane_layout)
        self._imsim_layout.addWidget(self._imsim_implane_layout)
        
        self.window = QMainWindow()
        self.window.setWindowTitle('Lenstronomy image simulator')
        self.window.setCentralWidget(self._imsim_layout)
        self.window.show()
        
        
        #-------------- Plot critical & caustic curves --------------#

        self.source_center_coordinates = self.LensModel.ray_shooting(0, 0, self.LensModel_kwargs)
        
        if self.previous_state_current_ROI is None or self.previous_state_current_ROI != ROI.getState() :
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
        
        self.critical_curve_plot = pg.PlotDataItem()
        self.critical_curve_plot.setPen( color=[255,0,255,255], width=4.0001 )
        self._lt_curve_coords_relative_broken = break_curves(fits_image.lt.lt_curve_coords_relative, distance_threshold=8.*fits_image.pix_deg_scale*3600) #sort_points(self.lt_curve_coords_relative, distance_threshold=1.0, angle_threshold=np.pi)
        x = (self._lt_curve_coords_relative_broken[0] - xr_topleft) / (fits_image.pix_deg_scale*3600)
        y = - (self._lt_curve_coords_relative_broken[1] - yr_topleft) / (fits_image.pix_deg_scale*3600)
        
        self.critical_curve_plot.setData(x, y)
        self.image_plane_plot.addItem(self.critical_curve_plot)
        
        self.filter_source = SourceFilter(self)
        self.source_plane_widget.installEventFilter(self.filter_source)
        
        self.filter_image = ImageFilter(self)
        self.image_plane_plot.installEventFilter(self.filter_image)
        
        
        #-------------- Plot source positions from multiple image catalog --------------#
        SX = []
        SY = []
        for mult in fits_image.lt.mult.cat :
            xr, yr = world_to_relative( mult['ra'], mult['dec'], (ra_center, dec_center) )
            sx, sy = self.LensModel.ray_shooting(xr, yr, self.LensModel_kwargs)
            SX.append(sx - self.source_center_coordinates[0])
            SY.append(sy - self.source_center_coordinates[1])
        
        source_coord_plot = pg.ScatterPlotItem()
        source_coord_plot.setData(SX, SY)
        self.source_plane_widget.addItem(source_coord_plot)
        
            
        
        
        
        #-------------- Lenstronomy definitions --------------#
        
        #-------------- PixelGrid --------------#
        transform_pix2angle = np.array([[1, 0], [0, 1]]) * fits_image.pix_deg_scale*3600
        
        npix = round(self._SquareOfInterest_side_arcsec / (fits_image.pix_deg_scale*3600))
        self.PixelGrid_kwargs = {'nx': npix,
                                        'ny': npix,  # number of pixels per axis
                                        'ra_at_xy_0': -offset_arcsec,
                                        'dec_at_xy_0': -offset_arcsec,
                                        'transform_pix2angle': transform_pix2angle} 
        self.PixelGrid = PixelGrid(**self.PixelGrid_kwargs)
        
        #-------------- DATA --------------#
        if False :
        #if not hasattr(self, 'filters_added') :
            if fits_image.filters is not None :
                names = np.array([name for name in fits_image.filters])
                #JWST_names = ['F070W','F090W','F115W','F140M','F150W','F150W2','F162M','F164N','F182M','F200W','F210M','F250M','F277W','F300M','F322W2','F323N','F335M','F356W','F360M','F405N','F410M','F430M','F444W','F460M','F466N','F470N','F480M']
                #mask = np.array([name in JWST_names for name in names])
                #names = names[mask]
                names = ['F200W']
                i = 0
                self._exposure_time = 0.
                data_shape = fits_image.filters[names[i]].image_data.shape
                self.filters_added = np.zeros(data_shape)
                same_sized_filters = []
                for name in names :
                    if fits_image.filters[name].image_data.shape == data_shape :
                        print(f'Adding filter {name}.')
                        self.filters_added += fits_image.filters[name].image_data #need correct units here
                        self._exposure_time += fits_image.filters[name].header['EXPTIME']
                        same_sized_filters.append(name)
                if fits_image.filters[names[i]].psf is not None :
                    psf_shape = fits_image.filters[names[i]].psf.data.shape
                    self.psf_added = np.zeros(psf_shape)
                    for name in same_sized_filters :
                        print(f'Adding psf {name}.')
                        self.psf_added += fits_image.filters[name].psf.data
        self.individual_filter = fits_image.filters['F150W']
        #-------------- ImageData --------------#
        if not hasattr(self, '_background_rms') :
            print('Calculating RMS...')
            self._background_rms = 0.23358145#np.std(self.individual_filter.image_data)
        self.ImageData_kwargs = {'image_data': self.individual_filter.image_data[round(y0)-npix:round(y0),round(x0):round(x0)+npix],
                                            'background_rms': self._background_rms,
                                            'exposure_time': self.individual_filter.header['EXPTIME'],
                                            'transform_pix2angle': transform_pix2angle,
                                            'ra_at_xy_0': -offset_arcsec,
                                            'dec_at_xy_0': -offset_arcsec}
        self.ImageData = ImageData(**self.ImageData_kwargs)
        
        #-------------- PSF --------------#
        if True :
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
    
        
        
        fits_image.lt.lt.set_field(initial_field)
        self.previous_state_current_ROI = self.ROI.getState()