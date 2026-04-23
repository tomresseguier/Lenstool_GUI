import copy
import numpy as np
from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.ImSim.image_linear_solve import ImageLinearFit
from lenstronomy.LensModel.lens_model_extensions import LensModelExtensions

from ...utils.utils_astro.utils_general import world_to_relative





def _fuse_curve_segments(x_segments, y_segments):
    x_fused = []
    y_fused = []
    for x_seg, y_seg in zip(x_segments, y_segments):
        x_fused.extend(np.asarray(x_seg).tolist() + [np.nan])
        y_fused.extend(np.asarray(y_seg).tolist() + [np.nan])
    return np.array(x_fused), np.array(y_fused)


def _update_curves(imsim):
    lens_model_extensions = LensModelExtensions(imsim.LensModel)
    xr_crit, yr_crit, xr_caustic, yr_caustic = lens_model_extensions.critical_curve_caustics(
        imsim.LensModel_kwargs,
        compute_window=imsim._SquareOfInterest_side_arcsec,
        grid_scale=getattr(imsim, '_curve_grid_scale', imsim._SquareOfInterest_side_arcsec / 100.0),
        center_x=0.0,
        center_y=0.0,
    )
    xr_crit_fused, yr_crit_fused = _fuse_curve_segments(xr_crit, yr_crit)
    xr_caustic_fused, yr_caustic_fused = _fuse_curve_segments(xr_caustic, yr_caustic)

    #imsim.source_center_coordinates = imsim.LensModel.ray_shooting(0.0, 0.0, imsim.LensModel_kwargs)
    x_caustic = xr_caustic_fused - imsim.source_center_coordinates[0]
    y_caustic = yr_caustic_fused - imsim.source_center_coordinates[1]
    imsim.caustic_plot.setData(x_caustic, y_caustic)

    pix_scale = imsim.fits_image.pix_deg_scale * 3600.0
    x_crit = (xr_crit_fused + imsim._SquareOfInterest_side_arcsec / 2.0) / pix_scale
    y_crit = (yr_crit_fused + imsim._SquareOfInterest_side_arcsec / 2.0) / pix_scale
    for cc in imsim.critical_curve_plots :
        cc.setData(x_crit, y_crit)

    current_lens_kwargs = copy.deepcopy(imsim.LensModel_kwargs)
    imsim._last_curve_lens_kwargs = current_lens_kwargs

    SX = []
    SY = []
    XR = []
    YR = []
    for mult in imsim.fits_image.lt.mult.cat :
        xr, yr = world_to_relative( mult['ra'], mult['dec'], imsim.center_world )
        XR.append(xr)
        YR.append(yr)
        sx, sy = imsim.LensModel.ray_shooting(xr, yr, imsim.LensModel_kwargs)
        SX.append(sx - imsim.source_center_coordinates[0])
        SY.append(sy - imsim.source_center_coordinates[1])
    imsim.source_coord_plot.setData(SX, SY)


def _kwargs_are_same(kwargs1, kwargs2) :
    if len(kwargs1) != len(kwargs2) :
        return False
    for kwargs_dict1, kwargs_dict2 in zip(kwargs1, kwargs2) :
        for key, value in kwargs_dict1.items() :
            if key not in kwargs_dict2 :
                return False
            if type(kwargs_dict2[key]) != type(value) :
                return False
            if type(kwargs_dict2[key]) == np.ndarray :
                if not np.array_equal(kwargs_dict2[key], value) :
                    return False
            else :
                if kwargs_dict2[key] != value :
                    return False
    return True


def simulate(imsim) :
    lm_dict_opt = imsim.lm_dict_opt

    if not _kwargs_are_same(imsim.LensModel_kwargs, imsim._last_curve_lens_kwargs) :
        _update_curves(imsim)
    
    ######### Plot the simulated image #########
    LightModel_source = LightModel(light_model_list=lm_dict_opt['models']['source_light_model_list'])
    # fixed_magnifications is weird wording but we want it here because we are in source plane
    PSModel_source = PointSource(point_source_type_list=lm_dict_opt['models']['point_source_model_list'], fixed_magnification_list=[True for src in lm_dict_opt['models']['point_source_model_list']])
    #PSModel_source = PointSource(point_source_type_list=lm_dict_opt['models']['point_source_model_list'], lens_model=imsim.LensModel, fixed_magnification_list=[True for src in lm_dict_opt['models']['point_source_model_list']])
    
    imsim.ImageModel = ImageModel( data_class=imsim.PixelGrid, psf_class=imsim.PSF,
                                    lens_model_class=imsim.LensModel,
                                    source_model_class=LightModel_source,
                                    point_source_class=PSModel_source,
                                    #lens_light_model_class=,
                                    kwargs_numerics=imsim.kwargs_numerics )
    
    imsim.ImageLinearFit = ImageLinearFit( data_class=imsim.ImageData, psf_class=imsim.PSF,
                                            lens_model_class=imsim.LensModel,
                                            source_model_class=LightModel_source,
                                            point_source_class=PSModel_source,
                                            #lens_light_model_class=,
                                            kwargs_numerics=imsim.kwargs_numerics )
    
    print('Start simulating image')
    imsim.simulated_image = imsim.ImageModel.image(kwargs_source=lm_dict_opt['kwargs']['kwargs_source'],
                                                    kwargs_ps=lm_dict_opt['kwargs']['kwargs_source_ps'],
                                                    #kwargs_lens_light=kwargs_light_lens,
                                                    kwargs_lens=imsim.LensModel_kwargs, unconvolved=False)

    # Preserve current linked image-plane view range when refreshing images.
    vb_ref = imsim.image_plane_observed.getView()
    x_range, y_range = vb_ref.viewRange()

    imsim.image_plane_simulated.setImage(imsim.simulated_image, autoRange=False)
    residual = imsim.image_data - imsim.simulated_image
    imsim.residual_image = residual
    imsim.image_plane_residual.setImage(residual, autoRange=False)
    
    # Re-apply the original range once updates are done.
    vb_ref.setXRange(x_range[0], x_range[1], padding=0)
    vb_ref.setYRange(y_range[0], y_range[1], padding=0)
    print('done')
    
    #solver = LensEquationSolver(imsim.lens_model)