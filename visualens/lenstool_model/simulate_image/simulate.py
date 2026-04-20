from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.ImSim.image_linear_solve import ImageLinearFit






def simulate(imsim) :
    lm_dict_opt = imsim.lm_dict_opt
    
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
    
    ######### Plot the critical curve #########
    x = (imsim._lt_curve_coords_relative_broken[0] - imsim._SquareOfInterest_xr_bottomleft) / (imsim.fits_image.pix_deg_scale*3600)
    y = imsim.simulated_image.shape[0] - (imsim._lt_curve_coords_relative_broken[1] - imsim._SquareOfInterest_yr_bottomleft) / (imsim.fits_image.pix_deg_scale*3600)
    for cc in getattr(imsim, 'critical_curve_plots', [imsim.critical_curve_plot]):
        cc.setData(x, y)
    # Re-apply the original range once updates are done.
    vb_ref.setXRange(x_range[0], x_range[1], padding=0)
    vb_ref.setYRange(y_range[0], y_range[1], padding=0)
    print('done')
    
    #solver = LensEquationSolver(imsim.lens_model)