import os
import numpy as np
import corner

from lenstronomy.Workflow.fitting_sequence import FittingSequence
from lenstronomy.Plots import chain_plot

from .simulate import simulate
from .utils import make_samples_dict
from .lenstronomy_model import lenstronomy_model





def is_mpi_running():
    # Common environment variables set by different MPI implementations
    mpi_vars = ['OMPI_COMM_WORLD_SIZE', 'PMI_SIZE', 'PMIX_SIZE']
    return any(var in os.environ for var in mpi_vars)

def optimize(imsim, fitting_kwargs_list=None) :
    print('Starting optimization')
    lm_dict_opt = imsim.lm_dict_opt
                    
    kwargs_lens_empty = [{} for _ in imsim.LensModel_list]
    lens_params = [kwargs_lens_empty, kwargs_lens_empty, imsim.LensModel_kwargs, kwargs_lens_empty, kwargs_lens_empty]
    source_params = [lm_dict_opt['kwargs_init']['kwargs_source'], lm_dict_opt['kwargs_sigma']['kwargs_source'], lm_dict_opt['kwargs_fixed']['kwargs_source'], lm_dict_opt['kwargs_lower']['kwargs_source'], lm_dict_opt['kwargs_upper']['kwargs_source']]
    ps_params = [lm_dict_opt['kwargs_init']['kwargs_source_ps'], lm_dict_opt['kwargs_sigma']['kwargs_source_ps'], lm_dict_opt['kwargs_fixed']['kwargs_source_ps'], lm_dict_opt['kwargs_lower']['kwargs_source_ps'], lm_dict_opt['kwargs_upper']['kwargs_source_ps']]
    
    imsim.optimization_params = {'lens_model': lens_params,
                                    #'lens_light_model': lens_light_params,
                                    'source_model': source_params,
                                    'point_source_model': ps_params}
    
    kwargs_likelihood = {'source_marg': False}
    single_band = [[imsim.ImageData_kwargs, imsim.PSF_kwargs, imsim.kwargs_numerics]]
    kwargs_data_joint = {'multi_band_list': single_band, 'multi_band_type': 'multi-linear'}
    kwargs_constraints = {'linear_solver': True} if imsim.use_linear_solver else {'linear_solver': False}
    
    kwargs_model = lm_dict_opt['models'].copy()
    kwargs_model['fixed_magnification_list'] = [ True for src in lm_dict_opt['models']['point_source_model_list'] ]
    
    fitting_seq = FittingSequence(kwargs_data_joint, kwargs_model, kwargs_constraints, kwargs_likelihood, imsim.optimization_params, mpi=is_mpi_running())
    
    if fitting_kwargs_list is None :
        fitting_kwargs_list = [['PSO', {'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100, 'threadCount': 16}],
                               ['MCMC', {'n_burn': 100, 'n_run': 100, 'n_walkers': 50, 'sigma_scale': .01, 'threadCount': 16}]]
        
    imsim.fitting_seq = fitting_seq
    imsim.chain_list = fitting_seq.fit_sequence(fitting_kwargs_list)
    imsim.result_kwargs = fitting_seq.best_fit()
    
    
    #for i in range(len(chain_list)):
    chain_plot.plot_chain_list(imsim.chain_list, 0)
    
    if False :
        sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = imsim.chain_list[1]
        print("number of non-linear parameters in the MCMC process: ", len(param_mcmc))
        print("parameters in order: ", param_mcmc)
        print("number of evaluations in the MCMC process: ", np.shape(samples_mcmc)[0])
        n_sample = len(samples_mcmc)
        print(n_sample)
        samples_mcmc_cut = samples_mcmc[int(n_sample*1/2.):]
        n, num_param = np.shape(samples_mcmc_cut)
        plot = corner.corner(samples_mcmc_cut[:,:], labels=param_mcmc[:], show_titles=True)
    
    
    
    lm_dict_opt['kwargs_opt'] = {'kwargs_source': imsim.result_kwargs['kwargs_source'], 'kwargs_source_ps': imsim.result_kwargs['kwargs_ps']}
        
    
    for i, samples in enumerate(imsim.chain_list) :
        if fitting_kwargs_list[i][0] == 'MCMC' :
            lm_dict_opt['kwargs_MCMC'] = make_samples_dict(lm_dict_opt, samples)
        #elif fitting_kwargs_list[i][0] == 'PSO' :
        #    lm_dict_opt['kwargs_PSO'] = None #format_samples(samples)
    
    
    if imsim.use_linear_solver :
        """ This does the linear solve necessary to find the amplitudes """
        imsim.ImageLinearFit.image_linear_solve(kwargs_lens=imsim.LensModel_kwargs,
                                                    kwargs_source=lm_dict_opt['kwargs_opt']['kwargs_source'],
                                                    kwargs_lens_light=None,
                                                    kwargs_ps=lm_dict_opt['kwargs_opt']['kwargs_source_ps'],
                                                    inv_bool=False)
    
    imsim.lm_optimized = lenstronomy_model(lm_dict_opt, imsim, chain_list=imsim.chain_list)
    imsim.lm_optimized.save('lm_last_optimization')
    #imsim.lm_optimized.plot()
    #imsim.lm_optimized.save('./last_optimized_lm.pkl')