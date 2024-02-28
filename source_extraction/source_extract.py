"""
source_extract.py from David Harvey's pyRRG
"""

#from pyextract import pysex
from source_extraction.pysex import pysex
from astropy.io import fits
from astropy.table import Table

import numpy as np
import os as os
from source_extraction.match_cat import run_match
from numpy.lib.recfunctions import append_fields as append_rec

module_dir = os.path.dirname(os.path.abspath(__file__))




def source_extract( image_path, weight_path=None, zero_point=None,
                    out_dir='PWD', outfile_name='SExtractor_cat.fits',
                    return_sources=True, rerun=False) :
    '''
    Given that source extration is a difficult thing for the photometry
    I will use the Mathilde mehtod, in Python to do a wrapper
    for source_extractor.
    
    The process is the same as in mathilde_setract.pro
    THe weight file is optional however should be used.
    zero_point is the zero point of the iamge, if not used, it
    will guess it from the header of the file and assume it is ACS

    This will run two runs of sextractor, cold and hot, where cold only finds
    big things and hot finds small things

    Then it combines the two catalogues
    
    '''
    
    if out_dir == 'PWD':
        out_dir = os.getcwd()
    out_path = out_dir + '/' + outfile_name
    if os.path.isfile(out_path) and not rerun :
        print('Previous SExtractor catalog found.')
        hot_sources = Table( fits.open(out_path)[1].data )
    else :
        config_dir = os.path.join(module_dir, 'SExtractor_config/from_pyRRG')
        check_sex_files(config_dir)
        
        header = fits.open( image_path )[0].header
        findPhot = np.array(['PHOTFLAM' in key for key in header.keys()])
        if np.all(findPhot == False) :
            header = fits.open( image_path )[1].header
        
        ### FIX THE PHOTOMETRY HERE!!! AND ADD WEIGHT FILE OPTION ###
        zero_point = 0. #acs_zero_point(header)
            
        conf_args = {'MAG_ZEROPOINT': zero_point,
                     'WEIGHT_TYPE': 'NONE',
                     'PARAMETERS_NAME': config_dir + '/rrg.param',
                     'STARNNW_NAME': config_dir + '/default.nnw',
                     'FILTER_NAME': config_dir + '/gauss_5.0_9x9.conv'}
        
        
        #F## COLD RUN ###
        cold_conf = config_dir + '/HFF_cold.param'
        cold_sources = pysex.run( image_path, \
                                  conf_file=cold_conf, \
                                  conf_args=conf_args, \
                                  param_file=config_dir+'/rrg.param')
        
        cold_sources = append_fits_field( cold_sources, 'RA', cold_sources['X_WORLD'])
        cold_sources = append_fits_field( cold_sources, 'DEC', cold_sources['Y_WORLD'])
        
        
    
        #Second hot 
        hot_conf = config_dir+'/HFF_hot.param'
        hot_sources = pysex.run( image_path, \
                                   conf_file=hot_conf, \
                                   conf_args=conf_args, \
                                   param_file=config_dir+'/rrg.param')
    
    
        hot_sources = append_fits_field( hot_sources, 'RA', hot_sources['X_WORLD'])
        hot_sources = append_fits_field( hot_sources, 'DEC', hot_sources['Y_WORLD'])
        
        #The NYMBER is a weird thing
        
        hot_sources['NUMBER'] = np.arange( len(hot_sources['NUMBER'])) +1
        cold_sources['NUMBER'] = np.arange( len(cold_sources['NUMBER'])) +1
        
        fits.writeto( 'cold_sources.fits', cold_sources, overwrite=True )
        fits.writeto( 'hot_sources.fits', hot_sources, overwrite=True )
    
    
        print('Matching cold and hot sources:')
        print( str(len(cold_sources)) + ' cold sources' )
        print( str(len(hot_sources)) + ' hot sources' )
        print( 'TOTAL: ' + str(len(cold_sources) + len(hot_sources)) + ' detections' )
        matched_sources = run_match( 'cold_sources.fits',
                                     'hot_sources.fits' )
        
        for iField in hot_sources.columns.names:
            hot_sources[iField][ matched_sources[1].data['NUMBER_2'] -1 ] = cold_sources[iField][ matched_sources[1].data['NUMBER_1'] - 1]
        
        print('Matching cold and hot sources:')
        print('Number of sources in matched catalog: ' + str(len(hot_sources)) + ' sources')
        
        #Need to retain numbering for bookeepin purposes
        hot_sources['NUMBER'] = np.arange(len(hot_sources))
        
        fits.writeto( out_dir + '/' + outfile_name, hot_sources, overwrite=True )
    if return_sources:
        return hot_sources
    
    
def acs_zero_point( header ):
    zpt = -2.5*np.log10(header['PHOTFLAM']) + header['PHOTZPT'] - 5.0*np.log10(header['PHOTPLAM'])+18.6921
    return zpt

def check_sex_files( config_dir ):
    #Check all the sex files to see if they exist
    if (not os.path.isfile(config_dir+'/HFF_hot.param')):
        raise ValueError('HFF_hot.param not found at ' + config_dir+'/HFF_hot.param')
    if (not os.path.isfile(config_dir+'/HFF_cold.param')):
        raise ValueError('HFF_cold.param not found at ' + config_dir+'/HFF_cold.param')
    if (not os.path.isfile(config_dir+'/rrg.param')):
        raise ValueError('rrg.param not found at ' + config_dir+'/rrg.param')


    
def append_fits_field( fits_array, name, array, format='D'):
    
    cols = [] 
    cols.append(
        fits.Column(name=name, format=format, array=array ))
                          
    orig_cols = fits_array.columns
    new_cols = fits.ColDefs(cols)
    new_fits = fits.BinTableHDU.from_columns(orig_cols + new_cols)
    return new_fits.data
