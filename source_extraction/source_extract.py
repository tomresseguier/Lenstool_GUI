"""
source_extract.py from David Harvey's pyRRG
"""

#from pyextract import pysex
from astropy.io import fits
from astropy.table import Table
from .pysex import pysex

import numpy as np
import os as os
import sys
from .match_cat import run_match
from numpy.lib.recfunctions import append_fields as append_rec

module_dir = os.path.dirname(os.path.abspath(__file__))




def source_extract( image_path, weight_path=None, pixel_scale=0.03, zero_point=None,
                    out_dir='PWD', outfile_name='SExtractor_cat.fits',
                    return_sources=True) :
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
    
    #for name in ['cold_sources.fits', 'hot_sources.fits', 'matched_A_B.fits'] :
    #    if os.path.isfile( os.path.join(os.getcwd(), name) ) :
    #        yesno = input("file " + "'" + name + "'" + " found. Do you want to rerun source extraction (takes a while)? [y][n]")
            #if yesno!='y' :
            #    sys.exit('Execution aborted by user.')
    
    if out_dir == 'PWD':
        out_dir = os.getcwd()
    out_path = out_dir + '/' + outfile_name
    #config_dir = os.path.join(module_dir, 'SExtractor_config/from_DC')
    config_dir = os.path.join(module_dir, 'SExtractor_config/from_pyRRG')
    
    if " " in config_dir :
        print('File in Drive, replacing:')
        print(config_dir)
        print('with:')
        config_dir = config_dir.replace('Mobile Documents/com~apple~CloudDocs', 'Mobile_Documents')
        print(config_dir)
    
    #check_sex_files(config_dir)
    
    header = fits.open( image_path )[0].header
    findPhot = np.array(['PHOTFLAM' in key for key in header.keys()])
    if np.all(findPhot == False) :
        header = fits.open( image_path )[1].header
    
    ### FIX THE PHOTOMETRY HERE!!! ###
    zero_point = acs_zero_point(header) #0.
        
    conf_args = {'PIXEL_SCALE': pixel_scale,
                 'MAG_ZEROPOINT': zero_point,
                 'WEIGHT_TYPE': 'NONE', #'MAP_WEIGHT'
                 'PARAMETERS_NAME': config_dir + '/rrg.param',
                 'STARNNW_NAME': config_dir + '/default.nnw',
                 'FILTER_NAME': config_dir + '/gauss_5.0_9x9.conv'}
    if weight_path is not None :
        conf_args['WEIGHT_IMAGE'] : weight_path
    
    print('conf_args:')
    print(conf_args)
    
    #F## COLD RUN ###
    cold_conf = config_dir + '/HFF_cold.param'
    
    cold_sources_path = os.path.join(os.getcwd(), 'cold_sources.fits')
    if os.path.isfile( cold_sources_path ) :
        yesno = input("file " + "'" + 'cold_sources.fits' + "'" + " found. Do you want to rerun source extraction (takes a while)? [y][n]")
        if yesno!='y' :
            with fits.open(cold_sources_path) as hdu :
                cold_sources = hdu[1].data
    else :
        cold_sources = pysex.run( image_path, \
                                  conf_file=cold_conf, \
                                  conf_args=conf_args, \
                                  param_file=config_dir+'/rrg.param')
        cold_sources = append_fits_field( cold_sources, 'RA', cold_sources['X_WORLD'])
        cold_sources = append_fits_field( cold_sources, 'DEC', cold_sources['Y_WORLD'])
    
    

    #Second hot 
    hot_conf = config_dir+'/HFF_hot.param'
    
    hot_sources_path = os.path.join(os.getcwd(), 'hot_sources.fits')
    if os.path.isfile( hot_sources_path ) :
        yesno = input("file " + "'" + 'hot_sources.fits' + "'" + " found. Do you want to rerun source extraction (takes a while)? [y][n]")
        if yesno!='y' :
            with fits.open(hot_sources_path) as hdu :
                hot_sources = hdu[1].data
    else :
        hot_sources = pysex.run( image_path, \
                                   conf_file=hot_conf, \
                                   conf_args=conf_args, \
                                   param_file=config_dir+'/rrg.param')
        hot_sources = append_fits_field( hot_sources, 'RA', hot_sources['X_WORLD'])
        hot_sources = append_fits_field( hot_sources, 'DEC', hot_sources['Y_WORLD'])
    
    #The NUMBER is a weird thing
    
    hot_sources['NUMBER'] = np.arange( len(hot_sources['NUMBER'])) +1
    cold_sources['NUMBER'] = np.arange( len(cold_sources['NUMBER'])) +1
    
    fits.writeto( 'cold_sources.fits', cold_sources, overwrite=True )
    fits.writeto( 'hot_sources.fits', hot_sources, overwrite=True )


    print('Matching cold and hot sources:')
    print( str(len(cold_sources)) + ' cold sources' )
    print( str(len(hot_sources)) + ' hot sources' )
    print( 'TOTAL: ' + str(len(cold_sources) + len(hot_sources)) + ' detections' )
    
    matched_sources_path = os.path.join(os.getcwd(), 'matched_A_B.fits')
    yesno = input("file " + "'" + 'matched_A_B.fits' + "'" + " found. Do you want to rerun matching (can take a while)? [y][n]")
    if yesno!='y' :
        with fits.open(matched_sources_path) as hdu :
            matched_sources = hdu[1].data
    else :
        matched_sources = run_match( 'cold_sources.fits',
                                     'hot_sources.fits' )
    
    for iField in hot_sources.columns.names:
        hot_sources[iField][ matched_sources[1].data['NUMBER_2'] -1 ] = cold_sources[iField][ matched_sources[1].data['NUMBER_1'] - 1]
    
    print('Matching cold and hot sources:')
    print('Number of sources in matched catalog: ' + str(len(hot_sources)) + ' sources')
    
    #Need to retain numbering for bookeepin purposes
    hot_sources['NUMBER'] = np.arange(len(hot_sources))
    
    fits.writeto( out_dir + '/' + outfile_name, hot_sources, overwrite=True )
    
    for name in ['cold_sources.fits', 'hot_sources.fits', 'matched_A_B.fits'] :
        os.rename( os.path.join(os.getcwd(), name), os.path.join(out_dir, name) )
        
    if return_sources:
        return Table(hot_sources)
    
    
def acs_zero_point( header ):
    if 'PHOTFLAM' in header and 'PHOTZPT' in header :
        zpt = -2.5*np.log10(header['PHOTFLAM']) + header['PHOTZPT'] - 5.0*np.log10(header['PHOTPLAM'])+18.6921
    else :
        zpt = 0.
    return zpt

def check_sex_files( config_dir ):
    #Check all the sex files to see if they exist
    if (not os.path.isfile(config_dir+'/DC.conf')):
        raise ValueError('DC.conf not found at ' + config_dir+'/DC.conf')
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




def _reg_path(filename):
    if isinstance(filename, str) and len(filename)>0 and filename[0] != os.path.sep:
        return os.path.abspath(filename)
    else:
        return filename


def source_extract_DIM( detection_paths, measurement_paths, pixel_scale=0.03, zero_point=None,
                        out_dir='PWD', outfile_name='SExtractor_cat_DIM.fits',
                        return_sources=True) :
    
    ref_image_path = detection_paths if type(detection_paths) is str else detection_paths[0]
    image_path = measurement_paths if type(measurement_paths) is str else measurement_paths[0]
    print('Measuring in ' + image_path, "\nDetecting in " + ref_image_path)
    
    if out_dir == 'PWD':
        out_dir = os.getcwd()
    out_path = out_dir + '/' + outfile_name
    config_dir = os.path.join(module_dir, 'SExtractor_config/from_DC')
    check_sex_files(config_dir)
    
    header = fits.open( image_path )[0].header
    #findPhot = np.array(['PHOTFLAM' in key for key in header.keys()])
    #if np.all(findPhot == False) :
    #    header = fits.open( image_path )[1].header
    
    ### FIX THE PHOTOMETRY HERE!!! AND ADD WEIGHT FILE OPTION ###
    zero_point = 0. #acs_zero_point(header)
        
    conf_args = {'PIXEL_SCALE': pixel_scale,
                 'MAG_ZEROPOINT': zero_point,
                 'PARAMETERS_NAME': config_dir + '/rrg.param',
                 'STARNNW_NAME': config_dir + '/default.nnw',
                 'FILTER_NAME': config_dir + '/gauss_2.0_5x5.conv'}
    if type(detection_paths) is list :
        conf_args['WEIGHT_TYPE'] = 'MAP_WEIGHT,MAP_WEIGHT'
        conf_args['WEIGHT_IMAGE'] = _reg_path(detection_paths[1]) + ',' + _reg_path(measurement_paths[1])
        
    conf = config_dir + '/DC.conf'
    sources = pysex.run( image=image_path, \
                         imageref=ref_image_path, \
                         conf_file=conf, \
                         conf_args=conf_args, \
                         param_file=config_dir+'/rrg.param')
    
    sources = append_fits_field( sources, 'RA', sources['X_WORLD'])
    sources = append_fits_field( sources, 'DEC', sources['Y_WORLD'])
    
    sources['NUMBER'] = np.arange(len(sources))
    
    fits.writeto( out_dir + '/' + outfile_name, sources, overwrite=True )
    
    if return_sources:
        return Table(sources)














