"""
match_cat.py from David Harvey's pyRRG
"""

import os
import numpy as np
from astropy.io import fits


module_dir = os.path.dirname(os.path.abspath(__file__))


def magnitude( flux, zpt, exptime, apcor ):
    return zpt - 2.5*np.log10(flux/exptime) - apcor

def run_match( cat_A, cat_B, search_rad=2 ):
    
    command_str ='stilts.sh tmatch2 in1="'+cat_A+'" in2="'+\
        cat_B+'" matcher=sky values1="RA DEC" values2="RA DEC" params="'\
        +str(search_rad)+'" out=matched_A_B.fits'
    
    ###########################################################################
    command_str = module_dir + '/stilts/' + command_str
    ###########################################################################
    
    os.system(command_str)

    matched_cat = fits.open('matched_A_B.fits')

    return matched_cat



