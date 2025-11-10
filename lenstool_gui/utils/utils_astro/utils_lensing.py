import numpy as np
import os
import sys
from astropy import units as u
from astropy import constants as astro_constants
#sys.path.append( os.path.join(os.path.expanduser("~"), 'Library/Mobile Documents/com~apple~CloudDocs/RESEARCH/PROCESS/') )
from lenstool_gui.utils.utils_astro.get_cosmology import get_cosmo
cosmo = get_cosmo()


def SigmaCrit(zl, zs) :
    #print("Using cosmology: ", cosmo)
    Dl  = cosmo.angular_diameter_distance_z1z2(0 , zl)
    Ds  = cosmo.angular_diameter_distance_z1z2(0 , zs)
    Dls = cosmo.angular_diameter_distance_z1z2(zl, zs)
    return ( astro_constants.c**2 / (4*np.pi*astro_constants.G) * Ds/Dl/Dls ).to(u.kg/u.m**2)
