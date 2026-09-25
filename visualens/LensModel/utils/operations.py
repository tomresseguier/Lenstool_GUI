from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np




def make_image_to_source(dx_map, dy_map, wcs) :    
    def transform_coords(ra, dec) :
        coord = SkyCoord(ra, dec, unit='deg')
        local_coord = WCS.world_to_pixel(wcs, coord)
        x, y = local_coord[0]*1., local_coord[1]*1.
        dx = dx_map[ int(y), int(x) ]
        dy = dy_map[ int(y), int(x) ]
        ra_source = ra + dx /3600 /np.cos( dec * np.pi/180 )
        dec_source = dec - dy /3600
        return ra_source, dec_source
    return transform_coords


def MakeFunctionFromMap(mmap, wcs) :    
    def get_magnification(ra, dec) :
        coord = SkyCoord(ra, dec, unit='deg')
        local_coord = WCS.world_to_pixel(wcs, coord)
        x, y = local_coord[0]*1., local_coord[1]*1.
        try :
            return mmap[ int(y), int(x) ]
        except :
            return np.nan
    return get_magnification


