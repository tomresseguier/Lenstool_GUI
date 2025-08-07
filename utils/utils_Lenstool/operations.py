from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import numpy as np




def make_image_to_source(dx_map, dy_map, wcs) :
    arcsec_per_pix = np.sqrt(wcs.pixel_scale_matrix[0, 0]**2+wcs.pixel_scale_matrix[0, 1]**2) * 3600
    
    def transform_coords(ra, dec) :
        coord = SkyCoord(ra, dec, unit='deg')
        local_coord = WCS.world_to_pixel(wcs, coord)
        x, y = local_coord[0]*1., local_coord[1]*1.
        
        dx = dx_map[ int(y), int(x) ]
        dy = dy_map[ int(y), int(x) ]
        x_source = x - dx/arcsec_per_pix
        y_source = y - dy/arcsec_per_pix
        
        source_coords = WCS.pixel_to_world(wcs, x_source, y_source)
        ra_source, dec_source = source_coords.ra.deg, source_coords.dec.deg
        return ra_source, dec_source
    
    return transform_coords


def make_magnification_function(mmap, wcs) :    
    def get_magnification(ra, dec) :
        coord = SkyCoord(ra, dec, unit='deg')
        local_coord = WCS.world_to_pixel(wcs, coord)
        x, y = local_coord[0]*1., local_coord[1]*1.
        return mmap[ int(y), int(x) ]
    return get_magnification


