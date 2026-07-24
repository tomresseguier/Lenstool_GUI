import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
#import time


def open_image(image_path) :
    # Manage the case where user passed a fits.hdu.hdulist.HDUList object
    if isinstance(image_path, fits.hdu.hdulist.HDUList) :
        hdus = image_path
    else :
        hdus = fits.open(image_path)

    # Make sure the hdus object is a fits.hdu.hdulist.HDUList object
    if isinstance(hdus, fits.hdu.hdulist.HDUList) :
        print(hdus.info())

        # Look at the image hdus only
        image_hdus = []
        for hdu in hdus :
            if hdu.is_image and isinstance(hdu.data, np.ndarray) :
                image_hdus.append(hdu)
        print(f"{len(image_hdus)} image hdus found.")
        
        if len(image_hdus)==0 :
            hdus.close()
            raise ValueError('No image found in FITS file.')
        
        # Additional check to make sure the image hdus are not empty
        true_image_hdus = []
        for i, image_hdu in enumerate(image_hdus):
            if image_hdu.data is not None and len(image_hdu.data.shape) >= 2:
                print(f"HDU {i} contains image data with shape {image_hdu.data.shape}")
                true_image_hdus.append(image_hdu)
        print(f"{len(true_image_hdus)} non empty images found.")
        
        
        def get_hdu_by_name(hdul, name):
            selected_hdus = []
            for hdu in hdul:
                header = hdu.header
                if header.get("EXTNAME") == name or header.get("FILETYPE") == name:
                    selected_hdus.append(hdu)
            return selected_hdus

        sci_hdus = get_hdu_by_name(true_image_hdus, 'SCI')
        wht_hdus = get_hdu_by_name(true_image_hdus, 'WHT')
        other_hdus = [hdu for hdu in true_image_hdus if (hdu not in sci_hdus and hdu not in wht_hdus)]
        #if len(sci_hdus) == 0 and len(wht_hdus) == 0 :
        #    print('No SCI or WHT HDU found.')

        # We use SCI HDU if found
        if len(sci_hdus) != 0 :
            if len(sci_hdus)==3 :
                x_sizes = [sci_hdu.data.shape[0] for sci_hdu in sci_hdus]
                y_sizes = [sci_hdu.data.shape[1] for sci_hdu in sci_hdus]
                if len(np.unique(x_sizes))==1 and len(np.unique(y_sizes))==1 :
                    print("Assuming RGB data.")
                data_red = sci_hdus[0].data
                data_green = sci_hdus[1].data
                data_blue = sci_hdus[2].data
                image = np.stack((data_red, data_green, data_blue), axis=2)
            elif len(sci_hdus)==1 :
                image = sci_hdus[0].data
            else :
                print("Multiple SCI HDUs found. Using first one.")
                image = sci_hdus[0].data
            ref_hdu = sci_hdus[0]
        # Secondly, we use other HDUs if three are found
        elif len(other_hdus) == 3 :
            x_sizes = [other_hdu.data.shape[0] for other_hdu in other_hdus]
            y_sizes = [other_hdu.data.shape[1] for other_hdu in other_hdus]
            if len(np.unique(x_sizes))==1 and len(np.unique(y_sizes))==1 :
                print("Assuming RGB data.")
            data_red = other_hdus[0].data
            data_green = other_hdus[1].data
            data_blue = other_hdus[2].data
            image = np.stack((data_red, data_green, data_blue), axis=2)
            ref_hdu = other_hdus[0]
        # Thirdly, we use the weight HDU if found
        elif len(wht_hdus) != 0 :
            ref_hdu = wht_hdus[0]
            image = ref_hdu.data
        elif len(other_hdus) != 0 :
            ref_hdu = other_hdus[0]
            image = ref_hdu.data
        else :
            print('No image found in FITS file.')
            hdus.close()
            return None

        wcs = WCS(ref_hdu.header)
        header = ref_hdu.header
        
        if 'ORIENTAT' in header :
            orientation = header['ORIENTAT']
        elif 'CD1_1' in header and 'CD2_2' in header :
            if 'CD1_2' in header and 'CD2_1' in header :
                cd = np.array([[header['CD1_1'], header['CD1_2']], [header['CD2_1'], header['CD2_2']]])
                #det = np.linalg.det(cd)
                #sign = np.sign(det)
                orientation = np.arctan2(cd[1,0], cd[1,1])
            else :
                orientation = 0.0
        elif 'PC1_1' in header and 'PC2_2' in header :
            if 'PC1_2' in header and 'PC2_1' in header :
                cd = np.array([[header['PC1_1'], header['PC1_2']], [header['PC2_1'], header['PC2_2']]])
                #det = np.linalg.det(cd)
                #sign = np.sign(det)
                orientation = np.arctan2(cd[1,0], cd[1,1])
            else :
                orientation = 0.0                    
        else :
            orientation = None
        orientation = np.rad2deg(orientation) if orientation is not None else None
        
        ### Finding the pixel scale ###
        #if 'CD1_1' in hdus[0].header.keys() :
        #    CD1_1 = hdus[0].header['CD1_1']
        #    CD1_2 = hdus[0].header['CD1_2']
        #    pix_deg_scale = np.sqrt(CD1_1**2+CD1_2**2)
        #elif 'CDELT1' in hdus[0].header.keys() :
        #    pix_deg_scale = abs(hdus[0].header['CDELT1'])
        #else :
        #    pix_deg_scale = input('Pixel scale not found in header. Please provide manually in degrees:')
        pix_deg_scale = np.sqrt(wcs.pixel_scale_matrix[0, 0]**2+wcs.pixel_scale_matrix[0, 1]**2)
        
        hdus.close()
        return image, pix_deg_scale, orientation, wcs, header
        
    else :
        print('Unable to extract image data from FITS file')
        hdus.close()
        return None




