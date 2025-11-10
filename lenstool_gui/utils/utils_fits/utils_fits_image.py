import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def open_image(image_path) :
    with fits.open(image_path) as hdus :
        if isinstance(hdus, fits.hdu.hdulist.HDUList) :
            image_hdus = []
            for hdu in hdus :
                if hdu.is_image and isinstance(hdu.data, np.ndarray) :
                    image_hdus.append(hdu)
            print(f"{len(image_hdus)} image hdus found.")
            
            if len(image_hdus)==0 :
                print('No image found in FITS file.')
            
            true_image_hdus = []
            for i, image_hdu in enumerate(image_hdus):
                if image_hdu.data is not None and len(image_hdu.data.shape) >= 2:
                    print(f"HDU {i} contains image data with shape {hdu.data.shape}")
                    true_image_hdus.append(image_hdu)
            print(f"{len(true_image_hdus)} non empty images found.")
            
            
            keyword = None
            if True in np.unique(['EXTNAME' in true_image_hdu.header for true_image_hdu in true_image_hdus]) :
                keyword = 'EXTNAME'
            elif True in np.unique(['FILETYPE' in true_image_hdu.header for true_image_hdu in true_image_hdus]) :
                keyword = 'FILETYPE'
            print("keyword is", keyword)
            if keyword is None :
                print("No 'SCI' extname found.")
                selected_hdus = true_image_hdus
            else :
                #sci_hdus = []
                selected_hdus = []
                wht_hdus = []
                for true_image_hdu in true_image_hdus :
                    if keyword in true_image_hdu.header :
                        if true_image_hdu.header[keyword]=='SCI' :
                            #sci_hdus.append(true_image_hdu)
                            selected_hdus.append(true_image_hdu)
                            print('Science image found.')
                        if true_image_hdu.header[keyword]=='WHT' :
                            wht_hdus.append(true_image_hdu)
                            print('Weight image found.')
                if len(selected_hdus)==0 and len(wht_hdus)==0 :
                    selected_hdus = true_image_hdus
            
            
            if len(selected_hdus)==3 :
                x_sizes = [selected_hdu.data.shape[0] for selected_hdu in selected_hdus]
                y_sizes = [selected_hdu.data.shape[1] for selected_hdu in selected_hdus]
                if len(np.unique(x_sizes))==1 and len(np.unique(y_sizes))==1 :
                    print("Assuming RGB data.")
                data_red = selected_hdus[0].data
                data_green = selected_hdus[1].data
                data_blue = selected_hdus[2].data
                image = np.dstack((data_red, data_green, data_blue))
            else :
                print("Using first hdu.")
                image = selected_hdus[0].data
            
            
            wcs = WCS(selected_hdus[0].header)
            header = selected_hdus[0].header
            
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
            
            return image, pix_deg_scale, orientation, wcs, header
        
        else :
            print('Unable to extract image data from FITS file')
            return None




