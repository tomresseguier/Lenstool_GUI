import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
from astropy.io import fits
from astropy.wcs import WCS


def plot_magnification(magnification, fig=None, ax=None, threshold=20, cmap='seismic', make_colorbar=True, absolute=False) :
    if fig is None :
        fig, ax = plt.subplots()
    if type(magnification)==str :
        with fits.open(magnification) as hdulist :
            magnification = hdulist[0].data
            magnification_wcs = WCS(hdulist[0].header)
    #else :
    #    magnification_wcs = magnification[1]
    #    print(magnification_wcs)
    #    magnification = magnification[0]
    
    magnification[ np.where(magnification>threshold) ] = threshold
    magnification[ np.where(magnification<-threshold) ] = -threshold
    if absolute :
        magnification = np.abs(magnification)
        cmap_alpha = cmap
    else :
        cmap_alpha = cmap + '_alpha'
        if cmap_alpha in plt.colormaps() :
            colormaps.unregister(cmap_alpha)
        ncolors = 256
        color_array = plt.get_cmap(cmap)(range(ncolors))
        split = int(ncolors/2)
        color_array[:split, 3] = np.linspace(1., 0., split)
        color_array[split:, 3] = np.linspace(0., 1., split)
        
        cmap_object = LinearSegmentedColormap.from_list(name=cmap_alpha, colors=color_array)
        colormaps.register(cmap_object)
    
    im = ax.imshow(magnification, cmap=cmap_alpha, origin='lower')#, transform=ax.get_transform(magnification_wcs)) #also cmap="PuOr" or cmap="BrBG"
    if make_colorbar :
        cb = plt.colorbar(im)
        cb.set_label('magnification')
    
    #with fits.open(mosaic_path) as hdu :
    #    header = hdu[0].header if 'NAXIS1' in hdu[0].header else hdu[1].header
    #ax.set_xlim([0,header['NAXIS1']])
    #ax.set_ylim([0,header['NAXIS2']])
    return fig, ax