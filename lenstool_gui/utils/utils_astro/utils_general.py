import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy import constants as astro_constants
from astropy.coordinates import SkyCoord
import stsynphot
import matplotlib.pyplot as plt
from .get_cosmology import get_cosmo
cosmo = get_cosmo()


#stsynphot.initialize()

def rearrange_points(x, y) :
    """
    This function takes in the coordinates of a set of points spread out in the
    field as two arrays and returns these coordinates sorted so that a plot of 
    the points will give a continuous line (in the output arrays, every 
    coordinate appears between its two closest neighbours)
    """
    x_bis = x.copy()
    y_bis = y.copy()
    
    new_x = np.zeros(len(x_bis)+1)
    new_y = np.zeros(len(y_bis)+1)
    new_x[0] = x_bis[0]
    new_y[0] = y_bis[0]
    new_x[-1] = x_bis[0]
    new_y[-1] = y_bis[0]
    
    idx = 0
    for i in range(len(x)-1) :     
        x_bis = np.delete(x_bis, idx)
        y_bis = np.delete(y_bis, idx)
        dist = np.sqrt( (x_bis - new_x[i])**2 + (y_bis - new_y[i])**2 )
        idx = np.argmin(dist)
        new_x[i+1] = x_bis[idx]
        new_y[i+1] = y_bis[idx]
    return new_x, new_y


def arcsec_to_kpc(arcsec, redshift) :
    separation_rad = arcsec * u.arcsec.to(u.rad)
    ang_diam_dist = cosmo.angular_diameter_distance(redshift)
    separation_kpc = separation_rad * ang_diam_dist.to(u.kpc) 
    return separation_kpc


def kpc_to_arcsec(separation_kpc, redshift) :
    ang_diam_dist = cosmo.angular_diameter_distance(redshift)
    separation_rad = separation_kpc * u.kpc / ang_diam_dist.to(u.kpc)
    arcsec = separation_rad * u.rad.to(u.arcsec)
    return arcsec


def plot_bandpass(band_name='acs,wfc,f814w'):
    bp = stsynphot.spectrum.ObsBandpass(band_name)
    plt.figure(figsize=(7, 4))
    plt.plot(bp.waveset.to('nm'), bp.throughput, color='blue', lw=2)
    plt.title(f'HST Bandpass: {band_name}', fontsize=14, fontweight='bold')
    plt.xlabel('Wavelength [nm]', fontsize=12)
    plt.ylabel('Throughput', fontsize=12)
    plt.ylim(0, 1.05)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def counts_to_muJy(electrons_per_second, photflam=None, photplam=None, photmjsr=None, pixel_area_arcsec2=None, band_name=None, return_magAB=False):
    """
    Converts counts (electrons/sec) to flux in µJy/arcsec^2 or AB magnitude.
    Use either PHOTFLAM and PHOTPLAM (or band_name), or PHOTMJSR and pixel_area_arcsec2.
    Parameters
    ----------
    electrons_per_second : float or numpy.ndarray
        Count rate in electrons per second.
    photflam : float, optional
        Inverse sensitivity (erg/s/cm^2/A per e-/s).
    photplam : float, optional
        Pivot wavelength in Angstrom. Can be derived from band_name.
    photmjsr : float, optional
        Conversion factor from e-/s to MJy/sr.
    pixel_area_arcsec2 : float, optional
        Pixel area in arcsec^2.
    band_name : str, optional
        Filter name string understood by ObsBandpass, e.g.
        'acs,wfc,f814w', 'wfc3,uvis1,f606w', 'wfc3,ir,f110w'.
        Used to calculate photplam if photplam is not provided.
    return_magAB : bool, optional
        If True, return AB magnitude instead of flux in µJy/arcsec^2. Default is False.
    Returns
    -------
    flux or magAB : float or numpy.ndarray
        Flux in µJy/arcsec^2 or AB magnitude.
    """
    if band_name is not None and photplam is None:
        photplam = get_photplam(band_name)

    if photflam is not None and photplam is not None:
        # Using PHOTFLAM
        f_lambda = electrons_per_second * photflam # in erg/s/cm^2/A
        f_nu = f_lambda * (photplam**2 / (astro_constants.c.to(u.AA/u.s).value)) # in erg/s/cm^2/Hz
        flux_jy_per_pixel = f_nu * 1e23
        
        if pixel_area_arcsec2 is None:
            raise ValueError("pixel_area_arcsec2 must be provided when using PHOTFLAM.")
        flux_uJy_per_arcsec2 = (flux_jy_per_pixel * 1e6) / pixel_area_arcsec2

        if return_magAB:
            magAB = -2.5 * np.log10(flux_jy_per_pixel / 3631)
            return magAB
        else:
            return flux_uJy_per_arcsec2

    elif photmjsr is not None and pixel_area_arcsec2 is not None:
        # Using PHOTMJSR
        surf_bright_mjy_per_sr = electrons_per_second * photmjsr
        flux_uJy_per_arcsec2 = surf_bright_mjy_per_sr * (u.MJy/u.sr).to(u.uJy/(u.arcsec**2))

        if return_magAB:
            pixel_area_sr = pixel_area_arcsec2 * (np.pi / (180.0 * 3600.0))**2
            flux_mjy = surf_bright_mjy_per_sr * pixel_area_sr
            flux_jy = flux_mjy * 1e6
            magAB = -2.5 * np.log10(flux_jy / 3631)
            return magAB
        else:
            return flux_uJy_per_arcsec2
    else:
        raise ValueError("Please provide either (photflam and photplam/band_name) or (photmjsr and pixel_area_arcsec2).")

def get_photplam(band_name='acs,wfc,f814w'):
    bp = stsynphot.spectrum.ObsBandpass(band_name)
    return bp.pivot().value


def flux_cgs_to_magAB(flux_cgs) :
    # Ensure flux is expressed as f_nu in cgs (erg cm^-2 s^-1 Hz^-1)
    if isinstance(flux_cgs, u.Quantity):
        f_nu_value = flux_cgs.to(u.erg / u.cm**2 / u.s / u.Hz).value
    else:
        f_nu_value = np.asarray(flux_cgs)
    magAB = -2.5 * np.log10(f_nu_value) - 48.60
    return magAB

def magAB_to_flux_cgs(magAB) :
    f_nu_value = 10**( -(np.asarray(magAB) + 48.60)/2.5 )
    return f_nu_value * (u.erg / u.cm**2 / u.s / u.Hz)



def muJy_to_magAB(flux_muJy) :
    magAB = -2.5 * np.log10(flux_muJy*1E-6/3631)
    return magAB

def magAB_to_muJy(magAB) :
    flux_jy = 10**( -(np.asarray(magAB) + 2.5 * np.log10(3631)) / 2.5 )
    flux_muJy = flux_jy * 1E6
    return flux_muJy

def apparent_to_absolute(app_mag, redshift) :
    distance = cosmo.luminosity_distance(redshift)
    abs_mag = app_mag - 5*np.log10( (distance.to(u.parsec) / (10*u.parsec)).value )
    return abs_mag


def pot_mass(a, s, sigma) :
    a = (a * u.kpc).to('m').value
    s = (s * u.kpc).to('m').value
    M_tot = 3/2 / astro_constants.G.value * (sigma*1E3)**2 * (s**2 - a**2) / s
    return M_tot / astro_constants.M_sun.value


def orientation_angle_diff(image_ref_path, image2_path) :
    with fits.open(image_ref_path) as hdu :
        wcs_ref = WCS(hdu[0].header)
    with fits.open(image2_path) as hdu :
        wcs2 = WCS(hdu[0].header)
    
    radec_origin = WCS.pixel_to_world(wcs2, 0, 0)
    radec_xoffset = WCS.pixel_to_world(wcs2, 100, 0)
    
    x_ref_origin, y_ref_origin = WCS.world_to_pixel(wcs_ref, radec_origin)
    x_xoffset, y_xoffset = WCS.world_to_pixel(wcs_ref, radec_xoffset)
    
    delta_x = x_xoffset - x_ref_origin
    delta_y = y_xoffset - y_ref_origin
    
    angle = np.arctan2(delta_y, delta_x)
    return angle


def world_to_relative(ra, dec, reference) :
    ref = SkyCoord(reference[0], reference[1], unit='deg')
    world_radec = SkyCoord(ra, dec, unit='deg')
    relative_coord = ( (ref.ra - world_radec.ra)*np.cos(world_radec.dec.rad), world_radec.dec - ref.dec )
    return relative_coord[0].arcsec, relative_coord[1].arcsec

def relative_to_world(xr, yr, reference) :
    ref = SkyCoord(reference[0], reference[1], unit='deg')
    dec = ref.dec.deg + yr*u.arcsec.to('deg')
    ra = ref.ra.deg - xr*u.arcsec.to('deg') / np.cos(dec*u.deg.to('rad'))
    return ra, dec


