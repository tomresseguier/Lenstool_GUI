import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.cosmology import WMAP9, FlatLambdaCDM
from astropy import constants as astro_constants
#import pysynphot
import matplotlib.pyplot as plt
from tqdm import tqdm

import sys
import os
#module_dir = os.path.dirname(os.path.abspath(__file__))
#main_dir = os.path.dirname(module_dir)
#sys.path.append(main_dir)
#sys.path.append(module_dir)

#os.chdir(main_dir)


DATA_dir = "../DATA/"
spt0615_image_path = DATA_dir + 'v7/color/spt0615_color_rgb.fits'



def mosaic_pixel_to_relative(x, y, reference=(93.9656102, -57.7801998)) :
    ref = SkyCoord(reference[0], reference[1], unit='deg')
    with fits.open(spt0615_image_path) as image :
        image_wcs = WCS(image[0])
    world_radec = WCS.pixel_to_world(image_wcs, x, y)
    #relative_coord = SkyCoord(world_radec.ra - ref.ra, world_radec.dec - ref.dec, unit='deg')
    relative_coord = ( (world_radec.ra - ref.ra)*np.cos(ref.dec.rad) , world_radec.dec - ref.dec )
    return -relative_coord[0].arcsec, relative_coord[1].arcsec


def world_to_relative(ra, dec, reference=(93.9656102, -57.7801998)) :
    ref = SkyCoord(reference[0], reference[1], unit='deg')
    world_radec = SkyCoord(ra, dec, unit='deg')
    relative_coord = ( (world_radec.ra - ref.ra)*np.cos(ref.dec.rad), world_radec.dec - ref.dec )
    return -relative_coord[0].arcsec, relative_coord[1].arcsec


def relative_to_mosaic_pixel(x, y, reference=(93.9656102, -57.7801998)) :
    arcsec_to_deg = 1/3600
    ref = SkyCoord(reference[0], reference[1], unit='deg')
    with fits.open(spt0615_image_path) as image :
        image_wcs = WCS(image[0])
    world_radec = (ref.ra.deg - x/np.cos(ref.dec.rad) * arcsec_to_deg, ref.dec.deg + y * arcsec_to_deg)
    world_coord = SkyCoord(world_radec[0], world_radec[1], unit='deg')
    pixel_xy = WCS.world_to_pixel( image_wcs, world_coord)
    if len(pixel_xy[0].shape)==0 :
        pixel_xy = (pixel_xy[0]*1., pixel_xy[1]*1.)
    return pixel_xy


def relative_to_world(x, y, reference=(93.9656102, -57.7801998)) :
    arcsec_to_deg = 1/3600
    ref = SkyCoord(reference[0], reference[1], unit='deg')
    world_radec = (ref.ra.deg - x/np.cos(ref.dec.rad) * arcsec_to_deg, ref.dec.deg + y * arcsec_to_deg)
    return world_radec

    
def sky_to_mosaic_pixel(ra, dec) :
    with fits.open(spt0615_image_path) as image :
        image_wcs = WCS(image[0])
    
    coord = SkyCoord(ra, dec, unit="deg")
    image_coord = WCS.world_to_pixel(image_wcs, coord)
    if len(image_coord[0].shape)==0 :
        image_coord = (image_coord[0]*1., image_coord[1]*1.)
    return image_coord


def mosaic_pixel_to_world(x, y, return_SkyCoord=False) :
    with fits.open(spt0615_image_path) as image :
        image_wcs = WCS(image[0])
    world_radec = WCS.pixel_to_world(image_wcs, x, y)
    if return_SkyCoord :
        to_return = world_radec
    else :
        to_return = (world_radec.ra.deg, world_radec.dec.deg)
    return to_return


def cat_len_displayer(catalogs_paths) :
    c=0
    for catalogs_path in catalogs_paths :
        cat = fits.open(catalogs_path)[1].data
        print(len(cat))
        c+=len(cat)
    print('total = ' + str(c))


def get_ref_coord(mass_model_dir) :
    para_file = open(mass_model_dir + "para.out", 'r')
    lines = para_file.readlines()
    runmode_index = np.where(['runmode' in line for line in lines])[0][0]
    i = runmode_index
    check = True
    while check :
        i+=1
        if 'reference' in lines[i] :
            check = False
            ref_ra = float(lines[i].split()[2])
            ref_dec = float(lines[i].split()[3])
        elif i - runmode_index > 100 :
            check = False
            print('runmode keyword not found')
    return ref_ra, ref_dec


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
    distance = WMAP9.comoving_distance(redshift)
    separation_kpc = separation_rad * distance.to(u.kpc) 
    return separation_kpc


def kpc_to_arcsec(separation_kpc, redshift) :
    distance = WMAP9.comoving_distance(redshift)
    separation_rad = separation_kpc * u.kpc / distance.to(u.kpc)
    arcsec = separation_rad * u.rad.to(u.arcsec)
    return arcsec


#def plot_band(band_name) :
#    band = pysynphot.ObsBandpass('acs,wfc1,' + band_name + ',mjd#59313')
#    fig, ax = plt.subplots()
#    ax.plot(band.binset, band(band.binset), 'b')
#    return ax


def counts_to_magAB(counts) :
    #band = pysynphot.ObsBandpass('acs,wfc1,f814w,mjd#59313')
    #spec = pysynphot.BlackBody(10000)
    #spec_norm = spec.renorm(1, 'counts', band)
    #obs = pysynphot.Observation(spec_norm, band)
    #AB_zero_point = obs.effstim('abmag')
    magAB = -2.5*np.log10(counts) + 25.934 #AB_zero_point
    return magAB


def magAB_to_flux_cgs(magAB) :
    flux_cgs = 10**( (-magAB + 48.60)/2.5 )
    return flux_cgs * (u.erg / u.cm**2 / u.second) # * u.Unit('erg cm^-2 s^-1')
    

def flux_cgs_to_magAB_HST(flux_cgs) :
    magAB = 48.60 - 2.5 * np.log10(flux_cgs)
    return magAB


def flux_muJy_to_magAB(flux_muJy) :
    magAB = -2.5 * np.log10(flux_muJy*1E-6/3631)
    return magAB


# This function does not give the correct luminosity
def flux_to_luminosity(flux, z=0.972) :
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    distance = cosmo.luminosity_distance(z)
    distance = cosmo.comoving_distance(z)
    
    luminosity = flux.cgs.value * 4 * np.pi * distance.cgs.value**2 * u.Unit('erg s^-1')
    return luminosity


def apparent_to_absolute(apparent_mag, z=0.972) :
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    distance = cosmo.luminosity_distance(z)
    absolute_magnitude = apparent_mag - 5*np.log10( (distance.to(u.parsec) / (10*u.parsec)).value )
    return absolute_magnitude


def abs_mag_to_luminosity(abs_mag) :
    abs_mag_sun = 4.83
    L = 10**(-0.4*(abs_mag-abs_mag_sun))                             #in solar luminosities
    #L = 10**(-0.4*(abs_mag-abs_mag_sun)) * astro_constants.L_sun    # in J/s
    return L


def angle_to_sky_distance(angle_arcsec, z=0.972):
    angle_rad = angle_arcsec * (1/3600) * (np.pi/180)
    cosmology = FlatLambdaCDM(H0=70, Om0=0.3)
    comoving_dist = cosmology.comoving_distance(z)
    distance_kpc = comoving_dist * angle_rad
    return distance_kpc.to('kpc')


def angular_diameter_distance(z) :
    cosmology = FlatLambdaCDM(H0=70, Om0=0.3)
    return cosmology.angular_diameter_distance(z)


def v_disp(galaxy, z=0.972) :
    abs_mag = apparent_to_absolute(galaxy['MAG_AUTO'], z)
    luminosity = abs_mag_to_luminosity(abs_mag)
    FJ_coef = 7.
    mass = FJ_coef * luminosity
    pix_arcsec_scale = 0.05
    radius_arcsec = galaxy['radius'] * pix_arcsec_scale
    radius_kpc = angle_to_sky_distance(radius_arcsec, z)
    v_disp = np.sqrt( (astro_constants.G * mass * astro_constants.M_sun) / radius_kpc.to('meter') ).to('km s^-1')
    return v_disp
    

def remove_png_margins(im) :
    x_shape = im.shape[0]
    y_shape = im.shape[1]
    
    threshold = 3000.
    
    i=0
    check = False
    while not check :
        i+=1
        check = not np.sum(im[0:i,:,:])>i*y_shape*4.-threshold
        if check :
            print(np.sum(im[0:i,:,:]), i*y_shape*4.)
    left_x_margin = i-1
    
    i=0
    check = False
    while not check :
        i+=1
        check = not np.sum(im[-i-1:-1,:,:])>i*y_shape*4.-threshold
    right_x_margin = -i-1
    
    i=0
    check = False
    while not check :
        i+=1
        check = not np.sum(im[:,0:i,:])>i*x_shape*4.-threshold
    bottom_y_margin = i-1
    
    i=0
    check = False
    while not check :
        i+=1
        check = not np.sum(im[:,-i-1:-1,:])>i*x_shape*4.-threshold
    top_y_margin = -i-1
    
    im_cropped = im[left_x_margin:right_x_margin, bottom_y_margin:top_y_margin, :]
    
    return im_cropped


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


def density_frac(mass, radius) :
    """
    Parameters
    ----------
    mass : in M_sun
    radius : in kpc
    """
    #cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    #rho_c = 8 * np.pi * astro_constants.G.value / (3 * ( H )**2)
    h=0.7
    rho_c = 2.7754E11 * h**2  * 1E-9
    Volume = 4./3 * np.pi * radius**3
    rho = mass / Volume
    return rho / rho_c








