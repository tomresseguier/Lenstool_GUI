from astropy.cosmology import WMAP9, Planck18, FlatLambdaCDM


def set_cosmo() :
    #cosmo = Planck18
    #cosmo = WMAP9
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    return cosmo
