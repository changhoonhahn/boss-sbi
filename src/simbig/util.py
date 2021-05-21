'''


module with some utility functions


'''
import numpy as np 
from astropy.stats import scott_bin_width
from scipy.interpolate import InterpolatedUnivariateSpline


def get_nofz(z, fsky, cosmo=None): 
    ''' calculate nbar(z) given redshift values and f_sky (sky coverage
    fraction)

    Parameters
    ----------
    z : array like
        array of redshift values 
    fsky : float 
        sky coverage fraction  
    cosmo : cosmology object 
        cosmology to calculate comoving volume of redshift bins 

    Returns
    -------
    number density at input redshifts: nbar(z) 

    Notes
    -----
    * based on nbdoykit implementation 

    '''
    # calculate nbar(z) for each galaxy 
    _, edges = scott_bin_width(z, return_bins=True)

    dig = np.searchsorted(edges, z, "right")
    N = np.bincount(dig, minlength=len(edges)+1)[1:-1]

    R_hi = cosmo.comoving_distance(edges[1:]) # Mpc/h
    R_lo = cosmo.comoving_distance(edges[:-1]) # Mpc/h

    dV = (4./3.) * np.pi * (R_hi**3 - R_lo**3) * fsky

    nofz = InterpolatedUnivariateSpline(0.5*(edges[1:] + edges[:-1]), N/dV, ext='const')
    
    return nofz
