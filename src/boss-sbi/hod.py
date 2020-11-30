'''

module for using halo occupation distribution (HOD) models to paint galaxies
onto halo catalogs constructed from N-body simulations 


'''
import numpy as np 
from nbodykit.hod import Zheng07Model, Leauthaud11Model, Hearin15Model



def hodGalaxies(halos, p_hod, seed=None): 
    ''' populate given halo catalog (halos) with galaxies
    based on HOD model with p_hod parameters 

    Parameters
    ----------
    p_hod : dict
        dictionary specifying the HOD parameters 
    '''
    # check HOD parameters
    if 'alpha' not in p_hod.keys(): 
        raise ValueError
    if 'logMmin' not in p_hod.keys(): 
        raise ValueError
    if 'logM1' not in p_hod.keys(): 
        raise ValueError
    if 'logM0' not in p_hod.keys(): 
        raise ValueError
    if 'sigma_logM' not in p_hod.keys(): 
        raise ValueError

    # populate using HOD
    hod = halos.populate(Zheng07Model, seed=seed, **p_hod) 
    
    # recalculate RSD velocity offset using stored rsd_factor attr from emanu.sims.data.hqHalos
    hod['VelocityOffset'] = hod['Velocity'] * hod.attrs['rsd_factor'] 
    return hod 

