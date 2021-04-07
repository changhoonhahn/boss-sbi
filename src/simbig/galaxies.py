'''

module for galaxy catalogs. Module includes functions to construct galaxy
catalogs using halo occupation distribution (HOD) models. HOD models paint
galaxies onto halo catalogs constructed from N-body simulations 



References
----------
* make_survey pipeline from BOSS: https://trac.sdss3.org/browser/repo/mockFactory/trunk/make_survey?order=name


'''
import os 
import numpy as np 
# --- nbodykit --- 
import nbodykit.lab as NBlab
from nbodykit.hod import Zheng07Model, Leauthaud11Model, Hearin15Model


def thetahod_lowz_ngc(): 
    ''' bestfit parameters of the lowz catalog from Table 2 of Manera et al.(2015)
    '''
    p_hod = {
            'logMmin': 13.20, 
            'sigma_logM': 0.62, 
            'logM0': 13.24, 
            'logM1': 14.32, 
            'alpha': 0.9
            }
    return p_hod 


def thetahod_lowz_sgc(): 
    ''' bestfit parameters of the lowz catalog from Table2 of Marnera et al.(2015) 

    Notes
    -----
    * Manera+(2015) actually uses a redshift dependent HOD. The HOD that's
        currently implemented is primarily for the 0.2 < z < 0.35 population,
        which has nbar~3x10^-4 h^3/Mpc^3
    '''
    p_hod = {
            'logMmin': 13.14, 
            'sigma_logM':0.55,
            'logM0': 13.43, 
            'logM1': 14.58, 
            'alpha': 0.93 
            }
    return p_hod 


def hodGalaxies(halos, p_hod, seed=None): 
    ''' populate given halo catalog (halos) with galaxies based on HOD model
    with p_hod parameters. Currently only supports the Zheng+(2007) model.

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


def BOSSGalaxies(sample='lowz-south'): 
    ''' Read in BOSS galaxy catalog. Data can be downloaded from 
    https://data.sdss.org/sas/dr12/boss/lss/
    

    Parameters
    ----------
    sample : string
        Specify the exact BOSS sample. For options see
        https://data.sdss.org/sas/dr12/boss/lss/
        (Default: 'cmass-north')

    Returns
    -------
    data : nbodykit.lab.FITSCatalog object
        BOSS galaxy catalog  
    '''
    if sample == 'cmass-north': 
        fgal = os.path.join(os.environ('BOSSSBI_DIR'), 'boss',
                'galaxy_DR12v5_CMASS_North.fits.gz') 
    elif sample == 'lowz-south': 
        fgal = os.path.join(
                os.path.dirname(os.path.realpath(__file__)), 'dat',
                'galaxy_DR12v5_LOWZ_South.fits.gz')
    else: 
        raise NotImplementedError

    data = NBlab.FITSCatalog(fgal)
    return data
