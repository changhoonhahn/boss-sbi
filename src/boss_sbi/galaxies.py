'''

module for galaxy catalogs. Module includes functions to construct galaxy
catalogs using halo occupation distribution (HOD) models. HOD models paint
galaxies onto halo catalogs constructed from N-body simulations 



References
----------
* make_survey pipeline from BOSS: https://trac.sdss3.org/browser/repo/mockFactory/trunk/make_survey?order=name


'''
import numpy as np 
# --- nbodykit --- 
import nbodykit.lab as NBlab
from nbodykit.hod import Zheng07Model, Leauthaud11Model, Hearin15Model



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



def BOSSGalaxies(sample='cmass-north'): 
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
    else: 
        raise NotImplementedError

    data = NBlab.FITSCatalog(fgal)
    return data
