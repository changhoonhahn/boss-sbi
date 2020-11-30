'''


module for reading in Quijote data 


'''
import os
import numpy as np 
import nbodykit.lab as NBlab
# --- emanu --- 
from . import readfof 
from . import readsnap as RS


quijote_zsnap_dict = {0.: 4, 0.5: 3, 1.:2, 2.: 1, 3.: 0}


def Nbody(): 
    '''
    '''
    return None 


def Halos(halo_folder, z=0.5, Om=None, Ol=None, h=None, Hz=None, Ob=0.049, ns=0.9624, s8=None, silent=True): 
    ''' read in Quijote halo catalog given the folder and snapshot # and store it as
    a nbodykit HaloCatalog object. The HaloCatalog object is convenient for 
    populating with galaxies and etc.


    Parameters
    ----------
    halo_folder : string
        directory that contains the halo catalogs e.g. on tiger it'd be
        something like: /projects/QUIJOTE/Halos/latin_hypercube/HR_0/

    Return 
    ------
    cat : nbodykit.lab.HaloCatalog 
        Quijote halo catalog  
    '''
    # if snapshot folder is not specified 
    # then all values have to be specified in kwargs
    assert all([tt is not None for tt in [Om, Ol, z, h, Hz]]) 

    # redshift snapshot 
    assert z in quijote_zsnap_dict.keys(), 'snapshots are available at z=0, 0.5, 1, 2, 3'
    snapnum = quijote_zsnap_dict[z]

    # this is a discrepant cosmology. it is not used for anything. but it's included for nbodykit 
    cosmo = NBlab.cosmology.Planck15.clone() 

    # read FOF catalog (~90.6 ms) 
    Fof = readfof.FoF_catalog(halo_folder, snapnum, read_IDs=False,
            long_ids=False, swap=False, SFR=False)
    group_data = {}  
    group_data['Length']    = Fof.GroupLen
    group_data['Position']  = Fof.GroupPos/1e3
    group_data['Velocity']  = Fof.GroupVel
    group_data['Mass']      = Fof.GroupMass*1e10
    # calculate velocity offset
    rsd_factor = (1.+z) / Hz
    group_data['VelocityOffset'] = group_data['Velocity'] * rsd_factor
    # save to ArryCatalog for consistency
    cat = NBlab.ArrayCatalog(group_data, BoxSize=np.array([1000., 1000., 1000.])) 
    cat = NBlab.HaloCatalog(cat, cosmo=cosmo, redshift=z, mdef='vir') 
    cat.attrs['h'] = h 
    cat.attrs['Hz'] = Hz 
    cat.attrs['rsd_factor'] = rsd_factor 
    return cat

