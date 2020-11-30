'''


module for reading in Quijote data 


'''
import os
import numpy as np 
import nbodykit.lab as NBlab
# --- emanu --- 
from . import readfof 
from . import readsnap as RS


def Nbody(): 
    '''
    '''
    return None 


def Halos(halo_folder, snap_folder, snapnum, Om=None, Ol=None, z=None, h=None, Hz=None, Ob=0.049, ns=0.9624, s8=None, silent=True): 
    ''' read in Quijote halo catalog given the folder and snapshot # and store it as
    a nbodykit HaloCatalog object. The HaloCatalog object is convenient for 
    populating with galaxies and etc.


    Parameters
    ----------
    halo_folder : string
        directory that contains the halo catalogs e.g. in my local directory it'd be 
        something like:
        /Users/ChangHoon/data/emanu/halos/hades/0.0eV/1

    snap_folder : string
        direcotry that contains the snapshot. 

    :param snapnum: 
        snapshot number 

        snapnum = 0 --> z=3
        snapnum = 1 --> z=2
        snapnum = 2 --> z=1
        snapnum = 3 --> z=0.5
        snapnum = 4 --> z=0


    Return 
    ------
    cat : nbodykit.lab.HaloCatalog 
        Quijote halo catalog  
    '''
    if snap_folder is not None: 
        # read in Gadget header (~65.1 microsec) 
        header = RG.header(os.path.join(snap_folder, 'snapdir_%s' % str(snapnum).zfill(3), 
            'snap_%s' % str(snapnum).zfill(3)))
        Om  = header.omega_m
        Ol  = header.omega_l
        z   = header.redshift
        h   = header.hubble 
        Hz  = header.Hubble #100.0 * np.sqrt(Om * (1.0 + z)**3 + Ol) # km/s/(Mpc/h)
    else: 
        # if snapshot folder is not specified 
        # then all values have to be specified in kwargs
        assert all([tt is not None for tt in [Om, Ol, z, h, Hz]]) 

    # this is a discrepant cosmology. it is not used for anything. but it's included for nbodykit 
    cosmo = NBlab.cosmology.Planck15.clone() 

    # read FOF catalog (~90.6 ms) 
    Fof = readfof.FoF_catalog(halo_folder, snapnum, read_IDs=False, long_ids=False, swap=False, SFR=False)
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

