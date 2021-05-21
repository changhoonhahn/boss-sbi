''' 

script to make accompanying random file for LOWZ-SGC forwardmodel


'''
import os 
import h5py 
import numpy as np 

import pymangle 
from simbig import forwardmodel as FM 

def random(): 
    ''' construct hdf5 for random catalog of RA and Dec from ransack output 
    '''
    # read ransack output 
    frand = os.path.join(os.environ['QUIJOTE_DIR'], 'chang', 'random_DR12v5_LOWZ_South.dat')
    rand_ra, rand_dec = np.loadtxt(frand, unpack=True, usecols=[0,1], skiprows=1) 

    # save to hdf5 file  
    f = h5py.File(frand.replace('.dat', '.hdf5'), 'w')
    f.create_dataset('ra', data=rand_ra) 
    f.create_dataset('dec', data=rand_dec) 
    f.close()
    return None 


def random_veto():
    ''' construct random catalog of RA, Dec with veto mask applied
    '''
    fvetos = [
            'badfield_mask_postprocess_pixs8.ply', 
            'badfield_mask_unphot_seeing_extinction_pixs8_dr12.ply',
            'allsky_bright_star_mask_pix.ply',
            'bright_object_mask_rykoff_pix.ply', 
            'centerpost_mask_dr12.ply', 
            'collision_priority_mask_dr12.ply']

    veto_dir = os.path.join(os.environ['QUIJOTE_DIR'], 'chang', 'simbig')
    for i, fveto in enumerate(fvetos): 
        if i == 0:
            frand = h5py.File(os.path.join(os.environ['QUIJOTE_DIR'], 'chang', 'random_DR12v5_LOWZ_South.hdf5'), 'r')
        else: 
            frand = h5py.File(os.path.join(os.environ['QUIJOTE_DIR'], 'chang', 'random_DR12v5_LOWZ_South.veto.hdf5'), 'r')
        rand_ra = frand['ra'][...]
        rand_dec = frand['dec'][...]
        frand.close() 

        veto = pymangle.Mangle(os.path.join(veto_dir, fveto))
        w_veto = veto.weight(rand_ra, rand_dec)
        in_veto = (w_veto > 0.)

        # save points outside of veto mask to hdf5 file  
        f = h5py.File(os.path.join(os.environ['QUIJOTE_DIR'], 'chang', 'random_DR12v5_LOWZ_South.veto.hdf5'), 'w')
        f.create_dataset('ra', data=rand_ra[~in_veto]) 
        f.create_dataset('dec', data=rand_dec[~in_veto]) 
        f.close()
    return None 


if __name__=="__main__": 
    #random()
    random_veto()
