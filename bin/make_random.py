''' 

script to make accompanying random file for LOWZ-SGC forwardmodel


'''
import os 
import h5py 
import numpy as np 

from boss_sbi import forwardmodel as FM 

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
    frand = h5py.File(os.path.join(os.environ['QUIJOTE_DIR'], 'chang', 'random_DR12v5_LOWZ_South.hdf5'), 'r')
    rand_ra = frand['ra'][...]
    rand_dec = frand['dec'][...]
    
    # identify points in veto mask 
    in_veto = FM.BOSS_veto(rand_ra, rand_dec)

    # save points outside of veto mask to hdf5 file  
    f = h5py.File(os.path.join(os.environ['QUIJOTE_DIR'], 'chang', 'random_DR12v5_LOWZ_South.veto.hdf5'), 'w')
    f.create_dataset('ra', data=rand_ra[~in_veto]) 
    f.create_dataset('dec', data=rand_dec[~in_veto]) 
    f.close()
    return None 


if __name__=="__main__": 
    #random()
    random_veto()
