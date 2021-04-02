'''

script to test the function `forwardmodel.BOSS_randoms`

'''
import os, time
import numpy as np 
from boss_sbi import util as UT
from boss_sbi.halos import Quijote_LHC_HR
from boss_sbi import galaxies as Galaxies
from boss_sbi import forwardmodel as FM 
# --- plotting --- 
import matplotlib as mpl
mpl.use('pdf')
import matplotlib.pyplot as plt
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['axes.xmargin'] = 1
mpl.rcParams['xtick.labelsize'] = 'x-large'
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['ytick.labelsize'] = 'x-large'
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['legend.frameon'] = False


# read in halo catalog 
t0 = time.time() 
halos = Quijote_LHC_HR(1, z=0.5)
print('halo readin takes %f sec' % ((time.time() - t0)))

# get LOWZ HOD parameters
theta_hod = Galaxies.thetahod_lowz_sgc()

# apply HOD 
t0 = time.time() 
hod = Galaxies.hodGalaxies(halos, theta_hod, seed=0) 
print('HOD takes %f sec' % ((time.time() - t0)))


for veto in [True, False]: 
    # apply forward model 
    t0 = time.time() 
    gals = FM.BOSS(hod, sample='lowz-south', seed=0, veto=veto, fiber_collision=False, silent=False)
    print('forward model takes %f sec' % ((time.time() - t0)))

    # get random 
    t0 = time.time() 
    rand = FM.BOSS_randoms(gals, sample='lowz-south', veto=veto) 
    print('randoms take %f sec' % ((time.time() - t0)))

    # compare footprint
    fig = plt.figure(figsize=(10,5))
    sub = fig.add_subplot(111)
    sub.scatter(np.array(rand['RA']), np.array(rand['DEC']), c='k', s=1, rasterized=True, label='Random') 
    sub.scatter(np.array(rand['RA'])-360., np.array(rand['DEC']), c='k', s=1, rasterized=True) 
    sub.scatter(np.array(gals['RA']), np.array(gals['DEC']), c='C0', s=1, rasterized=True, label='Galaxies') 
    sub.scatter(np.array(gals['RA'])-360, np.array(gals['DEC']), c='C0', s=1, rasterized=True) 
    sub.legend(loc='upper right', fontsize=15, handletextpad=0, markerscale=10)
    sub.set_xlabel('RA', fontsize=25) 
    sub.set_xlim(-50, 70) 
    sub.set_ylabel('Dec', fontsize=25) 
    if veto: 
        fig.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)), '_random_footprint.veto.png'), bbox_inches='tight') 
    else: 
        fig.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)), '_random_footprint.png'), bbox_inches='tight') 


    # compare redhsift distribtuion 
    fig = plt.figure(figsize=(6,6))
    sub = fig.add_subplot(111)
    _ = sub.hist(np.array(rand['Z']), density=True, range=(0.1, 0.4), bins=100, histtype='step')
    _ = sub.hist(np.array(gals['Z']), density=True, range=(0.1, 0.4), bins=100, histtype='step')
    sub.set_xlim(0.1, 0.4) 
    sub.set_xlabel('z', fontsize=25) 
    if veto: 
        fig.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)), '_random_zhist.veto.png'), bbox_inches='tight') 
    else: 
        fig.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)), '_random_zhist.png'), bbox_inches='tight') 


    # compare n(z)
    nbar_g = UT.get_nofz(np.array(gals['Z']), gals.attrs['fsky'], cosmo=gals.cosmo) 
    nbar_r = UT.get_nofz(np.array(rand['Z']), gals.attrs['fsky'], cosmo=gals.cosmo) 
    ng_nr = gals['Z'].shape[0] / rand['Z'].shape[0]
    
    fig = plt.figure(figsize=(6,6))
    sub = fig.add_subplot(111)
    sub.scatter(np.array(rand['Z']), ng_nr * nbar_r, c='k', s=1, rasterized=True, label=r'$N_g/N_r \times$ Random') 
    sub.scatter(np.array(gals['Z']), nbar_g, c='C0', s=1, rasterized=True, label='Galaxies') 
    sub.legend(loc='upper left', fontsize=15, handletextpad=0, markerscale=10)
    sub.set_xlim(0.1, 0.4) 
    sub.set_xlabel('z', fontsize=25) 
    sub.set_ylabel('$n(z)$', fontsize=25) 
    if veto: 
        fig.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)), '_random_nbarz.veto.png'), bbox_inches='tight') 
    else: 
        fig.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)), '_random_nbarz.png'), bbox_inches='tight') 

