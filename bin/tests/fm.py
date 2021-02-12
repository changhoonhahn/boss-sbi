#!/bin/env/python
'''

test boss_sbi.galaxies

'''
import os, time
import numpy as np 
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

# apply forward model 
t0 = time.time() 
gals = FM.BOSS(hod, sample='lowz-south', seed=0, silent=False)
print('forward model takes %f sec' % ((time.time() - t0)))

# read BOSS sample for comparison 
boss = Galaxies.BOSSGalaxies(sample='lowz-south') 
zlim = (np.array(boss['Z']) > 0.2) & (np.array(boss['Z']) < 0.37)

# compare footprint
fig = plt.figure(figsize=(10,5))
sub = fig.add_subplot(111)
sub.scatter(np.array(gals['RA']), np.array(gals['DEC']), c='C0', s=1, rasterized=True, label='Forward Model') 
sub.scatter(np.array(gals['RA'])-360, np.array(gals['DEC']), c='C0', s=1, rasterized=True) 
sub.scatter(np.array(boss['RA']), np.array(boss['DEC']), c='k', s=1, rasterized=True, label='LOWZ') 
sub.scatter(np.array(boss['RA'])-360., np.array(boss['DEC']), c='k', s=1, rasterized=True) 
sub.legend(loc='upper right', fontsize=15, handletextpad=0, markerscale=10)
sub.set_xlabel('RA', fontsize=25) 
sub.set_xlim(-50, 70) 
sub.set_ylabel('Dec', fontsize=25) 
fig.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)), '_fm_footprint.png'), bbox_inches='tight') 

# compare n(z) 
fig = plt.figure(figsize=(5,5))
sub = fig.add_subplot(111)
_ = sub.hist(np.array(boss['Z'])[zlim], color='k', histtype='step', density=True)
_ = sub.hist(np.array(gals['Z']), color='C0', histtype='step', density=True)
#sub.legend(loc='upper right', fontsize=15, handletextpad=0, markerscale=10)
sub.set_xlabel('redshift', fontsize=25) 
sub.set_xlim(0.15, 0.4) 
sub.set_ylabel('noramlized $n(z)$', fontsize=25) 
fig.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)), '_fm_nz.png'), bbox_inches='tight') 
