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

# apply forward model 
gals_nfc = FM.BOSS(hod, sample='lowz-south', seed=0, veto=True, fiber_collision=False, silent=False)
gals_fc = FM.BOSS(hod, sample='lowz-south', seed=0, veto=True, fiber_collision=True, silent=False)

# compare footprint
fig = plt.figure(figsize=(10,5))
sub = fig.add_subplot(111)
sub.scatter(np.array(gals_nfc['RA']), np.array(gals_nfc['DEC']), c='k', s=1, rasterized=True, label='no fiber coll.') 
sub.scatter(np.array(gals_nfc['RA'])-360., np.array(gals_nfc['DEC']), c='k', s=1, rasterized=True) 
sub.scatter(np.array(gals_fc['RA']), np.array(gals_fc['DEC']), c='C0', s=1, rasterized=True, label='fiber coll.') 
sub.scatter(np.array(gals_fc['RA'])-360, np.array(gals_fc['DEC']), c='C0', s=1, rasterized=True) 
sub.legend(loc='upper right', fontsize=15, handletextpad=0, markerscale=10)
sub.set_xlabel('RA', fontsize=25) 
sub.set_xlim(-50, 70) 
sub.set_ylabel('Dec', fontsize=25) 
fig.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)), '_fibercollisions.png'), bbox_inches='tight') 

# compare zoomed in footprint
fig = plt.figure(figsize=(6,6))
sub = fig.add_subplot(111)
sub.scatter(np.array(gals_nfc['RA']), np.array(gals_nfc['DEC']), c='k', s=2, rasterized=True, label='has fiber coll.') 
sub.scatter(np.array(gals_nfc['RA'])-360., np.array(gals_nfc['DEC']), c='k', s=2, rasterized=True) 
sub.scatter(np.array(gals_fc['RA']), np.array(gals_fc['DEC']), c='C1', s=3, rasterized=True, label='no fiber coll.') 
sub.scatter(np.array(gals_fc['RA'])-360, np.array(gals_fc['DEC']), c='C1', s=3, rasterized=True) 
sub.legend(loc='upper right', fontsize=15, handletextpad=0, markerscale=10)
sub.set_xlabel('RA', fontsize=25) 
sub.set_xlim(-2, 2) 
sub.set_ylabel('Dec', fontsize=25) 
sub.set_ylim(-2, 2) 
fig.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)), '_fibercollisions.zoomed.png'), bbox_inches='tight') 
