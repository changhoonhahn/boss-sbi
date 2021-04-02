#!/bin/env/python 
'''

script to test the power spectrum multipole calculation for survey geometry 


'''
import os, time
import numpy as np 
from boss_sbi.halos import Quijote_LHC_HR
from boss_sbi import galaxies as Galaxies
from boss_sbi import forwardmodel as FM 
from boss_sbi import obs as Obs 
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
gals = FM.BOSS(hod, sample='lowz-south', seed=0, veto=True, fiber_collision=False, silent=False)
print('forward model takes %f sec' % ((time.time() - t0)))

# get random 
t0 = time.time() 
rand = FM.BOSS_randoms(gals, sample='lowz-south', veto=True) 
print('randoms take %f sec' % ((time.time() - t0)))

# calculate power spectrum
t0 = time.time() 
k, p0k, p2k, p4k = Obs.Plk_survey(gals, rand, Ngrid=360, dk=0.005, P0=1e4, silent=False)
print('plk take %f sec' % ((time.time() - t0)))


# plot power spectrum 
fig = plt.figure(figsize=(5,5))
sub = fig.add_subplot(111)
sub.plot(k, p0k, c='k', label='monopole')
sub.plot(k, p2k, c='C0', label='quadrupole')
sub.plot(k, p4k, c='C1', label='hexadecapole')

sub.legend(loc='lower left', fontsize=10) 
sub.set_ylabel(r'$P_\ell(k)$', fontsize=25) 
sub.set_yscale('log') 
sub.set_ylim(1e3, 2e5)
sub.set_xlabel('$k$', fontsize=25) 
sub.set_xlim([3e-3, 1.]) 
sub.set_xscale('log') 
fig.savefig(os.path.join(os.path.dirname(os.path.realpath(__file__)), '_plk_survey.png'), bbox_inches='tight') 

