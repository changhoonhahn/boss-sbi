#!/bin/env/python
'''

test boss_sbi.galaxies

'''
import numpy as np 
from boss_sbi.halos import Quijote_LHC_HR
from boss_sbi import galaxies as Galaxies
from boss_sbi import forwardmodel as FM 

# read in halo catalog 
halos = Quijote_LHC_HR(1, z=0.5)

# get LOWZ HOD parameters
theta_hod = Galaxies.thetahod_lowz_sgc()
#theta_hod['logMmin'] = 13.
#theta_hod['logM1'] = 14.

# apply HOD 
gals = Galaxies.hodGalaxies(halos, theta_hod, seed=0) 

from boss_sbi.remap import Cuboid 
xyz = np.array(gals['Position']) / 1000.

C = Cuboid(u1=(1,1,0), u2=(0,1,0), u3=(0,0,1))

import time

t0 = time.time() 

xyz_t = np.empty(xyz.shape)
for i in range(xyz.shape[0]): 
    xyz_t[i,:] = C.Transform(xyz[i,0], xyz[i,1], xyz[i,2]) # transformed
print(xyz_t) 
print('%f sec' % ((time.time() - t0)))
xyz_t *= 1000.
print('%.1f < x < %.1f' % (xyz_t[:,0].min(), xyz_t[:,0].max()))
print('%.1f < y < %.1f' % (xyz_t[:,1].min(), xyz_t[:,1].max()))
print('%.1f < z < %.1f' % (xyz_t[:,2].min(), xyz_t[:,2].max()))
print()

xyz_t = np.dot(xyz_t, np.array([[0, -1, 0], [1, 0, 0,], [0, 0, 1]]))

print('%.1f < x < %.1f' % (xyz_t[:,0].min(), xyz_t[:,0].max()))
print('%.1f < y < %.1f' % (xyz_t[:,1].min(), xyz_t[:,1].max()))
print('%.1f < z < %.1f' % (xyz_t[:,2].min(), xyz_t[:,2].max()))

xyz_t += np.array([334.45, 738.4, -351.1])[None,:] # translate 
print('%.1f < x < %.1f' % (xyz_t[:,0].min(), xyz_t[:,0].max()))
print('%.1f < y < %.1f' % (xyz_t[:,1].min(), xyz_t[:,1].max()))
print('%.1f < z < %.1f' % (xyz_t[:,2].min(), xyz_t[:,2].max()))

import nbodykit.lab as NBlab
ra, dec, z = NBlab.transform.CartesianToSky(
        xyz_t, 
        gals.cosmo,
        velocity=gals['Velocity'], 
        observer=[0,0,0])

print(np.array(ra))
print(np.array(dec))
print(np.array(z))
    
gals['RA']  = ra
gals['DEC'] = dec 
gals['Z']   = z 


boss_poly = FM.BOSS_mask('lowz-south') 
in_footprint = FM.BOSS_angular(ra, dec, mask=boss_poly)
print('%i of %i in footprint' % (np.sum(in_footprint), len(in_footprint)))

radial_select = FM.BOSS_radial(np.array(z)[in_footprint], seed=0)
print('%i of %i in radial' % (np.sum(radial_select), len(in_footprint)))

rad_sel = np.zeros(len(in_footprint)).astype(bool)
rad_sel[np.arange(len(in_footprint))[in_footprint][radial_select]] = True 

print('%i of %i' % (np.sum(in_footprint & rad_sel), len(in_footprint)))

select_gals = gals[in_footprint & rad_sel]
print(select_gals['RA']) 



