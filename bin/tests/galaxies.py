#!/bin/env/python
'''

test boss_sbi.galaxies

'''
import numpy as np 
from boss_sbi.halos import Quijote_LHC_HR
from boss_sbi import galaxies as Galaxies

# read in halo catalog 
halos = Quijote_LHC_HR(1, z=0.5)

# get LOWZ HOD parameters
theta_hod = Galaxies.thetahod_lowz_ngc()

# apply HOD 
gals = Galaxies.hodGalaxies(halos, theta_hod, seed=0) 
print(np.array(gals['Position']))
