#!/bin/env/python
'''

test quijote halo readins

'''
import numpy as np 
from boss_sbi.halos import Quijote_LHC_HR

halos = Quijote_LHC_HR(1, z=0.5)
print(np.array(halos['Position']))
