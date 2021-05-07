import os
import numpy as np 
from simbig import halos as Halos

np.random.seed(918234) 

theta_x_pairs = []
for i in range(1000): 
    # read in halo catalog
    halos = Halos.Quijote_LHC_HR(i, z=0.5)

    # impose random halo mass limit as a proxy for baryonic effect 
    Mlim = np.random.uniform(12.5, 13.0)

    theta_cosmo = Halos.Quijote_LHC_cosmo(i)

    # observable: I'm goign to use Nhalo as a proxy for some observable 
    Nhalos = np.sum(np.array(halos['Mass']) > Mlim)
    
    # (parameter, data) pair
    theta_x = np.concatenate([theta_cosmo, [Mlim], [Nhalos]])
    theta_x_pairs.append(theta_x) 

np.save(os.path.join(os.environ['QUIJOTE_DIR'], 'chang', 'halo_cnf_data.npy'), np.array(theta_x_pairs))
