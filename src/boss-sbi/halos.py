'''

module to inferace with different N-body catalogs incl. Quijote 


'''
from .sims import quijote as Quijote 



def Quijote_LHC_HR(i, z=0.5): 
    ''' Read halo catalog from the high resolution Quijote LHC. 


    Parameters
    ---------- 
    i : int 
        ith realization of the Quijote LHC simulations 

    z : float
        redshift of the halo catalog. Quijote halo catalogs are available at
        z = 0, 0.5, 1., 2., and 3.

    Return
    ------
    cat : nbodykit.lab.HaloCatalog 
        Quijote HR LHC halo catalog  
    '''
    # directory that contains the Quijote LHC HR
    halo_folder = os.path.join(os.environ('QUJIOTE_DIR'),
            'Halos/latin_hypercube', 'HR_%i' % i)
    
    # look up cosmology of the LHC realization
    Om, Ol, h, Hz = Quijote_LHC_cosmo(i, z)
    
    # read halo catalog 
    halos = Quijote.Halos(halo_folder, zsnap, Om=Om, Ol=Ol)
    return halos



def Quijote_LHC_cosmo(i, z): 
    ''' cosmology look up for LHC realization i at redshift z 


    '''
    return Om, Ol, h, Hz

