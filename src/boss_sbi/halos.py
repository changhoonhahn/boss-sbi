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
    Om, Ob, h, ns, s8 = Quijote_LHC_cosmo(i)
    
    # read halo catalog 
    halos = Quijote.Halos(halo_folder, zsnap, Om=Om, Ob=Ob, h=h, ns=ns, s8=s8, Mnu=0.)
    return halos



def Quijote_LHC_cosmo(i): 
    ''' cosmology look up for LHC realization i at redshift z 

    '''
    fcosmo = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dat',
            'quijote_lhc_cosmo.txt')

    # Omega_m, Omega_l, h, ns, s8
    cosmo = np.loadtxt(fcosmo, unpack=True, usecols=range(5)) 

    Om = cosmo[0][i]
    Ob = cosmo[1][i]
    h  = cosmo[2][i]
    ns = cosmo[3][i]
    s8 = cosmo[4][i]

    return Om, Ob, h, ns, s8

