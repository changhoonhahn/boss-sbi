'''

module for measuring the observable of a galaxy catalog. Compile all the
different methods for measuring observables here.  


'''
import numpy as np
from . import util as UT 
# --- nbodykit --- 
import nbodykit.lab as nblab 
from nbodykit.algorithms.fftpower import FFTPower
from nbodykit.source.mesh.field import FieldMesh


def Plk_survey(galaxies, randoms, weights=None, Ngrid=360, dk=0.005, P0=1e4, silent=True):
    ''' Measure galaxy powerspectrum multipoles for a survey geometry using the
    `nbodykit`. This function uses the FKP estmiator for calculating the power 
    spectrum.
   

    Parameters
    ----------
    galaxies : GalaxyCatalog object

    randoms : nbodykit.ArrayCatalog object 

    weights : array_like, optional 
        weights for the galaxies

    Ngrid : int
        grid size for FFT 

    P0 : float 
        P0 value for FKP weights. (default: 1e4) 

    silent : boolean
        If True, the function will print out a bunch of stuff for sanity checks.

    
    Return
    ------
    k, p0k, p2k, p4k
        The power spectrum monopole, quadrupole, and hexadecapole 

    Notes
    -----
    * 05/21/2021: tested; nbar(z) calculation modified
    * 04/02/2021: implemented but not yet tested
    '''
    Ng = len(galaxies) # number of galaxies 
    Nr = len(randoms) # number of randoms
    
    # weights 
    if weights is None: 
        w_g = np.ones(Ng) 
    else: 
        w_g = weights
    w_r = np.ones(Nr) 
    if not silent: print('alpha = %f' % (np.sum(w_g)/np.sum(w_r)))

    # get nbar(z) for the galaxy and random samples
    ng_of_z = UT.get_nofz(np.array(galaxies['Z']), galaxies.attrs['fsky'], cosmo=galaxies.cosmo)
    nbar_g = ng_of_z(np.array(galaxies['Z']))
    nbar_r = ng_of_z(np.array(randoms['Z']))

    # calculate xyz positions
    pos_g = nblab.transform.SkyToCartesian(
            galaxies['RA'], galaxies['DEC'], galaxies['Z'], cosmo=galaxies.cosmo) 
    pos_r = nblab.transform.SkyToCartesian( 
            randoms['RA'], randoms['DEC'], randoms['Z'], cosmo=galaxies.cosmo) 
    
    _gals = nblab.ArrayCatalog({
        'Position': pos_g, 
        'NZ': nbar_g, 
        'WEIGHT': w_g, 
        'WEIGHT_FKP': 1./(1. + nbar_g * P0)
        })

    _rands = nblab.ArrayCatalog({ 
        'Position': pos_r, 
        'NZ': nbar_r,
        'WEIGHT': w_r,
        'WEIGHT_FKP': 1./(1. + nbar_r * P0)
    })

    fkp = nblab.FKPCatalog(_gals, _rands)
    mesh = fkp.to_mesh(Nmesh=Ngrid, nbar='NZ', fkp_weight='WEIGHT_FKP', comp_weight='WEIGHT', window='tsc')

    # compute the multipoles
    r = nblab.ConvolvedFFTPower(mesh, poles=[0,2,4], dk=dk, kmin=0.)
    
    k = r.poles['k'] 
    p0k = r.poles['power_0'].real - r.attrs['shotnoise']
    p2k = r.poles['power_2'].real 
    p4k = r.poles['power_4'].real 
    if not silent: 
        for key in r.attrs: print("   %s = %s" % (key, str(r.attrs[key])))

    return k, p0k, p2k, p4k


def Plk_box(galaxies, Lbox=1000., Ngrid=360, dk=0.005, LOS=[0,0,1]):
    ''' Measure galaxy powerspectrum multipoles for a galaxy sample in box using 
    `nbodykit`.

    Parameters
    ----------
    galaxies : GalaxyCatalog object

    
    Return
    ------
    k, p0k, p2k, p4k
        The power spectrum monopole, quadrupole, and hexadecapole 

    
    Notes 
    -----
    * CHH: I modified the code a bit based on https://nbodykit.readthedocs.io/en/latest/cookbook/fftpower.html

    '''
    # paint galaxies to mesh
    mesh = galaxies.to_mesh(window='tsc', Nmesh=Ngrid, BoxSize=Lbox, 
            compensated=True, position='Position')

    #compute the power spectrum moments using nbodykit 
    pk_moms = FFTPower(mesh, mode='2d', dk=dk, kmin=0., poles=[0,2,4], los=LOS)
    
    k = pk_moms.poles['k'] 
    # apply shot noise correction 
    p0k =  pk_moms.poles['power_0'].real - pk_moms.attrs['shotnoise']
    p2k =  pk_moms.poles['power_2'].real
    p4k =  pk_moms.poles['power_4'].real

    return k, p0k, p2k, p4k 


def B0k_survey(galaxies, randoms, weights=None, P0=1e4, Ngrid=360, Lbox=1400, step=3, Ncut=3, Nmax=40, fft='pyfftw', silent=True):
    ''' Measure the bispectrum monopole for a survey geometry using 
    the `pySpectrum` package.`pySpectrum` uses the Scoccimarro (2015)
    estimator for the bispectrum.
   

    Parameters
    ----------
    galaxies : GalaxyCatalog object
    
    randoms : nbodykit.ArrayCatalog object 

    weights : array_like, optional 
        weights for the galaxies
    
    P0 : float 
        P0 value for FKP weights. (default: 1e4) 

    Ngrid : int
        grid size for FFT 
    
    Lbox : float, optional 
        box size (default: 2600.)

    Nmax : int, optional 
        number of steps to include --- i.e. number of modes. (default: 40) 
    
    Ncut : int, optional 
        k minimum in units of fundamental mode (default: 3)

    step : int, optional 
        step size in units of fundamental mode (defualt: 3) 
    
    fft : string, optional
        specifies which fftw version to use. Options are 'pyfftw' and
        'fortran'. (default: 'pyfftw') 

    silent : boolean
        If True, the function will print out a bunch of stuff for sanity checks.
    
    Return
    ------
    k1, k2, k3, b123, q123
        the bispectrum B(k1, k2, k3) and reduced bispectrum Q(k1, k2, k3) 
    
    Notes
    -----
    * 05/27/2021: CH implemented 
    '''
    # import pyspectrum (see https://github.com/changhoonhahn/pySpectrum for details) 
    from pyspectrum import pyspectrum as pySpec

    Ng = len(galaxies) # number of galaxies 
    Nr = len(randoms) # number of randoms
    
    # weights 
    if weights is None: 
        w_g = np.ones(Ng) 
    else: 
        w_g = weights
    w_r = np.ones(Nr) 
    if not silent: print('alpha = %f' % (np.sum(w_g)/np.sum(w_r)))

    # get nbar(z) for the galaxy and random samples
    ng_of_z = UT.get_nofz(np.array(galaxies['Z']), galaxies.attrs['fsky'], cosmo=galaxies.cosmo)
    nbar_g = ng_of_z(np.array(galaxies['Z']))
    nbar_r = ng_of_z(np.array(randoms['Z']))

    # (RA, DEC, Z) for galaxies and random
    radecz_g = np.array([galaxies['RA'], galaxies['DEC'], galaxies['Z']]) 
    radecz_r = np.array([randoms['RA'], randoms['DEC'], randoms['Z']]) 

    # calculate bispectrum 
    bisp = pySpec.B0_survey(
            radecz_g, nbar_g, w=w_g,                        # galaxies
            radecz_r=radecz_r, nbar_r=nbar_r, w_r=w_r,      # randoms 
            P0_fkp=P0, 
            Ngrid=Ngrid,
            Lbox=Lbox, 
            step=step,
            Ncut=Ncut, 
            Nmax=Nmax, 
            cosmo=galaxies.cosmo,
            fft=fft, 
            nthreads=1, 
            silent=silent)

    k1 = bisp['meta']['kf'] * bisp['i_k1']
    k2 = bisp['meta']['kf'] * bisp['i_k2']
    k3 = bisp['meta']['kf'] * bisp['i_k3']

    b123 = bisp['b123']
    q123 = bisp['q123']

    return k1, k2, k3, b123, q123


def B0k_box(galaxies, Lbox=1400., Ngrid=360, step=3, Ncut=3, Nmax=40, fft='pyfftw', silent=True):
    ''' Measure galaxy bispectrum monopole for a periodic box using `pySpectrum`.

    Parameters
    ----------
    galaxies : GalaxyCatalog object
        galaxies in a periodic box. e.g. output from `simbig.galaxies.hodGalaxies`

    Ngrid : int
        grid size for FFT 
    
    Lbox : float, optional 
        box size (default: 2600.)

    Nmax : int, optional 
        number of steps to include --- i.e. number of modes. (default: 40) 
    
    Ncut : int, optional 
        k minimum in units of fundamental mode (default: 3)

    step : int, optional 
        step size in units of fundamental mode (defualt: 3) 
    
    fft : string, optional
        specifies which fftw version to use. Options are 'pyfftw' and
        'fortran'. (default: 'pyfftw') 

    silent : boolean
        If True, the function will print out a bunch of stuff for sanity checks.
    
    Return
    ------
    k1, k2, k3, b123, q123
        the bispectrum B(k1, k2, k3) and reduced bispectrum Q(k1, k2, k3) 
    
    Notes 
    -----
    * 05/27/2021: CHH, implemented function 
    '''
    # import pyspectrum (see https://github.com/changhoonhahn/pySpectrum for details) 
    from pyspectrum import pyspectrum as pySpec

    # x,y,z position of galaxies
    xyz = np.array(galaxies['Position']).T

    bisp = pySpec.Bk_periodic(xyz, Lbox=Lbox, Ngrid=Ngrid, step=step, Ncut=Ncut, Nmax=Nmax, fft=fft, nthreads=1, silent=silent)
    
    k1 = bisp['meta']['kf'] * bisp['i_k1']
    k2 = bisp['meta']['kf'] * bisp['i_k2']
    k3 = bisp['meta']['kf'] * bisp['i_k3']

    b123 = bisp['b123']
    q123 = bisp['q123']

    return k1, k2, k3, b123, q123



"""
    from skewspec import smoothing
    from skewspec.skew_spectrum import SkewSpectrum

    ############################################################
    # These are Pylians3 functions -> are we going to use them?
    import Pk_library as PKL
    import MAS_library as MASL
    import smoothing_library as SL
    ############################################################

    def Mk(galaxy_pos, Filter, R, p, ds, BoxSize, grid, MAS, threads):
        ''' Measure the marked spectrum using the `Pylians3` package  
        Input:
            galaxy_pos: (N,3) array
            FIlter:     'Top-Hat' or 'Gaussian'
            R:          parameter of the mark: scale to define local density
            p:          parameter of the mark
            ds:         parameter of the mark
            BoxSize
            grid:       scalar: size of the grid where we compute the density
            MAS:        'CIC'
            threads:    scalar
        Output:     
            Pk:         object with power spectrum: k = Pk.k3D
                                                    P0 = Pk.Pk[:,0]
                                                    P2 = Pk.Pk[:,1]
                                                    P4 = Pk.Pk[:,2]
        '''
        
        # calculate delta                                                                                                   
        delta = np.zeros((grid,grid,grid), dtype=np.float32)
        MASL.MA(galaxy_pos, delta, BoxSize, MAS)
        delta /= np.mean(delta, dtype=np.float64);  delta -= 1.0
        # smooth delta                                                                                                      
        W_k = SL.FT_filter(BoxSize, R, grid, Filter, threads)
        delta_smoothed = SL.field_smoothing(delta, W_k, threads)
        # marks                                                                                                             
        weight = np.zeros(galaxy_pos.shape[0], dtype=np.float32)
        MASL.CIC_interp(delta_smoothed, BoxSize, galaxy_pos, weight)
        mark = func_mark(weight,ds,p)
        delta_m = np.zeros((grid,grid,grid), dtype=np.float32)
        MASL.MA(galaxy_pos,delta_m,BoxSize,MAS,W=mark)
        delta_m /= np.mean(delta_m,dtype=np.float32);  delta_m -= 1.0
        # compute marked Pk                                                                                                 
        Pk = PKL.Pk(delta_m, BoxSize, axis, MAS, threads
        return Pk


    def Pk_skew(galaxies, Lbox=1000., Ngrid=360, Rsmooth):
        ''' Measure the redshift-space skew spectra using the `skewspec` package

        Parameters
        ----------
        galaxies : GalaxyCatalog object
        
        Return
        ------
        SkewSpec_arr: an array containing the values of k and the 14 skew spectra
        '''
       
        # Given an nbodykit CatalogSource object `cat' (e.g. containing a halo/galaxy
        # catalog), paint the overdensity delta on a 3D mesh using nbodykit.
        delta_mesh = FieldMesh(galaxies.to_mesh(Nmesh=Ngrid, BoxSize=Lbox,
            window='cic', interlaced=False, compensated=False).compute()-1)

        # Make a copy of the density and apply Gaussian smoothing
        delta_mesh_smoothed = FieldMesh(delta_mesh.compute(mode='real'))
        delta_mesh_smoothed = smoothing.GaussianSmoother(Rsmooth).apply_smoothing(
                              delta_mesh_smoothed)

        # Compute skew spectra
        
        # Set line of sight direction
        LOS = numpy.array([0,0,1])
       
        # Get list of the 14 quadratic fields S1-S14 derived from the tree-
        # level galaxy bispectrum in redshift space, see arXiv:2010.14267.
        # If redshift_space_spectra is False, get the 3 skew-spectra derived 
        # from the galaxy bispectrum in real space.
        skew_spectra = SkewSpectrum.get_list_of_standard_skew_spectra(
            LOS=LOS, redshift_space_spectra=True)
        
        # Compute skew spectrum and store in skew_spec.Pskew
        for skew_spec in skew_spectra:
            skew_spec.compute_from_mesh(mesh=delta_mesh_smoothed, 
                third_mesh=delta_mesh, power_kwargs=power_kwargs,
                store_key='default_key')
       
        # if we wanted to return an array of k and skew_spec values, we can return SkewSpec_arr
        mydtype = [('k', 'f8')]
        for ell in skew_spec.Pskew['default_key'].attrs['poles']:
            for skew_spec in skew_spectra:
                mydtype.append((skew_spec.name, 'f8'))
                arr = np.empty(shape=skew_spec.Pskew['default_key'].poles['k'].shape, dtype=mydtype)
                SkewSpec_arr['k'] = skew_spec.Pskew['default_key'].poles['k']
                for skew_spec in skew_spectra:
                    SkewSpec_arr[skew_spec.name] = skew_spec.Pskew['default_key'].poles['power_%d'%ell].real
        return SkewSpec_arr


    def func_mark(weight,delta_s,p):
        
        '''Functional form for the mark
        '''
        
        mark_m = ((delta_s+1.)/(delta_s+1.+weight+1e-6))**p
        return mark_m
"""
