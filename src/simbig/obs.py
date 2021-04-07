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


def Plk_survey(galaxies, randoms, Ngrid=360, dk=0.005, P0=1e4, silent=True):
    ''' Measure galaxy powerspectrum multipoles for a survey geometry using the
    `nbodykit`. This function uses the FKP estmiator for calculating the power 
    spectrum.
   

    Parameters
    ----------
    galaxies : GalaxyCatalog object

    randoms : nbodykit.ArrayCatalog object 

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
    * 04/02/2021: implemented but not yet tested
    '''
    # get nbar(z) for the galaxy and random samples
    nbar_g = UT.get_nofz(np.array(galaxies['Z']), galaxies.attrs['fsky'], cosmo=galaxies.cosmo)
    nbar_r = UT.get_nofz(np.array(randoms['Z']), galaxies.attrs['fsky'], cosmo=galaxies.cosmo)

    # calculate xyz positions
    pos_g = nblab.transform.SkyToCartesian(
            galaxies['RA'], galaxies['DEC'], galaxies['Z'], cosmo=galaxies.cosmo) 
    pos_r = nblab.transform.SkyToCartesian( 
            randoms['RA'], randoms['DEC'], randoms['Z'], cosmo=galaxies.cosmo) 

    Ng = pos_g.shape[0] # number of galaxies 
    Nr = pos_r.shape[0] # number of randoms

    # normalize nbar(z) for randoms 
    nbar_r *= Ng/Nr
    
    # weights 
    w_g = np.ones(Ng) 
    w_r = np.ones(Nr) 

    _gals = nblab.ArrayCatalog({
        'Position': pos_g, 
        'NZ': nbar_g, 
        'WEIGHT': w_g, 
        'WEIGHT_FKP': 1./(1.+nbar_g * P0)
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

    return k, p0k, p2k, p4k


def Plk_box(galaxies, Lbox=1000., Ngrid=360, k_bin_width=1):
    ''' Measure galaxy powerspectrum multipoles for a galaxy sample in box using 
    `nbodykit`.
   

    Parameters
    ----------
    galaxies : GalaxyCatalog object

    
    Return
    ------
    pk_moms_arr: an array containing values of k and power spectrum moments

    '''
    dk = 2.0 * np.pi / boxsize * k_bin_width
    kmin = 2.0 * np.pi / boxsize / 2.0

    # Set line of sight direction
    LOS = numpy.array([0,0,1])
    
    # paint galaxies to mesh
    delta_mesh = FieldMesh(galaxies.to_mesh(Nmesh=Ngrid, BoxSize=Lbox,
        window='cic', interlaced=False, compensated=False).compute()-1)

    #compute the power spectrum moments using nbodykit 
    if poles is None:
      poles = [0,2,4]
      pk_moms = FFTPower(first=delta_mesh,
                          second=None,
                          mode='2d',
                          dk=dk,
                          kmin=kmin,
                          poles=poles,
                          Nmu=5,
                          los=LOS)

    # if we wanted to return an array of k and pkl values, we can return PkMoms_arr
    for ell in pk_moms.attrs['poles']:
            mydtype = [('k', 'f8'), ('P', 'f8')]
            PkMoms_arr = np.empty(shape=pk_moms.poles['k'].shape, dtype=mydtype)
            PkMoms_arr['k'] = pk_moms.poles['k']
            PkMoms_arr['P'] = pk_moms.poles['power_%d'%ell].real
            
    return PkMoms_arr


def B0k(galaxies, Lbox=1000., Ngrid=360, step=3, Ncut=3, Nmax=40, fft='pyfftw'): 
    ''' Measure the bispectrum monopole using the `pySpectrum` package
   

    Parameters
    ----------
    galaxies : GalaxyCatalog object

    
    Return
    ------
    bisp : dictionary

    '''
    from pyspectrum import pyspectrum as pySpec
    bisp = pySpec.Bk_periodic(galaxies, Lbox=Lbox, Ngrid=Ngrid, step=step,
            Ncut=Ncut, Nmax=Nmax, fft='pyfftw')
    return bisp 


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
