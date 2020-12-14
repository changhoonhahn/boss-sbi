'''

module for measuring the observable of a galaxy catalog. Compile all the
different methods for measuring observables here.  


'''
from pyspectrum import pyspectrum as pySpec
############################################################
# These are Pylians3 functions -> are we going to use them?
import Pk_library as PKL
import MAS_library as MASL
import smoothing_library as SL
############################################################

def Plk():
    '''
    '''


def B0k(galaxies, Lbox=1000., Ngrid=360, step=3, Ncut=3, Nmax=40, fft='pyfftw'): 
    ''' Measure the bispectrum monopole using the `pySpectrum` package
   

    Parameters
    ----------
    galaxies : GalaxyCatalog object

    
    Return
    ------
    bisp : dictionary

    '''
    bisp = pySpec.Bk_periodic(galaxies, Lbox=Lbox, Ngrid=Ngrid, step=step,
            Ncut=Ncut, Nmax=Nmax, fft='pyfftw')
    return bisp 



def func_mark(weight,delta_s,p):
    
    '''Functional form for the mark
    '''
    
    mark_m = ((delta_s+1.)/(delta_s+1.+weight+1e-6))**p
    return mark_m


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
    Pk = PKL.Pk(delta_m, BoxSize, axis, MAS, threads)
    return Pk


def Pk_skew(galaxies,Ngrid,Lbox,Rsmooth):
    ''' Measure the redshift-space skew spectra using the `skewspec` package

    Parameters
    ----------
    galaxies : GalaxyCatalog object
    
    Return
    ------
    SkewSpec_arr : np array contaningthe values of k and skew-spectra
    '''
    ## AM: What is the format of input galaxy catalogue? I am assuming it is compatible with nbodykit CataligueSource object
    
    # Given an nbodykit CatalogSource object `cat' (e.g. containing a halo/galaxy
    # catalog), paint the overdensity delta on a 3D mesh using nbodykit.
    delta_mesh = FieldMesh(galaxies.to_mesh(Nmesh=Ngrid, BoxSize=Lbox,
        window='cic', interlaced=False, compensated=False).compute()-1)

    # Make a copy of the density and apply Gaussian smoothing
    delta_mesh_smoothed = FieldMesh(delta_mesh.compute(mode='real'))
    delta_mesh_smoothed = smoothing.GaussianSmoother(Rsmooth).apply_smoothing(
    delta_mesh_smoothed)

    # Compute skew spectra
    #Set line of sight direction
    LOS = numpy.array([0,0,1])
   
    # Get list of the 14 quadratic fields S1-S14 derived from the tree-
    # level galaxy bispectrum in redshift space, see arXiv:2010.14267.
    # If redshift_space_spectra is False, get the 3 skew-spectra derived 
    # from the galaxy bispectrum in real space.
    skew_spectra = SkewSpectrum.get_list_of_standard_skew_spectra(LOS=LOS, redshift_space_spectra=True)
    
    # Compute skew spectrum and store in skew_spec.Pskew
    for skew_spec in skew_spectra:
        skew_spec.compute_from_mesh(mesh=delta_mesh_smoothed,second_mesh=delta_mesh_smoothed,third_mesh=delta_mesh)
   
    ## AM: We can just return the Pskew object, which has the wavenumber and values of skew-spectra as its attribute
    ## otherwise, here is how we can return an array containing the values. What is the prefered format of the output?
    mydtype = [('k', 'f8')]
    for ell in skew_spec.Pskew['default_key'].attrs['poles']:
        for skew_spec in skew_spectra:
            mydtype.append((skew_spec.name, 'f8'))
            arr = np.empty(shape=skew_spec.Pskew['default_key'].poles['k'].shape, dtype=mydtype)
            SkewSpec_arr['k'] = skew_spec.Pskew['default_key'].poles['k']
            for skew_spec in skew_spectra:
                SkewSpec_arr[skew_spec.name] = skew_spec.Pskew['default_key'].poles['power_%d'%ell].real

   return SkewSpec_arr
