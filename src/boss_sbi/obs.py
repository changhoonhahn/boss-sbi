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

