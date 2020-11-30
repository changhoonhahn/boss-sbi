'''

module for measuring the observable of a galaxy catalog. Compile all the
different methods for measuring observables here.  


'''
from pyspectrum import pyspectrum as pySpec


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


