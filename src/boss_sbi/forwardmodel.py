''' 


module for forward modeling the BOSS survey: i.e. python version of
mksample 


'''
import os 
import pymangle 
from remap import Cuboid 


def BOSS(galaxies, sample='lowz-south'): 
    ''' Forward model the BOSS survey given a simulated galaxy catalog 
    '''
    assert samples == 'lowz-south', 'only LOWZ SGC has been implemented' 
    # use BoxRemap to transform the volume (https://arxiv.org/abs/1003.3178)
    C = Cuboid(u1=(1,1,0), u2=(0,1,0), u3=(0,0,1))
    
    #xyz = galaxies['Position'] 


    # apply redshift space distortion 

    # compute redshifts 

    # compute RA and Dec 

    # angular mask
    boss_poly = BOSS_mask(sample)
    in_footprint = BOSS_angular(Galaxies.ra, Galaxies.dec, mask=boss_poly)

    # radial mask
    in_nz = BOSS_radial(Galaxies.z)
    
    # apply masks 
    mask = in_footprint & in_nz

    # assign weights? 

    return Galaxies.apply_selection(mask)


def BOSS_mask(sample): 
    ''' read mangle polygon for specified sample 
    '''
    if sample == 'cmass-north': 
        boss_poly = pymangle.Mangle(
                os.path.join(os.environ('BOSSSBI_DIR'), 'boss',
                    'mask_DR12v5_CMASS_North.ply'))
    else: 
        boss_poly = pymangle.Mangle(
                os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dat',
                'boss_geometry_2014_05_28.ply')) 
    return boss_poly


def BOSS_radial(z): 
    ''' Downsample the redshifts to match the BOSS radial selection function 
    '''
    return None 


def BOSS_angular(ra, dec, mask=None): 
    ''' Given RA and Dec, check whether the galaxies are within the angular
    mask of BOSS
    '''
    inpoly = mask.contains(ra, dec)
    return inpoly 
