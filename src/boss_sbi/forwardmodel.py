''' 


module for forward modeling the BOSS survey: i.e. python version of
mksample 


'''
import os 
import numpy as np 
from .remap import Cuboid 
import nbodykit.lab as NBlab


def BOSS(galaxies, sample='lowz-south', seed=0): 
    ''' Forward model the BOSS survey given a simulated galaxy catalog 
    '''
    assert sample == 'lowz-south', 'only LOWZ SGC has been implemented' 
    assert np.all(galaxies.attrs['BoxSize'] == 1000.), 'only supported for 1Gpc/h cubic box'

    # use BoxRemap to transform the volume (https://arxiv.org/abs/1003.3178)
    # at the moment this takes about ~5sec --- but it can definitely be sped up.
    C = Cuboid(u1=(1,1,0), u2=(0,1,0), u3=(0,0,1))
    
    xyz = np.array(galaxies['Position']) / 1000.
    xyz_t = np.empty(xyz.shape)
    for i in range(xyz.shape[0]): 
        xyz_t[i,:] = C.Transform(xyz[i,0], xyz[i,1], xyz[i,2]) # transformed
    xyz_t *= 1000. 
    
    # rotate and translate BoxRemap-ed cuboid 
    xyz_t = np.dot(xyz_t, np.array([[0, -1, 0], [1, 0, 0,], [0, 0, 1]])) # rotate
    xyz_t += np.array([334.45, 738.4, -351.1])[None,:] # translate 
    
    # transform Cartesian to (RA, Dec, z) 
    ra, dec, z = NBlab.transform.CartesianToSky(
            xyz_t, 
            galaxies.cosmo,
            velocity=galaxies['Velocity'], 
            observer=[0,0,0])
    galaxies['RA']  = ra
    galaxies['DEC'] = dec 
    galaxies['Z']   = z 

    # angular mask
    boss_poly = BOSS_mask(sample)
    in_footprint = BOSS_angular(ra, dec, mask=boss_poly)

    # radial mask
    in_nz = BOSS_radial(z[in_footprint], sample=sample, seed=seed)
    in_radial_select = np.zeros(len(ra)).astype(bool) 
    in_radial_select[np.arange(len(ra))[in_footprint][in_nz]] = True
    
    select = in_footprint & in_radial_select

    return galaxies[select]


def BOSS_mask(sample): 
    ''' read mangle polygon for specified sample 
    '''
    import pymangle 
    if sample == 'lowz-south': 
        f_poly = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                'dat', 'mask_DR12v5_LOWZ_South.ply') 
    else: 
        raise NotImplementedError
    boss_poly = pymangle.Mangle(f_poly) 
    return boss_poly


def BOSS_angular(ra, dec, mask=None): 
    ''' Given RA and Dec, check whether the galaxies are within the angular
    mask of BOSS
    '''
    w = mask.weight(ra, dec)
    inpoly = (w > 0.) 
    return inpoly 


def BOSS_radial(z, sample='lowz-south', seed=0): 
    ''' Downsample the redshifts to match the BOSS radial selection function.
    This assumes that the sample consists of the same type of galaxies (i.e. 
    constant HOD), but selection effects randomly remove some of them 

    Notes
    -----
    * nbar file from https://data.sdss.org/sas/bosswork/boss/lss/DR12v5/
    '''
    if sample == 'lowz-south': 
        f_nbar = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                    'dat', 'nbar_DR12v5_LOWZ_South_om0p31_Pfkp10000.dat') 
        zmin, zmax = 0.2, 0.37 
    else: 
        raise NotImplementedError

    # zcen,zlow,zhigh,nbar,wfkp,shell_vol,total weighted gals
    zcen, zlow, zhigh, nbar, wfkp, shell_vol, tot_gal = np.loadtxt(f_nbar, 
            skiprows=2, unpack=True) 
    zedges = np.concatenate([zlow, [zhigh[-1]]])

    ngal_z, _ = np.histogram(np.array(z), bins=zedges)

    # fraction to downsample
    fdown_z = tot_gal/ngal_z.astype(float)

    # impose redshift limit 
    zlim = (z > zmin) & (z < zmax) 

    i_z = np.digitize(z, zedges)
    downsample = (np.random.rand(len(z)) < fdown_z[i_z])

    return zlim #& downsample 
