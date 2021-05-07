''' 


module for forward modeling the BOSS survey: i.e. python version of
mksample 


'''
import os, time
import h5py
import numpy as np 
from .remap import Cuboid 
import nbodykit.lab as NBlab

import pymangle 
from pydl.pydlutils.spheregroup import spherematch


def BOSS(galaxies, sample='lowz-south', seed=0, veto=True, fiber_collision=True, silent=True):
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
    if not silent: t0 = time.time() 
    boss_poly = BOSS_mask(sample)
    in_footprint = BOSS_angular(ra, dec, mask=boss_poly)
    if not silent: print('..applying angular mask takes %.f sec' % (time.time() - t0))

    # veto mask 
    if veto: 
        if not silent: t0 = time.time() 
        in_veto = BOSS_veto(ra, dec) 
        if not silent: print('..applying veto takes %.f sec' % (time.time() - t0))
        in_footprint = in_footprint & ~in_veto
    
    # radial mask
    if sample == 'lowz-south': 
        zmin, zmax = 0.2, 0.37 
        in_radial_select = (z > zmin) & (z < zmax) 
        if not silent: print('..applying radial selection')
    else: 
        raise NotImplementedError
    #if not silent: t0 = time.time() 
    #in_nz = BOSS_radial(z[in_footprint], sample=sample, seed=seed)
    #in_radial_select = np.zeros(len(ra)).astype(bool) 
    #in_radial_select[np.arange(len(ra))[in_footprint][in_nz]] = True
    #if not silent: print('..applying raidal takes %.f sec' % (time.time() - t0))

    select = in_footprint & in_radial_select

    if fiber_collision: # apply fiber collisions
        if not silent: t0 = time.time() 
        _fibcoll = BOSS_fibercollision(np.array(ra)[select], np.array(dec)[select])

        fibcoll = np.zeros(len(ra)).astype(bool) 
        fibcoll[np.arange(len(ra))[select][_fibcoll]] = True

        if not silent: print('..applying fiber collisions takes %.f sec' % (time.time() - t0))
    else: 
        fibcoll = np.zeros(len(ra)).astype(bool) 

    galaxies = galaxies[select & ~fibcoll]
    
    area = BOSS_area(sample=sample, veto=veto)
    fsky = area / (360.**2 / np.pi)
    if not silent: print("..footprint covers %.3f of sky" % fsky)
    galaxies.attrs['fsky'] = fsky 
    return galaxies


def BOSS_mask(sample): 
    ''' read mangle polygon for specified sample 
    '''
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


def BOSS_veto(ra, dec): 
    ''' given RA and Dec, find the objects that fall within one of the veto 
    masks of BOSS. At the moment it checks through the veto masks one by one.  
    '''
    in_veto = np.zeros(len(ra)).astype(bool) 
    fvetos = [
            'badfield_mask_postprocess_pixs8.ply', 
            'badfield_mask_unphot_seeing_extinction_pixs8_dr12.ply',
            'allsky_bright_star_mask_pix.ply',
            'bright_object_mask_rykoff_pix.ply', 
            'centerpost_mask_dr12.ply', 
            'collision_priority_mask_dr12.ply']

    veto_dir = os.path.join(os.environ['QUIJOTE_DIR'], 'chang', 'simbig')
    for fveto in fvetos: 
        veto = pymangle.Mangle(os.path.join(veto_dir, fveto))
        w_veto = veto.weight(ra, dec)
        in_veto = in_veto | (w_veto > 0.)
    return in_veto


def BOSS_fibercollision(ra, dec): 
    ''' apply BOSS fiber collisions 
    '''
    fib_angscale = 0.01722 # 62'' fiber collision angular scale 
    t0 = time.time() 
    m1, m2, d12 = spherematch(ra, dec, ra, dec, fib_angscale, maxmatch=2) 
    print('spherematch takes %f sec' % (time.time() - t0))

    notitself = (d12 > 0.0) 
    
    # only ~60% of galaxies within the angular scale are fiber collided 
    # since 40% are in overlapping regions with substantially lower 
    # fiber collision rates 
    notoverlap = (np.random.uniform(size=len(m1)) > 0.6)

    fibcollided = np.zeros(len(ra)).astype(bool)
    fibcollided[m1[notitself & notoverlap]] = True 
    return fibcollided 


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


def BOSS_area(sample='lowz-south', veto=True): 
    ''' get footprint area of BOSS sample. Currently everything is hardcoded


    Returns
    -------
    footprint of boss sample in sq. degrees


    Notes
    -----
    * veto fraction is hardcoded using values from https://data.sdss.org/sas/dr12/boss/lss/mksampleDR12/nbarutils.py

    '''
    if sample == 'lowz-south': 
        f_poly = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                'dat', 'mask_DR12v5_LOWZ_South.ply') 
    else: 
        raise NotImplementedError
    boss_poly = pymangle.Mangle(f_poly) 

    area = np.sum(boss_poly.areas * boss_poly.weights) # deg^2
    
    if veto: # remove veto mask areas 
        fveto = 0.5*(0.906614+0.906454) # from https://data.sdss.org/sas/dr12/boss/lss/mksampleDR12/nbarutils.py
        '''
        fvetos = [
                'badfield_mask_postprocess_pixs8.ply', 
                'badfield_mask_unphot_seeing_extinction_pixs8_dr12.ply',
                'allsky_bright_star_mask_pix.ply',
                'bright_object_mask_rykoff_pix.ply', 
                'centerpost_mask_dr12.ply', 
                'collision_priority_mask_dr12.ply']
        for fveto in fvetos: 
            veto = pymangle.Mangle(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dat', fveto))
            area -= veto.area 
        '''
        area *= fveto 
    return area 


def BOSS_randoms(boss_gals, sample='lowz-south', veto=True): 
    ''' given forward modeled galaxy catalog (output from the BOSS function) construct 
    accompanying random catalog. 

    Parameters 
    ----------
    boss_gals : catalog object
        catalog output from `forwardmodel.BOSS` function 
    sample : string
        specify which BOSS sample. currently only supports 'lowz-south'
    veto : boolean
        if True, veto mask is applied on the random 

    Returns
    -------
    nbodykit.ArrayCatalog with RA, DEC, and Z of random catalog

    '''
    from astropy.stats import scott_bin_width
    if sample not in ['lowz-south']: 
        raise NotImplementedError('%s not yet supported' % sample)

    if not veto: # without veto mask 
        frand = os.path.join(os.environ['QUIJOTE_DIR'], 'chang', 
                'random_DR12v5_LOWZ_South.hdf5')
    else: # with veto mask
        frand = os.path.join(os.environ['QUIJOTE_DIR'], 'chang', 
                'random_DR12v5_LOWZ_South.veto.hdf5')
    
    # read RA and Dec values 
    rand = h5py.File(frand, 'r') 
    rand_ra     = rand['ra'][...]
    rand_dec    = rand['dec'][...]

    # generate redshifts that match input galaxy redshift distribution. This 
    # implementation is similiar to the implementation in nbodykit 
    w, bins = scott_bin_width(np.array(boss_gals['Z']), return_bins=True)
    hist, edges = np.histogram(np.array(boss_gals['Z']), bins=bins) 
    cutoffs = np.cumsum(hist) / np.sum(hist)

    prng = np.random.uniform(size=len(rand_ra))
    rand_z = edges[:-1][cutoffs.searchsorted(prng)] + w * np.random.uniform(size=len(rand_ra))

    # radial mask
    if sample == 'lowz-south': 
        zmin, zmax = 0.2, 0.37 
        in_radial_select = (rand_z > zmin) & (rand_z < zmax) 
    else: 
        raise NotImplementedError
        
    return NBlab.ArrayCatalog({'RA': rand_ra[in_radial_select], 'DEC': rand_dec[in_radial_select], 'Z': rand_z[in_radial_select]})


def Box_RSD(cat, LOS=[0,0,1], Lbox=1000.):
    ''' Given a halo/galaxy catalog in a periodic box, apply redshift space 
    distortion specified LOS along LOS

    Parameters
    ----------
    cat : CatalogBase
        nbodykit.Catalog object
    LOS : array_like 
        3 element list specifying the direction of the line-of-sight
    Lbox : float
        box size in Mpc/h
    '''
    pos = np.array(cat['Position']) + np.array(cat['VelocityOffset']) * LOS

    # impose periodic boundary conditions for particles outside the box
    i_rsd = np.arange(3)[np.array(LOS).astype(bool)][0]
    rsd_pos = pos[:,i_rsd] % Lbox
    pos[:,i_rsd] = np.array(rsd_pos)
    return pos
