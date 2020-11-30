''' 


module for forward modeling the BOSS survey: i.e. python version of
mksample 


'''
import os 
import pymangle 

# BOSS Mangle Polygon 
boss_poly = pymangle.Mangle(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dat',
        'boss_geometry_2014_05_28.ply')
        ) 



def BOSS(galaxies): 
    ''' Forward model the BOSS survey given GalaxyCatalog object 
    '''

    # angular mask
    in_footprint = BOSS_angular(Galaxies.ra, Galaxies.dec)

    # radial mask
    in_nz = BOSS_radial(Galaxies.z)
    
    # apply masks 
    mask = in_footprint & in_nz

    return Galaxies.apply_selection(mask)


def BOSS_radial(z): 
    ''' Downsample the redshifts to match the BOSS radial selection function 
    '''
    return None 


def BOSS_angular(ra, dec): 
    ''' Given RA and Dec, check whether the galaxies are within the angular
    mask of BOSS
    '''
    inpoly = boss_poly.contains(ra, dec)
    return inpoly 


'''
    below is the original BOSS DR12 mksample code https://data.sdss.org/sas/dr12/boss/lss/mksampleDR12/
    import re
    import os
    import get_collate_targets
    import match_sample
    import numpy as np
    import sys
    import mkweightfiles
    import make_catalog_z
    import time
    import fitsio
    import maskutils
    import fibercollidezfail

    def photokey(darr):
      """
      SDSS magic photoid is some combo of these things, so they can be used as a unique key for each object.
      """
      return (darr['RUN'], darr['RERUN'], darr['CAMCOL'], darr['FIELD'], darr['ID'])

    def cleanval(val):
      """
      Strip off comments and spaces.
      """
      myval = ''
      myval = myval + val  ## this should function to copy val rather than point to it(?)
      if re.search('#',val):
        myval = val.split('#')[0]
      myval = myval.strip()  ## strips spaces, tabs, returns, etc.
      return myval

    def getlist(valstr,defaultlist):
      try:
        valstr = valstr.split('[')[1].split(']')[0]
        mylist = []
        for vv in valstr.split(','):
          if type(defaultlist[0]) is str:
            mylist.append(str(cleanval(vv)))
          if type(defaultlist[0]) is int:
            mylist.append(int(cleanval(vv)))
          if type(defaultlist[0]) is float:
            mylist.append(float(cleanval(vv)))
      except:
        print 'input not a list!!  aborting'
        print valstr
        sys.exit(1)

      return mylist

    def runtests(runparams,sample,fkeep):
      """
      This function knows all the file names and reports on sanity of each one.
      Examines one sample at a time.
      fkeep is the fraction of randoms retained after all the veto masks.
      for cmass dr12v4, fkeep = [0.93332873, 0.90672223] (grep fkeep OUT*dr12v4)
      for lowz dr12v4, fkeep = [0.93325963, 0.90661381]
      """
      trimmedf = runparams['data_dir'] + "trimmed-collate-%s-%s.fits" % (sample, runparams['runtag'])
      dd = fitsio.read(trimmedf)
      print 'total number of initial targets',len(dd)

      ## check mask stuff.
      m1 = maskutils.mask(runparams['geom_file'])
      m1.maskstats(verbose=1)
      if sample == 'cmass':
        b1 = maskutils.mask(runparams['data_dir']+'binarymask-dr12v4.fits')
      else:
        b1 = maskutils.mask(runparams['data_dir']+'binarymask-dr12v4-lowzchunkcut.fits')
      b1.maskstats(verbose=1)

      mSfname = runparams['data_dir']+"mask-%s-%s-%s-%s.fits" % (sample,runparams['runtag'],'S',runparams['cattypetag'])
      mdat = []
      for fkeepval, NS in zip(fkeep,['N','S']):
        mfname = runparams['data_dir']+"mask-%s-%s-%s-%s.fits" % (sample,runparams['runtag'],NS,runparams['cattypetag'])
        mm = maskutils.mask(mfname)
        print 'maskstats ',NS
        mdat.append(mm.maskstats(verbose=1))

      degperstr = (180./np.pi)**2
      print 'total N/S mask area excluding vetos',mdat[0][0,0]*fkeep[0]*degperstr,mdat[1][1,0]*fkeep[1]*degperstr
      print 'total N/S effective area excluding vetos',mdat[0][0,1]*fkeep[0]*degperstr,mdat[1][1,1]*fkeep[1]*degperstr

      for NS in ['N','S']:
        collided_file = runparams['data_dir'] + "collided-%s-%s-%s-%s.fits" % (sample,runparams['runtag'],NS,runparams['cattypetag'])
        ## keep only galaxies inside the binary mask.
        dd = fitsio.read(collided_file)
        xx = np.where(b1.weight[dd['ISECT']] > 0.0001)[0]
        dd = dd[xx]
        print 'working on',collided_file
        print 'imatch:',np.histogram(dd['IMATCH'], bins=np.arange(-0.5,10.5,1.0))
        zfail = np.where(dd['IMATCH'] == 5)[0]
        icol = np.where(dd['IMATCH'] == 3)[0]
        print 'total zfail correction matches?',dd['WEIGHT_CP'][zfail].sum(),(dd['WEIGHT_NOZ'][:] - 1.0).sum()
        print 'total close pair correction matches?',len(icol), (dd['WEIGHT_CP'][:]-1.0).sum()
        if sample == 'lowz':
          assert (np.fabs(dd['WEIGHT_NOZ'][:] - np.floor(dd['WEIGHT_NOZ'][:])) < 2.0e-6).all()
          print 'lowz sample has integer weight_noz values'

      ## some stats on galaxy
      for NS in ['N','S']:
        catalog_base = runparams['data_dir'] + "%s-%s-%s-%s" % (sample,runparams['runtag'],NS,runparams['cattypetag'])
        print 'working on ',catalog_base
        dd = fitsio.read(catalog_base+'.dat.fits')
        print 'total gals/cp weights',len(dd),(dd['WEIGHT_NOZ'] + dd['WEIGHT_CP'] - 1.0).sum()
        print 'so shot noise goes up by ',len(dd)/float((dd['WEIGHT_NOZ'] + dd['WEIGHT_CP'] - 1.0).sum())
        assert (np.fabs(dd['WEIGHT_SYSTOT'] - dd['WEIGHT_SEEING']*dd['WEIGHT_STAR']) < 2.0e-5).all()
        print 'sys weights mean/std:'
        for ff in ['WEIGHT_SEEING','WEIGHT_STAR','WEIGHT_SYSTOT']:
          print ff.split('_')[1], dd[ff].mean(), dd[ff].std()
        rr = fitsio.read(catalog_base+'.ran.fits')
        ## did zindx work?
        try:
          assert (np.fabs(rr['Z'] - dd['Z'][rr['ZINDX']]) < 2.0e-5).all()
        except:
          print 'still failing ZINDX, need to set that up after the fact.'

        ## make sure our seeing and extinction cuts in the mask propagated!
        di = np.where((dd['EB_MINUS_V'] > 0.15) | (dd['PSF_FWHM'][:,1] > 2.3) | (dd['PSF_FWHM'][:,2] > 2.1) | (dd['PSF_FWHM'][:,3] > 2.0))[0]
        ri = np.where((rr['EB_MINUS_V'] > 0.15) | (rr['PSF_FWHM'][:,1] > 2.3) | (rr['PSF_FWHM'][:,2] > 2.1) | (rr['PSF_FWHM'][:,3] > 2.0))[0]
        print 'these should be really close to 0 if we got the badfield mask correct',len(di),len(ri),len(di)/float(len(dd)), len(ri)/float(len(rr))



      return 0



    def read_mksample_params(pfname,logfp):

      ## set up defaults for all parameters.
      runparams = {}

      #############################################
      ## generic parameters.
      runparams['runtag'] = 'dr12test'   ## append to catalog file names to identify catalog run.
      runparams['runtagCMASS'] = 'dr12vblah' ## the tag for the LOWZ component of the combined catalog.
      runparams['nrancats'] = 1          ## number of random catalogs to generate.
      runparams['RoDfac'] = 52.0          ## raw number of randoms to generate compared to number of data, where this number of randoms is before veto (so increase by ~5% from final value)
      runparams['RoDfacsys'] = 10.5          ## raw number of randoms to generate compared to number of data, where this number of randoms is before veto (so increase by ~5% from final value)
      runparams['samplerngseeds'] = [578708690, 2072652973, 2069565457,1147349261, 1125647807, 844666025, 2099787806]
      runparams['catalog_type'] = 0
      runparams['cattypetag'] = 'Anderson'
      runparams['zminsys'] = 0.43
      runparams['zmaxsys'] = 0.7
      runparams['legacynbropt'] = 1   ## set to 1 to allow legacy objects to be nearest neighbors of collided BOSS gals.
      runparams['dolegacycollide'] = 0
      runparams['ifib2zfailweight'] = 0
      runparams['sysseewgtfile'] = None
      runparams['sysstwgtfile'] = None
      runparams['nbar_new_omegam'] = 0.274
      runparams['nbar_new_tag'] = None
      runparams['P0'] = 20000.0
      #############################################
      ## which catalogs to generate (for each sample, choose N, S, or both.
      runparams['cmassopt'] = [0,0]          ## generate cmass catalog
      runparams['lowzopt'] = [0,0]           ## generate lowz catalog
      runparams['lowzearly2opt'] = [0,0]      ## generate lowz early catalog
      runparams['lowzearly3opt'] = [0,0]      ## generate lowz early catalog
      runparams['cmasslowzopt'] = [0,0]
      runparams['cmasslowzearly2opt'] = [0,0]
      runparams['cmasslowzearly3opt'] = [0,0]
      #############################################
      ## tasks to complete in this run.  Useful if you only need to rerun certain portions of the code.
      runparams['rungetcollate'] = 1     ## rerun getcollate, necessary when there's a new collate file and/or when veto masks change.
      runparams['runbinarymask'] = 1
      runparams['runspallmatch'] = 1
      runparams['runfiberczfail'] = 1
      runparams['runcombinestep'] = 0
      runparams['runmakemask'] = 1
      runparams['makesysrandoms'] = 1
      runparams['getsysrandoms'] = 1
      runparams['runcatradec'] = 1
      runparams['runcatz'] = 1
      runparams['fillranBOSSim'] = 2
      runparams['nbar_newcosmo'] = 0
      runparams['nbar_reassign'] = 0
      runparams['dolowztrim'] = 0
      runparams['cptodir'] = None  ## give the directory if you want to automatically copy to SAS dir.
      runparams['docombineN'] = None ## list the three runtags you want to include.
      #############################################
      ## files and paths.
      ## two choices: a full path (beginning with '/') or we assume the location should be prefaced by input_dir
      ## must have 'file' in the parameter name to do this check.
      runparams['input_dir'] = 'inputfiles/'
      runparams['data_dir'] = '../data/'
      runparams['collatefile'] = 'bosstile-final-collated-boss2-boss38.fits'
      runparams['collatefilephoto'] = 'bosstile-final-collated-boss2-boss38-galimaging.fits'
      runparams['known_file'] = 'spAll_primary_good.fits'
      runparams['spall_file'] = 'spAll-v5_7_0.fits'
      runparams['platelist_file'] = 'platelist.fits'
      runparams['geom_file'] = '/home/bareid/boss/bosslss/trunk/geometry/boss_geometry_2014_05_28.fits'
      runparams['collided_file_for_targcat'] = None
      ## add more here.

      mykeys = runparams.keys()

      ifpp = open(pfname,'r')
      for line in ifpp:
        ## skip comment lines.
        if re.match('^#',line):
          continue
        if not re.search('=',line):
          print "Can't interpret this line, aborting!!"
          sys.exit(1)
        if not len(line.split('=')) == 2:
          print 'only one equals per line!'
          sys.exit(1)

        mykey = cleanval(line.split('=')[0])
        myval = cleanval(line.split('=')[1])
        if not (mykey in runparams.keys()):
          print 'this key value not found, aborting!', mykey
          print 'this is the list of keys you have coded:'
          print runparams.keys()
          sys.exit(1)

        ## use default values for typing.
        if type(runparams[mykey]) is str:
          runparams[mykey] = str(myval)
        if type(runparams[mykey]) is int:
          runparams[mykey] = int(myval)
        if type(runparams[mykey]) is float:
          runparams[mykey] = float(myval)
        if type(runparams[mykey]) is list:
          runparams[mykey] = getlist(myval,runparams[mykey])

        ## additionally check for the optional weight files, whose default is None and can't be typed.
        if mykey == 'sysseewgtfile' or mykey == 'sysstwgtfile' or mykey == 'collided_file_for_targcat' or mykey == 'nbar_new_tag' or mykey == 'cptodir':
          runparams[mykey] = str(myval)
        if mykey == 'docombineN':
          runparams[mykey] = getlist(myval,['dr12vX','dr12vX','dr12vX'])

        print 'setting key/val:',mykey,runparams[mykey]

      ## do check on filenames; append data_dir/ input_dir as necessary.
      for mykey, myval in runparams.iteritems():
        if re.search('file',mykey):
          if myval is None: continue
          if mykey == 'sysseewgtfile' or mykey == 'sysstwgtfile': continue
          if re.match('^/',myval):
            print 'leaving',myval,'for key',mykey,'as is'
          else:
            if mykey == 'collided_file_for_targcat':
              print 'appended data_dir to key',mykey
              runparams[mykey] = str(runparams['data_dir']) + str(myval)
              print runparams[mykey]

            else:
              print 'appended input_dir to key',mykey
              runparams[mykey] = str(runparams['input_dir']) + str(myval)
              print runparams[mykey]

      print 'Here is the dictionary of inputs to mksample'
      ## print values to logfile.
      for mykey, myval in runparams.iteritems():
        logfp.write('%s = %s\n' % (str(mykey),str(myval)))

      ## check that ifib2 weights are on for CMASS.
      if (np.array(runparams['cmassopt']) > 0).any() and (runparams['runfiberczfail'] == 1):
        assert runparams['ifib2zfailweight'] == 1

      ## check that collided_base stays None if not doing cattype 5.
      if runparams['catalog_type'] != 5:
        assert runparams['collided_file_for_targcat'] is None

      return runparams

    def cleanup(data_dir,dryrun=1):
      blahlist = ['nbar*-sys.dat', 'nbar*ran[0-9][0-9].dat','*-ran[0-9][0-9]*dat.fits','*sys-full.dat.fits']
      rmstr = 'rm '
      for blah in blahlist:
        mystr = 'ls -l %s/%s' % (data_dir,blah)
        os.system(mystr)
        rmstr = rmstr + ' %s/%s' % (data_dir,blah)
      print 'run this',rmstr
      #os.system(rmstr)

    def writesubheader(ofps,jobname,jobtimehours):
      ofps.write('#PBS -N %s\n' % jobname)
      ofps.write('#PBS -l nodes=1:ppn=1,walltime=%s:00:00\n' % (jobtimehours))
      ofps.write('#PBS -o %s.$PBS_JOBID.out\n' % (jobname))
      ofps.write('#PBS -e %s.$PBS_JOBID.err\n' % (jobname))
      ofps.write('#PBS -V\n')
      ofps.write('#\ncd $PBS_O_WORKDIR\n')


    if __name__ == '__main__':

      qsublist = []

      ## BOSS fiber collision scale.
      fbrad = 62./3600.

      logfp = open(sys.argv[1]+'.log','a')
      logfp.write('*'*50+'\n')
      logfp.write('Starting a run at time '+time.strftime("%c")+'\n')
      runparams = read_mksample_params(sys.argv[1],logfp)

      idl_catalog_type = runparams['catalog_type']
      if idl_catalog_type == 3:
        idl_catalog_type = 0
      if idl_catalog_type == 4:
        idl_catalog_type = 1
      if idl_catalog_type == 5:
        idl_catalog_type = 2


      runmatrix = np.zeros([7,2],dtype='int')  ## first index is sample in cmass, lowz, lowzearly2, lowzearly3,cmass+lowz, cmass+lowzearly2,cmass+lowzearly3 order.  Second index is N/S.
      runmatrix[:,:] = 1 ## default is to run everything! make this input from a text file later.

      runmatrix[0,:] = np.array(runparams['cmassopt'])  # this should be a list.
      runmatrix[1,:] = np.array(runparams['lowzopt'])
      runmatrix[2,:] = np.array(runparams['lowzearly2opt'])
      runmatrix[3,:] = np.array(runparams['lowzearly3opt'])
      runmatrix[4,:] = np.array(runparams['cmasslowzopt'])
      runmatrix[5,:] = np.array(runparams['cmasslowzearly2opt'])
      runmatrix[6,:] = np.array(runparams['cmasslowzearly3opt'])
      logfp.write('run matrix!\n')
      logfp.write('%s\n' % str(runmatrix))

      ## start hard coded choices.
      ##
      ## use "runsysrandoms" instead.
      sampletag = ['cmass','lowz','lowzearly2','lowzearly3-6','cmasslowz','cmasslowzearly2','cmasslowzearly3-6']
      combinedcatlist = ['cmasslowz','cmasslowzearly2','cmasslowzearly3-6']
      combinedlowzchunkcutlist = [1,2,3]
      lowzchunkcutlist = [0,1,2,3,1,2,3]  ## cut out lowz chunk < 7, but keep them for hte lowzearly catalog.
      NStaglist = ['N','S']
      NSintlist = [0,1]  ## input to match_sample corresponding to N (0) or S (1).
      mincomp = 0.7
      minzcomp = 0.8
      imagingdropt = 1  ## set to 1 for BOSS [uses BOSS targeting geometry and resolve evolution], 0 for eBOSS [final DR8 imaging everywhere in the sky]
      ## end hard coded choices.

      assert runparams['catalog_type'] == 0 or runparams['catalog_type'] == 3 or \
             runparams['catalog_type'] == 4 or runparams['catalog_type'] == 5
             # more coding to add other options!

      rdtargtype = np.dtype([('RA', '>f8'), ('DEC', '>f8'),('IPOLY', '>i4'), ('ISECT', '>i4')])

      ## trimmed output file names.
      trimmedflist = []  ## same order as sampletag.
      for sample in sampletag:
        trimmedflist.append(runparams['data_dir'] + "trimmed-collate-%s-%s.fits" % (sample, runparams['runtag']))

      logfp.write('this is trimmedflist: %s\n' % (str(trimmedflist)))

      for si in range(len(sampletag)):
        ## if you are NOT running N and not running S for a particular sample, don't output it from get_collate, so versions don't get overwritten by mistake.
        ## get_collate_targets only writes input files that are not None.
        if (runmatrix[si,:] == 0).all():
          print 'setting to None:',trimmedflist[si]
          trimmedflist[si] = None

      ## run get_collate_targets again if necessary.
      if runparams['rungetcollate'] == 1:
        get_collate_targets.get_collate_targets(collatefile=runparams['collatefile'],\
            collatefilephoto=runparams['collatefilephoto'], \
            cmass_output = trimmedflist[0], \
            lowz_output = trimmedflist[1], \
            lowzearly2_output = trimmedflist[2], \
            lowzearly3_output = trimmedflist[3])

        ## fill in imaging information for systematics.
        ## This depends on BOSSTILELIST_DIR
        ## currently has all the resolve version paths hard-coded for riemann.
        DorRopt = 0
        for trimmedf in trimmedflist:
          if trimmedf is None: continue

          mystr = """idl << EOF\n.com fill_BOSS_imaging\nfill_BOSS_imaging,"%s",%d,%d\nEOF\n""" % (trimmedf.split('.fits')[0],imagingdropt, DorRopt)
          logfp.write('running fill_BOSS_imaging with the following command:\n')
          logfp.write('%s\n' % (mystr))
          t1 = time.time()
          os.system(mystr)
          t2 = time.time()
          logfp.write('time in seconds to run fill_BOSS_imaging for %s: %s\n' % (trimmedf,str(t2-t1)))

      ## run for cattypes 3,4,5
      if runparams['catalog_type'] >= 3 and runparams['catalog_type'] <=4 and runparams['runbinarymask'] == 1:  ## new 'Reid' algorithm for mask and fiber collisions.
        maskutils.makebinary(platelistfits = runparams['platelist_file'], geomfitsfname = runparams['geom_file'], binarymaskfitsout = runparams['data_dir'] + "binarymask-%s.fits" % (runparams['runtag']), splitNSopt = 1, lowzchunkcutopt = 0)
        maskutils.makebinary(platelistfits = runparams['platelist_file'], geomfitsfname = runparams['geom_file'], binarymaskfitsout = runparams['data_dir'] + "binarymask-%s-lowzchunkcut.fits" % (runparams['runtag']), splitNSopt = 1, lowzchunkcutopt = 1)
        maskutils.makebinary(platelistfits = runparams['platelist_file'], geomfitsfname = runparams['geom_file'], binarymaskfitsout = runparams['data_dir'] + "binarymask-%s-lowzearlychunk2targ.fits" % (runparams['runtag']), splitNSopt = 1, lowzchunkcutopt = 2)
        maskutils.makebinary(platelistfits = runparams['platelist_file'], geomfitsfname = runparams['geom_file'], binarymaskfitsout = runparams['data_dir'] + "binarymask-%s-lowzearlychunk3-6targ.fits" % (runparams['runtag']), splitNSopt = 1, lowzchunkcutopt = 3)
        maskutils.makebinary(platelistfits = runparams['platelist_file'], geomfitsfname = runparams['geom_file'], binarymaskfitsout = runparams['data_dir'] + "binarymask-%s-lowzearlychunk2targ-chunk2.fits" % (runparams['runtag']), splitNSopt = 1, lowzchunkcutopt = 4)
        maskutils.makebinary(platelistfits = runparams['platelist_file'], geomfitsfname = runparams['geom_file'], binarymaskfitsout = runparams['data_dir'] + "binarymask-%s-lowzearlychunk3-6targ-chunk3-6.fits" % (runparams['runtag']), splitNSopt = 1, lowzchunkcutopt = 5)



      for si, sample,trimmedf,lowzchunkcut in zip(range(len(sampletag)),sampletag,trimmedflist,lowzchunkcutlist):
        ## seed random number generator for this sample.
        np.random.seed(seed=runparams['samplerngseeds'][si])
        ## plus 1 is to make separate randoms for Ashley's code.
        ## first index is for mkweights, rest are for nrancats.
        ## this will make the results reproducible no matter what nrancats is.
        ransackseedlist = np.random.randint(0,np.iinfo(np.int32).max,[runparams['nrancats']+1,2])
        zseedlist = np.random.randint(0,np.iinfo(np.int32).max,[runparams['nrancats']+1,2])
        ## random seeds for the random fiber collision choices (legacy close pairs)
        fbseedlist = np.random.randint(0,np.iinfo(np.int32).max,[2])

        logfp.write('%s %s %d, seeded random number generator with %d\n' % (sample,trimmedf, lowzchunkcut,runparams['samplerngseeds'][si]))
        logfp.write('ransack seeds: %s\n' % (ransackseedlist))
        logfp.write('zseeds: %s\n' % (zseedlist))

        for nsi in range(len(NStaglist)):
          NStag = NStaglist[nsi]
          NSint = NSintlist[nsi]
          if runmatrix[si,nsi] == 0: continue
          match_file = "match-%s-%s-%s.fits" % (sample,runparams['runtag'],NStag)
          match_output_file = runparams['data_dir'] + match_file ## python and idl have different conventions with data_dir assumption.
          collided_file = "collided-%s-%s-%s-%s.fits" % (sample,runparams['runtag'],NStag,runparams['cattypetag'])
          mask_base = "mask-%s-%s-%s-%s" % (sample,runparams['runtag'],NStag,runparams['cattypetag'])
          catalog_base = "%s-%s-%s-%s" % (sample,runparams['runtag'],NStag,runparams['cattypetag'])
          nbar_base = "nbar-"+catalog_base
          catalog_tag_mkweights = "%s-%s" % (sample,runparams['runtag'])
          catalog_app_mkweights = "-%s-sys.dat.fits" % (runparams['cattypetag'])
          catalog_appr_mkweights = "-%s-sys.ran.fits" % (runparams['cattypetag'])
          catalog_base_mkweights = "%s-%s-%s-%s-sys" % (sample,runparams['runtag'],NStag,runparams['cattypetag'])
          nbar_base_mkweights = "nbar-"+catalog_base_mkweights
          cpmaskpairfits = "maskcp-%s-%s-%s-%s.fits" % (sample,runparams['runtag'],NStag,runparams['cattypetag'])
          qpmmaskfits = "maskqpm-%s-%s-%s-%s.fits" % (sample,runparams['runtag'],NStag,runparams['cattypetag'])
          if lowzchunkcut == 0:
            binarymaskfits = runparams['data_dir'] + "binarymask-%s-%s.fits" % (runparams['runtag'],NStag)
          elif lowzchunkcut == 1:
            binarymaskfits = runparams['data_dir'] + "binarymask-%s-lowzchunkcut-%s.fits" % (runparams['runtag'],NStag)
            ## for cmasslowz catalog, trim cmass back to lowzchunkcut region as well.
            binarymaskfitstrim = runparams['data_dir'] + "binarymask-%s-lowzchunkcut-%s.fits" % (runparams['runtag'],NStag)
          elif lowzchunkcut == 2:
            binarymaskfits = runparams['data_dir'] + "binarymask-%s-lowzearlychunk2targ-%s.fits" % (runparams['runtag'],NStag)
            binarymaskfitstrim = runparams['data_dir'] + "binarymask-%s-lowzearlychunk2targ-chunk2-%s.fits" % (runparams['runtag'],NStag)

          elif lowzchunkcut == 3:
            binarymaskfits = runparams['data_dir'] + "binarymask-%s-lowzearlychunk3-6targ-%s.fits" % (runparams['runtag'],NStag)
            binarymaskfitstrim = runparams['data_dir'] + "binarymask-%s-lowzearlychunk3-6targ-chunk3-6-%s.fits" % (runparams['runtag'],NStag)
          else:
            print 'lowzchunkcut not understood!',lowzchunkcut
            sys.exit(1)

          ## beth is here.
          if sample in combinedcatlist:
            assert runparams['rungetcollate'] == 0
            assert runparams['runbinarymask'] == 0
            assert runparams['runspallmatch'] == 0
            assert runparams['runfiberczfail'] == 0
            assert runparams['runmakemask'] == 0
            assert runparams['makesysrandoms'] == 0
            assert runparams['getsysrandoms'] == 0
          if runparams['runcombinestep'] == 1:
            assert sample in combinedcatlist
            ## read in CMASS, assign systematic weights
            binarymaskfitsCMASS = runparams['data_dir'] + "binarymask-%s-%s.fits" % (runparams['runtag'],NStag)
            collided_file_CMASS = runparams['data_dir'] + "collided-cmass-%s-%s-%s.fits" % (runparams['runtagCMASS'],NStag,runparams['cattypetag'])
            seewgtfcmass = 'seefitscmass-%s-%s.dat' % (runparams['runtagCMASS'],NStag)
            stfcmass = 'nstlinfitscmass-%s128.dat' % (runparams['runtagCMASS'])
            ddC = mkweightfiles.fillfitsfromparams(collided_file_CMASS,stfcmass,seewgtfcmass,sample='cmass',writeopt=0)
            ## restrict to within cmass binary mask, so we'll get the same sys renormalization for all cmass samples.
            bmaskCMASS = maskutils.mask(binarymaskfitsCMASS)
            ## restrict CMASS to binary mask.
            xtmp = np.where(bmaskCMASS.weight[ddC['ISECT']] > 0.001)[0]
            ddC = ddC[xtmp]
            if 0==1: ## testing!  do we match weights in cmass catalogs?
              targ_objid = {}
              for i in range(len(ddC)):
                targ_objid[(ddC['RUN'][i], ddC['RERUN'][i], ddC['CAMCOL'][i], ddC['FIELD'][i], ddC['ID'][i])] = i
              ctmp = fitsio.read(runparams['data_dir'] + "%s-%s-%s-%s.dat.fits" % ('cmass',runparams['runtagCMASS'],NStag,runparams['cattypetag']))
              cntchk = 0
              for i in range(len(ctmp)):
                photoid = photokey(ctmp[i])
                if photoid in targ_objid:
                  for mytmp in ['WEIGHT_STAR','WEIGHT_SEEING', 'WEIGHT_SYSTOT']:
                    assert np.fabs(ctmp[mytmp][i] - ddC[mytmp][targ_objid[photoid]]) < 2.0e-6
                  cntchk += 1
              ## passed this check!
              ## >> checked cmass weights are correct 607515
              print 'checked cmass weights are correct',cntchk
              del ctmp

            ## renormalize systematic weights.
            ## this is not exactly correct because some imatch = 2 objects will be tossed because of completeness; but close enough for these purposes.
            xkeep = np.where((ddC['IMATCH'] == 1) | (ddC['IMATCH'] == 2))[0]
            mymean = (ddC['WEIGHT_SYSTOT'][xkeep]*(ddC['WEIGHT_NOZ'][xkeep] + ddC['WEIGHT_CP'][xkeep] - 1.0)).sum()/((ddC['WEIGHT_NOZ'][xkeep] + ddC['WEIGHT_CP'][xkeep] - 1.0)).sum()
            ddC['WEIGHT_SYSTOT'] = ddC['WEIGHT_SYSTOT']/mymean
            mymeanchk = (ddC['WEIGHT_SYSTOT'][xkeep]*(ddC['WEIGHT_NOZ'][xkeep] + ddC['WEIGHT_CP'][xkeep] - 1.0)).sum()/((ddC['WEIGHT_NOZ'][xkeep] + ddC['WEIGHT_CP'][xkeep] - 1.0)).sum()
            print 'mean sys',mymean,'chk',mymeanchk
            assert np.fabs(mymeanchk - 1.0) < 2.0e-4

            bmask = maskutils.mask(binarymaskfits) ## this will have LOWZ footprint.

            ## restrict CMASS to LOWZ binary mask.
            xtmp = np.where(bmask.weight[ddC['ISECT']] > 0.001)[0]
            ddC = ddC[xtmp]
            targ_objid_CMASS = {}
            for i in range(len(ddC)):
              targ_objid_CMASS[(ddC['RUN'][i], ddC['RERUN'][i], ddC['CAMCOL'][i], ddC['FIELD'][i], ddC['ID'][i])] = i

            ## read in LOWZ catalog.
            lowztag = sample.strip('cmass')
            match_file_LOWZ =  runparams['data_dir'] + "match-%s-%s-%s.fits" % (lowztag,runparams['runtag'],NStag)
            collided_file_tot =  runparams['data_dir'] + "collided-%s-%s-%s-%s.fits" % (sample,runparams['runtag'],NStag,runparams['cattypetag'])
            ddL = fitsio.read(match_file_LOWZ)
            ## restrict LOWZ to binary mask.
            xtmp = np.where(bmask.weight[ddL['ISECT']] > 0.001)[0]
            ddL = ddL[xtmp]

            xkeep = []
            dups = 0
            for i in range(len(ddL)):
              photoid = photokey(ddL[i])
              if photoid in targ_objid_CMASS:
                icmass = targ_objid_CMASS[photoid]
                assert photokey(ddC[icmass]) == photoid
                assert ddL['ISECT'][i] == ddC['ISECT'][icmass]
                dups += 1
              else:
                xkeep.append(i)
            xkeep = np.array(xkeep)
            print 'lowz dups:',dups,len(ddL), len(ddL[xkeep])
            ddL = ddL[xkeep]

            ## rerun lowz fiber collisions with duplicates removed.
            fibercollidezfail.fibercollide(ddL,binarymaskfits,legacyopt=runparams['legacynbropt'])
            fitsio.write(collided_file_tot,ddL,clobber=True)

            if runparams['catalog_type'] == 3 and runparams['dolegacycollide'] == 1:
              fibercollidezfail.legacycollide(ddL,maskcpfits = runparams['data_dir'] + cpmaskpairfits, fbseed = fbseed, fbrad = fbrad)
              ltmp = np.where(ddL['IMATCH'] == 8)[0]
              print 'legacy corrected this many fiber collision galaxies:',len(ltmp)

            tb1 = time.time()
            fibercollidezfail.wgtzfail(ddL,zfailopt=1,ifib2reweightopt=0)
            tb2 = time.time()
            print 'time to run zfail:',tb2-tb1
            print 'final histogram after fibercollidezfail:'
            h, xtmp = np.histogram(ddL['IMATCH'],bins=np.arange(-0.5,10.5,1.0))
            print h
            fitsio.write(runparams['data_dir']+collided_file_tot,ddL,clobber=True)

            if lowztag == 'lowzearly3-6':
              assert NStag == 'N'
              ## assign LOWZ systematic weights, if lowzearly3-6.
              seewgtflowz = 'seefits%s-%s-%s.dat' % (lowztag,runparams['runtag'],NStag)
              stflowz = None
              ddL = mkweightfiles.fillfitsfromparams(collided_file_tot,stflowz,seewgtflowz,sample=lowztag,writeopt=0)
              if 0==1:
                targ_objid = {}
                for i in range(len(ddL)):
                  targ_objid[(ddL['RUN'][i], ddL['RERUN'][i], ddL['CAMCOL'][i], ddL['FIELD'][i], ddL['ID'][i])] = i
                ltmp = fitsio.read(runparams['data_dir'] + "%s-%s-%s-%s.dat.fits" % (lowztag,runparams['runtag'],NStag,runparams['cattypetag']))
                cntchk = 0
                for i in range(len(ltmp)):
                  photoid = photokey(ltmp[i])
                  if photoid in targ_objid:
                    for mytmp in ['WEIGHT_STAR','WEIGHT_SEEING', 'WEIGHT_SYSTOT']:
                      assert np.fabs(ltmp[mytmp][i] - ddL[mytmp][targ_objid[photoid]]) < 2.0e-6
                    cntchk += 1
                print 'checked lowz weights are correct',cntchk
                ## passed test!
                #checked lowz weights are correct 144708
                del ltmp

              ## renormalize systematic weights across binary mask.
              ## this is not exactly correct because some imatch = 2 objects will be tossed because of completeness; but close enough for these purposes.
              xkeep = np.where((ddL['IMATCH'] == 1) | (ddL['IMATCH'] == 2))[0]
              mymean = (ddL['WEIGHT_SYSTOT'][xkeep]*(ddL['WEIGHT_NOZ'][xkeep] + ddL['WEIGHT_CP'][xkeep] - 1.0)).sum()/((ddL['WEIGHT_NOZ'][xkeep] + ddL['WEIGHT_CP'][xkeep] - 1.0)).sum()
              ddL['WEIGHT_SYSTOT'] = ddL['WEIGHT_SYSTOT']/mymean
              mymeanchk = (ddL['WEIGHT_SYSTOT'][xkeep]*(ddL['WEIGHT_NOZ'][xkeep] + ddL['WEIGHT_CP'][xkeep] - 1.0)).sum()/((ddL['WEIGHT_NOZ'][xkeep] + ddL['WEIGHT_CP'][xkeep] - 1.0)).sum()
              print 'mean sys',mymean,'chk',mymeanchk
              assert np.fabs(mymeanchk - 1.0) < 2.0e-4

            ## no sys weights for lowz
            else:
              pass

            ## concatenate catalogs.
            ddCL = np.zeros(len(ddC)+len(ddL),dtype=ddC.dtype)
            ddCL[:len(ddC)] = ddC
            ddCL[len(ddC):] = ddL
            fitsio.write(collided_file_tot,ddCL,clobber=True)
            ## then run fibercollidezfail.calc_cp_comp()
            dummy1, dummy2 = fibercollidezfail.calc_cp_comp(ddCL,binarymaskfits,cpmaskqpmfits = runparams['data_dir'] + qpmmaskfits, cpmaskpairfits = runparams['data_dir'] + cpmaskpairfits,combinedopt=1)
            ## maskutils.makemaskcomp()
            mcomp = maskutils.makemaskcomp(runparams['data_dir']+collided_file_tot,binarymaskfits,runparams['data_dir']+mask_base,zfailopt = -1, verboseopt = 1, combinedopt = 1) ## don't toss sectors based on redshift completeness.
            maskutils.maskfitstoply(runparams['data_dir']+mask_base+'.fits',runparams['geom_file'])


          if runparams['catalog_type'] == 5:
            ##
            assert runparams['collided_file_for_targcat'] is not None
            if not re.search('targ',runparams['cattypetag']):
              print 'mksample convention is to have "targ" in the cattypetag for target catalog generation.  Aborting!'
              sys.exit(1)
            if not re.search(runparams['runtag'],runparams['collided_file_for_targcat']):
              print 'collided file should have the same runtag as target file!',runparams['runtag'],runparams['collided_file_for_targcat']
              sys.exit(1)
            if not re.search('-'+NStag+'-',runparams['collided_file_for_targcat']):
              print 'N/S misaligned collided file! ',runparams['collided_file_for_targcat']
              sys.exit(1)

            for ff in ['rungetcollate','runbinarymask','runspallmatch','runfiberczfail','runmakemask','makesysrandoms','getsysrandoms','runcatz']:
              if runparams[ff] != 0:
                print 'catalog_type = 5 not currently set up for this task.  Aborting!',ff
                sys.exit(1)
            if runparams['runcatradec'] != 1:
              print 'Nothing to do!'
              sys.exit(1)

            ## this is so clunky
            cc = fitsio.read(runparams['collided_file_for_targcat'])  ## want this to have the same data structure as all other catalogs.
            bmask = maskutils.mask(binarymaskfits)

            ## trim the galaxy catalog by the binary mask.
            xx = np.where(bmask.weight[cc['ISECT']] > 0.001)[0]
            print 'masked collided file by binary mask.  retained',len(xx),'of',len(cc)
            cc = cc[xx]
            myngal = len(cc)
            fitsio.write(runparams['data_dir']+catalog_base+'.dat.fits',cc,clobber=True)
            del cc
            ## fill in systematics parameters if this is cmass.
            if runparams['sysseewgtfile'] is not None and runparams['sysstwgtfile'] is not None:
              assert si == 0
              mkweightfiles.fillfitsfromparams(runparams['data_dir'] + catalog_base + '.dat.fits', runparams['sysstwgtfile'], runparams['sysseewgtfile'],sample=sample)
            else:
              print 'you need to add systematics for cmass!'
              assert si != 0
            for rcati in range(runparams['nrancats']):
              ransackseed = ransackseedlist[rcati+1,nsi]
              Nran = int(np.floor(myngal*runparams['RoDfac']))
              ra, dec, ipoly = maskutils.runransack(binarymaskfits,runparams['geom_file'],ransackseed=ransackseed, Nran=Nran)
              rr = np.recarray(len(ra),dtype=rdtargtype)
              rr['RA'] = ra
              rr['DEC'] = dec
              tmpfits = "tmp.%s.fits" % (str(ransackseed))
              fitsio.write(tmpfits,rr,clobber=True)
              ## free memory
              del ra
              del dec
              del ipoly
              del rr
              #os.system("""idl -e 'vetofile,"tmp.fits","tmp.fits"'""")
              #print 'yo the string',"""idl -e 'vetofile,""" + """"%s","%s"'""" % (tmpfits,tmpfits)
              os.system("""idl -e 'vetofile,""" + """"%s","%s"'""" % (tmpfits,tmpfits))
              rr = fitsio.read(tmpfits)
              mystr = 'rm %s' % (tmpfits)
              os.system(mystr)
              catalog_base_ii = catalog_base
              if rcati > 0:
                catalog_base_ii = catalog_base + '-ran%02d' % (rcati)
              fitsio.write(runparams['data_dir']+catalog_base_ii+'.ran.fits',rr,clobber=True)



          if runparams['runspallmatch'] == 1:
            match_sample.match_sample(target_file = trimmedf, known_file = runparams['known_file'], NSopt = NSint, LRopt = 0, spall_file = runparams['spall_file'], output_file = match_output_file)

          if runparams['catalog_type'] == 0 and (runparams['runfiberczfail'] == 1 or runparams['runmakemask'] == 1):
            if lowzchunkcut >= 2:
              print 'this catalog type not supported for lowzchunkcut = 2!  Abort.'
              print 'need to modify make_mask accordingly.'
              sys.exit(1)
            ## make a giant IDL string to call fiber collisions and make_mask code.
            myidlstr = """idl << EOF\n.com D_c\n.com fiberc_zfail\n.com make_mask\n.com make_catalog_radec\n.com fill_BOSS_imaging\n"""
            ## now add the rest of the procedures.
            ## fiberc_zfail
            if runparams['runfiberczfail'] == 1:
              myidlstr = myidlstr + """fiberc_zfail,"%s","%s"\n""" % (match_file,collided_file)
            if runparams['runmakemask'] == 1:
              myidlstr = myidlstr + """make_mask,"%s","%s",%s,%s,%d,%d\n""" % (collided_file, mask_base, str(mincomp),str(minzcomp),idl_catalog_type,lowzchunkcut)
            myidlstr = myidlstr + "EOF\n"
            logfp.write('running in IDL:\n%s\n' % (myidlstr))
            os.system(myidlstr)

          if runparams['catalog_type'] >= 3 and (runparams['runfiberczfail'] == 1 or runparams['runmakemask'] == 1):
            assert (runparams['runfiberczfail'] == 1 and runparams['runmakemask'] == 1)  ## they are intertwined!!
            ## only writing these options for now!
            assert runparams['catalog_type'] >= 3 or runparams['catalog_type'] == 4

            fbseed = fbseedlist[nsi]
            if runparams['runfiberczfail'] == 1 and runparams['runmakemask'] == 1:
              dd = fitsio.read(match_output_file)
              if runparams['catalog_type'] == 3:  #only do fiber collisions for type 3.
                ## this has lowzchunkcut applied if you're working on lowz sample.
                fibercollidezfail.fibercollide(dd,binarymaskfits,legacyopt=runparams['legacynbropt'])
                ## see main block of fibercollidezfail.py for generating tables in the catalog paper with fiber collision stats.
                dummy1, dummy2 = fibercollidezfail.calc_cp_comp(dd,binarymaskfits,cpmaskqpmfits = runparams['data_dir'] + qpmmaskfits, cpmaskpairfits = runparams['data_dir'] + cpmaskpairfits)

              ## make mask before legacy collisions happen; write out collided file first!!!
              fitsio.write(runparams['data_dir']+collided_file,dd,clobber=True)
              mcomp = maskutils.makemaskcomp(runparams['data_dir']+collided_file,binarymaskfits,runparams['data_dir']+mask_base,zfailopt = -1, verboseopt = 1) ## don't toss sectors based on redshift completeness.
              maskutils.maskfitstoply(runparams['data_dir']+mask_base+'.fits',runparams['geom_file'])
              ## now run legacy and zfail
              if runparams['catalog_type'] == 3 and runparams['dolegacycollide'] == 1:
                fibercollidezfail.legacycollide(dd,maskcpfits = runparams['data_dir'] + cpmaskpairfits, fbseed = fbseed, fbrad = fbrad)
                ltmp = np.where(dd['IMATCH'] == 8)[0]
                print 'legacy corrected this many fiber collision galaxies:',len(ltmp)

              tb1 = time.time()
              fibercollidezfail.wgtzfail(dd,zfailopt=1,ifib2reweightopt=runparams['ifib2zfailweight'])
              tb2 = time.time()
              print 'time to run zfail:',tb2-tb1
              print 'final histogram after fibercollidezfail:'
              h, xtmp = np.histogram(dd['IMATCH'],bins=np.arange(-0.5,10.5,1.0))
              print h
              fitsio.write(runparams['data_dir']+collided_file,dd,clobber=True)

          ## run systematics with a small set of randoms for both N and S.
          if runparams['makesysrandoms']:
            ransackseed = ransackseedlist[0,nsi]
            zseed = zseedlist[0,nsi]
            print 'running ',catalog_base_mkweights,'for sys weights with random seeds',ransackseed,zseed
            myidlstr = """idl << EOF\n.com D_c\n.com fiberc_zfail\n.com make_mask\n.com make_catalog_radec\n.com fill_BOSS_imaging\n"""
            myidlstr = myidlstr + """make_catalog_radec,"%s","%s","%s","%s", %s, %d, %s, %s\n""" % (collided_file,mask_base,catalog_base_mkweights,nbar_base_mkweights, str(runparams['RoDfacsys']),idl_catalog_type, str(ransackseed), str(zseed))
            DorRoptran = 1
            myidlstr = myidlstr + """fill_BOSS_imaging,"%s",%d,%d\n""" % (catalog_base_mkweights+".ran",imagingdropt, DorRoptran)
            myidlstr = myidlstr + "EOF\n"
            t1 = time.time()
            os.system(myidlstr)
            t2 = time.time()
            print 'time to make_catalog_radec and fill_BOSS_imaging for sys randoms for ',catalog_base_mkweights,':',t2-t1

        ## go outside N/S loop.
        if runparams['getsysrandoms'] and si==0 and runmatrix[si,:].sum() > 0:
          ## for no, sys randoms are only for CMASS.
          assert runmatrix[si,0] == 1 and runmatrix[si,1] == 1  ## need to run N/S at the same time.
          mkweightfiles.mkweightfiles(file=catalog_tag_mkweights,ranopt=1,plotopt=0,input_app=catalog_app_mkweights,input_appr=catalog_appr_mkweights,input_zmin=runparams['zminsys'],input_zmax=runparams['zmaxsys'])
          ## all the weight files get tagged with catalog_tag_weights, e.g., cmass-dr12v1-N (appendix not added, like -Anderson-sys.  So then I think I just run with a different appendix to fill systematic weights later, and ranopt = 0
        if runparams['getsysrandoms'] and si==3 and runmatrix[si,:].sum() > 0:  #lowzearly3-6
          assert runparams['zminsys'] == 0.2 and runparams['zmaxsys'] == 0.5  ## email from Ashley Nov 6 2014
          print 'hi beth'
          mkweightfiles.mkweightfilesLZ36(file=catalog_tag_mkweights,ranopt=1,plotopt=0,input_app=catalog_app_mkweights,input_appr=catalog_appr_mkweights,input_zmin=runparams['zminsys'],input_zmax=runparams['zmaxsys'])


        ## now the systematic weights exist, we can run through and assign them.
        for nsi in range(len(NStaglist)):
          NStag = NStaglist[nsi]
          NSint = NSintlist[nsi]
          if runmatrix[si,nsi] == 0: continue

          collided_file = "collided-%s-%s-%s-%s.fits" % (sample,runparams['runtag'],NStag,runparams['cattypetag'])
          mask_base = "mask-%s-%s-%s-%s" % (sample,runparams['runtag'],NStag,runparams['cattypetag'])
          catalog_base = "%s-%s-%s-%s" % (sample,runparams['runtag'],NStag,runparams['cattypetag'])
          nbar_base = "nbar-"+catalog_base
          catalog_tag_mkweights = "%s-%s" % (sample,runparams['runtag'])
          catalog_app_mkweights = "-%s-sys.dat.fits" % (runparams['cattypetag'])
          catalog_appr_mkweights = "-%s-sys.ran.fits" % (runparams['cattypetag'])
          catalog_base_mkweights = "%s-%s-%s-%s-sys" % (sample,runparams['runtag'],NStag,runparams['cattypetag'])
          nbar_base_mkweights = "nbar-"+catalog_base_mkweights

          ## define the nbar_base_new for use below, even if nbar_newcosmo != 1.
          if runparams['nbar_new_tag'] is not None:
            nbar_base_new = nbar_base + '-' + runparams['nbar_new_tag']
          else:
            nbar_base_new = nbar_base + '-newtest'
          if runparams['nbar_newcosmo'] == 1:
            ## in principle can input a different P0; in practice use the default P0
            make_catalog_z.nbar_newcosmo(nbar_base,nbar_base_new,runparams['nbar_new_omegam'],P0=runparams['P0'])

          for rcati in range(runparams['nrancats']):
            ransackseed = ransackseedlist[rcati+1,nsi]
            zseed = zseedlist[rcati+1,nsi]
            catalog_base_ii = catalog_base
            nbar_base_ii = nbar_base
            mkweights_app_ii = '-%s.dat.fits' % (runparams['cattypetag'])

            ## only give an extra name if its later than the first catalog.
            if rcati > 0:
              catalog_base_ii = catalog_base + '-ran%02d' % (rcati)
              nbar_base_ii = nbar_base + '-ran%02d' % (rcati)
              mkweights_app_ii = '-%s-ran%02d.dat.fits' % (runparams['cattypetag'],rcati)

            if runparams['runcatradec'] == 1:
              if runparams['catalog_type'] == 5: continue
              logfp.write('running %s with random seeds %s %s\n' % (catalog_base_ii,str(ransackseed),str(zseed)))

              myidlstr = """idl << EOF\n.com D_c\n.com fiberc_zfail\n.com make_mask\n.com make_catalog_radec\n.com fill_BOSS_imaging\n"""
              myidlstr = myidlstr + """make_catalog_radec,"%s","%s","%s","%s", %s, %d, %s, %s\n""" % (collided_file,mask_base,catalog_base_ii,nbar_base_ii, str(runparams['RoDfac']),idl_catalog_type, str(ransackseed), str(zseed))
              ## do not run this now, push to the end.
              #DorRoptran = 1
              #myidlstr = myidlstr + """fill_BOSS_imaging,"%s",%d,%d\n""" % (catalog_base_ii+".ran",imagingdropt, DorRoptran)
              myidlstr = myidlstr + "EOF\n"
              logfp.write('%s\n' % (myidlstr))
              t1 = time.time()
              os.system(myidlstr)
              t2 = time.time()
              logfp.write('time to make_catalog_radec for randoms for %s: %s\n' % (catalog_base_ii,t2-t1))

            if runparams['runcatz'] == 1:
              ## need systematics weights filled in for redshift distribution of randoms to be correct.
              if si==0 or si == 3:  ## cmass or lowzearly3-6 only for now!
                ## standard way, fill from the text file in the local directory.
                if si == 0 and (runparams['sysseewgtfile'] is None or runparams['sysstwgtfile'] is None) or \
                   si == 3 and runparams['sysseewgtfile'] is None:
                  mkweightfiles.fillfits(v=catalog_tag_mkweights,NS=NStag,input_app=mkweights_app_ii,sample=sample)
                  ## make sure it worked!
                  ddtmp1 = fitsio.read(runparams['data_dir'] + catalog_base_mkweights+'.dat.fits')
                  ddtmp2 = fitsio.read(runparams['data_dir'] + catalog_base_ii + '.dat.fits')
                  print runparams['data_dir'] + catalog_base_mkweights+'.dat.fits'
                  print runparams['data_dir'] + catalog_base_ii + '.dat.fits'
                  assert (ddtmp1['WEIGHT_SYSTOT'] == ddtmp2['WEIGHT_SYSTOT']).all()
                  logfp.write('passed systematics check, is std nonzero? %s\n' % (str(ddtmp2['WEIGHT_SYSTOT'].std())))
                  del ddtmp1
                  del ddtmp2
                else: # fill hte systematic weights from the parameters files.
                  mkweightfiles.fillfitsfromparams(runparams['data_dir'] + catalog_base_ii + '.dat.fits', runparams['sysstwgtfile'], runparams['sysseewgtfile'],sample=sample)

              ## in principle can input a different P0; in practice use the default P0
              make_catalog_z.make_catalog_z(collided_file, mask_base, catalog_base_ii, nbar_base_ii, runparams['RoDfac'], idl_catalog_type, np.uint64(ransackseed), np.uint64(zseed),P0=runparams['P0'])


            if runparams['nbar_reassign'] == 1:
              ## in practice, make sure we ran nbar file to make sure P0 is consistent between the nbar file and the input parameter to reassign.
              assert runparams['nbar_newcosmo'] == 1
              make_catalog_z.reassign_nbar_fkp(catalog_base_ii,nbar_base_new, P0=runparams['P0'])

            ## trim lowz2 and lowz3-6, save new catalogs.
            ## also trim combined catalogs.
            if runparams['dolowztrim'] == 1 and \
              (sample == 'lowzearly2' or sample == 'lowzearly3-6' or \
               sample == 'cmasslowz' or sample == 'cmasslowzearly2' or \
               sample == 'cmasslowzearly3-6'):
              fbase = runparams['data_dir'] + catalog_base_ii
              print 'trimming',fbase,'using',binarymaskfitstrim
              b1 = maskutils.mask(binarymaskfitstrim)
              for fend in ['.dat.fits','.ran.fits']:
                fendout = '_trimmed' + fend
                ee = fitsio.read(fbase + fend)
                xx = np.where(b1.weight[ee['ISECT']] > 0.0001)[0]
                ee = ee[xx]
                fitsio.write(fbase + fendout,ee,clobber=True)


      #sampletag = ['cmass','lowz','lowzearly2','lowzearly3-6']

            if runparams['cptodir'] is not None:
              if rcati == 0:
                ## move nbar, mask, main catalog, first random catalog
                flist = [nbar_base + '.dat',mask_base+'.fits',mask_base+'.ply',\
                         catalog_base_ii + '.dat.fits', catalog_base_ii + '.ran.fits']
                if sample == 'lowzearly3-6' or sample == 'lowzearly2' or \
                   sample == 'cmasslowz' or sample == 'cmasslowzearly2' or \
                   sample == 'cmasslowzearly3-6':  ## add trimmed files.
                  flist.append(catalog_base_ii + '_trimmed.dat.fits')
                  flist.append(catalog_base_ii + '_trimmed.ran.fits')
              else:  ## just copy random file.
                flist = [catalog_base_ii + '.ran.fits']
                if sample == 'lowzearly3-6' or sample == 'lowzearly2' or \
                   sample == 'cmasslowz' or sample == 'cmasslowzearly2' or \
                   sample == 'cmasslowzearly3-6':  ## add trimmed files.
                  flist.append(catalog_base_ii + '_trimmed.ran.fits')

              for ff in flist:
                mycmd = 'cp %s %s' % (runparams['data_dir']+ff,runparams['cptodir'])
                print mycmd
                os.system(mycmd)

            if runparams['fillranBOSSim'] == 1 or runparams['fillranBOSSim'] == 2 or runparams['fillranBOSSim'] == 3:
              DorRoptran = 1
              myidlstr = """idl << EOF\n.com D_c\n.com fiberc_zfail\n.com make_mask\n.com make_catalog_radec\n.com fill_BOSS_imaging\n"""
              myidlstr = myidlstr + """fill_BOSS_imaging,"%s",%d,%d\n""" % (catalog_base_ii+".ran",imagingdropt, DorRoptran)
              myidlstr = myidlstr + "EOF\n"
              logfp.write('%s\n' % (myidlstr))
              if runparams['fillranBOSSim'] == 1: ## run fill from here.
                t1 = time.time()
                os.system(myidlstr)
                t2 = time.time()
                logfp.write('time to fill_BOSS_imaging for sys randoms for %s: %s\n' % (catalog_base_ii,str(t2-t1)))
              else:  ## just write riemann sub files to run them in parallel.

                sfname = catalog_base_ii + '.sub'
                qsublist.append(sfname)
                ofps = open(sfname,'w')
                jobname = 'fillBOSS%03d' % (int(np.floor(np.random.random()*999)))
                writesubheader(ofps,jobname,jobtimehours=96)
                ofps.write(myidlstr)
                ofps.close()


      if runparams['docombineN'] is not None:
        NStag = NStaglist[0]
        NSint = NSintlist[0]
        combinetag = 'cmasslowztot'
        mtotf = runparams['data_dir'] + "mask-%s-%s-%s-%s" % (combinetag,runparams['runtag'],NStag,runparams['cattypetag'])+'.fits'
        catalog_base_tot = "%s-%s-%s-%s" % (combinetag,runparams['runtag'],NStag,runparams['cattypetag'])
        cptotf = runparams['data_dir'] + "maskcp-%s-%s-%s-%s.fits" % (combinetag,runparams['runtag'],NStag,runparams['cattypetag'])


        ## list of catalogs/randoms to combine.
        ## bigrlist is a list of rcat lists.
        bigclist = [[],[]]
        cnamelist = [runparams['data_dir'] + catalog_base_tot + '.dat.fits',\
                     runparams['data_dir'] + catalog_base_tot + '.ran.fits']

        for rcati in range(1,runparams['nrancats']):
          bigclist.append([])
          cnamelist.append( runparams['data_dir'] + catalog_base_tot + '-ran%02d' % (rcati) + '.ran.fits')

        #combinedcatlist = ['cmasslowz','cmasslowzearly2','cmasslowzearly3-6']
        for si, sample, lowzchunkcut, runtag in zip(range(len(combinedcatlist)), combinedcatlist, combinedlowzchunkcutlist,runparams['docombineN']):
          mask_base = "mask-%s-%s-%s-%s" % (sample,runtag,NStag,runparams['cattypetag'])
          catalog_base = "%s-%s-%s-%s" % (sample,runtag,NStag,runparams['cattypetag'])
          cpmaskpairfits = "maskcp-%s-%s-%s-%s.fits" % (sample,runtag,NStag,runparams['cattypetag'])
          bigclist[0].append(runparams['data_dir'] + catalog_base + '_trimmed.dat.fits')

          for rcati in range(runparams['nrancats']):
            catalog_base_ii = catalog_base
            if rcati > 0:
              catalog_base_ii = catalog_base + '-ran%02d' % (rcati)
            bigclist[rcati+1].append(runparams['data_dir'] + catalog_base_ii + '_trimmed.ran.fits')

          if lowzchunkcut == 0:
            binarymaskfits = runparams['data_dir'] + "binarymask-%s-%s.fits" % (runtag,NStag)
          elif lowzchunkcut == 1:
            binarymaskfits = runparams['data_dir'] + "binarymask-%s-lowzchunkcut-%s.fits" % (runtag,NStag)
            ## for cmasslowz catalog, trim cmass back to lowzchunkcut region as well.
            binarymaskfitstrim = runparams['data_dir'] + "binarymask-%s-lowzchunkcut-%s.fits" % (runtag,NStag)
          elif lowzchunkcut == 2:
            binarymaskfits = runparams['data_dir'] + "binarymask-%s-lowzearlychunk2targ-%s.fits" % (runtag,NStag)
            binarymaskfitstrim = runparams['data_dir'] + "binarymask-%s-lowzearlychunk2targ-chunk2-%s.fits" % (runtag,NStag)

          elif lowzchunkcut == 3:
            binarymaskfits = runparams['data_dir'] + "binarymask-%s-lowzearlychunk3-6targ-%s.fits" % (runtag,NStag)
            binarymaskfitstrim = runparams['data_dir'] + "binarymask-%s-lowzearlychunk3-6targ-chunk3-6-%s.fits" % (runtag,NStag)
          else:
            print 'lowzchunkcut not understood!',lowzchunkcut
            sys.exit(1)

          if si==0:
            mtot = maskutils.mask(runparams['data_dir'] + mask_base+'.fits')
            cptot = maskutils.mask(runparams['data_dir'] + cpmaskpairfits)
          else:
            mtmp = maskutils.mask(runparams['data_dir'] + mask_base+'.fits')
            cptmp = maskutils.mask(runparams['data_dir'] + cpmaskpairfits)
            btrim = maskutils.mask(binarymaskfitstrim)
            ## fill in both sectors and polygons appropriately.
            xsec = np.where(btrim.weight > 0.001)[0]
            xply = np.where(btrim.mf['WEIGHT'] > 0.001)[0]
            mtot.weight[xsec] = mtmp.weight[xsec]
            mtot.mf['WEIGHT'][xply] = mtmp.mf['WEIGHT'][xply]
            assert (mtot.weight[btrim.mf['SECTOR'][xply]] == mtot.mf['WEIGHT'][xply]).all()
            cptot.weight[xsec] = cptmp.weight[xsec]
            cptot.mf['WEIGHT'][xply] = cptmp.mf['WEIGHT'][xply]
            assert (cptot.weight[btrim.mf['SECTOR'][xply]] == cptot.mf['WEIGHT'][xply]).all()

        print 'begin cp check'
        for nn in range(1,5):
          xtmp = np.where(cptot.ntiles == nn)
          print 'unique cp values for ntiles = ',nn,':',np.unique(cptot.weight[xtmp])
        print 'end cp check'

        ## write masks.
        fitsio.write(mtotf,mtot.mf,clobber=True)
        fitsio.write(cptotf,cptot.mf,clobber=True)

        ## make polygon files.
        maskutils.maskfitstoply(mtotf,runparams['geom_file'])
        maskutils.maskfitstoply(cptotf,runparams['geom_file'])

        ## all of these have data_dir prepended already.
        flistsas = [mtotf, mtotf.split('.fits')[0]+'.ply',\
                    cptotf, cptotf.split('.fits')[0]+'.ply']

        ## output concatenated files.
        ## note that ZINDX points to any of the 3 untrimmed galaxy catalogs (depending on which chunk the random is located in).
        ## no obvious way to get around this since the random redshifts can be drawn from objects outside the trimmed region.
        for cattotname, cflist in zip(cnamelist,bigclist):
          #print 'slicing these files into this one'
          #print bigclist
          #print cattotname
          clist = []
          for ctmp in cflist:
            clist.append(fitsio.read(ctmp))
          fitsio.write(cattotname,np.concatenate(tuple(clist)),clobber=True)
          flistsas.append(cattotname)

        ## copy to SAS
        if runparams['cptodir'] is not None:
          for ff in flistsas:
            mycmd = 'cp %s %s' % (ff,runparams['cptodir'])
            print mycmd
            os.system(mycmd)

          ## these are wrong, fixed manually dr12v4 --> dr12v5
          for ff in flistsas:
            ## make symbolic links for SAS.
            newfname = '-S-'.join(ff.split('-N-')).strip(runparams['data_dir'])
            oldfname = 'cmasslowz'.join(newfname.split(combinetag))
            mycmd = 'ln -s %s %s' % (oldfname,newfname)
            print mycmd
            #os.system(mycmd)


      if runparams['fillranBOSSim'] == 3:
        for ss in qsublist:
          mystr = 'qsub %s' % (ss)
          os.system(mystr)
          print mystr




      print 'finished!'
      logfp.close()

'''
