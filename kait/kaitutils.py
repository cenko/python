import pyraf
from pyraf import iraf
import copy, os, shutil, glob, sys, string, re, math, operator, time
import pyfits
from types import *
from mx.DateTime import *

from iqpkg import *
import ephem

# Necessary packages
iraf.images()
iraf.immatch()
iraf.imfilter()
iraf.noao()
iraf.imred()
iraf.ccdred()
iraf.digiphot()
iraf.apphot()

yes=iraf.yes
no=iraf.no
INDEF=iraf.INDEF
hedit=iraf.hedit
imgets=iraf.imgets
imcombine=iraf.imcombine

pyrafdir="python/pyraf/"
pyrafdir_key='PYRAFPARS'

if os.environ.has_key(pyrafdir_key):
    pardir=os.environ[pyrafdir_key]
else:
    pardir=os.environ['HOME']+'/'+pyrafdir

if pardir[-1] != '/':
    pardir += '/'

globclob=yes
globver=yes

# Kait Specific Values
CD11 = CD22 = -2.22e-4
CD12 = CD21 = 0
PIXSCALE = 0.8
AP_GAIN = 3.6
AP_READN = 15.6
FL_GAIN = 4.5
FL_READN = 12.0
TFL = DateTime(2007,5,1)
CRSFX = 'crm'
MASKSFX = 'mask'
SIGMA = 1.5
SATVAL = 50000.0
FWHM = 2.0
REFMAG = 'R2MAG'
ZPTKEY = 'ZEROPT'
ZPUKEY = 'ZEROPTU'

########################################################################

def kaitproc(inlist, clobber=globclob, verbose=globver):

    infiles = iraffiles(inlist)

    # Big loop
    for image in infiles:
    
        # Fix bad header keywords
        root,ext = image.split('.')
        fimg = pyfits.open(image)
        fimg.verify('fix')
        image2 = '%s.fits' % root
        print fimg[0].header
        fimg.writeto(image2)
        #iraf.imcopy(image, image2)

        # Update the WCS keywords
        update_head(image2, ['CD1_1', 'CD2_2', 'CD1_2', 'CD2_1', 'PIXSCALE'], 
                    [CD11, CD22, CD12, CD21, PIXSCALE])
        update_head(image2, ['WCSDIM','LTM1_1','LTM2_2','WAT0_001','WAT1_001',
                    'WAT2_001'], [2,1.0,1.0,'system=image',
                    'wtype=tan axtype=ra','wtype=tan axtype=dec'])
        delete_head(image2, ['CROTA1','CROTA2','CDELT1','CDELT2','CSOURCE'])
    
        # Update RA and Dec to J2000.0
        [ra0, dec0, epoch] = get_head(image2, ['RA','DEC','EPOCH'])
        point = ephem.readdb('Point,f|M|F7,%s,%s,0.0,%f' % (ra0,dec0,epoch))
        point.compute(epoch='2000.0')
        coo=astrocoords(str(point.a_ra),str(point.a_dec))
        update_head(image2, ['RA','DEC','EPOCH','CRVAL1','CRVAL2'], 
                    [coo.sxg()[0],coo.sxg()[1],2000.0,coo.radeg(),coo.dcdeg()])

        # Need Date to determine gain and readnoise (camera switch in May 2007)
        dateobs=get_head(image2,"DATE-OBS")
        tobs=DateTime(int(dateobs[6:]),int(dateobs[3:5]),int(dateobs[:2]))
        if tobs-TFL > 0:
            gain=FL_GAIN
            readn=FL_READN
        else:
            gain=AP_GAIN
            readn=AP_READN
        update_head(image2,['GAIN','READN'],[gain,readn])

        # Clean cosmic rays
        #tmpimg=iraf.mktemp("iqcr")+".fits"
        #check_exist("%s_%s.fits" % (root,CRSFX),"w",clobber)
        #iraf.lacos_im(image2,tmpimg,"%s_%s.fits" % (root,CRSFX),
                      #gain=gain,readn=readn,skyval=0.0,sigclip=4.5,
                      #sigfrac=0.5,objlim=2.0,niter=1)
        #iraf.imdel(image2,verify=no,go_ahead=yes)
        #iraf.imcopy(tmpimg,image2,verbose=no)
        #iraf.imdel(tmpimg,verify=no,go_ahead=yes) 

        # Detect Objects
        iraf.iqobjs(image2,SIGMA,SATVAL,skyval="0.0",masksfx=MASKSFX,
                    wtimage="",fwhm=FWHM,pix=PIXSCALE,aperture=2*FWHM/PIXSCALE,
                    gain=gain,minlim=no,clobber=yes,verbose=no)
        
        # Fit WCS
        iraf.iqwcs(image2,objkey='OBJECT',rakey='RA',deckey='DEC',
                   pixscl=PIXSCALE,pixtol=0.05,starfile='!STARFILE',nstar=40,
                   catalog='web',ubhost="localhost",diffuse=yes,
                   clobber=yes,verbose=verbose)


        # Calculate zero-point
        if get_head(image2,'IQWCS'):
            object=get_head(image2,'OBJECT')
            iraf.iqzeropt(image2,REFMAG,starfile="!STARFILE",
                          catalog=object.replace(" ","_")+".cat",pixtol=3.0,
                          useflags=yes,maxnum=50,method="mean",rejout=1,
                          fencelim=0.50,sigma=1.5,maxfrac=0.25,zptkey=ZPTKEY,
                          zpukey=ZPUKEY,clobber=yes,verbose=verbose)

    print 'Exiting Successfully'
    return
  
#############################################

def kaitsub(refimage,inlist,kernel=9,satval=50000.0,dspace=0,dbkg=0,
            nstamps=2,clobber=globclob):

    infiles = iraffiles(inlist)

    # RA and Dec of center pixel in reference image
    [nx, ny] = get_head(refimage, ['NAXIS1', 'NAXIS2'])
    xcenpix = (nx / 2.0)
    ycenpix = (ny / 2.0)
    [[ra, dec]] = impix2wcs(refimage, xcenpix, ycenpix)

    # Big loop
    for image in infiles:

        # Separate image rootname from extension
        root,ext = image.split('.')
 
        # Make sure WCS fit was successful
        if not check_head(image, 'IQWCS'):
            print 'Skipping image %s: IQWCS not successful' % image
            infiles.remove(image)
            continue
        
        # Remap image using objects detected on both frames and
        # create stamp list
        stars=Starlist(get_head(image,'STARFILE'))
        refstars=Starlist(get_head(refimage,'STARFILE'))
        refstars.pix2wcs(refimage)
        refstars.wcs2pix(image)
        match,refmatch=stars.match(refstars,useflags=yes,tol=10.0)
        nstars=len(match)
        if not (nstars>2):
            print 'Could not find star matches between reference and %s' % image
            infiles.remove(image)
            continue
        refmatch.pix2wcs(image)
        refmatch.wcs2pix(refimage)
        matchfile=open('%s.match' % root, 'w')
        stampfile=open('%s.stamps' % root, 'w')
        for i in range(len(match)):
            matchfile.write('%10.3f%10.3f%10.3f%10.3f\n' % (refmatch[i].xval,
                            refmatch[i].yval,match[i].xval,match[i].yval))
            stampfile.write('%10.3f%10.3f\n' % (refmatch[i].xval,
                            refmatch[i].yval))
        matchfile.close()
        stampfile.close()
        check_exist('%s.geodb' % root, 'w', clobber=clobber)
        iraf.geomap('%s.match' % root,'%s.geodb' % root,1.0,nx,1.0,ny,
                    fitgeom="general",verbose=no,interactive=no)
        check_exist('%s.shift.fits' % root, 'w', clobber=clobber)
        iraf.geotran(image,'%s.shift' % root,'%s.geodb' % root,
                     '%s.match' % root,geometry="geometric",
                     boundary="constant",verbose=no)
 
        # Run the subtraction
        check_exist('%s.sub.fits' % root, 'w', clobber=clobber)
        if os.path.exists("stamps.lis"):
            cmd = '$REDUCTION/hotpants -inim %s.shift.fits -tmplim %s -outim %s.sub.fits -tu %.2f -tg %.2f -tr %.2f -iu %.2f -ig %.2f -ir %.2f -r %.2f -ssf stamps.lis -afssc 0 -n t -ko %i -bgo %i -savexy %s.st -nsx %i -nsy %i' % (root, refimage, root, satval, AP_GAIN, AP_READN, satval, FL_GAIN, FL_READN, kernel, dspace, dbkg, root, nstamps, nstamps)
        else:
            cmd = '$REDUCTION/hotpants -inim %s.shift.fits -tmplim %s -outim %s.sub.fits -tu %.2f -tg %.2f -tr %.2f -iu %.2f -ig %.2f -ir %.2f -r %.2f -ssf %s.stamps -afssc 0 -n t -ko %i -bgo %i -savexy %s.st -nsx %i -nsy %i' % (root, refimage, root, satval, AP_GAIN, AP_READN, satval, FL_GAIN, FL_READN, kernel, root, dspace, dbkg, root, nstamps, nstamps)
            #cmd = '$REDUCTION/pois -k %i %i -s %.2f %.2f -S %s.stamps %s %s.shift.fits %s.sub.fits' % (kernel,kernel,satval,satval,root,refimage,root,root)
        fcmd = os.popen(cmd, 'r')
        sublines = fcmd.readlines()
        fcmd.close()

    # Return
    print 'Exiting successfully'
    return

############################################################

def kaitphot(refimage, inlist, sncoo, aperture=10.0, dx=20, dy=20):

    infiles = iraffiles(inlist)

    # Find location of SN
    [[snra, sndec]] = impix2wcs(refimage, sncoo[0], sncoo[1])
    sn=[Star(name="New_SN",xval=sncoo[0],yval=sncoo[1])]
    snstars=Starlist(stars=sn)
    print "Coordinates of new SN: %s %s (J2000.0)" % (snra, sndec)

    # Grab reference rootname
    refroot,null = refimage.split('.')
    
    # Big loop
    for image in infiles:

        # Image root
        root,null = image.split('.')

        # Grab zeropt and zeropt error from image header
        zp,zpu = get_head(image,[ZPTKEY,ZPUKEY])

        # Detect objects in difference image
        iraf.iqobjs('%s.sub.fits' % root,SIGMA,SATVAL,skyval="0.0",
                    zeropt=float(zp),masksfx=MASKSFX,wtimage="",fwhm=FWHM,
                    pix=PIXSCALE,aperture=2*FWHM/PIXSCALE,gain=FL_GAIN,
                    minlim=no,clobber=yes,verbose=no)

        # See if object was detected
        stars=Starlist('%s.sub.fits.stars' % root)
        match,snmatch=stars.match(snstars,useflags=yes)
        
        # If so, use SN magnitude and statistical error
        if len(match):
            snmag = match[0].mag
            snerr = match[0].magu
   
        # Otherwise calculate upper limit
        else:
            statsec = '[%i:%i,%i:%i]' % (int(sncoo[0]-dx/2), int(sncoo[0]+dx/2),
                                         int(sncoo[1]-dy/2), int(sncoo[1]+dy/2))
            iraf.iterstat('%s.sub%s' % (root, statsec), nsigrej=5, maxiter=5,
                          prin=no,verbose=no)
            sigma=float(iraf.iterstat.sigma)
            snmag = - 2.5 * log10(3 * sigma * sqrt(pi * power(aperture,2))) + \
                    float(zp)
            snerr = 99.999

        # Want UT of Observation
        [dateobs,ut] = get_head(image,['DATE-OBS','UT'])
        tobs=DateTime(int(dateobs[6:]),int(dateobs[3:5]),int(dateobs[:2]),
                      int(ut[:2]),int(ut[3:5]),int(ut[6:]))
        
        # Return SN Mag / UL
        print '%s\t%s %.3f %10.3f%10.3f%10.3f%10.3f' % (root,tobs.Format("%b"),
              (((tobs.second/3600.0)+tobs.minute)/60.0+tobs.hour)/24.0+
              tobs.day,snmag,snerr,float(zp),float(zpu))
   
    print 'Exiting Successfully'
    return

###############################################################################

def kaitlc(infiles, refstars, ot, t0, filters, outfile):

    for image in infiles:
        refstars.wcs2pix(image)
        root,null = image.split('.')
        xyfile=open('%s.xy' % root, 'w')
        for star in refstars:
            xyfile.write('%10.3f%10.3f\n' % (star.xval, star.yval))
        xyfile.close()
        iraf.apphot.phot(image,coords='%s.xy' % root,output='%s.mag' % root,verbose=no)
        ot.wcs2pix(image)
        coofile=open('%s.coo' % root, 'w')
        for star in ot:
           coofile.write('%10.3f%10.3f\n' % (star.xval, star.yval))
        coofile.close()
        iraf.apphot.phot(image,coords='%s.coo' % root,output='%s.grb' % root,verbose=no)
        [date,ut,exptime,filter]=get_head(image,['DATE-OBS','UT','EXPTIME','FILTERS'])
        day=int(date[:2])
        month=int(date[3:5])
        year=int(date[6:])
        hour=int(ut[:2])
        minute=int(ut[3:5])
        second=float(ut[6:])
        t=DateTime(year,month,day,hour,minute,second)+exptime/24.0/3600.0/2.0
        dt=(t-t0).seconds
        grb=Starlist('%s.grb' % root)
        if len(grb):
           mag=grb[0].mag
           smagu=grb[0].magu
        else:
           mag=99.999
           smagu=99.999
        stars=Starlist('%s.mag' % root)
        refstars.set_mag(filters[filter])
        zp,zpmagu=stars.zeropt(refstars,method='mean',rejout=1,sigma=1.5,fencelim=0.50)
        if not (mag==99.999):
           mag+=zp
           magu=sqrt(pow(smagu,2)+pow(zpmagu,2))
        else:
           magu=99.999
        outfile.write("%s\t%s\t%10.2f%10.2f\t%10.3f%10.3f%10.3f%10.3f\n" % (t,filter,dt,exptime,mag,magu,smagu,zpmagu))

###############################################################################

def kaitpsflc(infiles, refstars, ot, t0, filters, outfile):

    for image in infiles:
        refstars.wcs2pix(image)
        ot.wcs2pix(image)
        root,null=image.split('.')
        [date,ut,exptime,filter]=get_head(image,['DATE-OBS','UT','EXPTIME','FILTERS'])
        day=int(date[:2])
        month=int(date[3:5])
        year=int(date[6:])
        hour=int(ut[:2])
        minute=int(ut[3:5])
        second=float(ut[6:])
        t=DateTime(year,month,day,hour,minute,second)+exptime/24.0/3600.0/2.0
        dt=(t-t0).seconds
        alsfile=open('%s.als.1' % root)
        lines=alsfile.readlines()
        stars=[]
        i=0
        while i<len(lines):
            if lines[i].startswith('#'):
                i+=1
            else:
                temp1=lines[i].split()
                star=Star(name='PSF-%02i' % int(temp1[0]),
                          xval=float(temp1[1]), yval=float(temp1[2]),
                          mag=float(temp1[3]),magu=float(temp1[4]))
                stars.append(star)
                i+=2
        psfstars=Starlist(stars=stars)
        a,b=psfstars.match(ot)
        if len(a):
            mag=a[0].mag
            smagu=a[0].magu
        else:
            mag=99.999
            smagu=99.999
        refstars.set_mag(filters[filter])
        zp,zpmagu=psfstars.zeropt(refstars,method='mean',rejout=1,sigma=1.5,fencelim=0.50)
        if not (mag==99.999):
           mag+=zp
           magu=sqrt(pow(smagu,2)+pow(zpmagu,2))
        outfile.write("%s\t%s\t%10.2f%10.2f\t%10.3f%10.3f%10.3f%10.3f\n" % (t,filter,dt,exptime,mag,magu,smagu,zpmagu))

    return

###############################################################################

def gettmid(filelist):

   inlis=iraffiles('@'+filelist)
   im1=inlis[0]
   [date,ut,exptime]=get_head(im1,['DATE-OBS','UT','EXPTIME'])
   day=int(date[:2])
   month=int(date[3:5])
   year=int(date[6:])
   hour=int(ut[:2])
   minute=int(ut[3:5])
   second=float(ut[6:])
   t1=DateTime(year,month,day,hour,minute,second)+exptime/24.0/3600.0/2.0
   counter=0
   for i in range(1,len(inlis)):
      imi = inlis[i]
      [date2,ut2,exptime2] = get_head(imi, ['DATE-OBS','UT','EXPTIME'])
      day=int(date2[:2])
      month=int(date2[3:5])
      year=int(date2[6:])
      hour=int(ut2[:2])
      minute=int(ut2[3:5])
      second=float(ut2[6:])
      t2=DateTime(year,month,day,hour,minute,second)+exptime/24.0/3600.0/2.0
      counter+=(t2-t1).days
   midt=t1+(counter/len(inlis))
   return midt

############################################################################

def kait_disp():

    inlis=iraffiles('*_001.sub.fit')
    for image in inlis:
        root=image[:len(image)-12]
        iraf.display("%s_001.fit" % root, 1)
        if os.path.exists('%s.z.can' % root):
            iraf.tvmark(1,'%s.z.can' % root, radii=10)
        if os.path.exists('%s.mp.cor' % root):
            iraf.tvmark(1,'%s.mp.cor' % root, radii=10)
        iraf.display('%s.fit' % root,2)
        iraf.display('%s_001.sub.fit' % root, 3)
        if os.path.exists('%s.z.cans' % root):
            iraf.tvmark(3, '%s.z.cans' % root, radii=10)
        iraf.imexam()

###########################################################################

def kait_anet(imfile):
   
   inf = open(imfile)
   lines = inf.readlines()
   for line in lines:
      image = line.rstrip()
      #root='ccd.'+image.split('.')[1]+'.0'
      root = image.split('.')[0]
      anetcmd = "anet.py %s > %s.wcs" % (image, root)
      print anetcmd
      success = 0
      while not success:
         os.system(anetcmd)
         outf = open("%s.wcs" % root)
         outlines = outf.readlines()
         outf.close()
         if re.search("Solved", outlines[len(outlines)-2]):
            update_head(image, 'IQWCS', 1)
            success=1
         elif re.search("Failed", outlines[len(outlines)-2]) or \
              re.search("Running", outlines[len(outlines)-2]):
            update_head(image, 'IQWCS', 0)
            success=1
         else:
            success=0
