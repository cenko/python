
#############################################################
# $Id$
#############################################################

import pyraf
from pyraf import iraf
import copy, os, shutil, glob, sys, string, re, math, operator
import pyfits
from types import *
#from mx.DateTime import *
from subprocess import PIPE, Popen
import numpy as np

# Access to the iqutils
from iqutils import *
import iqpkg

# Necessary packages
iraf.images()
iraf.immatch()
iraf.noao()
iraf.imred()
iraf.ccdred()
iraf.digiphot()
iraf.daophot()
iraf.images()
iraf.imcoords()
iraf.stsdas()
iraf.hst_calib()
iraf.nicmos()

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

PSFEXDIR = "/Users/scenko/python/ppqlc"

######################################################################

def psfphot(inlist, ra, dec, reffilt, interact, fwhm, readnoise, gain, 
            threshold,refimage=None,starfile=None,maxnpsf=5, 
            clobber=globclob,verbose=globver,skykey='SKYBKG',
            filtkey='FILTER',pixtol=3.0):


    """ perform PSF-based photometry on a single target star (SN?) at RA, Dec and  
        also on a set of comparison stars, using daophot.  simultaneously 
        perform aperture photometry on all the comparison stars (after 
        subtracting off contributions from neighbors) to enable absolute 
        photometry by comparison to aperture photometry of standard stars 
        observed in other fields """

    # Defaults / constants
    psfmult=5.0         #standard factor (multiplied by fwhm to get psfradius)
    psfmultsmall=3.0    #similar to psfmult, adjusted for nstar and substar

    # Necessary package
    iraf.imutil()

    # Parse inputs
    infiles=iraffiles(inlist)

    # Which file is reffilt?  call it refimage
    if refimage==None:
        for image in infiles:
            if check_head(image, filtkey):
                try:
                    imgfilt = get_head(image, filtkey)
                    if imgfilt == reffilt:
                        refimage = image
                        break
                except:
                    pass
            
    if not refimage:
        print "BAD USER!  No image corresponds to the filter: %s" % reffilt
        return
    else:
        refroot='s'+refimage.split('.')[0]

    #first make sure to add back in background of sky
    iraf.iqsubsky(inlist, sub=no, skykey=skykey)

    #put reference image first on list
    infiles.remove(refimage)
    infiles.insert(0,refimage)

    #setup for keywords
    if gain == "!GAIN":
        try: gainval = float(get_head(image, gain))
        except:
            print "Bad header keyword for gain."
    else:
        gainval = float(gain)

    if readnoise == "!READNOISE":
        try: readval = float(get_head(image, readnoise))
        except:
            print "Bad header keyword for readnoise."
    else:
        readval = float(readnoise)

    # Process each file in turn
    for image in infiles:

        # Check that the image is there
        check_exist(image,"r")

        # Grab image root name
        root=image.split('.')[0]

        # Map image to reference image
        if not (image==refimage):
            [nx,ny]=get_head(image,['NAXIS1','NAXIS2'])
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
            for i in range(len(match)):
                matchfile.write('%10.3f%10.3f%10.3f%10.3f\n' % 
                               (refmatch[i].xval,refmatch[i].yval,
                                match[i].xval,match[i].yval))
            matchfile.close()
            check_exist('%s.geodb' % root, 'w', clobber=clobber)
            iraf.geomap('%s.match' % root,'%s.geodb' % root,1.0,nx,1.0,ny,
                        verbose=no,interactive=no)
            check_exist('s%s.fits' % root, 'w', clobber=clobber)
            iraf.geotran(image,'s%s' % root,'%s.geodb' % root,
                         '%s.match' % root,geometry="geometric",
                         boundary="constant",verbose=no)
        else:
            iraf.imcopy(image,'s%s' % root)
        root='s%s' % root
 
        #get sky level and calculate sigma
        #if check_head(image, skykey):
        #    try:
        #        sky=float(get_head(image, skykey))
        #    except:
        #        print "No sky levels in header."

        #sigma= (((sky * gainval) + readval**2)**.5) / gainval        
        iraf.iterstat(image)
        
        # Saturation level
        if not check_head(image, "SATURATE"):
        	saturate = 60000.0
        else:
        	saturate = get_head(image, "SATURATE")
        	        
        # Update datapars and daopars
        iraf.datapars.fwhmpsf=fwhm
        iraf.datapars.sigma=iraf.iterstat.sigma
        iraf.datapars.datamin=iraf.iterstat.median-10*iraf.iterstat.sigma
        iraf.datapars.datamax=0.90*saturate
        iraf.datapars.readnoise=readval
        iraf.datapars.epadu=gainval
        iraf.datapars.filter=filtkey
        iraf.daopars.psfrad=psfmult*fwhm
        iraf.daopars.fitrad=fwhm
        iraf.daopars.function="gauss,moffat15,moffat25,lorentz,penny1"

        #find stars in image unless a starlist is given
        if image==refimage and starfile==None:
            iraf.daophot.daofind(root,'refimage.coo.1',threshold=threshold,verify=no,
                         verbose=verbose)
        elif image==refimage:
            shutil.copy(starfile,'refimage.coo.1')

        #initial photometry
        iraf.daophot.phot(root,'refimage.coo.1','default',aperture=fwhm,verify=no,
                  verbose=verbose)

        #select stars for psf the first time
        refstarsfile = "refimage.pst.1"
        if image == refimage:
            iraf.pstselect(root,'default',refstarsfile,maxnpsf,
                           interactive=yes,verify=no,verbose=verbose)

        #fit the psf
        iraf.psf(root,'default',refstarsfile,'default','default','default',
                 interactive=interact,verify=no,verbose=verbose)

        #identify neighboring/interfering stars to selected stars
        groupingfile = root+".psg.1"
        iraf.nstar(root,groupingfile,'default','default','default',
                   psfrad= psfmultsmall * fwhm,verify=no,verbose=verbose)

        #subtract out neighboring stars from image
        iraf.substar(root,'default',refstarsfile,'default','default',
                     psfrad=psfmultsmall*fwhm,verify=no,verbose=verbose)

        #repeat psf to get better psf model
        #IRAF's interactive version usually crashes
        subtractedimage = root+".sub.1"
        iraf.psf(subtractedimage,root+".nst.1",refstarsfile,'%s.psf.2' % root,
                 '%s.pst.2' % root,'%s.psg.2' % root,interactive=interact,
                 verify=no,verbose=verbose)

        #Need to make sure SN was detected by daofind
        stars=Starlist('%s.mag.1' % root)
        SN=Star(name='SN',radeg=ra,dcdeg=dec,fwhm=2.0,fwhmw=2.0)
        SNlis=Starlist(stars=[SN])
        SNlis.wcs2pix(image)
        if (len(stars.match(SNlis)[0])==0):
            #No match - need to add to daofind file
            print "No match!"
            coofile=open('refimage.coo.1', 'a+')
            coofile.write('%10.3f%10.3f%9.3f%8.3f%13.3f%12.3f%8i\n' % (SNlis[0].xval, SNlis[0].yval,99.999,0.500,0.000,0.000,999))
            coofile.close()    

        #repeat aperture photometry to get good comparisons to standard fields
        iraf.daophot.phot(root,'refimage.coo.1','default',aperture=psfmult*fwhm,
                  verify=no,verbose=verbose)

        # allstar run
        iraf.allstar(root,'default','default','default','default','default',
                     verify=no,verbose=verbose)

######################################################################

def als2reg(image, number):

    """Convert an IRAF .als file to an iqutils Starlist""" 

    root = image[:-5]
    inf = open("%s.als.%i" % (root, number))
    lines = inf.readlines()
    
    stars = []
    i=0
    while i < len(lines):
        
        line = lines[i]
        if re.match("#", line):
            i+=1
        else:
            temp = line.strip().split()
            newstar = Star(name="IRAF-%i" % int(temp[0]), xval=float(temp[1]),
                           yval=float(temp[2]), mag=float(temp[3]),
                           magu=float(temp[4]))
            stars.append(newstar)
            i+=2

    newstars=Starlist(stars=stars)
    return newstars

##########################################################################

def make_ascii(inlist, ot, refstars, outfile, filter, pixtol=3.0):

    inlis = iraffiles(inlist)
    outf = open(outfile, "w")
    refstars.set_mag("%sMAG" % filter)
    
    for image in inlis:
        stars = als2reg("%s" % image, 1)
        refstars.wcs2pix("%s" % image)
        zp, zpu = stars.zeropt(refstars, method="mean", rejout=0,tol=pixtol)
        c,d = stars.match(refstars,tol=pixtol)
        ot.wcs2pix("%s" % image)
        a, b = stars.match(ot, tol=pixtol)
        zpu = zpu / sqrt(len(c))
        if (len(a) > 0):
            mag = a[0].mag + zp
            dmag = a[0].magu 
        else:
            mag = 99.99
            dmag = 99.99
        utshut,exptime = get_head(image, ["UTSHUT", "EXPTIME"])
        t1 = DateTime(int(utshut[:4]), int(utshut[5:7]), int(utshut[8:10]), 
                       int(utshut[11:13]), int(utshut[14:16]), 
                       float(utshut[17:]))
        #t1 = DateTime(2012,1,1,0,0,0)
        #utshut = "Null"
        #exptime = 1.0
        outf.write("%s\t%10.5f%10.2f%10.3f%10.3f%10.3f\n" % (utshut, t1.mjd, exptime, mag, dmag, zpu))

    outf.close()

##########################################################################

def psfex_phot(inlis, stamps, filt, ot, outfile):

    '''PSF Photometry using the PSFEx package'''
	
    stamps.set_mag("%sMAG" % filt)
    outf = open(outfile, "w")

    for image in inlis:

        root = image.rstrip(".fits")

        # Convert stamps file to pixel coordinates
        stamps.wcs2pix(image)
        outf2 = open("%s.stamps" % root, "w")
        for star in stamps:
            outf2.write("%10.3f%10.3f\n" % (star.xval, star.yval))
        outf2.close()

        # Run SExtractor to get vignettes
        scmd = "sex -c %s/prepsfex.sex " % PSFEXDIR + \
               "-CATALOG_NAME %s.pcat " % root + \
               "-ASSOC_NAME %s.stamps %s" % (root, image)
        junk = Popen(scmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        # Create PSF from vignettes
        pcmd = "psfex -c %s/psfex.conf %s.pcat" % (PSFEXDIR, root)
        junk = Popen(pcmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        # Now do PSF-photometry
        scmd = "sex -c %s/psfex.sex -PSF_NMAX 1 " % PSFEXDIR + \
               "-PSF_NAME %s.psf " % root + \
               "-CATALOG_NAME %s.cat %s" % (root, image)
        junk = Popen(scmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        # Read in resulting photometry file
        pfile = np.loadtxt("%s.cat" % root)
        stars = []
        for line in pfile:
            star = Star(name="SEX_%03i" % line[0], xval=line[1], yval=line[2], 
                        mag=line[5], magu=line[6])
            stars.append(star)

        # Zeropoint
        stars = Starlist(stars=stars)
        zp, zpu = stars.zeropt(stamps,method="mean",rejout=0)

        # OT 
        ot.wcs2pix(image)
        a, b = stars.match(ot)
        if len(a)==1:
            otmag = a[0].mag+zp; otmagu = a[0].magu
        else:
            otmag = otmagu = 99.0

        [mjd, exptime] = get_head(image, ["OBSMJD", "EXPTIME"])
        outf.write("%15.7f%8.2f%10.3f%10.3f%10.3f\n" % (float(mjd), exptime, 
                                                        otmag, otmagu, zpu))

    outf.close()
              
##########################################################################
	

# location of parameter files

_parfile=pardir + "psfphot.par"
t=iraf.IrafTaskFactory(taskname="psfphot",
                       value=_parfile,function=psfphot)

_parfile=pardir + "iqcals.par"

