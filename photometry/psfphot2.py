
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

def psfphot(image, clobber=globclob, verbose=globver, pixtol=3.0,
            maxnpsf=5, interact=yes):

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

    # Detect stars
    iqpkg.iqobjs(image, 3.0, 50000.0, wtimage="", skyval="!MEDSKY")

    root = image[:-5]
    [gain, rnoise, fwhm] = get_head(image, ["GAIN", "READNOI", "SEEPIX"])
    fwhm = float(fwhm); rnoise = float(rnoise)

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
    iraf.datapars.datamax=70000.0
    iraf.datapars.readnoise=rnoise
    iraf.datapars.epadu=gain   
    iraf.daopars.psfrad=psfmult*fwhm
    iraf.daopars.fitrad=fwhm
    iraf.daopars.function="gauss,moffat15,moffat25,lorentz,penny1"

    # coo file
    stars = Starlist("%s.stars" % image)
    outf = open("%s.coo.1" % image[:-5], "w")
    for star in stars:
        outf.write("%10.3f%10.3f\n" % (star.xval, star.yval))
    outf.close()

    #initial photometry
    iraf.daophot.phot(root,'default','default',aperture=fwhm,verify=no,
                      verbose=verbose)

    iraf.datapars.datamax=30000.0
    iraf.pstselect(root,'default','default',maxnpsf,interactive=yes,
                   verify=no,verbose=verbose)

    iraf.psf(root,'default','default','default','default','default',
             interactive=interact,verify=no,verbose=verbose)

    iraf.allstar(root,'default','default','default','default','default',
                 verify=no,verbose=verbose)

    iraf.iterstat("%s.sub.fits" % root)

    iraf.datapars.sigma=iraf.iterstat.sigma
    iraf.datapars.datamin=iraf.iterstat.median-10*iraf.iterstat.sigma

    iraf.datapars.datamax=70000.0
    iraf.daophot.phot("%s.sub.fits" % root, "SN.coo", 'default', 'default',
                      aperture=fwhm, verify=no, verbose=verbose)

    iraf.datapars.datamax=30000.0
    iraf.daopars.fitrad=fwhm*2.0
    iraf.allstar("%s.sub.fits" % root, 'default', "%s.psf.1.fits" % root, 
                 'default', 'default', 'default', verify=no, verbose=no)

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

