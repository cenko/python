
# Global python packages
import pyraf
from pyraf import iraf
import os, shutil, math
import pyfits
from types import *
import time

# Local python packages
from iqutils import *

# IRAF modules
iraf.images()
iraf.noao()
iraf.imred()
iraf.ccdred()
iraf.specred()

# IRAF variables
yes=iraf.yes
no=iraf.no
INDEF=iraf.INDEF
hedit=iraf.hedit
imgets=iraf.imgets
imcombine=iraf.imcombine

# PyRAF setup
pyrafdir="python/pyraf/"
pyrafdir_key='PYRAFPARS'

if os.environ.has_key(pyrafdir_key):
    pardir=os.environ[pyrafdir_key]
else:
    pardir=os.environ['HOME']+'/'+pyrafdir

if pardir[-1] != '/':
    pardir += '/'

# Global variables
REDBIAS='[1201:1231,*]'
REDTRIM='[1:1200,41:200]'
REDGAIN=3.0
REDRDNOISE=12.5
BLUEBIAS1='[2052:2080,*]'
BLUEBIAS2='[2082:2110,*]'
BLUETRIM1='[1:1024,31:300]'
BLUETRIM2='[1025:2048,31:300]'
BLUEGAIN=1.2
BLUERDNOISE=3.7


############################################################################

def redbias(images, biassec=REDBIAS, trimsec=REDTRIM):

    '''Subtract overscan and trim red frames'''

    for image in images:
        iraf.ccdproc(image, ccdtype='', noproc=no, fixpix=no, overscan=yes, 
                     trim=yes, zerocor=no, darkcor=no, flatcor=no, illumcor=no,
                     fringecor=no, readcor=no, scancor=no, biassec=biassec, 
                     trimsec=trimsec)
        update_head(image, 'DISPAXIS', 1)

    return

#############################################################################

def bluebias(images, biassec1=BLUEBIAS1, trimsec1=BLUETRIM1, 
             biassec2=BLUEBIAS2, trimsec2=BLUETRIM2):

    '''Subtract overscan and trim blue frames'''
     
    for image in images:

        root,ext=image.split('.')
        iraf.ccdproc(image, output='%s_1' % root, ccdtype='', noproc=no, 
                     fixpix=no, overscan=yes, trim=yes, zerocor=no, 
                     darkcor=no, flatcor=no, illumcor=no, fringecor=no, 
                     readcor=no, scancor=no, biassec=biassec1, 
                     trimsec=trimsec1)
        iraf.ccdproc(image, output='%s_2' % root, ccdtype='', noproc=no, 
                     fixpix=no, overscan=yes, trim=yes, zerocor=no, 
                     darkcor=no, flatcor=no, illumcor=no, fringecor=no, 
                     readcor=no, scancor=no, biassec=biassec2, 
                     trimsec=trimsec2)
        iraf.imjoin('%s_1,%s_2' % (root, root), 'j%s' % image, 1)
        update_head('j%s' % image, 'DISPAXIS', 1)
        update_head('j%s' % image, ['CCDSEC', 'DATASEC'], ['[1:2048,1:270]', 
                    '[1:2048,1:270]'])

    return

############################################################################

def make_flat(images, outflat, gain=1.0, rdnoise=0.0, xwindow=50,
              ywindow=50, hmin=0, hmax=65535, lowclip=0.7, highclip=1.3):

    '''Construct flat field from individual frames'''
    
    flatimages=','.join(images)
    iraf.flatcombine(flatimages, output='flat1', combine='median', 
                     reject='avsigclip', ccdtype='', process=no, subsets=no,
                     delete=no, clobber=no, scale='median', lsigma=3.0,
                     hsigma=3.0, gain=gain, rdnoise=rdnoise)
    iraf.fmedian('flat1', 'flat2', xwindow, ywindow, hmin=hmin, hmax=hmax)
    iraf.imarith('flat1', '/',  'flat2', outflat)
    iraf.imreplace(outflat, 1.0, lower=INDEF, upper=lowclip)
    iraf.imreplace(outflat, 1.0, lower=highclip, upper=INDEF)

    return

#############################################################################

def reference_arc(image, output, reference, 
                  coordlist='home$linelists/licklinelist.dat'):

    '''Construct reference wavelength solution from arc lamps'''

    iraf.apall(image, output=output, references=reference, interactive=no,
               find=no, recenter=no, resize=no, edit=no, trace=no, fittrace=no,
               extract=yes, extras=yes, review=no, background='none')
    iraf.identify(output, coordlist=coordlist)

    return

############################################################################

def red_standard(image, arcs, flats, object=None, biassec=REDBIAS, 
                 trimsec=REDTRIM, outflat='Flat-Red.fits', gain=REDGAIN,
                 rdnoise=REDRDNOISE, arc='Arc-Red.fits', 
                 caldir='home$standards/'):

    '''Reduce and calibrate standard star observation with red CCD'''

    
    # Bias subtract everything first
    redbias(image, biassec=biassec, trimsec=trimsec)
    redbias(arcs, biassec=biassec, trimsec=trimsec)
    redbias(flats, biassec=biassec, trimsec=trimsec)

    # Create and apply flat-field
    make_flat(flats, outflat, gain=gain, rdnoise=rdnoise)
    iraf.ccdproc(image[0], ccdtype='', noproc=no, fixpix=no, overscan=no, 
                 trim=no, zerocor=no, darkcor=no, flatcor=yes, illumcor=no,
                 fringecor=no, readcor=no, scancor=no, flat=outflat)
    arcimages=','.join(arcs)
    iraf.ccdproc(arcs, ccdtype='', noproc=no, fixpix=no, overscan=no, 
                 trim=no, zerocor=no, darkcor=no, flatcor=yes, illumcor=no,
                 fringecor=no, readcor=no, scancor=no, flat=outflat)

    # Extract spectrum of standard
    if object==None:
        object=get_head(image, 'OBJECT')
    iraf.apall(image[0], output=object, references='', interactive=yes, 
               find=yes, recenter=yes, resize=yes, edit=yes, trace=yes, 
               fittrace=yes, extract=yes, extras=yes, review=no, 
               background='fit', weights='variance', pfit='fit1d', 
               readnoise=rdnoise, gain=gain)

    # Extract arc and fit wavelength solution
    iraf.imarith(arcs[0], '+', arcs[1], 'Arc-Sum.fits')
    reference_arc('Arc-Sum.fits', arc, image[0])

    # Apply wavelength solution to standard
    iraf.refspec(object, references=arc, sort="", group="", override=yes,
                 confirm=no, assign=yes)
    iraf.dispcor(object, '%s.w' % object, confirm=no, listonly=no)

    # Remove absorption features and smooth
    iraf.splot('%s.w' % object, 1, 1)
    iraf.gauss('temp1[*,1,1]', '%s.smooth' % object, 3.0)
    iraf.sarith('%s.w' % object, '/', '%s.smooth' % object, '%s.s' % object)

    # Create and apply telluric correction
    iraf.sarith('%s.w' % object, '/', 'temp1', 'telluric')
    iraf.splot('telluric', 1, 1)
    iraf.telluric('%s.s' % object, '%s.t' % object, 'telluric', xcorr=no,
                  tweakrms=no, interactive=no, sample='6850:6950,7575:7700')

    # Define bandpasses for standard star calculation
    iraf.standard('%s.t' % object, '%s.std' % object, 
                  extinction='onedstds$kpnoextinct.dat', caldir=caldir,
                  observatory='Lick', interact=yes, star_name=object, 
                  airmass='', exptime='')

    # Determine sensitivity function
    iraf.sensfunc('%s.std' % object, '%s.sens' % object, 
                  extinction='onedstds$kpnoextinct.dat', 
                  newextinction='extinct.dat', observatory='Lick',
                  function='legendre', order=3, interactive=yes)

    return

############################################################################

def blue_standard(image, arcs, flats, object=None, biassec1=BLUEBIAS1,
                  trimsec1=BLUETRIM1, biassec2=BLUEBIAS2, 
                  trimsec2=BLUETRIM2, outflat='Flat-Blue.fits', gain=BLUEGAIN,
                  rdnoise=BLUERDNOISE, arc='Arc-Blue.fits', 
                  caldir='home$standards/'):

    '''Reduce and calibrate standard star observation with blue CCD'''
  
    # Bias subtract everything first
    bluebias(image, biassec1=biassec1, trimsec1=trimsec1, biassec2=biassec2,
             trimsec2=trimsec2)
    bluebias(arcs, biassec1=biassec1, trimsec1=trimsec1, biassec2=biassec2,
             trimsec2=trimsec2)
    bluebias(flats, biassec1=biassec1, trimsec1=trimsec1, biassec2=biassec2,
             trimsec2=trimsec2)

    # Create and apply flat-field
    for i in range(len(flats)):
        flats[i]='j%s' % flats[i]
    make_flat(flats, outflat, gain=gain, rdnoise=rdnoise)
    iraf.ccdproc('j%s' % image[0], ccdtype='', noproc=no, fixpix=no, 
                 overscan=no, trim=no, zerocor=no, darkcor=no, flatcor=yes, 
                 illumcor=no, fringecor=no, readcor=no, scancor=no, 
                 flat=outflat)
    iraf.ccdproc('j%s' % arcs[0], ccdtype='', noproc=no, fixpix=no, 
                 overscan=no, trim=no, zerocor=no, darkcor=no, flatcor=yes, 
                 illumcor=no, fringecor=no, readcor=no, scancor=no, 
                 flat=outflat)

    # Extract spectrum of standard
    if object==None:
        object=get_head(image, 'OBJECT')
    iraf.apall('j%s' % image[0], output=object, references='', interactive=yes, 
               find=yes, recenter=yes, resize=yes, edit=yes, trace=yes, 
               fittrace=yes, extract=yes, extras=yes, review=no, 
               background='fit', weights='variance', pfit='fit1d', 
               readnoise=rdnoise, gain=gain)

    # Extract arc and fit wavelength solution
    reference_arc('j%s' % arcs[0], arc, 'j%s' % image[0])

    # Apply wavelength solution to standard
    iraf.refspec(object, references=arc, sort="", group="", override=yes,
                 confirm=no, assign=yes)
    iraf.dispcor(object, '%s.w' % object, confirm=no, listonly=no)

    # Remove absorption features and smooth
    iraf.splot('%s.w' % object, 1, 1)
    iraf.gauss('temp1[*,1,1]', '%s.smooth' % object, 3.0)
    iraf.sarith('%s.w' % object, '/', '%s.smooth' % object, '%s.s' % object)

    # Define bandpasses for standard star calculation
    iraf.standard('%s.s' % object, '%s.std' % object, 
                  extinction='onedstds$kpnoextinct.dat', caldir=caldir,
                  observatory='Lick', interact=yes, star_name=object, 
                  airmass='', exptime='')

    # Determine sensitivity function
    iraf.sensfunc('%s.std' % object, '%s.sens' % object, 
                  extinction='onedstds$kpnoextinct.dat', 
                  newextinction='extinct.dat', observatory='Lick',
                  function='legendre', order=3, interactive=yes)

    return

#############################################################################

def red_science(image, flats, spath, object=None, arc='Arc-Red.fits', 
                smooth='BD284211.smooth.fits', telluric='telluric.fits',
                sens='BD284211.sens', biassec=REDBIAS, trimsec=REDTRIM, 
                outflat='Flat-Red.fits', gain=REDGAIN, rdnoise=REDRDNOISE):

    '''Full reduction of KAST science spectra on red CCD'''

    # Bias subtract everything first
    redbias(image, biassec=biassec, trimsec=trimsec)
    redbias(flats, biassec=biassec, trimsec=trimsec)

    # Create and apply flat-field
    make_flat(flats, outflat, gain=gain, rdnoise=rdnoise)
    iraf.ccdproc(image[0], ccdtype='', noproc=no, fixpix=no, overscan=no, 
                 trim=no, zerocor=no, darkcor=no, flatcor=yes, illumcor=no,
                 fringecor=no, readcor=no, scancor=no, flat=outflat)

    # Cosmic ray rejection
    #iraf.lacos_spec(image[0], 'c%s' % image[0], 'cm%s' % image[0],
    #                gain=gain, readn=rdnoise, xorder=9, yorder=3,
    #                sigclip=4.5, sigfrac=0.5, objlim=1.0, niter=3)

    # Extract spectrum
    if object==None:
        object=get_head(image, 'OBJECT')
    iraf.apall('%s' % image[0], output=object, references='', interactive=yes, 
               find=yes, recenter=yes, resize=yes, edit=yes, trace=yes, 
               fittrace=yes, extract=yes, extras=yes, review=no, 
               background='fit', weights='variance', pfit='fit1d', 
               readnoise=rdnoise, gain=gain)

    # Apply wavelength solution to standard
    shutil.copy('%s/%s' % (spath, arc), '.')
    shutil.copy('%s/database/id%s' % (spath, arc.rstrip('.fits')), 'database')
    iraf.refspec(object, references=arc, sort="", group="", override=yes,
                 confirm=no, assign=yes)
    iraf.dispcor(object, '%s.w' % object, confirm=no, listonly=no)

    # Smooth
    shutil.copy('%s/%s' % (spath, smooth), '.')
    iraf.sarith('%s.w' % object, '/', smooth, '%s.s' % object)

    # Create and apply telluric correction
    shutil.copy('%s/%s' % (spath, telluric), '.')
    iraf.telluric('%s.s' % object, '%s.t' % object, telluric, xcorr=yes,
                  tweakrms=yes, interactive=yes, sample='6850:6950,7575:7700')

    # Flux calibration
    shutil.copy('%s/%s.0001.fits' % (spath, sens), '.')
    iraf.calibrate('%s.t' % object, '%s.f' % object, extinct=yes, flux=yes,
                   extinction='onedstds$kpnoextinct.dat', 
                   observatory='Lick', sensitivity=sens, airmass='',
                   exptime='')

    return

############################################################################

def blue_science(image, spath, object=None, flat='Flat-Blue.fits', 
                 arc='Arc-Blue.fits', smooth='BD284211.smooth.fits',
                 sens='BD284211.sens', biassec1=BLUEBIAS1, trimsec1=BLUETRIM1,
                 biassec2=BLUEBIAS2, trimsec2=BLUETRIM2, gain=BLUEGAIN, 
                 rdnoise=BLUERDNOISE):

    '''Full reduction of KAST science spectra on blue CCD'''

    # Bias subtract everything first
    bluebias(image, biassec1=biassec1, trimsec1=trimsec1, biassec2=biassec2,
             trimsec2=trimsec2)

    # Apply flat-field
    shutil.copy('%s/%s' % (spath, flat), '.')
    iraf.ccdproc('j%s' % image[0], ccdtype='', noproc=no, fixpix=no, 
                 overscan=no, trim=no, zerocor=no, darkcor=no, flatcor=yes, 
                 illumcor=no, fringecor=no, readcor=no, scancor=no, 
                 flat=flat)

    # Cosmic ray rejection
   #iraf.lacos_spec('j%s' % image[0], 'cj%s' % image[0], 'cmj%s' % image[0],
   #                gain=gain, readn=rdnoise, xorder=9, yorder=3,
   #                sigclip=4.5, sigfrac=0.5, objlim=1.0, niter=3)

    # Extract spectrum
    if object==None:
        object=get_head(image, 'OBJECT')
    iraf.apall('j%s' % image[0], output=object, references='', 
               interactive=yes, find=yes, recenter=yes, resize=yes, edit=yes, 
               trace=yes, fittrace=yes, extract=yes, extras=yes, review=no, 
               background='fit', weights='variance', pfit='fit1d', 
               readnoise=rdnoise, gain=gain)

    # Apply wavelength solution to standard
    shutil.copy('%s/%s' % (spath, arc), '.')
    shutil.copy('%s/database/id%s' % (spath, arc.rstrip('.fits')), 'database')
    iraf.refspec(object, references=arc, sort="", group="", override=yes,
                 confirm=no, assign=yes)
    iraf.dispcor(object, '%s.w' % object, confirm=no, listonly=no)

    # Smooth
    shutil.copy('%s/%s' % (spath, smooth), '.')
    iraf.sarith('%s.w' % object, '/', smooth, '%s.s' % object)

    # Flux calibration
    shutil.copy('%s/%s.0001.fits' % (spath, sens), '.')
    iraf.calibrate('%s.s' % object, '%s.f' % object, extinct=yes, flux=yes,
                   extinction='onedstds$kpnoextinct.dat', 
                   observatory='Lick', sensitivity=sens, airmass='',
                   exptime='')

    return

###########################################################################

def spec_combine(red, blue, cname, matchw1=5400, matchw2=5500, outw1=3500, 
                 outw2=10000):

    '''Combines blue and red spectra into single composite'''

    iraf.scombine("%s,%s" % (red, blue), cname, combine="average",
                  reject="avsigclip", w1=outw1, w2=outw2, dw=INDEF,
                  nw=INDEF, scale="median", zero="none", weight="none",
                  sample="%i:%i" % (matchw1, matchw2))

    #iraf.sarith(red, '/', blue, 'temp2', w1=matchw1, w2=matchw2)
    #iraf.iterstat('temp2[*,1,1]', nsigrej=5, maxiter=3, prin=no, verbose=no)
    #mean=float(iraf.iterstat.mean)
    #iraf.sarith(red, '/', mean, 'temp3')
    #iraf.scopy('temp3', 'red', w1=matchw1, w2=outw2)
    #iraf.scopy(blue, 'blue', w1=outw1, w2=matchw2)
    #iraf.scombine('red,blue', cname, w1=outw1, w2=outw2)

    return

###########################################################################

def do_all_science(object, blueimage, bluepath, bluestd, redimage, redflats, redpath,
                   redstd, btrimsec1=BLUETRIM1, btrimsec2=BLUETRIM2, rtrimsec=REDTRIM):

    if not os.path.exists("combine"):
        os.mkdir("combine")

    pyraf.iraffunctions.chdir('blue')
    blue_science(blueimage, '../%s' % bluepath, object=object, smooth="%s.smooth.fits" %
                 bluestd, sens="%s.sens" % bluestd, trimsec1=btrimsec1, 
                 trimsec2=btrimsec2)
    iraf.scopy("%s.f" % object, "../combine/blue", w1=3500, w2=5500)

    pyraf.iraffunctions.chdir("../red")
    red_science(redimage, redflats, "../%s" % redpath, object=object, 
                smooth="%s.smooth.fits" % redstd, sens="%s.sens" % redstd,
                trimsec=rtrimsec)
    iraf.scopy("%s.f" % object, "../combine/red", w1=5400, w2=10000)

    pyraf.iraffunctions.chdir("../combine")
    spec_combine("red", "blue", object)

    pyraf.iraffunctions.chdir("../")
    return





