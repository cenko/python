#!/usr/bin/env python

######################################################################
# $Id: saofocus.py,v 1.1 2004/07/29 19:45:39 derekfox Exp derekfox $ #
######################################################################

version='1.1'

# Global packages
import pyraf
from pyraf import iraf
import copy, os, shutil, glob, sys, string, re, math
import getopt
import pyfits
import numarray
import time
from types import *

# Local packages
from iqutils import *
import iqpkg
import add_wcs

# Necessary IRAF packages
iraf.imred()
iraf.ccdred()
iraf.images()
iraf.imutil()

# Shortcuts
yes=iraf.yes
no=iraf.no
INDEF=iraf.INDEF
hedit=iraf.hedit
imgets=iraf.imgets

# Pyraf parameters, where to find them
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

######################################################################

def saofocus(inpat,prepfile=yes,endhow="never",endwhen="",
             clobber=globclob,verbose=globver):

    """ derive best focus from an SAOFOCUS run on P60 """

    # Defaults
    twait=30
    reduced={}
    bpmname="BPM.pl"
    saokey="IMGTYPE"
    saore="SAOFOCUS"
    fseqkey="SAOFVALS"
    sigma=2.0
    satval=50000.0
    masksfx="mask"
    pix=0.3787 # arcsec per pixel
    mindiff=2  # minimum tolerable diff b/ best and worst seeing (pix)

    # Parse end-condition "time"
    if endhow=="time":
        re1=re.search("^(\d+):(\d+)",endwhen)
        if re1:
            tnow=time.gmtime()
            # Force UT calculation w/ ignoring of DST issues
            reftime=time.mktime([tnow[0],tnow[1],tnow[2],
                                 int(re1.group(1)),int(re1.group(2)),0,
                                 tnow[6],tnow[7],0]) - time.timezone
            if reftime<time.time():
                reftime+=86400
            if verbose:
                print "Running until %s" % \
                      time.strftime("%d %b %H:%M",time.gmtime(reftime))
        else:
            print "Failed to parse %s as UT time" % endwhen
            print "Running until stopped..."
            endhow="never"

    ##################################################

    # Big Loop
    done=no
    while not done:

        # Parse inputs
        allfiles=glob.glob(inpat)

        newfiles=[]
        for image in allfiles:
            if not reduced.has_key(image):
                # Interested only in SAOFOCUS images
                saoval=get_head(image,saokey)
                if re.search(saore,saoval,re.I):
                    reduced[image]=no
                    newfiles.append(image)
                else:
                    # Skip
                    reduced[image]=yes

        for image in newfiles:

            if verbose:
                print "Processing new SAOFOCUS image %s" % image

            # Only attempt reduction once for each image
            reduced[image]=yes

            # Demosaic
            ccdproc=iraf.ccdproc
            ccdproc.overscan=yes
            ccdproc.trim=yes
            ccdproc.zerocor=no
            iraf.iqmosaic(image,outpfx="",biassec="!BIASSEC",
                          trimsec="!TRIMSEC",joinkeys="DETSEC",
                          splitkeys="AMPSEC,CCDSEC,TRIM,OVERSCAN",
                          clobber=yes,verbose=verbose)
            # Fix the WCS that IRAF has now clobbered
            iraf.add_wcs(image,instrument="p60new")

            # Get some important keywords
            [raoff,decoff,saonseq,saoshift] = \
                   get_head(image,['RAOFF','DECOFF','NUMFOC',
                                   'SAOSHIFT'])
            # Type conversions
            saonseq=int(saonseq)
            saoshift=float(saoshift)

            # No bias-subtraction or flat-fielding

            # Subtract stair-step background if requested
            if prepfile:
                qtmp=iraf.mktemp("saofoc")+".fits"
                iraf.imcopy(image,qtmp,verbose=no)
                for iamp in [-1,1]:
                    for istep in range(saonseq):
                        qsec=iraf.mktemp("saofsc")+".fits"
                        y1=1024+iamp*(1+istep*saoshift)+max(0,iamp)
                        y2=y1+iamp*(saoshift-1)
                        if istep==0:
                            y1-=iamp
                        elif istep==saonseq-1:
                            y2=1+max(0,iamp)*2047
                        ymin=min(y1,y2)
                        ymax=max(y1,y2)
                        statsec='[*,%d:%d]' % (ymin,ymax)
                        iraf.iterstat(image+statsec,nsigrej=5,
                                      maxiter=10,prin=no,verbose=no)
                        medreg=float(iraf.iterstat.median)
                        iraf.imarith(image+statsec,'-',medreg,
                                     qsec,verbose=no,noact=no)
                        iraf.imcopy(qsec,qtmp+statsec,verbose=no)
                        iraf.imdel(qsec,verify=no,go_ahead=yes)
                iraf.imdel(image,verify=no,go_ahead=yes)
                iraf.imrename(qtmp,image,verbose=no)
            
            # No renaming of the file

            # No Bad Pixel Masking or rudimentary WCS

            # Object-detection
            doobj=1
            if check_head(image,'STARFILE'):
                if os.path.exists(get_head(image,'STARFILE')):
                    doobj=0
            if doobj:
                iraf.iqobjs(image,sigma,satval,skyval="0.0",
                            masksfx=masksfx,wtimage="none",minlim=no,
                            clobber=yes,verbose=verbose)

            # Establish our "region of interest" for pixels
            pixrough=800
            pixminx=1024-float(decoff)/pix-pixrough/2
            pixmaxx=pixminx+pixrough
            pixmaxy=1023
            pixminy=pixmaxy-float(raoff)/pix-(1+saonseq)*saoshift-pixrough/2

            # The array of focus positions (from low Y to high Y)
            farrkeys=[]
            for i in range(saonseq):
                farrkeys.append("SAOFOC%02d" % (i+1))
            farr=get_head(image,farrkeys)

            # Starlist from Sextractor
            starfile=get_head(image,'STARFILE')
            stars=Starlist(starfile)
            stars.mag_sort()

            # Find some of the stars in the sequence
            xvals=[]
            yvals=[]
            good1=[]
            for star in stars:
                # We're going in order from bright to faint
                if star.xval>pixminx and star.xval<pixmaxx and \
                   star.yval<pixmaxy and star.yval>pixminy and \
                   star.elong<3:
                    good1.append(star)
                    xvals.append(star.xval)
                    yvals.append(star.yval)
                    if len(good1)>=1.2*saonseq:
                        break

            # Sanity checking
            if len(xvals)<1:
                print "No stars found in search region"
                update_head(image,'SAOFOCUS',0,
                            "saofocus processing failed for this image")
                continue

            # Guess for x is median of the xvals of the goodstars
            xctr=median(xvals,pick=1)

            # First guess for y is median of the goodstars that fall
            # within +/- (pixfine) pix of this X value
            pixfine=8
            yval2=[]
            good2=[]
            for star in good1:
                if abs(star.xval-xctr)<pixfine:
                    good2.append(star)
                    yval2.append(star.yval)

            # Sanity checking
            if len(yval2)<1:
                print "No stars found in refined search region"
                update_head(image,'SAOFOCUS',0,
                            "saofocus processing failed for this image")
                continue

            yctr=median(yval2,pick=1)

            # This will be our reference star
            usestar=good2[yval2.index(yctr)]
            [refx,refy]=[usestar.xval,usestar.yval]
            
            # Identify additional stars, working out from the
            # reference star
            usestars=[usestar]

            # Search increasing & decreasing in Y
            # (we use a bigger box in Y, and correct for the last star
            # position, because of potential for tracking errors)
            for isgn in [+1,-1]:
                shifty=refy
                for i in range(1,saonseq):
                    shifty += isgn*saoshift
                    tgtstars=stars.starsinbox([[refx-pixfine,refx+pixfine],
                                               [shifty-1.5*pixfine,
                                                shifty+1.5*pixfine]])
                    if len(tgtstars)>0:
                        minmag=99.99
                        # Pick the brightest in the box
                        for tstar in tgtstars:
                            if tstar.mag<minmag:
                                svstar=tstar
                                minmag=tstar.mag
                        usestars.append(svstar)
                        shifty=svstar.yval
                    else:
                        # Quit searching when there are no stars in box
                        break

            # Exclude faint stars if we have too many
            # (could also drop stars until we reach the number, but
            # then we would want to check that all the y-values are in
            # strict succession...)
            magrough=2.0
            if len(usestars)>saonseq:
                usemags=[]
                for star in usestars:
                    usemags.append(star.mag)
                minmag=min(usemags)
                use2=[]
                for star in usestars:
                    if star.mag<minmag+magrough:
                        use2.append(star)
                if len(use2)==saonseq:
                    usestars=use2

            # Check for the right number
            if len(usestars)!=saonseq:
                print "Failed to get right number of target stars"
                update_head(image,'SAOFOCUS',0,
                            "saofocus processing failed for this image")
                continue

            # Construct the array of x- and y-values
            xval3=[]
            yval3=[]
            for star in usestars:
                xval3.append(star.xval)
                yval3.append(star.yval)
            # Choose a single x-value for the run
            xctr2="%.3f" % median(xval3)
            # Increasing y-value means increasing index in farr
            yval3.sort()
            ylis=','.join(map(str,yval3))
            flis=','.join(map(str,farr))

            # Run iqfocus
            iraf.iqfocus(image,xctr2,ylis,flis,method="acorr",
                         boxsize=saoshift,focuskey="BESTFOC",
                         seekey="SEEING",fseekey="SAOSEE",
                         update=yes,clobber=yes,verbose=verbose)

            # Check for successful iqfocus run
            if not check_head(image,'BESTFOC'):
                print "iqfocus run failed to find good focus"
                update_head(image,'SAOFOCUS',0,
                            "saofocus processing failed for this image")
                continue

            # Grab focus positions and seeing values
            bestfoc,saosee=get_head(image,['BESTFOC','SAOSEE'])

            # Some information
            if verbose:
                print "Focus settings:  "+flis
                print "Seeing values:   "+saosee
                print "Best focus:      %.3f" % bestfoc

            # Check for best-focus "near edge" condition
            for i in range(saonseq):
                if ("%.3f" % bestfoc)==("%.3f" % float(farr[i])):
                    ibest=i
                    break
            if (saonseq>3 and (ibest==0 or ibest==saonseq-1)) or \
               (saonseq>5 and (ibest==1 or ibest==saonseq-2)):
                update_head(image,'SAOEDGE',1,
                            "Best seeing is near edge of focus run")

            # Check for insufficient variation across the focus run
            maxsee=0
            minsee=0.5*saoshift
            saoseevals=eval(saosee)
            for seeval in saoseevals:
                minsee=min(minsee,seeval)
                maxsee=max(maxsee,seeval)
            if (maxsee-minsee) < mindiff:
                print "Failed to observe significant variation "+\
                      "across focus run"
                update_head(image,'SAOFOCUS',0,
                            "saofocus processing failed for this image")
                continue

            # Convert seeing value to arcsec
            seepix=get_head(image,'SEEPIX')
            update_head(image,'BESTSEE',seepix,
                        "Best seeing in pixels for this series")
            try:
                seeasec=float(seepix)*pix
                update_head(image,'SEEING',seeasec,
                            "Best seeing in arcsec for this series")
            except:
                print "Failed to convert seeing to arcsec for %s" % image

            # Successful processing
            update_head(image,'SAOFOCUS',1,
                        "saofocus processing succeeded for this image")

        ##############################

        # Test end conditions
        if endhow=="never":
            done=no
        elif endhow=="once":
            done=yes
        elif endhow=="time":
            if time.time()>reftime:
                done=yes
        
        # Wait a little while
        if not done:
            time.sleep(twait)

######################################################################

# location of parameter files
_parfile=pardir + "saofocus.par"
t=iraf.IrafTaskFactory(taskname="saofocus",
                       value=_parfile,function=saofocus)

######################################################################
######################################################################

def main():

    # Command line
    try:
        opts, args = getopt.getopt(sys.argv[1:], 
                                   "honw:",
                                   ["help","once","noprep","when="])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    endhow="never"
    endwhen=""
    focmode=""
    prepfile=1

    for opt, val in opts:
        # Usage info only
        if opt in ("-h", "--help"):
            usage()
            sys.exit(1)
        # End of run specification: Specified time
        elif opt in ("-w", "--when"):
            endhow="time"
            endwhen=val
        # End of run specification: Run once
        elif opt in ("-o", "--once"):
            endhow="once"
        # Prep file by subtracting stairstep background?
        elif opt in ("-n", "--noprep"):
            prepfile=0
        else:
            sys.stderr.write("Unmatched option %s\n" % opt)

    inpat=args[0]

    print "Running saofocus with inpat=%s" % inpat
    saofocus(inpat,prepfile=prepfile,endhow=endhow,endwhen=endwhen)
    sys.exit(0)

######################################################################

def usage():

    (xdir,xname)=os.path.split(sys.argv[0])
    print "Usage:  %s [-o | -w HH:MM] [-n] \"filepattern\"" % xname
    print
    print "   -o          single processing loop"
    print "   -w HH:MM    run until UT time HH:MM"
    print "   -n          'noprep' run without stairstep back-sub"

##############################

# Running as executable
if __name__=='__main__':

    main()
