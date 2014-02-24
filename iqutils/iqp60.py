#!/usr/bin/env python

######################################################################
# $Id: iqp60.py,v 2.4 2009/06/15 05:38:23 p60 Exp $ #
######################################################################

version='1.10'

# Global packages
import pyraf
from pyraf import iraf
import copy, os, shutil, glob, sys, string, re, math
import getopt
import pyfits
import time
from types import *

# Local packages
import add_wcs
import iqpkg 
from iqutils import *

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

def p60status(statfile,image,scistat):

    if os.path.exists(statfile):
        # Open for appending
        fstat=open(statfile,'a')
    else:
        # Open for writing
        fstat=open(statfile,'w')

    # Extract time only from image name
    imgtime=image[0:14]

    for key in scistat.keys():
        lineout="%14s  %12s  %s" % (imgtime,key,str(scistat[key]))
        fstat.write(lineout+"\n")

    fstat.close()

####################

def p60parsebin(binstring):

    bin=[1,1]
    re1=re.search("(\d+) +(\d+)",binstring)
    if re1:
        try:
            bin=[int(re1.group(1)),int(re1.group(2))]
        except:
            print "Failed to parse binning string %s... setting to default" % \
                  binstring
    else:
        print "Failed to parse binning string %s... setting to default" % \
              binstring

    return bin

####################

def p60parseroi(detsec):

    # Default region of interest
    roi=[1,2048,1,2048]

    re1=re.search("\[(\d+):(\d+),(\d+):(\d+)\]",detsec)
    if re1:
        try:
            llc=[int(re1.group(1)),int(re1.group(3))]
            urc=[int(re1.group(2)),int(re1.group(4))]
            roi=[llc[0],llc[1],
                 urc[0]-llc[0]+1,urc[1]-llc[1]+1]
        except:
            print "Failed to parse ROI string %s... setting to default" % \
                  detsec
    else:
        print "Failed to parse ROI string %s... setting to default" % \
              detsec

    return roi

####################

def p60roiflat(fullflat,outflat,bin,roi):

    """p60roiflat(fullflat,outflat,bin,roi):
       converts full-chip flatfield fullflat into binning/roi
       flatfield outflat by extracting & summing"""

    ff=pyfits.open(fullflat)

    lenx=roi[2]
    cutx=lenx/bin[0]
    leny=roi[3]
    cuty=leny/bin[1]

    # First extract the ROI
    dfull=ff[0].data
    fshape=dfull.shape
    if leny<=fshape[0] and lenx<=fshape[1]:
        # Note IRAF convention:  Indexing from (1,1)
        #    + Python convention:  Last named index is dropped
        dnew=dfull[roi[1]-1:roi[1]-1+roi[3],roi[0]-1:roi[0]-1+roi[2]]
    else:
        print "ROI is larger than the full frame I am working with"
        dnew=dfull.copy()

    # Then average over the binned pixels
    if bin[0]>1 or bin[1]>1:
        # Average pixel values for binning
        dcut=numarray.resize(dnew,[cuty,cutx])
        dcut[:,:]=0.0
        ncmb=bin[0]*bin[1]
        # Use array slices to average over the binsize
        for i in range(0,bin[0]):
            for j in range(0,bin[1]):
                dcut[:,:] += dnew[j::bin[1],i::bin[0]]
        dcut /= float(ncmb)
        ff[0].data=dcut
    else:
        # No averaging to be done
        ff[0].data=dnew

    # Write the new file
    ff.writeto(outflat)
    ff.close()

####################

def p60roibpm(fullbpm,outbpm,bin,roi):

    """p60roibpm(fullbpm,outbpm,bin,roi):
       converts full-chip bad-pixel mask fullbpm into binning/roi
       bad-pixel mask outbpm by extracting & summing

       Note:  Works only with FITS files, not .pl files."""
    
    ff=pyfits.open(fullbpm)

    lenx=roi[2]
    cutx=lenx/bin[0]
    leny=roi[3]
    cuty=leny/bin[1]

    # First extract the ROI
    dfull=ff[0].data
    fshape=dfull.shape
    if leny<=fshape[0] and lenx<=fshape[1]:
        # Note IRAF convention:  Indexing from (1,1)
        #    + Python convention:  Last named index is dropped
        dnew=dfull[roi[1]-1:roi[1]-1+roi[3],roi[0]-1:roi[0]-1+roi[2]]
    else:
        print "ROI is larger than the full frame I am working with"
        dnew=dfull.copy()

    # Then sum over the binned pixels
    if bin[0]>1 or bin[1]>1:
        dcut=numarray.resize(dnew,[cuty,cutx])
        dcut[:,:]=0
        # Use array slices to average over the binsize
        for i in range(0,bin[0]):
            for j in range(0,bin[1]):
                dcut[:,:] += dnew[j::bin[1],i::bin[0]]
        # Clip the values to the range [0,1]
        ff[0].data=numarray.clip(dcut,0,1)
    else:
        # No averaging to be done
        ff[0].data=dnew

    # Write the new file
    ff.writeto(outbpm)
    ff.close()

##################################################

def iqp60(inpat,endhow="never",endwhen="",focmode="all",
          logfile="",clobber=globclob,verbose=globver):

    """ intelligent processing of P60 data files """

    # Defaults
    twait=30
    reduced={}
    filtkey="FILTER"
    biaskey="IMGTYPE"
    biasre="BIAS"
    biasroot="Bias"
    biaspfx="b"
    flatkey="IMGTYPE"
    flatre="DOMEFLAT"
    flatpre="Flat-"
    flatpfx="f"
    bpm0root="BPM0"
    bpmroot="BPM"
    binkey="CCDSUM"
    detseckey="DETSEC"
    binroikey="BINROI"
    binroidef="A"
    nfocus=10
    focuskey="IMGTYPE"
    focusre="^FOCUS"
    fseqkey="FOCUSSEQ"
    skipkey="IMGTYPE"
    skipre="SAOFOCUS"
    statsec=""
    reqkeys=['IMGTYPE','OBJECT','FILTER','BINROI',
             'CCDSUM','DETSEC','ROISEC']
    sigma=3.0
    minraw=0.0
    maxraw=65535.0
    satval=50000.0
    masksfx="mask"
    rawpix=0.3787
    zptkey="ZEROPT"
    zpukey="ZEROPTU"
    # For now:  
    #catmags={'B':'BMAG','V':'BMAG','g':'BMAG',
    #         'R':'RMAG','r':'RMAG','I':'IMAG',
    #         'i':'IMAG','ip':'IMAG','z':'IMAG',
    #         'zp':'IMAG','zl':'IMAG','zs':'IMAG'}
    ## In the future, once Sloan mags are written to the ubcone output:
    catmags={'B':'BMAG','V':'VMAG','g':'GPMAG',
             'R':'RMAG','r':'RPMAG','I':'IMAG',
             'i':'IMAG','ip':'IPMAG','z':'ZPMAG',
             'zp':'ZPMAG','zl':'ZPMAG','zs':'ZPMAG',
             'upr':'UPMAG','gpr':'GPMAG','rpr':'RPMAG',
             'ipr':'IPMAG','zpr':'ZPMAG'}
    exptmkey="EXPTIME"
    exptmstd=60.0
    minzpts={}
    minskys={}
    fwhm0=1.5
    gain=2.3
    readn=7.5
    crsfx="crm"
    statfile="Status_Science.txt"

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

    # Parse focmode setting
    fonly=0
    if focmode=="nofocus":
        skipre="FOCUS"
    elif focmode=="focus":
        fonly=1
        skipre="SCIENCE"

    # Parse logfile setting
    if len(logfile)>0:
        check_exist(logfile,'w',yes)
        try:
            log=open(logfile,'w')
            sys.stdout=log
        except:
            print "Failed to open logfile %s for writing" % logfile

    # Setup ccdproc options
    ccdproc=iraf.ccdred.ccdproc
    ccdproc.overscan=yes
    ccdproc.trim=yes
    ccdproc.zerocor=no

    ccdproc.interactive=no
    ccdproc.function="chebyshev"
    ccdproc.order=3
    ccdproc.sample="*"
    ccdproc.naverage=3
    ccdproc.niterate=2
    ccdproc.low_reject=3.0
    ccdproc.high_reject=3.0
    ccdproc.grow=0.0

    # Demosaic all the files
    allfiles=glob.glob(inpat)
    alllist=','.join(allfiles)
    iraf.iqmosaic(alllist,outpfx="",biassec="!BIASSEC",
                  trimsec="!TRIMSEC",joinkeys="CCDSIZE,DETSEC,ROISEC",
                  splitkeys="AMPSEC,CCDSEC,TRIM,OVERSCAN,CCDPROC",
                  clobber=yes,verbose=verbose)
    iraf.add_wcs(alllist,instrument="p60new",binkey=binkey)

    # Set ccdproc params for bias-subtraction only
    ccdproc.overscan=no
    ccdproc.trim=no
    ccdproc.zerocor=yes

    # Determine if we need to make calibration files
    dobias=no
    doflats=no
    
    for image in allfiles:

        # Criteria for bias image
        biasval,flatval=get_head(image,[biaskey,flatkey])
        if re.search(biasre,biasval,re.I):
            dobias=yes
            reduced[image]=yes
        elif re.search(flatre,flatval,re.I):
            doflats=yes
            reduced[image]=yes

    # Calibrations
    if (dobias or doflats) and not fonly:
        iraf.iqcals(alllist,dobias=dobias,biasproc=no,
                    biaskey=biaskey,biasre=biasre,biasroot=biasroot,
                    doflats=doflats,flatproc=yes,flatkey=flatkey,
                    flatre=flatre,filtkey=filtkey,flatpre=flatpre,
                    flatscale="mode",statsec="",normflat=yes,
                    dobpm=no,bpmroot=bpmroot,classkey=binroikey,
                    classdef=binroidef,mosaic=no,
                    clobber=clobber,verbose=verbose)

    # Delete binroi-specific BPMs and flatfields
    oldbpms=glob.glob(bpmroot+'-'+'*.pl')
    for oldbpm in oldbpms:
        iraf.imdel(oldbpm,verify=no,go_ahead=yes)
    oldflats=glob.glob(flatpre+'*-*.fits')
    for oldflat in oldflats:
        iraf.imdel(oldflat,verify=no,go_ahead=yes)

    # Clip all flatfields and update BPM as we go
    bpm0name=bpm0root+'.pl'
    bpmname=bpmroot+'.pl'
    if os.path.exists(bpmname):
        check_exist(bpm0name,'w',clobber=yes)
        iraf.imrename(bpmname,bpm0name,verbose=no)
        allflats=glob.glob(flatpre+'*.fits')
        for flatimg in allflats:
            iraf.iqclip(flatimg,lthresh=0.75,hthresh=1.25,bookend=no,
                        replace=1.0,maskin=bpm0name,maskout=bpmname,
                        maskval=1,clobber=yes,verbose=no)
            iraf.imdel(bpm0name,verify=no,go_ahead=yes)
            iraf.imrename(bpmname,bpm0name,verbose=no)
            update_head(flatimg,"CLIPFLAT",1,"Has flatfield been clipped?")
            update_head(flatimg,"CLIPREPL",1.0,
                        "Replacement value for clipped pixels")
        iraf.imrename(bpm0name,bpmname,verbose=no)

    ##################################################

    # Lists of focusfiles
    allfocus=[]
    focusfiles=[]
    # Should additional calibrations be taken later on...
    biasfiles=[]
    flatfiles=[]

    # Big Loop
    done=no
    while not done:

        # Parse inputs
        allfiles=glob.glob(inpat)

        newfiles=[]
        for image in allfiles:
            if not reduced.has_key(image):
                # Exclude Bias & Flats
                if re.search(biasroot,image,re.I) or \
                   re.search(flatpre,image,re.I):
                    reduced[image]=yes
                    continue
                # Exclude calibration files (for now)
                biasval,flatval=get_head(image,[biaskey,flatkey])
                if re.search(biasre,biasval,re.I):
                    biasfiles.append(image)
                    reduced[image]=yes
                elif re.search(flatre,flatval,re.I):
                    flatfiles.append(image)
                    reduced[image]=yes
                # Exclude "skippable" files
                skipval=get_head(image,skipkey)
                if re.search(skipre,skipval,re.I):
                    reduced[image]=yes
                # Queue file for processing
                else:
                    reduced[image]=no
                    newfiles.append(image)
            elif not reduced[image]:
                # Attempted reduction before, but failed?
                newfiles.append(image)

        if newfiles:

            # Brief pause to avoid file conflicts (?)
            time.sleep(3.0)

            # Demosaic
            newlist=','.join(newfiles)
            ccdproc.overscan=yes
            ccdproc.trim=yes
            ccdproc.zerocor=no
            iraf.iqmosaic(newlist,outpfx="",biassec="!BIASSEC",
                          trimsec="!TRIMSEC",joinkeys="DETSEC",
                          splitkeys="AMPSEC,CCDSEC,TRIM,OVERSCAN",
                          clobber=yes,verbose=verbose)

            iraf.add_wcs(newlist,instrument="p60new",binkey=binkey)

        # Reset ccdproc params for bias-subtraction only
        ccdproc.overscan=no
        ccdproc.biassec=""
        ccdproc.trim=no
        ccdproc.trimsec=""
        ccdproc.zerocor=yes
        ccdproc.flatcor=no

        ##############################

        # Process new images
        for image in newfiles:

            # Required keyword check
            if not check_head(image,reqkeys):
                print "Image %s is missing required keywords" % image
                continue

            # Attempt processing only once per image
            reduced[image]=yes
            if verbose:
                print "Reducing new image %s" % image

            # Track useful status information
            scistat={}

            # Prepare for processing
            pix=rawpix
            binroi=binroidef
            biasname=biasroot+'.fits'
            bpmname=bpmroot+'.pl'
            filt=get_head(image,filtkey)
            flatname=flatpre+filt+".fits"

            # Calibration files for BINROI setting
            if check_head(image,binroikey):
                binroi=get_head(image,binroikey)
                if binroi != binroidef:
                    # Get the right calibration files
                    fullbias=biasname
                    biasname=biasroot+'-'+binroi+'.fits'
                    fullbpm =bpmroot+'.pl'
                    bpmname =bpmroot+'-'+binroi+'.pl'
                    fullflat=flatname
                    flatname=flatpre+filt+'-'+binroi+".fits"
                    # Details of this binning/roi combo
                    [binval,detsec]=get_head(image,[binkey,detseckey])
                    binning=p60parsebin(binval)
                    roi=p60parseroi(detsec)
                    # Pixel scale (no anisotropic binning, sorry)
                    pix=binning[0]*rawpix
                    # Check for existance of bias; make it if we have to
                    if not os.path.exists(biasname):
                        print ("Failed to find bias for %s=%s; "+ \
                               "constructing one from full-frame bias") % \
                              (binroikey,binroi)
                        p60roiflat(fullbias,biasname,binning,roi)
                    # Produce the flatfield if necessary
                    if not os.path.exists(flatname):
                        p60roiflat(fullflat,flatname,binning,roi)
                    # Produce the bad-pixel mask if necessary
                    if not os.path.exists(bpmname):
                        # Get BPM into FITS format
                        bpm1=iraf.mktemp("iqbpm")+".fits"
                        iraf.imcopy(fullbpm,bpm1,verbose=no)
                        # Adjust BPM for binning/ROI
                        bpm2=iraf.mktemp("iqbpm")+".fits"
                        p60roibpm(bpm1,bpm2,binning,roi)
                        # Convert BPM to PL format
                        iraf.imcopy(bpm2,bpmname,verbose=no)
                        # Delete temporary files
                        iraf.imdel(bpm1,verify=no,go_ahead=yes)
                        iraf.imdel(bpm2,verify=no,go_ahead=yes)

            # Bias subtraction
            check_exist(biasname,"r")
            ccdproc.zero=biasname
            image1=biaspfx+image
            check_exist(image1,'w',clobber)
            iraf.ccdproc(image,output=image1,Stdout=1)

            # Flatfielding
            check_exist(flatname,"r")
            iraf.iqflatten(image1,flatname,outpfx=flatpfx,
                           normflat=no,clipflat=no,
                           statsec=statsec,subsky=yes,
                           vignflat=no,clobber=yes,verbose=no)
            imagef=flatpfx+image1

            # Clip (bookend) flattened & sky-subtracted image
            skybkg=float(get_head(image1,'SKYBKG'))
            iraf.iqclip(image1,lthresh=minraw-skybkg,hthresh=maxraw-skybkg,
                        bookend=yes,maskin="",maskout="",clobber=yes,
                        verbose=no)

            # Rename the file, if everything is well-behaved
            reim=re.search('^%s%s(.+[^r])r?\.fits' % (flatpfx,biaspfx),
                           imagef,re.I)
            if reim:
                image2=reim.group(1)+'p.fits'
                check_exist(image2,"w",yes)
                iraf.imrename(imagef,image2,verbose=verbose)
            else:
                image2=imagef

            # Set bad pixels to zero
            update_head(image2,'BPM',bpmname,"Bad Pixel Mask")
            iraf.iqmask(image2,mask='!BPM',method='constant',value='0.0',
                        clobber=yes,verbose=no)

            # Defringing happens later, in bunches

            # Cosmic Ray Rejection
            tmpimg=iraf.mktemp("iqcr")+".fits"
            check_exist("%s_%s.fits" % (image2[:15],crsfx),"w",clobber)
            iraf.lacos_im(image2,tmpimg,"%s_%s.fits" % (image2[:15],crsfx),
                          gain=gain,readn=readn,skyval=skybkg,sigclip=4.5,
                          sigfrac=0.5,objlim=1.0,niter=1)
            iraf.imdel(image2,verify=no,go_ahead=yes)
            iraf.imcopy(tmpimg,image2,verbose=no)
            iraf.imdel(tmpimg,verify=no,go_ahead=yes)

            # Object-detection
            if scistat.has_key("SEEING_%s" % file):
                try:
                    fwhm=float(scistat["SEEING_%s" % file])
                except:
                    fwhm=fwhm0
            else:
                fwhm=fwhm0
            iraf.iqobjs(image2,sigma,satval,skyval="!SKYBKG",masksfx=masksfx,
                        wtimage="",fwhm=fwhm,pix=rawpix,aperture=2*fwhm/rawpix,
                        gain=gain,minlim=no,clobber=yes,verbose=no)

            # Track some status information
            nstars=get_head(image2,'NSTARS')
            scistat['NSTARS']=nstars
            if minskys.has_key(filt):
                minsky=minskys[filt]
                scistat["SKYBKG_%s" % filt]=skybkg/minsky
            
            # Is it a focus frame?
            focusval=get_head(image2,focuskey)

            # Focus frame actions
            if re.search(focusre,focusval,re.I):

                focusfiles.append(image2)
                # Find its place in the sequence
                focusseq=get_head(image2,fseqkey)
                ref=re.search("(\d+)of(\d+)",focusseq,re.I)
                if ref:
                    ifocus=int(ref.group(1))
                    nfocus=int(ref.group(2))

            # Non-focus frame actions
            else:

                # Make a new keyword for mosaic objects
                objkey='OBJECT'
                if check_head(image2,'MSCSIZE'):
                    mscsize=get_head(image2,'MSCSIZE')
                    mscdims=mscsize.split(',')
                    if int(mscdims[0])>1 or int(mscdims[1])>1:
                        object=get_head(image2,'OBJECT')
                        mscposn=get_head(image2,'MSCPOSN')
                        objmsc=object+'-'+mscposn
                        update_head(image2,'OBJMSC',objmsc,
                                    'Object Name/Mosaic Position')
                        objkey='OBJMSC'

                # Refine WCS (also calculates seeing in arcsec)
                iraf.iqwcs(image2,objkey=objkey,rakey='RA',
                           deckey='DEC',pixtol=0.05,starfile='!STARFILE',
                           nstar=40,catalog='web',ubhost="localhost",
                           diffuse=yes,clobber=yes,verbose=verbose,
                           nstarmax=40)

                # Retry if unsuccessful
                if not check_head(image2,'IQWCS'):
                    iraf.iqwcs(image2,objkey=objkey,rakey='RA',
                               deckey='DEC',pixtol=0.05,starfile="!STARFILE",
                               nstar=40,catalog='web',ubhost="localhost",
                               diffuse=yes,clobber=yes,verbose=verbose,
                               nstarmax=1000)

                # Calculate further quantities if WCS was successful
                extinct='INDEF'
                if check_head(image2,'IQWCS'):

                    scistat['WCSGOOD'] = 1
                    scistat["SEEING_%s" % filt]=get_head(image2,'SEEING') 

                    [rapt,dcpt,lenx,leny]=get_head(image2,
                                            ['RA','DEC','NAXIS1','NAXIS2'])
                    [[ractr,dcctr]]=impix2wcs(image2,0.5*(1+float(lenx)),
                                            0.5*(1+float(leny)))
                    coopt=astrocoords(rapt,dcpt)
                    cooctr=astrocoords(ractr,dcctr)
                    [offra,offdc]=radecdiff(coopt.radeg(),coopt.dcdeg(),
                                            cooctr.radeg(),cooctr.dcdeg())
                    update_head(image2,'PTOFFSET',"%.2f %.2f" % \
                                (offra,offdc),
                                "Offset from nominal pointing (arcsec)")

                    scistat['PTOFFSET']="%.2f,%.2f" % (offra,offdc)

                    # Get zero-point against SDSS/NOMAD/USNO-B 
                    if catmags.has_key(filt):
                        
                        if check_head(image2,'OBJMSC'):
                            object=get_head(image2,'OBJMSC')
                        else:
                            object=get_head(image2,'OBJECT')
                        zpcat=""

                        # See if attempt already made for SDSS
                        if not (os.path.exists("%s.sdss" % object)):
                            getsdss(image2, "%s.sdss" % object)

                        # If successful, set to SDSS
                        sdssfile=open("%s.sdss" % object)
                        if (len(sdssfile.readlines())>2):
                            catalog="%s.sdss" % object
                            zpcat="SDSS"

                        # If SDSS unsuccessful, update USNO-B 
                        if not zpcat:
                            update_usnob("%s.cat" % object, 
                                         outfile="%s.reg" % object)
                            catalog="%s.reg" % object
                            zpcat="USNO-B"
                       
                        # Update image header
                        update_head(image2,"ZPCAT",zpcat)

                        # Run the zeropoint calculation
                        iraf.iqzeropt(image2,catmags[filt],
                                      starfile="!STARFILE",catalog=catalog,
                                      pixtol=3.0,useflags=yes,maxnum=50,
                                      method="mean",rejout=1,fencelim=0.50,
                                      sigma=1.5,maxfrac=0.25,
                                      zptkey=zptkey,zpukey=zpukey,
                                      clobber=yes,verbose=verbose)

                        # Use that to get an extinction estimate
                        zeropt,zeroptu,exptime=get_head(image2,[zptkey,
                                                        zpukey,exptmkey])
                        zeropt=float(zeropt)
                        zeroptu=float(zeroptu)
                        exptime=float(exptime)

                        scistat["ZEROPT_%s" % filt] = \
                            zeropt - 2.5*math.log10(exptime/exptmstd)
                        
                        if minzpts.has_key(filt):
                            extinct = zeropt - \
                                      2.5*math.log10(exptime/exptmstd) - \
                                      minzpts[filt]
                            scistat["EXTINCT_%s" % filt]=extinct

                # Fix the SEEING keyword if WCS fit was not successful
                else:
                    scistat['WCSGOOD'] = 0

                    try:
                        seepix=float(get_head(image2,'SEEPIX'))
                        seeing= "%.3f" % (pix*seepix)
                    except:
                        seeing='INDEF'

                    update_head(image2,'SEEING',seeing,
                                "Estimated seeing in arcsec or INDEF")

            # Write Archive Keywords
            update_head(image2,'PROCESSD',1,
                        "Image has been processed by iqp60")
            update_head(image2,'PROCVER',version,
                        "Version number of iqp60 used")
            update_head(image2,'EXTINCT',extinct,
                        "Estimated mags extinction or INDEF")
            update_head(image2,'PROCPROB','None',
                        "Problems encountered in iqp60 processing")

            # Update status file
            p60status(statfile,image2,scistat)

            # Clean up
            check_exist(image1,"w",yes)

        ##############################

        # If a set of focus images have been processed, find best focus
        if len(focusfiles)>=nfocus:

            print "%d focus images have been processed" % len(focusfiles)

            # Pick a starfile
            focusfiles.sort()
            medfile=focusfiles[len(focusfiles)/2]
            starfile=get_head(medfile,'STARFILE')

            # Update seeing values in all images based on good
            # stars in the starfile.   
            ffiles=','.join(focusfiles)
            iraf.iqseeing(ffiles,stars=starfile,
                          seekey="SEEING",method="acorr",useflags=yes,
                          skipstars=0,usestars=10,boxsize=64,
                          strictbox=no,imgsec="[940:1980,64:1980]",
                          update=yes,clobber=yes,verbose=no)

            # Grab focus positions and seeing values
            focx=[]
            focy=[]
            for image in focusfiles:
                [fpos,seep]=list2float(get_head(image,['FOCUSPOS','SEEING']))
                fpos=float(fpos)
                seep=float(seep)
                focx.append(fpos)
                pix=rawpix
                if check_head(image,binkey):
                    binval=get_head(image,binkey)
                    binels=binval.split()
                    binning=[int(binels[0]),int(binels[1])]
                    # Let's hope we don't ever do asymmetric binning
                    pix=binning[0]*rawpix
                focy.append(pix*seep)
                update_head(image,'SEEING',pix*seep,
                            "Estimated seeing in arcsec")

            # Choose best-seeing as minimum of these
            bestsee=min(focy)
            bestfoc=focx[focy.index(bestsee)]

            # This was too much verbosity for Brad
            if verbose and 0:
                print "Focus settings:  "+liststr(focx,"%.2f")
                print "Seeing values:   "+liststr(focy,"%.1f")
                print
                print "Best focus:  %.2f" % bestfoc

            # Expunge the focusfiles list
            allfocus.append(focusfiles)
            focusfiles=[]

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
_parfile=pardir + "iqp60.par"
t=iraf.IrafTaskFactory(taskname="iqp60",
                       value=_parfile,function=iqp60)

######################################################################
######################################################################

def main():

    # Command line
    try:
        opts, args = getopt.getopt(sys.argv[1:], 
                                   "how:f:l:",
                                   ["help","once","when=","focmode=",
                                    "log"])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    endhow="never"
    endwhen=""
    focmode="all"
    logfile=""

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
        # Focus mode
        elif opt in ("-f", "--focmode"):
            focmode=val
        elif opt in ("-l", "--log"):
            logfile=val
        else:
            sys.stderr.write("Unmatched option %s\n" % opt)

    inpat=args[0]

    print "Running iqp60 with inpat=%s" % inpat
    iqp60(inpat,endhow=endhow,endwhen=endwhen,focmode=focmode,
          logfile=logfile)
    sys.exit(0)

######################################################################

def usage():

    (xdir,xname)=os.path.split(sys.argv[0])
    print "Usage:  %s [-o | -w HH:MM] [-f focmode] \"filepattern\"" % xname
    print
    print "   -o          single processing loop"
    print "   -w HH:MM    run until UT time HH:MM"
    print "   -f focmode  'nofocus' to skip all focus images"
    print "               'focus' to process only focus images"
    print "               'all' (default) to process all images"
    print "   -l logfile  log operations to logfile"

##############################

# Running as executable
if __name__=='__main__':

    main()
