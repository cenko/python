#!/usr/bin/env python

######################################################################
# $Id: iqnirc.py,v 1.2 2005/02/24 19:39:16 derekfox Exp derekfox $
######################################################################

version='1.2'

# Global packages
import pyraf
from pyraf import iraf
import copy, os, shutil, glob, sys, string, re, math
import getopt
import pyfits
#import numarray
import numpy.numarray as numarray
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
iraf.imfit()
iraf.proto()
iraf.stsdas()
iraf.stsdas.hst_calib()
iraf.stsdas.hst_calib.nicmos()

# Shortcuts
yes=iraf.yes
no=iraf.no
INDEF=iraf.INDEF
hedit=iraf.hedit
imgets=iraf.imgets

# Pyraf parameters: Where to find them
pyrafdir="python/pyraf/"
pyrafdir_key='PYRAFPARS'

if os.environ.has_key(pyrafdir_key):
    pardir=os.environ[pyrafdir_key]
else:
    pardir=os.environ['HOME']+'/'+pyrafdir

if pardir[-1] != '/':
    pardir += '/'

# Some local defaults
def_inpat="s?????.fits"
def_lampon="lowlamp"
def_lampoff="lampsoff"

globclob=yes
globver=yes

######################################################################

# Utility routines go here...

def qombo_split(qombo):

    hyph=qombo.rfind('-')
    if hyph<0:
        print "Improperly formed object-filter qombo '%s'" % qombo
        return

    obj,filt=qombo[0:hyph],qombo[hyph+1:]

    return [obj,filt]

####################

def nirc_filter(file):

    filtk1="FWONAME"
    filtk2="FWINAME"
    open="clear"
    blocked="BLOCK"

    try:
        [filt1,filt2]=get_head(file,[filtk1,filtk2])
    except:
        print "Error accessing filter keywords for %s" % file
        return ""

    if filt1==open:
        #re1=re.search('^([^_]+)_',filt2)
        #filter=re1.group(1)
        filter=filt2
    elif filt2==open:
        #re1=re.search('^([^_]+)_',filt1)
        #filter=re1.group(1)
        filter=filt1
    else:
        filter=blocked

    return filter

##################################################

def nirc_sky(image,bpm):

    usebpm=0
    if len(bpm)>0:
        re1=re.search("^\!(\w+)",bpm)
        if re1:
            bpmkey=re1.group(1)
            bpmfile=get_head(image,bpmkey)
        else:
            bpmfile=bpm
        if len(bpmfile)>0 and os.path.exists(bpmfile):
            usebpm=1

    if usebpm:
        mimstat=iraf.mimstatistics(image,imasks='!BPM',omasks="",
                     fields='image,npix,mean,stddev,min,max,mode',
                     lower='INDEF',upper='INDEF',nclip=4,
                     lsigma=5.0,usigma=5.0,binwidth=0.1,
                     format=no,Stdout=1,Stderr=1)
        mimels=mimstat[0].split()
        skybkg=None

        if len(mimels)==6:
            try:
                skybkg=float(mimels[5])
            except:
                pass

        if not skybkg:
            usebpm=0
            
    if not usebpm:
        iraf.iterstat(image,nsigrej=5,maxiter=10,
                      prin=no,verbose=no)
        skybkg=float(iraf.iterstat.median)

    return skybkg

##################################################

def nirc_flatscale(infiles,inflat,outflat,clipfrac=0.03,
                   clobber=globclob):

    # "infiles" is a list of filenames
    inlist=','.join(infiles)
    bpmfile=get_head(infiles[0],'BPM')

    # Make a median combination of input files
    imcmb1=iraf.mktemp("iqnircc")+".fits"
    iraf.imcombine(inlist,imcmb1,combine="median",reject="none",
                   project=no,outtype="real",outlimits="",offsets="",
                   masktype="none",blank=0.0,scale="!SKYBKG",zero="none",
                   weight="none",statsec="",lthreshold="INDEF",
                   hthreshold="INDEF",nkeep=1)
    [naxis1,naxis2]=get_head(imcmb1,['NAXIS1','NAXIS2'])
    npixall=float(naxis1)*float(naxis2)

    # Calculate sky background & divide, subtract 1
    skybkg=nirc_sky(imcmb1,bpmfile)
    imcmb2=iraf.mktemp("iqnircc")+".fits"
    iraf.imarith(imcmb1,'/',skybkg,imcmb2,verbose=no,noact=no)
    iraf.imdel(imcmb1,verify=no,go_ahead=yes)
    iraf.imarith(imcmb2,'-',1.0,imcmb1,verbose=no,noact=no)

    # Surface fit to median image
    imsurf1=iraf.mktemp("iqnircs")+".fits"
    iraf.imsurfit(imcmb1,imsurf1,xorder=6,yorder=6,type_output="fit",
                  function="chebyshev",cross_terms=yes,
                  xmedian=21,ymedian=21,median_percent=50.0,
                  lower=0.0,upper=0.0,ngrow=0,niter=0,
                  regions="all",rows="*",columns="*",border=50,
                  sections="",circle="",div_min="INDEF")
    # Corresponding bad pixel mask
    imbpm1=iraf.mktemp("iqnircb")+".pl"
    iraf.imexpr("abs(x)>%.3f ? 0 : 1" % clipfrac,imbpm1,x=imsurf1)

    # Subtract 1.0 from flatfield
    imflat1=iraf.mktemp("iqnircf")+".fits"
    iraf.imarith(inflat,'-',1.0,imflat1,verbose=no,noact=no)

    # Surface fit to the flatfield
    imsurf2=iraf.mktemp("iqnircs")+".fits"
    iraf.imsurfit(imflat1,imsurf2,xorder=6,yorder=6,type_output="fit",
                  function="chebyshev",cross_terms=yes,
                  xmedian=21,ymedian=21,median_percent=50.0,
                  lower=0.0,upper=0.0,ngrow=0,niter=0,
                  regions="all",rows="*",columns="*",border=50,
                  sections="",circle="",div_min="INDEF")
    # Corresponding bad pixel mask
    imbpm2=iraf.mktemp("iqnircb")+".pl"
    iraf.imexpr("abs(x)>%.3f ? 0 : 1" % clipfrac,imbpm2,x=imsurf2)

    # Combine bad pixel masks for median + flat
    imbpm3=iraf.mktemp("iqnircb")+".pl"
    iraf.imexpr("(x>0 || y>0) ? 1 : 0",imbpm3,x=imbpm1,y=imbpm2)

    # Calculate the ratio image
    imratio=iraf.mktemp("iqnircr")+".fits"
    iraf.imexpr("z>0 ? 0 : x/y",imratio,x=imsurf1,y=imsurf2,z=imbpm3)

    # Mimstat on the ratio image
    mimstat=iraf.mimstatistics(imratio,imasks=imbpm3,omasks="",
                     fields='image,npix,mean,stddev,min,max,mode',
                     lower='INDEF',upper='INDEF',nclip=4,
                     lsigma=5.0,usigma=5.0,binwidth=0.1,
                     format=no,Stdout=1,Stderr=1)
    mimels=mimstat[0].split()
    npix=float(mimels[1])
    xmult=float(mimels[2])

    # Check that a reasonable number of pixels have made the grade
    check_exist(outflat,'w',clobber=clobber)
    if npix<0.05*npixall:
        print "Less than 5% of pixels passed the cut... preserving flatfield"
        iraf.imcopy(inflat,outflat,verbose=no)
        xmult=1.0
    else:
        # Create the final flatfield image
        iraf.imexpr("x*%.3f + 1" % xmult,outflat,x=imflat1)

    # Update header keywords
    update_head(outflat,'RESCALE',1,"Flatfield has been rescaled")
    update_head(outflat,'ORIGFLAT',inflat,
                "Input flatfield name (before rescaling)")
    update_head(outflat,'XMULT',xmult,
                "Multiplied ORIGFLAT by this factor to rescale")

    # Clean up
    iraf.imdel(imcmb1,verify=no,go_ahead=yes)
    iraf.imdel(imcmb2,verify=no,go_ahead=yes)
    iraf.imdel(imsurf1,verify=no,go_ahead=yes)
    iraf.imdel(imsurf2,verify=no,go_ahead=yes)
    iraf.imdel(imbpm1,verify=no,go_ahead=yes)
    iraf.imdel(imbpm2,verify=no,go_ahead=yes)
    iraf.imdel(imbpm3,verify=no,go_ahead=yes)
    iraf.imdel(imflat1,verify=no,go_ahead=yes)
    iraf.imdel(imratio,verify=no,go_ahead=yes)

##################################################

class nircred:

    def __init__(self,inpat="s?????.fits",final=yes,
                 skipkey="OBJECT",skipre="FOCUS",statsec="",
                 objkey="OBJECT",filtset="FILTER",
                 sigma=2.0,satval=90000.0,pix=0.2461,masksfx="mask",
                 lampkey="OBJECT",lampon="lowlamp",
                 lampoff="lampsoff",lampfkey="LAMPFILT",
                 flatkey="OBJECT",flatre="Dome",flatpre="Flat-",
                 flatlow=0.2,flathigh=5.0,
                 onofflist=['Ks','K','J','H'],bpmroot="BPM",
                 bpmffilt="On-Ks,Off-Ks,On-H,On-J",
                 flatpfx="f",fringepre="Fringe-",coaddkey="FRMCOADD",
                 exptmkey="TINT",defrgpfx="s",frgamp=0.50,
                 clipreg="[1:1024,1:1024]",clpsfx="-cut",
                 medsfx="-med",bpmsfx="-bpm",crmsfx="-crm",
                 exppfx="E",xrgpfx="X",xshext="xsh",
                 catmags={'J':'JMAG','H':'HMAG','Ks':'KMAG','K':'KMAG'},
                 logfile="",clobber=globclob,verbose=globver):

        self.inpat=inpat
        self.final=final
        self.skipkey=skipkey
        self.skipre=skipre
        self.statsec=statsec
        self.objkey=objkey
        self.filtset=filtset
        self.sigma=sigma
        self.satval=satval
        self.pix=pix
        self.masksfx=masksfx
        self.lampkey=lampkey
        self.lampon=lampon
        self.lampoff=lampoff
        self.lampfkey=lampfkey
        self.flatkey=flatkey
        self.flatre=flatre
        self.flatpre=flatpre
        self.flatlow=flatlow
        self.flathigh=flathigh
        self.onofflist=onofflist
        self.bpmroot=bpmroot
        self.bpmffilt=bpmffilt
        self.flatpfx=flatpfx
        self.fringepre=fringepre
        self.coaddkey=coaddkey
        self.exptmkey=exptmkey
        self.defrgpfx=defrgpfx
        self.frgamp=frgamp
        self.clipreg=clipreg
        self.clpsfx=clpsfx
        self.medsfx=medsfx
        self.bpmsfx=bpmsfx
        self.crmsfx=crmsfx
        self.exppfx=exppfx
        self.xrgpfx=xrgpfx
        self.xshext=xshext
        self.catmags=catmags
        self.logfile=logfile
        self.clobber=clobber
        self.verbose=verbose

        # Parse logfile setting
        if len(logfile)>0:
            check_exist(logfile,'w',yes)
            try:
                self.log=open(logfile,'w')
                sys.stdout=self.log
                sys.stderr=self.log
            except:
                print "Failed to open logfile %s for writing" % logfile
        else:
            self.log=None

    ##############################

    def science_images(self,qombos=[],pfx=''):

        """ Returns list of relevant science images """

        science=[]
        allfiles=glob.glob(self.inpat)

        for image in allfiles:

            # Exclude flats by filename
            if re.search(self.flatpre,image,re.I):
                continue
            # Get image characteristics
            [obj,flatval,skipval]=get_head(image,[self.objkey,
                                           self.flatkey,self.skipkey])
            filt=nirc_filter(image)
            # Exclude Biases
            if re.search('BLOCK',filt,re.I):
                continue
            # Exclude calibration files (for now)
            if re.search(self.flatre,flatval,re.I):
                continue
            # Exclude "skippable" files
            if re.search(self.skipre,skipval,re.I):
                continue

            if len(qombos)>0:
                # Only interested in some images
                qombo=obj+'-'+filt
                if qombo in qombos:
                    science.append(image)
            else:
                # Interested in all images
                science.append(image)

        # Only interested in associated files that have the right prefix
        if len(pfx)>0:
            newsci=[]
            for image in science:
                if os.path.exists(pfx+image):
                    newsci.append(pfx+image)
            science=newsci

        # The end
        return science

    ####################

    def all_qombos(self):

        qombos=[]
        science=self.science_images()

        for image in science:
            obj=get_head(image,self.objkey)
            filt=nirc_filter(image)
            qombo=obj+'-'+filt
            if qombo not in qombos:
                qombos.append(qombo)

        return qombos
    ####################

    def all_objects(self):

        objs=[]
        science=self.science_images()

        for image in science:
            obj=get_head(image,self.objkey)
            if obj not in objs:
                objs.append(obj)

        return objs

    ####################

    def all_filters(self):

        filts=[]
        science=self.science_images()

        for image in science:
            filt=nirc_filter(image)
            if filt not in filts:
                filts.append(filt)

        return filts

    ####################

    def check_calib(self,filts):

        # Convert string input to list
        if type(filts) is StringType:
            filt=filts
            filts=[filt]

        # Check for BPM
        ok=os.path.exists(self.bpmroot+".pl")

        # Check for Flatfield images
        if ok:
            for filt in filts:
                qflat=self.flatpre+filt+'.fits'
                if not os.path.exists(qflat):
                    ok=0
                    break

        # Calibrations okay?
        return ok

    ##############################

    def calib(self):

        dolamps=no
        lampfilts=[]
        lampfiles=[]
        lamponfiles={}
        lampofffiles={}

        allfiles=glob.glob(self.inpat)

        # Look for flatfield images
        for image in allfiles:

            # Criteria for flatfield images
            lampval=get_head(image,self.lampkey)
            if re.search(self.lampon,lampval,re.I) or \
               re.search(self.lampoff,lampval,re.I):
                dolamps=yes
                filter=nirc_filter(image)
                if filter not in lampfilts:
                    lampfilts.append(filter)
                if re.search(self.lampon,lampval,re.I):
                    update_head(image,self.lampfkey,'On-'+filter,
                                "Lamp status + Filter setting")
                    lampfiles.append(image)
                    if lamponfiles.has_key(filter):
                        lamponfiles[filter].append(image)
                    else:
                        lamponfiles[filter]=[image]
                elif re.search(self.lampoff,lampval,re.I):
                    update_head(image,self.lampfkey,'Off-'+filter,
                                "Lamp status + Filter setting")
                    lampfiles.append(image)
                    if lampofffiles.has_key(filter):
                        lampofffiles[filter].append(image)
                    else:
                        lampofffiles[filter]=[image]

        lamplist=','.join(lampfiles)

        # Use flatfield images to make our calibration files
        if dolamps:

            # Make lamp-on and lamp-off flats & bad-pixel mask
            iraf.iqcals(lamplist,dobias=no,biasproc=no,
                        doflats=yes,flatproc=no,flatkey=self.flatkey,
                        flatre=self.flatre,filtkey=self.lampfkey,
                        flatpre=self.flatpre,flatscale="median",
                        statsec=self.statsec,normflat=no,
                        dobpm=yes,bpmmethod='flatratio',
                        bpmroot=self.bpmroot,
                        bpmffilt=self.bpmffilt,classkey="",
                        mosaic=no,clobber=self.clobber,
                        verbose=self.verbose)

            # Subtract lamp-off flat from lamp-on flat, if necessary
            for filter in lamponfiles.keys():

                onflat=self.flatpre+'On-'+filter+'.fits'
                offflat=self.flatpre+'Off-'+filter+'.fits'
                outflat=self.flatpre+filter+'.fits'
                check_exist(outflat,'w',self.clobber)

                if filter in self.onofflist:

                    if not lampofffiles.has_key(filter):
                        print "Need lamp-off flats for filter '%s'" % filter
                        continue

                    iraf.imarith(onflat,'-',offflat,outflat,
                                 verbose=self.verbose,noact=no)

                else:

                    # For some filters "lampsoff" counts are negligible
                    iraf.imcopy(onflat,outflat,verbose=yes)

                # Normalize flatfield
                medflat=nirc_sky(outflat,"")
                iraf.imarith(outflat,'/',medflat,outflat,verbose=no,
                             noact=no)

                # Update header keywords
                update_head(outflat,[self.filtset,'BPM','NORMFLAT','MEDFLAT'],
                            [filter,self.bpmroot+'.pl',1,medflat],
                            ["WIRC filter","Bad pixel mask",
                             "Has flatfield been normalized?",
                             "Median value of original flatfield"])

    ##############################
            
    def preproc(self,qombos):

        if type(qombos) is StringType:
            qombo=qombos
            qombos=[qombo]

        allfiles=self.science_images(qombos)

        for image in allfiles:

            if self.verbose:
                print "Preprocessing image %s" % image

            # Add BPM as a header keyword
            update_head(image,'BPM',self.bpmroot+'.pl',
                        "Bad Pixel Mask")

            # Calculate sky background
            skybkg=nirc_sky(image,'!BPM')
            update_head(image,'SKYBKG',skybkg,
                        'Value of sky background')

            # Rudimentary WCS
            iraf.add_wcs(image,instrument="nirc")

            # Done with Step #2
            update_head(image,'IQWRCSTP',2,
                        "Stage of IQWIRC processing")

    ##############################

    def flatten(self,qombos,flat=None,flatscale=no):

        if type(qombos) is StringType:
            qombo=qombos
            qombos=[qombo]

        # Loop over input list of qombos
        for qombo in qombos:

            # Check for presence of calibration files
            obj,filt=qombo_split(qombo)
            if not self.check_calib(filt):
                self.calib()

            # Corresponding files
            offiles=self.science_images(qombo)
            if len(offiles)==0:
                continue

            # Make a text ".lis" file to record this
            putlines(qombo+'.lis',offiles,clobber=yes)

            # Rescale flatfield if requested
            inflat=self.flatpre+filt+'.fits'
            outflat=self.flatpre+qombo+'.fits'
            check_exist(outflat,"w",clobber=self.clobber)

            if flatscale:
                nirc_flatscale(offiles,inflat,outflat,clipfrac=0.02,
                               clobber=yes)
            else:
                iraf.imcopy(inflat,outflat,verbose=no)
            check_exist(outflat,"r")

            # Clip flatfield and update BPM
            iraf.iqclip(outflat,lthresh=self.flatlow,
                        hthresh=self.flathigh,bookend=no,
                        replace=1.0,maskin='!BPM',maskval=1,
                        maskout=self.bpmsfx)

            # Rename the updated BPM
            newbpm1=get_head(outflat,'BPM')
            newbpm2=newbpm1.replace(self.flatpre,'')
            check_exist(newbpm2,'w',clobber=self.clobber)
            iraf.imrename(newbpm1,newbpm2,verbose=no)
            update_head(outflat,'BPM',newbpm2)

            # Further preprocessing
            for image in offiles:

                # Keyword updates
                update_head(image,[self.filtset,'BPM'],
                                  [filt,newbpm2])

                # Flatfielding
                iraf.iqflatten(image,outflat,outpfx=self.flatpfx,
                               normflat=no,statsec=self.statsec,
                               subsky=yes,vignflat=no,clipflat=no,
                               clobber=yes,verbose=no)

                # Set bad pixels to zero
                image2=self.flatpfx+image
                iraf.iqmask(image2,mask='!BPM',method='constant',value=0.0,
                            clobber=yes,verbose=no)

                # Write Useful Keywords
                update_head(image2,'PROCESSD',1,
                            "Image has been processed by iqnirc")
                update_head(image2,'PROCVER',version,
                            "Version number of iqnirc used")
                update_head(image2,'ZEROPT','INDEF',
                            "Zero-point relative to 2MASS or INDEF")
                update_head(image2,'PROCPROB','None',
                            "Problems encountered in iqnirc processing")

            # List of flatfielded images
            offiles2=pfx_list(offiles,self.flatpfx)
            oflist2=','.join(offiles2)

            # Object-detection
            iraf.iqobjs(oflist2,self.sigma,self.satval,
                        skyval="!SKYBKG",masksfx=self.masksfx,
                        wtimage="",minlim=no,
                        clobber=yes,verbose=no)

            # Attempt WCS refinement
            iraf.iqwcs(oflist2,objkey=self.objkey,rakey='RA',
                       deckey='DEC',pixscl=self.pix,pixtol=0.05,
                       starfile='!STARFILE',catalog='ir',
                       diffuse=yes,clobber=yes,verbose=self.verbose)

            # Done with Step #3
            update_head(offiles2,'IQWRCSTP',3,
                        "Stage of IQWIRC processing")

    ##############################

    def mkfringe(self,qombos):

        if type(qombos) is StringType:
            qombo=qombos
            qombos=[qombo]

        # Loop over input list of qombos
        for qombo in qombos:

            # Corresponding files
            offiles2=self.science_images(qombo,pfx=self.flatpfx)
            if len(offiles2)==0:
                continue
            oflist2=','.join(offiles2)
            offrg=self.fringepre+qombo+'.fits'

            # Construct fringe image
            iraf.iqfringe(oflist2,offrg,combine="average",
                          fringeamp=self.frgamp,
                          skykey="SKYBKG",subkey="SKYSUB",bpmkey="MASKNAME",
                          masktype="goodvalue",maskvalue=0,
                          clobber=self.clobber,verbose=no)

            # Done with Step #4
            update_head(offiles2,'IQWRCSTP',4,
                        "Stage of IQWIRC processing")

    ##############################

    def defringe(self,qombos):

        if type(qombos) is StringType:
            qombo=qombos
            qombos=[qombo]

        # Loop over input list of qombos
        for qombo in qombos:

            # Corresponding files
            offiles2=self.science_images(qombo,pfx=self.flatpfx)
            if len(offiles2)==0:
                continue
            oflist2=','.join(offiles2)

            # Check for existence of fringe image
            offrg=self.fringepre+qombo+'.fits'
            if not os.path.exists(offrg):
                self.mkfringe(qombo)

            # Defringe images
            iraf.iqdefringe(oflist2,offrg,outpfx=self.defrgpfx,
                            skykey="SKYBKG",minlim=no,
                            clobber=self.clobber,verbose=no)

            # List of defringed images
            offiles3=pfx_list(offiles2,self.defrgpfx)
            oflist3=','.join(offiles3)

            # Clipping of defringed images
            clipreg=self.clipreg
            if len(clipreg)>0:

                re1=re.search('(\d+):(\d+),(\d+):(\d+)',clipreg)

                if re1:
                    szx=int(re1.group(2))-int(re1.group(1))+1
                    szy=int(re1.group(4))-int(re1.group(3))+1
                
                    for image in offiles3:

                        # Check that it's not already clipped
                        naxis1,naxis2=get_head(image,['NAXIS1','NAXIS2'])
                        if naxis1<=szx and naxis2<=szy:
                            continue

                        # Clip the defringed images
                        iraf.imcopy(image+clipreg,image,verbose=no)

                        # Clip the corresponding BPM
                        oldbpm=get_head(image,'BPM')
                        newbpm=oldbpm.replace('.pl',self.clpsfx+'.pl')
                        if os.path.exists(oldbpm) and \
                               not os.path.exists(newbpm):
                            iraf.imcopy(oldbpm+clipreg,newbpm,verbose=yes)

                        # Update the image headers
                        update_head(image,'BPM',newbpm)

                else:
                    print "Failed to parse clipping region... not clipping"

            ####################

            # Attempt WCS refinement
            for image in offiles3:

                # Skip images with good WCS
                if check_head(image,'IQWCS'):
                    if int(get_head(image,'IQWCS')):
                        continue

                # Object-detection
                iraf.iqobjs(image,self.sigma,self.satval,skyval="!SKYBKG",
                            masksfx=self.masksfx,wtimage="",minlim=no,
                            clobber=yes,verbose=no)

                # WCS refinement
                iraf.iqwcs(image,objkey=self.objkey,rakey='RA',
                           deckey='DEC',pixscl=self.pix,pixtol=0.05,
                           starfile='!STARFILE',catalog='ir',
                           diffuse=yes,clobber=yes,verbose=self.verbose)

            # Figure out which have good WCS
            wcsdict=keysplit(offiles3,'IQWCS')

            # Interpolate WCS onto files without good WCS so far
            if wcsdict.has_key(""):
                wcsbadfiles=wcsdict[""]
                if len(wcsbadfiles)==len(offiles3):
                    print "No images with good WCS in this set"
                else:
                    wcsgoodfiles=wcsdict[1]
                    wcsinterp(wcsbadfiles,wcsgoodfiles)

            # Done with Step #5
            update_head(offiles3,'IQWRCSTP',5,
                        "Stage of IQWIRC processing")

    ##############################

    def coadd(self,qombos):

        if type(qombos) is StringType:
            qombo=qombos
            qombos=[qombo]

        # Are we doing a final combination or median-only?
        if self.final:
            medonly=no
        else:
            medonly=yes
        
        # Loop over input list of qombos
        for qombo in qombos:

            # Corresponding files
            ffpfx=self.defrgpfx+self.flatpfx
            offiles3=self.science_images(qombo,pfx=ffpfx)
            if len(offiles3)==0:
                continue
            oflist3=','.join(offiles3)

            # Perform the coaddition
            output=qombo+'.fits'
            iraf.iqcoadd(oflist3,output,medsfx=self.medsfx,
                         xshext=self.xshext,bpm="!BPM",exppfx=self.exppfx,
                         xrgpfx=self.xrgpfx,bpmsfx=self.bpmsfx,
                         crmsfx=self.crmsfx,expkey=self.exptmkey,
                         window=15,crmeth="sigma",crthresh=7.0,
                         overres=2,useold=yes,medonly=medonly,
                         trimout=yes,trimmax=no,subpixel=yes,
                         doclean=yes,clobber=self.clobber,verbose=no)

            # Done with Step #6
            update_head(output,'IQWRCSTP',6,
                        "Stage of IQWIRC processing")

    ##############################

    def zeropt(self,qombos):

        if type(qombos) is StringType:
            qombo=qombos
            qombos=[qombo]

        # Loop over input list of qombos
        for qombo in qombos:

            # Do the split
            obj,filt=qombo_split(qombo)

            # Corresponding files
            image=qombo+'.fits'
            if not os.path.exists(image):
                continue
            wtimage=image.replace('.fits','.weight.fits')
            if not os.path.exists(wtimage):
                wtimage=""

            # Object-detection
            iraf.iqobjs(image,self.sigma,self.satval,skyval="!SKYBKG",
                        masksfx=self.masksfx,wtimage=wtimage,wtcut=0,
                        minlim=no,clobber=yes,verbose=no)
            starfile=get_head(image,'STARFILE')

            # Set pointer to NOMAD Starlist
            catfile=obj+'.cat'
            update_head(image,'CATALOG',catfile,
                        "Reference catalog for astrometry & photometry")

            # Estimate seeing using catalog stars
            stars=Starlist(starfile)
            catstars=Starlist(catfile)
            catstars.wcs2pix(image)
            match,catmatch=stars.match(catstars,tol=3.0,useflags=yes,
                                       image=image,maxnum=50)
            if len(match)>2:
                fwhmpix=median(match.fwhms())
                seeing=fwhmpix*self.pix
            else:
                seeing="INDEF"
            update_head(image,'SEEING',seeing,
                        "Estimated seeing in arcsec or INDEF")

            # Calculate zero-point, if possible
            if self.catmags.has_key(filt):
                iraf.iqzeropt(image,self.catmags[filt],
                              starfile="!STARFILE",catalog="!CATALOG",
                              pixtol=3.0,useflags=yes,maxnum=50,
                              zptkey="ZEROPT",zpukey="ZEROPTU",
                              clobber=yes,verbose=self.verbose)
            elif self.verbose:
                print "No zero-point estimation for filter '%s'" % filt

            # Magnitude of supernova?
            #########################

            # Done with Step #7
            update_head(image,'IQWRCSTP',7,
                        "Stage of IQWIRC processing")

######################################################################
######################################################################

def iqnirc(inpat="s?????.fits",final=yes,
           forcecals=no,lampon="lowlamp",
           lampoff="lampsoff",
           logfile="",objlist="",filtlist="",objflist="",
           steplist="",clobber=globclob,verbose=globver):

    """ intelligent processing of WIRC data """

    # Define the nircred object
    wred=nircred(inpat=inpat,final=final,lampon=lampon,lampoff=lampoff,
                 logfile=logfile,clobber=clobber,verbose=verbose)

    ##############################

    # Parse the input object, filter, and objfilt lists

    objects=[]
    if len(objlist)>0:
        objects=objlist.split(',')

    filters=[]
    if len(filtlist)>0:
        filters=filtlist.split(',')

    objfilts=[]
    if len(objflist)>0:
        objfilts=objflist.split(',')

    steps=range(2,8)
    if len(steplist)>0:
        steps=[]
        stepstr=steplist.split(',')
        for step in stepstr:
            steps.append(int(step))

    ####################

    # Full set of available qombos
    all_qombos=wred.all_qombos()
    qombos=[]

    if len(objects)==0 and len(filters)==0 and len(objfilts)==0:
        # Do everything
        qombos=all_qombos

    else:

        if len(objects)>0 or len(filters)>0:
            if len(objects)>0 and len(filters)==0:
                # All filters for these objects
                all_filters=wred.all_filters()
                for obj in objects:
                    for filt in all_filters:
                        qombo=obj+'-'+filt
                        if qombo not in qombos and qombo in all_qombos:
                            qombos.append(qombo)
            elif len(filters)>0 and len(objects)==0:
                # All objects in these filters
                all_objects=wred.all_objects()
                for filt in filters:
                    for obj in all_objects:
                        qombo=obj+'-'+filt
                        if qombo not in qombos and qombo in all_qombos:
                            qombos.append(qombo)
            else:
                # Each of these objects in each of these filters
                for obj in objects:
                    for filt in filters:
                        qombo=obj+'-'+filt
                        if qombo not in qombos and qombo in all_qombos:
                            qombos.append(qombo)

        if len(objfilts)>0:
            # Add specified qombos to the list
            for qombo in objfilts:
                if qombo not in qombos and qombo in all_qombos:
                    qombos.append(qombo)

    ##############################

    # Check for calibration files
    if not forcecals:
        allfilt=[]
        for qombo in qombos:
            obj,filt=qombo_split(qombo)
            if filt not in allfilt:
                allfilt.append(filt)
        calok=wred.check_calib(allfilt)

    # Calibrations before all
    if forcecals or not calok:
        wred.calib()

    # Process images object-filter by object-filter
    for qombo in qombos:
        if 2 in steps:
            wred.preproc(qombo)
        if 3 in steps:
            wred.flatten(qombo)
        if 4 in steps:
            wred.mkfringe(qombo)
        if 5 in steps:
            wred.defringe(qombo)
        if 6 in steps:
            wred.coadd(qombo)
        if 7 in steps:
            wred.zeropt(qombo)

######################################################################

# location of parameter files
_parfile=pardir + "iqnirc.par"
t=iraf.IrafTaskFactory(taskname="iqnirc",
                       value=_parfile,function=iqnirc)

######################################################################

def main():

    # Command line
    try:
        opts, args = getopt.getopt(sys.argv[1:], 
                          "hcdg:F:f:l:o:r:q:s:",
                         ["help","cals","draft","glob=","flaton=","flatoff=",
                          "log=","object=","filter=","qombo=","step="])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    # Pipeline chain
    procsteps={'calib':1,
               'preproc':2,
               'flatten':3,
               'mkfringe':4,
               'defringe':5,
               'coadd':6,
               'zeropt':7}

    # Defaults
    forcecals=no
    final=yes
    inpat=def_inpat
    lampon=def_lampon
    lampoff=def_lampoff
    logfile=""
    objects=[]
    filters=[]
    objfilts=[]
    steps=[]

    for opt, val in opts:
        # Usage info only
        if opt in ("-h", "--help"):
            usage()
            sys.exit(1)
        # New glob pattern
        elif opt in ("-g", "--glob"):
            inpat=val
        # Lamp-on object keyword
        elif opt in ("-F", "--flaton"):
            lampon=val
        # Lamp-on object keyword
        elif opt in ("-f", "--flatoff"):
            lampoff=val
        # Logfile spec
        elif opt in ("-l", "--log"):
            logfile=val
        # Redo/clobber calibrations
        elif opt in ("-c", "--cals"):
            forcecals=yes
        # "Draft" coaddition only
        elif opt in ("-d", "--draft"):
            final=no
        # Specify object
        elif opt in ("-o", "--object"):
            objs=val.split(',')
            for obj in objs:
                objects.append(obj)
        # Specify filter
        elif opt in ("-r", "--filter"):
            flts=val.split(',')
            for flt in flts:
                filters.append(flt)
        # Specify object-filter
        elif opt in ("-q", "--qombo"):
            qmbs=val.split(',')
            for qmb in qmbs:
                objfilts.append(qmb)
        # Specify processing steps
        elif opt in ("-s", "--step"):
            stps=val.split(',')
            for stp in stps:
                if re.search('\d',stp):
                    # Numerical steps
                    if re.search('\d+-\d+',stp):
                        # Inclusive range of steps
                        re1=re.search('(\d+)-(\d+)',stp)
                        for i in xrange(int(re1.group(1)),int(re1.group(2))+1):
                            steps.append(i)
                    else:
                        # Single step
                        steps.append(int(stp))
                else:
                    # Alpha steps
                    if procsteps.has_key(stp.lower()):
                        steps.append(procsteps[stp.lower()])
                    else:
                        print "Unrecognized processing step '%s' ignored" % \
                              stp
        else:
            sys.stderr.write("Unmatched option %s\n" % opt)

    # Run the program
    objlist=','.join(objects)
    filtlist=','.join(filters)
    objflist=','.join(objfilts)

    # Convert list of steps to string list
    ssteps=[]
    if len(steps)>0:
        for step in steps:
            ssteps.append("%d" % step)
    steplist=','.join(ssteps)
    
    iqnirc(inpat=inpat,final=final,forcecals=forcecals,
           lampon=lampon,lampoff=lampoff,
           logfile=logfile,objlist=objlist,filtlist=filtlist,
           objflist=objflist,steplist=steplist)

    sys.exit(0)

######################################################################

def usage():

    (xdir,xname)=os.path.split(sys.argv[0])
    print ("Usage:  %s [-h] [-c] [-d] [-g glob] [-F lampon] [-f lampsoff] "+
           "[-l logfile] [-o object] [-r filter] [-q obj-filt] "+
           "[-s step1,step2]") % xname
    print "   -h              print help (this text)"
    print "   -c              force recreation of all calibration files [no]"
    print "   -d              perform 'draft' coaddition only (median) [no]"
    print "   -g glob         change the glob pattern for raw images [%s]" %\
          def_inpat
    print "   -F lampon       object keyword for lamp-on flats [%s]" %\
          def_lampon
    print "   -f lampoff      object keyword for lamp-off flats [%s]" %\
          def_lampoff
    print "   -l logfile      log operations to logfile"
    print "   -o object       object name to process [all]"
    print "   -r filter       filter name to process [all]"
    print "   -q obj-filt     object-filter combination to process [all]"
    print "   -s step1,step2  processing steps as names or integers (1-7)"
    print
    print "Processing steps:"
    print "    1 - calib      Create calibration files"
    print "    2 - preproc    Preprocess raw images"
    print "    3 - flatten    Flatten images, BPM mask, WCS/1"
    print "    4 - mkfringe   Create fringe image"
    print "    5 - defringe   Apply fringe image, WCS/2"
    print "    6 - coadd      Coadd images"
    print "    7 - zeropt     Derive zero-pt by comparison to 2MASS"
    print

##############################

# Running as executable
if __name__=='__main__':

    main()
