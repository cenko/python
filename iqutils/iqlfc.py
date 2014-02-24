
import pyraf
from pyraf import iraf
import copy, os, shutil, glob, sys, string, re, math
import pyfits
import numarray
import time
from types import *

import add_wcs
import iqflatten
import iqobjs

yes=iraf.yes
no=iraf.no
INDEF=iraf.INDEF
hedit=iraf.hedit
imgets=iraf.imgets

pardir="/home/afterglow/derekfox/python/pyraf/"

globclob=yes
globver=yes

######################################################################

def iqlfc(inpat,clobber=globclob,verbose=globver):

    """ intelligent processing of LFC imaging data """

    # Necessary packages
    if not iraf.mscred._loaded:
        iraf.mscred()

    # Defaults
    twait=30
    reduced={}
    filtkey="FILTER"
    trimsec="[1:2048,1:2048]"
    biasname="Bias.fits"
    biaspfx="b"
    flatpre="Flat-"
    flatpfx="f"
    statsec=""
    sigma=2.0
    satval=25000.0
    masksfx="mask"
    pix=0.378

    # Perform the LFC assembly
    allfiles=glob.glob(inpat)

    for file in allfiles:
        if not file.endswith('.fits'):
            os.rename(file,file+".fits")

    assembled={}
    allnames=[]
    for file in allfiles:
        if file.endswith('.fits'):
            (root,num,chip,ext)=file.split('.')
        else:
            (root,num,chip)=file.split('.')
        if not assembled.has_key(num):
            name="%s.%s" % (root,num)
            iraf.lfcassemble(name)
            assembled[num]=yes
            allnames.append(name)

def dummy(dum1,dum2):

    # Setup ccdproc options
    ccdproc=iraf.mscred.ccdproc
    ccdproc.ccdtype=""
    ccdproc.noproc=no
    ccdproc.fixpix=no
    ccdproc.oversca=no
    ccdproc.trim=yes
    ccdproc.zerocor=yes
    ccdproc.flatcor=no
    ccdproc.illumco=no
    ccdproc.fringec=no
    ccdproc.readaxi="line"
    ccdproc.trimsec=trimsec
    ccdproc.zero=biasname

    # Basic sanity checks
    check_exist(biasname,"r")
    
    # Big Loop
    while(1):

        # Parse inputs
        allfiles=glob.glob(inpat)

        newfiles=[]
        for image in allfiles:
            if not reduced.has_key(image):
                newfiles.append(image)

        for image in newfiles:

            if verbose:
                print "Reducing new image %s" % image

            # Bias subtraction
            image1=biaspfx+image
            s1=ccdproc(image,output=image1,Stdout=1)

            # Flatfielding
            filt=get_head(image1,filtkey)
            flatname=flatpre+filt+".fits"
            check_exist(flatname,"r")
            iraf.iqflatten(image1,flatname,outpfx=flatpfx,
                           normflat=yes,statsec=statsec,vignflat=no,
                           clobber=yes,verbose=no)
            image2=flatpfx+image1

            # Rudimentary WCS
            iraf.add_wcs(image2,instrument="p60ccd")

            # Object-detection
            iraf.iqobjs(image2,sigma,satval,masksfx=masksfx,
                        wtimage="none",clobber=yes,verbose=no)

            # Refine WCS
            ############

            # Clean up
            check_exist(image1,"w")

            # Done with processing
            reduced[image]=yes
            
        # Wait a little while
        time.sleep(twait)

######################################################################
######################################################################

def check_exist(filename, status, clobber=globclob):

    """ check_exist(filename, status, clobber=yes)
    checks to see if filename exists
    if status==r, must exist, otherwise prints error and exits
    if status==w, if exists and clobber=no then prints error and exits
    else deletes
    """     

    if (status == "r"):
        # check to see if it exists for reading
        # (i.e. must be present)
        if (not (os.path.exists(filename))):
            print "Couldn't open input file: %s" % filename
            sys.exit(1)
    else:
        # check to see if it exists for writing
        # (i.e. must not exist or clobber=yes)
        if (os.path.exists(filename)):
            if (clobber):
                os.remove(filename)
            else:
                print "File %s already exists and clobber=no" % filename

######################################################################

def getlines(filename,blanks=no):

    """ lines=getlines(filename,blanks)
    opens file filename
    reads lines and strips off whitespace
    if (blanks==1), will return fully blank lines
    otherwise excises
    """

    infile=open(filename)
    lines=infile.readlines()
    infile.close()
    lines2=[]

    if (blanks==1):
        return lines
    
    for i in range(0,len(lines)):
        if (re.search("\S+",lines[i])):
            lines2.append(lines[i].rstrip())

    return lines2

######################################################################

def check_head(file,key,extn=0):

    check_exist(file,"r")

    try:
        fimg=pyfits.open(file)
        if extn>len(fimg):
            print "Requested extension [%d] does not exist" % extn
            out=no

        head=fimg[extn].header
        out=head.has_key(key)
        fimg.close()

    except:
        print "Error reading header of %s" % file
        out=no

    return out

######################################################################

def get_head(file,keys,extn=0):

    """ reads one or more header keywords from a FITS file
        using PyFITS """

    vals=[]
    check_exist(file,"r")

    try:
        fimg=pyfits.open(file)
        if extn>len(fimg):
            print "Requested extension [%d] does not exist" % extn
            vals=[""]*len(keys)
            return vals

        head=fimg[extn].header
        if type(keys)==StringType:
            key=keys
            if head.has_key(key):
                vals.append(head[key])
            else:
                print "Error reading keyword %s from %s" % (key,file)
                vals.append("")
        elif type(keys)==ListType:
            for key in keys:
                if head.has_key(key):
                    vals.append(head[key])
                else:
                    print "Error reading keyword %s from %s" % (key,file)
                    vals.append("")
        else:
            print "Bad variable type for keys"
        fimg.close()
    except:
        print "Error reading header of %s" % file
        vals=[""]*len(keys)

    if len(vals)==1:
        vals=vals[0]

    return vals

######################################################################

def update_head(file,key,value,comment=""):

    """ updates a single header keyword using pyfits
    """

    check_exist(file,"r")
    try:
        inf=pyfits.open(file,"update")
        if (comment != ""):
            inf[0].header.update(key,value,comment=comment)
        else:
            inf[0].header.update(key,value)
        inf.close()
    except:
        print "Error updating header of %s" % file
      
######################################################################

# location of parameter files
_parfile=pardir + "iqlfc.par"
t=iraf.IrafTaskFactory(taskname="iqlfc",
                       value=_parfile,function=iqlfc)

