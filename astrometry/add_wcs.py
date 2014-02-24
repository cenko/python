#!/usr/bin/env python

""" use python/pyfits to add WCS information to FITS headers
    can operate as either an executable or as an IRAF task"""


import copy, os, shutil, glob, sys, string, re, types
import pyfits
import math

from types import *

yes=0
no=1

try:
    import pyraf
    from pyraf import iraf
    useiraf=1
except:
    useiraf=0

##############################

# Pyraf parameters: Where to find them
pyrafdir="python/pyraf/"
pyrafdir_key='PYRAFPARS'

if os.environ.has_key(pyrafdir_key):
    pardir=os.environ[pyrafdir_key]
else:
    pardir=os.environ['HOME']+'/'+pyrafdir

if pardir[-1] != '/':
    pardir += '/'

instlist=['nirc','lrisr','esi','p60ccd','p60new','c40ir','cosmicR',
          'p60ir','emmi','c40ccd','wirc','panic','nirc2n','nirc2m',
          'nirc2w', 'c100ccd', 'lrisb', 'lfc', 'pairitel', 'lfosc',
          'nickel']

######################################################################

def iraffiles(files,nfiles=0):

    if type(files) is not StringType:
        print "Input filelist is not a string"
        exit(1)

    fout=[]
    fmult=files.split(",")
    for fcand in fmult:
        re1=re.search("^(.+//)?@(.+)(//.+)?$",fcand)
        re2=re.search("[\*\?]",fcand)
        if re1:
            # Using the IRAF "@file.lis" convention
            flist=re1.group(2)
            if os.path.exists(flist):
                fflist=getlines(flist)
                for fmem in fflist:
                    if re1.group(1):
                        fmem=re1.group(1)[:-2]+fmem
                    if re1.group(3):
                        fmem=fmem+re1.group(3)[2:]
                    if (fitsfile(fmem)!=""):
                        fout.append(fitsfile(fmem))
        elif re2:
            # Using UNIX wildcards
            flist=glob.glob(fcand)
            for fmem in flist:
                if (fitsfile(fmem)!=""):
                    fout.append(fitsfile(fmem))
        else:
            # Just plain filenames (?)
            if fitsfile(fcand)!="":
                fout.append(fitsfile(fcand))
            
    return fout

######################################################################

def fitsfile(file):
    outfile=""

    if file.endswith('.fits') or file.endswith('.imh'):
        if os.path.exists(file):
            outfile=file
        else:
            print "Can't find requested file %s" % file
    else:
        if os.path.exists(file+'.fits'):
            outfile=file+'.fits'
        elif os.path.exists(file):
            outfile=file
        elif os.path.exists(file+'.imh'):
            outfile=file+'.imh'
        else:
            print "Can't find requested file %s or variants" % file
        
    return outfile

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

def add_wcs(inlist, instrument, binkey="binning"):

    """add_wcs(inlist, instrument)
    adds header WCS information to FITS file with format defined by instrument
    uses pyfits
    """

    # Parse input list
    infiles=iraffiles(inlist)

    # Some defaults
    pazero=0.0

    # Definitions
    # 
    #   pa         Position angle of +Y axis, East of North, degrees
    #   pixscale   Arcsec/pixel
    #   sign       Does pixel orientation give on-sky (1), or mirror (-1)?
    #   raindeg    RA is in degrees (1) or hours (0)
    
    # do basic setup
    if (instrument == "emmi"):
        rakey="RA"
        deckey="DEC"
        rotkey=""
        pixscale=0.167
        raindeg=1
    elif (instrument == 'pairitel'):
        rakey="RA"
        deckey="DEC"
        rotkey="CROTA2"
        pixscale=2.044
        raindeg=1
    elif (instrument == 'nickel'):
        rakey="RA"
        deckey="DEC"
        rotkey="TUB"
        pixscale=0.368
        raindeg=0
    elif (instrument == 'lfosc'):
        rakey = "RA"
        deckey = "DEC"
        rotkey = ""
        pixscale = 1.083
        raindeg=0
        pazero=270
    elif (instrument == 'c100ccd'):
        rakey = "RA"
        deckey = "DEC"
        rotkey="CASSPOS"
        pixscale = 0.259
        raindeg = 0
    elif (instrument == "c40ccd"):
        rakey="RA"
        deckey="DEC"
        rotkey=""
        pixscale=0.435
        raindeg=0
    elif (instrument == "panic"):
        rakey="RA"
        deckey="DEC"
        rotkey="ROTANG"
        pixscale=0.125
        pazero=226.504
        raindeg=1    
    elif (instrument == "nirc"):
        rakey="RA"
        deckey="DEC"
        rotkey="ROTPOSN"
        pixscale=0.15
        pazero=181.25
        raindeg=1
    elif (instrument == "nirc2n"):
        rakey="RA"
        deckey="DEC"
        rotkey="ROTPOSN"
        pixscale=0.009942
        pazero=0.7
        raindeg=1    
    elif (instrument == "nirc2m"):
        rakey="RA"
        deckey="DEC"
        rotkey="ROTPOSN"
        pixscale=0.019829
        pazero=0.7
        raindeg=1    
    elif (instrument == "nirc2w"):
        rakey="RA"
        deckey="DEC"
        rotkey="ROTPOSN"
        pixscale=0.039686
        pazero=0.7
        raindeg=1    
    elif (instrument == "lfc"):
        rakey="RA"
        deckey="DEC"
        rotkey="ROTANGLE"
        pixscale=0.36
        pazero=0.0
        raindeg=0
    elif (instrument == "lrisr"):
        rakey="RA"
        deckey="DEC"
        rotkey="ROTPOSN"
        pixscale=0.135
        pazero=180.0
        raindeg=0
    elif (instrument == "lrisb"):
        rakey="RA"
        deckey="DEC"
        rotkey="ROTPOSN"
        pixscale=0.135
        pazero=180.0
        #pazero=0.0
        raindeg=0
    elif (instrument == "esi"):
    	rakey="RA"
	deckey="DEC"
	rotkey="INSTANGL"
	pixscale=0.153
        raindeg=0
    elif (instrument == "p60ccd"):
    	rakey="RA"
	deckey="DEC"
	rotkey=""
	pixscale=0.378
        raindeg=0
    elif (instrument == "p60new"):
    	rakey="RA"
	deckey="DEC"
	rotkey=""
	pixscale=0.3787
        #pazero=-2.76 # rotation for pre-Sept 2004
        pazero=0.74
        raindeg=0
    elif (instrument == "c40ir"):
    	rakey="RA"
	deckey="DEC"
	rotkey=""
	pixscale=0.6
        raindeg=0
    elif (instrument == "cosmicR"):
    	rakey="RA"
	deckey="DEC"
	rotkey=""
	pixscale=0.4
        raindeg=0
    elif (instrument == "p60ir"):
    	rakey="RA_OBS"
	deckey="DEC_OBS"
	rotkey=""
	pixscale=0.62
        raindeg=1
    elif (instrument == "wirc"):
    	rakey="RA"
	deckey="DEC"
	rotkey=""
	pixscale=0.249
        raindeg=0
    elif (instrument == "jcam0"):
    	rakey="RA"
	deckey="DEC"
	rotkey=""
	pixscale=0.37
        raindeg=0
    elif (instrument == "jcam1"):
    	rakey="RA"
	deckey="DEC"
	rotkey=""
	pixscale=0.37
        raindeg=0
    else:
        print "Unrecognized instrument: %s" % instrument
        sys.exit(1)

    # Big loop
    for file in infiles:
    
        try:
            # open image
            fimg=pyfits.open(file,"update")
            # get header
            pa=0
            hdr=fimg[0].header
            ra0=hdr.get(rakey)
            dec0=hdr.get(deckey)
            if (str(dec0).find(",") >= 0):
                dec0=re.sub(r',','',dec0)

            if (str(ra0).find(':') >= 0):
                ra=str(ra0).split(':')
                ra0=float(ra[0])+(float(ra[1])+float(ra[2])/60.0)/60.0
            if (str(dec0).find(':') >= 0):
                dec=str(dec0).split(':')
                sign=1
                if re.search("-\d",dec[0]):
                    sign=-1
                    dec[0]=abs(float(dec[0]))
                dec0=float(dec[0])+(float(dec[1])+float(dec[2])/60.0)/60.0
                dec0*=sign        
            if (len(rotkey) > 0):
                pa=hdr[rotkey]
            shape=fimg[0].data.shape
            x0=(shape[0]+1)/2
            y0=(shape[1]+1)/2
            if (not raindeg):
                ra0=15*float(ra0)
            # check for binning
            if hdr.has_key(binkey):
                # binning keyword syntax:  "nx ny"
                binning=hdr.get(binkey)
                bin2=binning.split(' ')
                if len(bin2)==2:
                    try:
                        if int(bin2[0])>1 or int(bin2[1])>1:
                            pixmult=int(bin2[0])
                            if bin2[0]!=bin2[1]:
                                print "Can't handle anisotropic binning, sorry"
                            pixscale *= pixmult
                    except:
                        print "Failure reading binning keyword %s from image %s" % \
                              (binkey,file)
            # update basic info
            hdr.update("PIXSCALE",pixscale)
            hdr.update("PIXSCAL1",pixscale)
            hdr.update("PIXSCAL2",pixscale)
            # According to the FITS WCS standard, we should not write
            # both the CD matrix and the CDELT keywords...
            hdr.update("CDELT1",pixscale/3600)
            hdr.update("CDELT2",pixscale/3600)
            hdr.update("CRVAL1",ra0)
            hdr.update("CRVAL2",dec0)
            hdr.update("CRPIX1",x0)
            hdr.update("CRPIX2",y0)
            hdr.update("CTYPE1","RA---TAN")
            hdr.update("CTYPE2","DEC--TAN")
            hdr.update("WCSDIM",2)
            hdr.update("WAT0_001","system=image")
            hdr.update("WAT1_001","wtype=tan axtype=ra")
            hdr.update("WAT2_001","wtype=tan axtype=dec")
            hdr.update("LTM1_1",1.0)
            hdr.update("LTM2_2",1.0)
        except:
            print "Error updating %s" % file

        # figure out matrix
        sign=1
        if (instrument == "nirc"):
            pa=pa-pazero
            pa=pa*3.141592654/180
            pa*=-1
            cd1_1=-sign*pixscale*math.cos(pa)/3600
            cd2_2=sign*pixscale*math.cos(pa)/3600
            cd1_2=-sign*pixscale*math.sin(pa)/3600
            cd2_1=-sign*pixscale*math.sin(pa)/3600
        elif (instrument == "pairitel"):
            pa=float(pa)*math.pi/180.0
            cd1_1=-sign*pixscale*math.cos(pa)/3600
            cd2_2=sign*pixscale*math.cos(pa)/3600
            cd1_2=-sign*pixscale*math.sin(pa)/3600
            cd2_1=-sign*pixscale*math.sin(pa)/3600
        elif (instrument == "lfosc"):
            pa=pa-pazero
            pa=pa*math.pi/180.0
            cd1_1=-sign*pixscale*math.cos(pa)/3600
            cd2_2=sign*pixscale*math.cos(pa)/3600
            cd1_2=-sign*pixscale*math.sin(pa)/3600
            cd2_1=sign*pixscale*math.sin(pa)/3600
        elif (instrument == "nickel"):
            pa=pa*math.pi/180.0
            cd1_1=-sign*pixscale*math.cos(pa)/3600
            cd2_2=-sign*pixscale*math.cos(pa)/3600
            cd1_2=-sign*pixscale*math.sin(pa)/3600
            cd2_1=sign*pixscale*math.sin(pa)/3600
        elif (instrument.startswith("nirc2")):
            pa=pa-pazero
            pa=pa*3.141592654/180
            pa*=-1
            cd1_1=-sign*pixscale*math.cos(pa)/3600
            cd2_2=sign*pixscale*math.cos(pa)/3600
            cd1_2=-sign*pixscale*math.sin(pa)/3600
            cd2_1=-sign*pixscale*math.sin(pa)/3600
        elif (instrument == "panic"):
            pa=pa-pazero
            pa=pa*3.141592654/180
            cd1_1=sign*pixscale*math.cos(pa)/3600
            cd2_2=sign*pixscale*math.cos(pa)/3600
            cd1_2=sign*pixscale*math.sin(pa)/3600
            cd2_1=-sign*pixscale*math.sin(pa)/3600
        elif (instrument == "lfc"):
            pa=pa-pazero
            pa=pa*3.141592654/180
            pa*=-1
            cd1_1=-sign*pixscale*math.cos(pa)/3600
            cd2_2=sign*pixscale*math.cos(pa)/3600
            cd1_2=-sign*pixscale*math.sin(pa)/3600
            cd2_1=-sign*pixscale*math.sin(pa)/3600
        elif (instrument == "lrisr"):
            pa=pa-pazero
            pa=pa*3.141592654/180
            pa*=-1
            cd1_1=-sign*pixscale*math.cos(pa)/3600
            cd2_2=-sign*pixscale*math.cos(pa)/3600
            cd1_2=sign*pixscale*math.sin(pa)/3600
            cd2_1=-sign*pixscale*math.sin(pa)/3600
        elif (instrument == "lrisb"):
            pa=pa-pazero
            pa=pa*3.141592654/180
            pa*=-1
            cd1_1=-sign*pixscale*math.cos(pa)/3600
            cd2_2=-sign*pixscale*math.cos(pa)/3600
            cd1_2=sign*pixscale*math.sin(pa)/3600
            cd2_1=-sign*pixscale*math.sin(pa)/3600
        elif (instrument == "c100ccd"):
            pa=pa-pazero
            pa=pa*3.141592654/180
            cd1_1= sign*pixscale*math.cos(pa)/3600
            cd2_2=sign*pixscale*math.cos(pa)/3600
            cd1_2= sign*pixscale*math.sin(pa)/3600
            cd2_1= sign*pixscale*math.sin(pa)/3600
        elif (instrument == "c40ccd"):
            cd1_1=-sign*pixscale/3600
            cd2_2=-sign*pixscale/3600
            cd1_2=-0*sign*pixscale/3600
            cd2_1=-0*sign*pixscale/3600
        elif (instrument == "p60ccd"):
            cd1_1=sign*pixscale/3600
            cd2_2=sign*pixscale/3600
            cd1_2=-0*sign*pixscale/3600
            cd2_1=-0*sign*pixscale/3600
        elif (instrument == "p60new"):
            pa=pa-pazero
            pa *= 3.141592654/180
            cd1_1=-sign*pixscale*math.sin(pa)/3600
            cd2_2=+sign*pixscale*math.sin(pa)/3600
            cd1_2=+sign*pixscale*math.cos(pa)/3600
            cd2_1=+sign*pixscale*math.cos(pa)/3600
        elif (instrument == "c40ir"):
            cd1_1=-sign*pixscale/3600
            cd2_2=sign*pixscale/3600
            cd1_2=-0*sign*pixscale/3600
            cd2_1=-0*sign*pixscale/3600
        elif (instrument == "cosmicR"):
            cd1_1=-sign*pixscale/3600
            cd2_2=-sign*pixscale/3600
            cd1_2=-0*sign*pixscale/3600
            cd2_1=-0*sign*pixscale/3600
        elif (instrument == "esi"):
            pa=pa*3.141592654/180
            cd1_1=-sign*pixscale*math.cos(pa)/3600
            cd2_2=-sign*pixscale*math.cos(pa)/3600
            cd1_2=sign*pixscale*math.sin(pa)/3600
            cd2_1=-sign*pixscale*math.sin(pa)/3600
        elif (instrument == "p60ir"):
            cd1_1=-sign*pixscale/3600
            cd2_2=sign*pixscale/3600
            cd1_2=-0*sign*pixscale/3600
            cd2_1=-0*sign*pixscale/3600
        elif (instrument == "wirc"):
            cd1_1=-1*sign*pixscale/3600
            cd2_2=1*sign*pixscale/3600
            cd1_2=0*sign*pixscale/3600
            cd2_1=0*sign*pixscale/3600
        elif (instrument == "jcam0"):
            cd1_1=0*sign*pixscale/3600
            cd2_2=0*sign*pixscale/3600
            cd1_2=1*sign*pixscale/3600
            cd2_1=-1*sign*pixscale/3600
        elif (instrument == "jcam1"):
            cd1_1=0*sign*pixscale/3600
            cd2_2=0*sign*pixscale/3600
            cd1_2=1*sign*pixscale/3600
            cd2_1=-1*sign*pixscale/3600
        elif (instrument == "emmi"):
            cd1_1=sign*pixscale/3600
            cd2_2=-sign*pixscale/3600
            cd1_2=-0*sign*pixscale/3600
            cd2_1=-0*sign*pixscale/3600


        try:
            # update matrix
            hdr.update("CD1_1",cd1_1)
            hdr.update("CD1_2",cd1_2)
            hdr.update("CD2_1",cd2_1)
            hdr.update("CD2_2",cd2_2)
            fimg.flush()
            fimg.close()
            print "Image %s updated" % file
        except:
            print "Error updating %s" % file

######################################################################

def usage():
    (xdir,xname)=os.path.split(sys.argv[0])
    print "Usage:  %s [-i instrument] <filename(s)>" % xname
    print "\tInstrument is one of: %s" % " | ".join(instlist)

######################################################################

def main():
    
    files=[]
    instrument=""

    if (len(sys.argv)==1):
        usage()
        sys.exit(1)

    i=1
    while (i<len(sys.argv)):
        arg=sys.argv[i]
        isarg=0
        if (arg.find("-i") > -1):
            instrument=sys.argv[i+1]
            i+=1
            isarg=1
        if (arg.find("-h") > -1):
            usage()
            sys.exit(1)
        if (not isarg):
            if (len(files) > 0):
                files+=",%s" % arg
            else:
                files=arg
        i+=1

    filenames=files.split(",")

    for file in filenames:
        add_wcs(file,instrument)

######################################################################
# Running as executable
if __name__=='__main__':
    main()
else:
    if useiraf:
        _parfile=pardir + "add_wcs.par"
        t=iraf.IrafTaskFactory(taskname="add_wcs",value=_parfile,
                               function=add_wcs)

######################################################################
