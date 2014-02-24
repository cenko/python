
#############################################################
# $Id$
#############################################################

import pyraf
from pyraf import iraf
import copy, os, shutil, glob, sys, string, re, math, operator
import pyfits
from types import *

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
iraf.photcal()
iraf.images()
iraf.imcoords()

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

######################################################################

def calibphot(inlist, fwhm,
              rdnoise='2.6', gain='2.5',
              interact=yes, config='',
              objkey="OBJECT", airkey="AIRMASS",
              filtkey='FILTER', skykey='SKYBKG', phosfx="pho",
              possfx="pos", 
              clobber=globclob, verbose=globver):

    """ perform aperture photometry on a standard field to fit the
    transformation equations.  these transformation equations will
    then be applied to comparison stars' aperture photometries
    outputted by psfphot.py, in order to find the true magnitudes of
    the comparison stars.  this will then provide and apply the
    corrections to the magnitude of a single target star (SN?). """

    # Defaults / constants
    alpha=list('abcdefghijklmnopqrstuvwxyz')
    psfmult=5     #standard factor (multiplied by fwhm to get psfradius)
    imgprop={}
    imglist={}
    allobjs={}
    allfilt={}

    # Necessary package
    iraf.imutil()

    # Parse gain/readnoise inputs
    gainkey=None
    re1=re.search('^!(.+)',gain)
    re2=re.search('^[\d\.]+',gain)
    if re1:
        gainkey=re1.group(1)
    elif re2:
        try:
            gainval=float(gain)
        except:
            print "Error converting gain"
            return
    else:
        print "Error finding gain"
        return
    
    rdnskey=None
    re1=re.search('^!(.+)',rdnoise)
    re2=re.search('^[\d\.]+',rdnoise)
    if re1:
        rdnskey=re1.group(1)
    elif re2:
        try:
            readval=float(rdnoise)
        except:
            print "Error converting readnoise"
            return
    else:
        print "Error finding readnoise"
        return
    
    # Parse inputs
    infiles=iraffiles(inlist)

    #first make sure to add back in background of sky
    iraf.iqsubsky(inlist, sub=no, skykey=skykey)

    # Process each file in turn
    for image in infiles:

        # Check that the image is there
        check_exist(image,"r")

        # Grab all the useful keywords
        [sky,objname,filter,airmass] = \
              get_head(image,[skykey,objkey,filtkey,airkey])

        # Get gain & readnoise, if necessary
        if gainkey:
            gainval=get_head(image,gainkey)
        if rdnskey:
            readval=get_head(image,rdnskey)

        # Calculate theoretical sigma of sky
        if len(str(sky))>0:
            sigma= (((sky * gainval) + readval**2)**.5) / gainval
        else:
            print "Failed to retrieve sky value from image header"
            return

        # check standard field name
        if len(objname)==0:
            print "Failed to retrieve object name from image header"
            return

        # Generate posfile and phofile names
        posfile=objname+'.'+possfx
        phofile=objname+'.'+phosfx
        
        # check filter name
        if len(filter)==0:
            print "Failed to retrieve filter name from image header"
            return

        # Check airmass value
        if len(airmass)>0:
            try:
                airval=float(airmass)
            except:
                print "Failed to convert airmass value to float"
                return
        else:
            print "Failed to retrieve airmass value from image header"
            return

        # Save keyword values for this image
        imgprop[image]=[gainval,readval,sky,objname,filter,airmass]

        # Add dictionary entries for the current image
        listkey=objname+'-'+filter
        if imglist.has_key(listkey):
            imglist[listkey].append(image)
        else:
            imglist[listkey]=[image]

        # Keep track of all filters/object fields
        allobjs[objname]=1
        allfilt[filter]=1

        #add in extra range of aperture radii, leave in background?
        fwhmpsf=float(fwhm)
        ap1 = fwhmpsf
        ap15 = 1.5 * fwhmpsf
        ap2 = 2 * fwhmpsf
        ap25 = 2.5 * fwhmpsf
        ap3 = 3 * fwhmpsf
        innersky = 3.5 * fwhmpsf

        #do photometry on all stars in field
        check_exist(image+".mag.std.1", "w", clobber=yes)
        iraf.phot(image, output = image+".mag.std.1", coords=posfile,
                  sigma=sigma, datamin = sky - 3.5*sigma, apertures = ap25,
                  dannulus="10", annulus=innersky, verify="no")
    
    ##############################

    #make imsets file
    imsetsfile="stdfield.dat"
    imsetstuff = open(imsetsfile,"w")
    for objname in allobjs.keys():
        maxnum=0
        for key in imglist.keys():
            if re.search(objname,key):
                maxnum=max(maxnum,len(imglist[key]))
        for i in range(maxnum):
            imset=objname+alpha[i]+" :"
            for key in imglist.keys():
                if re.search(objname,key):
                    if i<len(imglist[key]):
                        imset+=" "+imglist[key][i]
            imsetstuff.write(imset+"\n")
    imsetstuff.close()

    #make observations file (specify aperture number?)
    obsfile = "standobs"
    check_exist(obsfile, "w", clobber=yes)
    iraf.mknobsfile(photfiles = "*.mag.std.1",
                    idfilters = " ".join(allfilt.keys()),
                    imsets = imsetsfile, observations = obsfile,
                    tolerance="500")
    
    #rename stars in catalog file so they agree, print good stars to new file
    catpholines = getlines(phofile)
    s=1
    newcatpholines = [catpholines[0]]
    lastline = catpholines[-1]
    lastlinedata = lastline.split()
    digits = len(lastlinedata[0])+1
    for line in catpholines:
        if re.search("vary?", line):
            continue
        ess=str(s)
        idlength=len(objname)+len(ess)+1
        numspace= digits - idlength + 1
        newline = " "*numspace+objname+"-"+ess+line[digits:]
        s=s+1
        if not re.search("99.999",newline):
            newcatpholines.append(newline)
    newcatdata = "\n".join(newcatpholines)
    newphofile = phofile+".2"
    newphostuff = open(phofile+".2", "w")
    newphostuff.write(newcatdata)
    newphostuff.close()
    
    #rename stars in standobs to be field-# instead of fielda-#, fieldb-#, etc
    #what is field?
    obsdatalines=getlines(obsfile)
    for i in range(len(obsdatalines)):
        line=obsdatalines[i]
        for objname in allobjs.keys():
            if re.search(objname,line):
                index=line.index(objname)
                lline=list(line)
                lline[index:index+len(objname)+1]=list(" "+objname)
                line="".join(lline)
                obsdatalines[i]=line
                break
    newobsdata = "\n".join(obsdatalines)
    fixobsfile = "fix"+obsfile
    fixobsstuff = open(fixobsfile, "w")
    fixobsstuff.write(newobsdata)
    fixobsstuff.close()
    
    #if no config file exists, make one
    #first, make a file that has stetson catalog formats in it
    if not config:
        newconfig="fcat.dat"
        check_exist(newconfig, "w", clobber=yes)
        catformat=open(newconfig, "w")
        catinfo = ["#Declare the catalog variables", " ", "catalog", " ",
                   " B\t2", "Berror\t3", "V\t6", "Verror\t7",
                   "R\t10", "Rerror\t11", "I\t14", "Ierror\t15"]
        catdata = "\n".join(catinfo)
        catformat.write(catdata)
        catformat.close()
        config="cat.cfg"

        transfile="ftrans.dat"
        check_exist(transfile, "w", clobber=yes)
        transstuff=open(transfile, "w")
        transinfo = ["transformation", " ", "fit v1 = 0, v2 = 0.17, v3= 0", "fit b1 = 0, b2 = 0.35, b3= 0", "fit r1 = 0, r2 = 0.08, r3= 0", " ", "const v4 = 0", "const b4 = 0", "const r4 = 0", " ", "Vfit : V = mV + v1 + v2*XV + v3*(mB-mV) + v4*(mB-mV)*XV", "Bfit : B = mB + b1 + b2*XB + b3*(mB-mV) + b4*(mB-mV)*XB", "Rfit : R = mR + r1 + r2*XR + r3*(mV-mR) + r4*(mV-mR)*XR"]
        transdata = "\n".join(transinfo)+"\n"
        transstuff.write(transdata)
        transstuff.close()

        check_exist(config, "w", clobber=yes)
        iraf.mkconfig(config, catalog = newconfig,
                      observations = "f"+obsfile+".dat",
                      transform=transfile, check=no, edit=no)

    #fit parameters of transformation equations (specify aperture in ansfile?)
    ansfile = objname+".ans.1"
    check_exist(ansfile, "w", clobber=yes)
    iraf.fitparams(observations = fixobsfile, catalogs = newphofile,
                   config=config, parameters = ansfile, interactive = interact)
    
    #apply to standard stars to check to make sure fits are okay?
    calibfile = objname+".calib.1"
    check_exist(calibfile, "w", clobber=yes)
    iraf.evalfit(observations = fixobsfile, config = config,
                 parameters = ansfile, calib = calibfile)
        

######################################################################

# location of parameter files

_parfile=pardir + "calibphot.par"
t=iraf.IrafTaskFactory(taskname="calibphot",
                       value=_parfile,function=calibphot)


