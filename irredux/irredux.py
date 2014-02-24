
"""IR reduction
DLK 2002-11-11"""


import pyraf
from pyraf import iraf
import copy, os, shutil, glob, sys, string, re, types, math
import pyfits
#import numarray
import numpy

yes=iraf.yes
no=iraf.no
INDEF=iraf.INDEF
hedit=iraf.hedit
imcombine=iraf.imcombine
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

# stsdas.toolbox.imgtools (for imcalc)
iraf.stsdas()
iraf.toolbox()
iraf.imgtools()
# images.immatch (for geomap/geotran)
iraf.images()
iraf.immatch()
# imred.ccdred and imred.crutil
iraf.imred()
iraf.ccdred()
iraf.crutil()
# digiphot and digiphot.apphot
iraf.digiphot()
iraf.apphot()

globclob=yes
globver=yes
globlog=""

######################################################################

def make_bpm(inlist, outname, bpmkey, modlist, expkey, coaddkey, combine, subsect, lthreshold, hthreshold, clean, clobber=globclob, verbose=globver):
    """Bad Pixel Mask Creation"""

    check_exist(inlist,"r",clobber)
    
    # check for FITS or PL extension
    extn=os.path.splitext(outname)[1]
    if (extn != ".pl" and extn != ".fits" and extn != ".FITS"):
        outname+=".pl"
    check_exist(outname,"w",clobber)

    # get scale for combination
    # need a file for scales
    check_exist(outname+".scales","w",yes)
    infile=open(inlist,'r')
    sclfile=open(outname+".scales",'w')
    files=infile.readlines()
    infile.close()
    n=len(files)
    for i in range(0,n):
        file=files[i]
        file=file.strip()
        check_exist(file,"r",clobber)
        # get exposure & coadd
        imgets(file,expkey)
        exp=float(imgets.value)
        imgets(file,coaddkey)
        coadd=float(imgets.value)
        scale=1.0/(coadd*exp)
        sclfile.write("%f\n" % scale)
    sclfile.close()

    s1="comb.fits"
    s2="sigma.fits"
    s3="rej.fits"
    check_exist(s1,"w",yes)
    check_exist(s2,"w",yes)
    check_exist(s3,"w",yes)
    check_exist("sigma_norm.fits","w",yes)
    check_exist("sigma_temp.fits","w",yes)

    # setup imcombine
    iraf.unlearn(imcombine)
    i=imcombine
    i.sigma=s2
    i.combine=combine
    i.reject="none"
    i.scale="@" + outname + ".scales"
    # combine
    s=imcombine("@" + inlist, s1,Stdout=1)
    if (verbose):
        print "\n".join(s)

    # get stats
    iraf.iterstat(s2+subsect,nsigrej=5,maxiter=10,prin=no,verbose=verbose)
    mn=float(iraf.iterstat.mean)
    sig=float(iraf.iterstat.sigma)
    hthresh=mn+hthreshold*sig
    lthresh=mn-lthreshold*sig

    # produce normalized output
    iraf.imarith(s2,"-",mn,"sigma_temp.fits")
    iraf.imarith("sigma_temp.fits","/",sig,"sigma_norm.fits")
    os.remove("sigma_temp.fits")

    # find rejected bit
    s="if im1 .ge. %f .or. im1 .le. %f then 1 else 0" % (hthresh,lthresh)
    iraf.imcalc(s2,s3,s,verbose=verbose)
    iraf.imcopy(s3,outname,verbose=verbose)

    hedit(outname,"MEAN",mn,add=yes,verify=no,show=verbose,update=yes)
    hedit(outname,"SIGMA",sig,add=yes,verify=no,show=verbose,update=yes)
    hedit(outname,"LTHRESH",lthreshold,add=yes,verify=no,show=verbose,update=yes)
    hedit(outname,"HTRESH",hthreshold,add=yes,verify=no,show=verbose,update=yes)

    # update headers
    if (len(bpmkey) > 0):
        try:
            check_exist(modlist,"r",clobber)
            hedit("@"+modlist,bpmkey,outname,add=yes,verify=no,show=verbose,update=yes)
        except:
            s2=s2
    
    os.remove(outname+".scales")
    if (clean):
        os.remove(s1)
        os.remove(s2)
        os.remove(s3)
        os.remove("sigma_norm.fits")

    print "BPM File %s created" % outname

######################################################################
def make_dark(objfile,darkfile,outprefix,expkey,coaddkey,combine,reject,lthreshold,hthreshold,nlow,nhigh,nkeep,mclip,lsigma,hsigma,rdnoise,gain,snoise,sigscale,pclip,grow,clobber=globclob,verbose=globver):
    """Dark creation"""

    check_exist(objfile,"r",clobber)
    check_exist(darkfile,"r",clobber)

    outlist=objfile + ".dark"
    check_exist(outlist,"w",clobber)
    # read in list of object files
    # also get tint & coadd value for each
    params={}
    files=getlines(objfile)
    if (verbose):
        print "File\t\ttin\ttcoadd"
    for i in range(0,len(files)):
        s=files[i]
        imgets(s,expkey)
        tint=(int(10*float(imgets.value)))/10.0
        imgets(s,coaddkey)
        coadd=int(imgets.value)
        s2=",%s" % s
        s3="%s_%s" % (tint,coadd)
        try:
            params[s3]=params[s3]+s2
        except:
            params[s3]=s            
        if (verbose):
            print "%s\t%3.1f\t%d" % (s,tint,coadd)    

    dparams={}
    dfiles=getlines(darkfile)
    for i in range(0,len(dfiles)):
        s=dfiles[i]
        imgets(s,expkey)
        tint=(int(10*float(imgets.value)))/10.0
        imgets(s,coaddkey)
        coadd=int(imgets.value)
        s2=",%s" % s
        s3="%s_%s" % (tint,coadd)
        try:
            dparams[s3]=dparams[s3]+s2
        except:
            dparams[s3]=s           

    if (verbose):
        print "Read in %d object files and %d dark files" % (len(files),len(dfiles))

    outfile=open(outlist,"w")
    for param in params.keys():
        try:
            gooddfiles=dparams.get(param).split(",")
        except:
            print "No appropriate dark files found for %s" % param
            sys.exit(2)
            
        (tint,coadd)=param.split("_")
        if (verbose):
            print "For tint=%s and coadd=%s there are %d dark frames to use" % (tint,coadd,len(gooddfiles))
        if (len(gooddfiles) > 0):
            outname=outprefix + "_" + param + ".fits"
            check_exist(outname,"w",clobber)
            iraf.unlearn(imcombine)
            imcombine.combine=combine
            imcombine.reject=reject
            imcombine.scale="none"
            imcombine.zero="none"
            imcombine.weight="none"
            imcombine.lthreshold=lthreshold
            imcombine.hthreshold=hthreshold
            imcombine.nlow=nlow
            imcombine.nhigh=nhigh
            imcombine.mclip=mclip
            imcombine.lsigma=lsigma
            imcombine.rdnoise=rdnoise
            imcombine.gain=gain
            imcombine.snoise=snoise
            imcombine.sigscale=sigscale
            imcombine.pclip=pclip
            imcombine.grow=grow
            imcombine.rejmask=""
            imcombine.bpmasks=""
            imcombine.sigma=""
            sout=imcombine(",".join(gooddfiles),outname,Stdout=1)
            if (verbose):
                print "\n".join(sout)
            hedit(outname,"title","Dark Frame" + param,add=yes,verify=no,show=no)
            for s in params[param].split(","):
                outfile.write("%s %s\n" % (s,outname))
    outfile.close()
    if (verbose):
        print "List of files in %s" % outlist
        print "Dark images created!"

######################################################################
def dark_subtr(inlist,outpfx="d",clobber=globclob,verbose=globver):
    """ dark subtraction"""

    check_exist(inlist,"r")

    outpfx=outpfx.strip()
    if (len(outpfx) == 0):
        print "Warning: this will overwrite old files"
    outlist=inlist+ ".sub"
    check_exist(outlist,"w",clobber)    

    files=getlines(inlist)
    outfile=open(outlist,"w")
    for s in files:
        s2=s.split()
        if (len(s2)==2):
            (image,dark)=s.split()
        else:
            (image,dark,junk)=s.split()
        check_exist(image,"r")
        check_exist(dark,"r")
        if (verbose):
            print "Processing %s with %s -> %s" % (image,dark,outpfx + image)
        # get processing status
        imgets(image,"DARKSUB",Stdout=1,Stderr=1)
        if (int(imgets.value) != 0):
            print "Warning: %s has already been dark subtracted" % image
        else:
            if (len(outpfx) != 0):
                check_exist(outpfx+image,"w",clobber)
            iraf.imarith(image,"-",dark,outpfx+image)
            update_head(outpfx+image,"DARKSUB",1,"Has dark-subtraction occured?")
            update_head(outpfx+image,"DARKNAME",dark,"Name of dark frame")
            update_head(outpfx+image,"NAME0",image,"Original name of file")
            outfile.write("%s%s\n" % (outpfx,image))
    outfile.close()
    if (verbose):
        print "List of output written to %s" % outlist
        print "Dark subtraction finished!"

######################################################################
def make_flat(inlist,window,donorm,prefix,filtkey,expkey,coaddkey,bpmkey,seqkey,combine,reject,skymethod,scaletype,zero,statsec,maskscale,maskcomb,minseq,maxseq,lthreshold,hthreshold,nlow,nhigh,nkeep,mclip,lsigma,hsigma,rdnoise,gain,snoise,sigscale,pclip,grow,clobber=globclob,verbose=globver):
    """ flatfield creation"""

    check_exist(inlist,"r")
    outlist=inlist + ".flat"
    check_exist(outlist,"w",clobber)
    outfile=open(outlist,"w")

    # get input files
    files=getlines(inlist)
    if (verbose):
        print "File\t\tFilt\t\tTint\t\tCoadd"
    params={}
    paramstouse={}
    medians={}
    invmed={}
    inseq={}
    for file in files:
        check_exist(file,"r")
        s=file
        imgets(s,expkey)
        tint=(int(10*float(imgets.value)))/10.0
        imgets(s,coaddkey)
        coadd=int(imgets.value)
        imgets(s,filtkey)
        filt=imgets.value
        if (len(seqkey) > 0):
            # check sequence
            imgets(s,seqkey)
            seq=float(imgets.value)
            if (seq >= minseq and seq <= maxseq):
                inseq[file]=yes
            else:
                inseq[file]=no
        else:
            inseq[file]=no                        
        s2="%s_%s_%s" % (tint,coadd,filt)
        s3=",%s" % s
        try:
            params[s2]=params[s2]+s3
        except:
            params[s2]=s
        try:
            paramstouse[s2]=paramstouse[s2]+s3
        except:
            paramstouse[s2]=s           
        if (verbose):
            print "%s\t%s\t%s\t%s" % (file,filt,tint,coadd)
        if (verbose):
            print "Computing initial sky value for %s" % file

        if (maskcomb):
            # need to define BPM for imcombine later
            # make sure mask exists
            imgets(file,bpmkey,Stdout=1,Stderr=1)
            if (len(imgets.value) == 0):
                print "No mask defined for image %s" % file
                sys.exit(1)
            if (bpmkey != "BPM"):
                # if BPM isn't defined through BPM keyword, 
		# imcombine won't understand.
		# so define it
                update_head(file,"BPM",iraf.imgets.value,"Current BPM")
            medians[file]=1.0*float(iraf.mimstatistics(images=file,imasks=imgets.value,fields="midpt",nclip=10,lsigma=5,usigma=5,format=no,Stdout=1)[0])
                           
        else:
            iraf.iterstat(file,nsigrej=5,maxiter=10,prin=no,verbose=no,lower=INDEF,upper=INDEF)
            medians[file]=float(iraf.iterstat.median)


        # then store in header
        invmed[file]=1.0/medians[file]
        update_head(file,"SKYMED",medians[file],"Value of SKY median")
        update_head(file,"MEDSCALE",invmed[file],"Inverse of SKY median, for scaling")

    if (verbose):
        print "Read %d files" % len(files)

    # go through combinations
    iparam=0
    for param in params.keys():
        iparam+=1
        (tint,coadd,filt)=param.split("_")
        goodfiles=params[param].split(",")
        filestouse=paramstouse[param].split(",")
        if (verbose):
            print "For filter=%s tint=%s coadd=%s there are %d files" % (filt,tint,coadd,len(filestouse))
        # begin making files
        if (window == 0 or window > len(filestouse)):
            window_use=len(filestouse)
        else:
            window_use=window
        # make sure it's odd
        if (window_use % 2 == 0):
            window_use+=1
        flat_num=0
        for l in range(0,len(goodfiles)):
            flat_num+=1
            # use from l+(-nleft,...,0,...,nright)
            nleft=int(window_use/2)
            nright=int(window_use/2)
            # handle edges
            #if ((l-nleft) < 0):
                #nleft=l
                #nright=window_use-l-1
            #if ((l+nright) >= len(goodfiles)):
                #nleft+=((nright+l)-len(goodfiles)+1)
                #nright=(len(goodfiles)-l)-1
            #if (l==nleft):
            #    nleft-=1
            #i1=l-nleft
            #i2=l+nright+1
            new=[]
            for m in range(l+1,l+nright+1):
                 new.append(goodfiles[m%len(goodfiles)])
            for m in range(l-1,l-nleft-1,-1):
                 new.append(goodfiles[m%len(goodfiles)])
            files=",".join(new)
            # We've selected the files, now use them  
            #files=",".join(goodfiles[i1:i2])

            # now make flats
            flatname="%s_%s_%s.fits" % (prefix,iparam,flat_num)
            check_exist(flatname,"w",clobber)
            #combine
            imcombine.combine=combine
            imcombine.reject=reject
            imcombine.bpmasks=""
            imcombine.expmasks=""
            imcombine.sigma=""
            imcombine.scale="none"
            imcombine.weight="none"
            imcombine.scale="!MEDSCALE"
            imcombine.weight="!MEDSCALE"
            imcombine.zero=zero
            imcombine.lthreshold=lthreshold
            imcombine.hthreshold=hthreshold
            imcombine.nlow=nlow
            imcombine.nhigh=nhigh
            imcombine.nkeep=nkeep
            imcombine.mclip=mclip
            imcombine.lsigma=lsigma
            imcombine.hsigma=hsigma
            imcombine.rdnoise=rdnoise
            imcombine.gain=gain
            imcombine.snoise=snoise
            imcombine.sigscale=sigscale
            imcombine.pclip=pclip
            imcombine.grow=grow
            imcombine.maskval=0
            if (maskcomb):
                imcombine.masktype="goodvalue"			                        
            else:
                imcombine.masktype="none"
            check_exist(flatname+".exp.pl","w",clobber)
            check_exist(flatname+".rej.pl","w",clobber)
            # need to see if rejecting too many points
            if (reject == "minmax"):
                if (imcombine.nhigh+imcombine.nlow > window_use):
                    while (imcombine.nhigh+imcombine.nlow > window_use):
                        imcombine.nlow-=1
                        imcombine.nhigh-=1
                    if (verbose):
                        print "Changed nhigh=%d nlow=%d" % (imcombine.nhigh,imcombine.nlow)
            s=imcombine(input=files,output=flatname,rejmask="",nrejmask=flatname+".rej.pl",expmasks=flatname+".exp",Stdout=1)
            if (verbose):
                print "\n".join(s)
            check_exist(flatname,"r")
            # divide by sky, to keep consistent
            med=medians[goodfiles[l]]
            if (donorm):
                iraf.imarith(flatname,"/",med,flatname+ "2.fits",verbose=verbose)
                os.remove(flatname)
                os.rename(flatname+"2.fits",flatname)
            update_head(flatname,"EXPMAP",flatname+".exp.pl","Location of exposure map")

            if (verbose):
                print "Images combined"
                print "Computing stats with %s" % skymethod
            scaleuse=scaletype
            if (skymethod == "imstat"):
                # the old way
                if (not maskscale):
                    s=iraf.imstat(flatname+statsec,fields=scaleuse,format=no,Stdout=1)
                else:
                    print "Cannot use imstat with maskscale=yes"
                    sys.exit(1)
                scaleval=float(s)
            if (skymethod == "iterstat"):
                # the new way
                if (not maskscale):
                    iraf.iterstat(image=flatname+statsec,nsigrej=5,maxiter=10,prin=no,verbose=no,lower=INDEF,upper=INDEF)
                    if (scaleuse == "mean"):
                        scaleuseval=float(iraf.iterstat.mean)
                    if (scaleuse == "median"):
                        scaleval=float(iraf.iterstat.median)
                    if (scaleuse == "mode"):
                        scaleval=float(iraf.iterstat.valmode)

                else:
                    imgets(goodfiles[l],bpmkey,Stdout=1,Stderr=1)
                    if (len(imgets.value) > 0):
                        mask=imgets.value
                    else:
                        mask=""

                    if (len(mask) == 0 or mask == 0):
                        print "Mask does not exist for image %s" % goodfiles[l]
                        sys.exit(1)
                    check_exist("_"+mask,"w",yes)
                    if (scaleuse == "median"):
                        scaleuse="midpt"
                    s=iraf.mimstatistics(flatname+statsec,imask=mask,fields=scaleuse,lower=INDEF,upper=INDEF,lsigma=3,usigma=3,format=no,Stdout=1)
                    if (scaleuse == "midpt"):
                        scaleuse="median"

                    scaleval=float(s[0])
            if (skymethod == "sky"):
                # another new way
                # this automatically uses the mask
                if (scaleuse=="median"):
                    print "Warning: sky doesn't compute median - using mean"
                    scaleuse="mean"
                    
                iraf.imgets(flatname,bpmkey,Stdout=1,Stderr=1)
                if (len(imgets.value) > 0 and maskscale):
                    mask=imgets.value
                else:
                    mask=""		
                    iraf.sky.lower=-200
                    iraf.sky.upper=65000
                    iraf.sky.subsky=no
                    iraf.sky.width=8
                    iraf.sky.stat=scaleuse
                    iraf.sky.skyname=""
                    iraf.sky(input=flatname,verbose=verbose)
                    scaleval=float(iraf.sky.skyvalue)

            if (skymethod == "none"):
                scaleval=1.0
            if (donorm):
                iraf.imarith(operand1=flatname,op="/",operand2=scaleval,result=flatname,title="Flatfield for "+goodfiles[l],verbose=verbose)
            else:
                med=medians[goodfiles[l]]
                iraf.imarith(operand1=flatname,op="/",operand2=scaleval/med,result=flatname,title="Flatfield for "+goodfiles[l],verbose=verbose)
                
            outfile.write("%s %s %s\n" % (goodfiles[l],flatname,medians[goodfiles[l]]))
            update_head(flatname,"FORIMG",goodfiles[l],"File flatfield is for")
            update_head(flatname,"SKYBKG",scaleval,"Value of SKY background")
    outfile.close()
    if (verbose):
        print "List of output written to %s" % outlist
        print "Flat creation done!"
    
######################################################################
def flat_divide(inlist, outpfx="f", clobber=globclob, verbose=globver):
    """ flatfield division"""

    check_exist(inlist,"r")
    outlist=inlist+".div"
    check_exist(outlist,"w",clobber)
    outfile=open(outlist,"w")
    if (len(outpfx) == 0):
        print "Warning: this will overwrite old files"

    data=getlines(inlist)
    for datum in data:
        (image,flat,sky)=datum.split()
        check_exist(image,"r")
        check_exist(flat,"r")
        if (verbose):
            print "Processing %s with %s sky=%s" % (image,flat,sky)
        # get processing status
	imgets(image,"FLATDIV",Stdout=1,Stderr=1)
        if (int(imgets.value) != 0):
            print "Warning: %s has already been flatfielded!" % image
        else:
            if (len(outpfx) != 0):
                check_exist(outpfx+image,"w",clobber)
            iraf.imarith(operand1=image,op="/",operand2=flat,result=outpfx+image,verbose=verbose)
	    iraf.imarith(operand1=outpfx+image,op="-",operand2=sky,result=outpfx+image,verbose=verbose)
            
	    update_head(outpfx+image,"FLATDIV",1,"Has flat-field division occured?")
	    update_head(outpfx+image,"FLATNAME",flat,"Name of flatfield")
	    update_head(outpfx+image,"SKYBKG",sky,"Value of SKY background")
	    update_head(outpfx+image,"NAME1",image,"Previous name of file")
            outfile.write("%s\n" % (outpfx+image))
    outfile.close()
    if (verbose):
        print "List of output written to %s" % outlist
        print "Flatfielding finished!"
        
######################################################################
def irfixpix(inlist, outpfx="b", fake=no, bpmkey="BPM",clobber=globclob,verbose=globver):
    check_exist(inlist,"r")
    """ bad-pixel fixing"""

    batchnum=20

    outlist=inlist+".fix"
    check_exist(outlist,"w",clobber)
    outfile=open(outlist,"w")
    if (len(outpfx) == 0):
        print "Warning: this will overwrite old files"

    files=getlines(inlist)
    # go through all files, prepare lists for crfix
    inpfiles=""
    maskfiles=""
    for file in files:
        check_exist(file,"r")
        # check to see if has already been done
        imgets(file,"BPMFIX",Stdout=1,Stderr=1)
        if (len(imgets.value) > 0 and int(imgets.value) > 0):
            print "Warning: %s has already been fixed.  Skipping." % file
            continue
        imgets(file,bpmkey,Stdout=1,Stderr=1)
        bpmfile=imgets.value
        if (len(bpmfile) <= 1):
            print "No BPM file for image %s" % file
            sys.exit(1)
        check_exist(bpmfile,"r")
        if (len(outpfx) > 0):
            check_exist(outpfx+file,"w",clobber)
            iraf.imcopy(file,outpfx+file,verbose=verbose)
        if (len(inpfiles) > 0):
            inpfiles+=","
        if (len(maskfiles) > 0):
            maskfiles+=","
        inpfiles+=outpfx+file
        maskfiles+=bpmfile
        outfile.write("%s\n" % (outpfx+file))
        if (verbose):
            print "Process %s with %s -> %s" % (file,bpmfile,outpfx+file)
        update_head(outpfx+file,"BPMFIX",bpmfile,"File containing fixed BPM")
        update_head(outpfx+file,"NAME2",file,"Previous name of file")

    if (len(inpfiles) > 0 and not fake):
        if (verbose):
            print "Running crfix"

        # make sure package is loaded
        if not iraf.imred._loaded:
            iraf.imred()
        if not iraf.imred.crutil._loaded:
            iraf.imred.crutil()

        if (len(inpfiles.split(",")) > batchnum):
            # need to split it up into batches for processing
            j=int(len(inpfiles.split(","))/batchnum)+1
            for i in range(0,int(len(inpfiles.split(","))/batchnum)+1):
                inpf=",".join(inpfiles.split(",")[batchnum*i:batchnum*(i+1)])
                maskf=",".join(maskfiles.split(",")[batchnum*i:batchnum*(i+1)])
                if (verbose):
                    print "\tProcessing batch %d/%d" % (i+1,j)
                iraf.crfix(input=inpf,output=inpf,crmask=maskf)
        else:
            iraf.crfix(input=inpfiles,output=inpfiles,crmask=maskfiles)

    if (verbose):
        print "List of output written to %s" % outlist
        print "BPM-fixing finished!"  

######################################################################

def ircrzap(inlist,outpfx="c",fake=no,bpmkey="BPM0",bpmthresh=5,bpmupdate=no,clobber=globclob,verbose=globver):
    """ cosmic-ray zapping"""

    batchnum=20

    check_exist(inlist,"r")
    outlist=inlist+".zap"
    check_exist(outlist,"w",clobber)
    outfile=open(outlist,"w")
    if (len(outpfx) == 0):
        print "Warning: this will overwrite old files"

    files=getlines(inlist)
    inpfiles=""
    outfiles=""
    crmaskfiles=""
    for image in files:
        check_exist(image,"r")
        # check to see if has already been done
        imgets(image,"CRZAP",Stdout=1,Stderr=1)
        if (len(imgets.value) > 0 and int(imgets.value) > 0):
            print "Warning: %s has already been zapped.  Skipping." % file
            continue
        # get static bpm
        imgets(image,bpmkey,Stdout=1,Stderr=1)
        bpmfile=imgets.value
        if (len(bpmfile) <= 1):
            print "No BPM file for image %s" % file
            sys.exit(1)
        check_exist(bpmfile,"r")
        if (len(outpfx) > 0):
            check_exist(outpfx+image,"w",clobber)
        outfile.write("%s\n" % (outpfx+image))
        if (len(inpfiles) > 0):
            inpfiles+=","
        if (len(outfiles) > 0):
            outfiles+=","
        if (len(crmaskfiles) > 0):
            crmaskfiles+=","
        rootname=os.path.splitext(os.path.basename(image))[0]
        crmask="crmask_"+rootname+".pl"
        check_exist(crmask,"w",clobber)
        inpfiles+=image
        outfiles+=outpfx+image
        crmaskfiles+=crmask
        if (verbose):
            print "Zap %s -> %s mask=%s" % (image,outpfx+image,crmask)
    if (len(inpfiles) > 0 and not fake):
        if (verbose):
            print "Running craverage"
        # make sure package is loaded
        try:
            s=iraf.lpar(iraf.crfix,Stdout=1)            
        except:
            iraf.imred.crutil()
            
        craverage=iraf.craverage
        craverage.average=""
        craverage.sigma=""
        craverage.navg=5
        craverage.nrej=0
        craverage.nbkg=5
        craverage.nsig=25
        craverage.var0=0
        craverage.var1=0
        craverage.var2=0
        craverage.crval=1
        craverage.lcrsig=10
        craverage.hcrsig=5
        craverage.crgrow=0
        craverage.objval=0

        if (len(inpfiles.split(",")) > batchnum):
            # need to split it up into batches for processing
            j=int(len(inpfiles.split(","))/batchnum)+1
            for i in range(0,int(len(inpfiles.split(","))/batchnum)+1):
                inpf=",".join(inpfiles.split(",")[batchnum*i:batchnum*(i+1)])
                outf=",".join(outfiles.split(",")[batchnum*i:batchnum*(i+1)])
                crf=",".join(crmaskfiles.split(",")[batchnum*i:batchnum*(i+1)])
                if (verbose):
                    print "\tProcessing batch %d/%d" % (i+1,j)
                craverage(input=inpf,output=outf,crmask=crf)
        else:
            craverage(input=inpfiles,output=outfiles,crmask=crmaskfiles)
        
        # update headers
        if (verbose):
            print("Updating headers")
        files=outfiles.split(",")
        masks=crmaskfiles.split(",")
        for i in range(0,len(files)):
            update_head(files[i],"CRMASK",masks[i],"CR mask")
        if (verbose):
            print "Done with zapping!"
        outfile.close()
    if (bpmupdate):
        check_exist("crsum.pl","w",clobber)
        # all files must have the same BPM
	# get from last file
        check_exist("temp_sum.pl","w",yes)

        imsum=iraf.imsum
        imsum.option="sum"
        imsum.low_rej=0
	imsum.high_rej=0
	# sum up all individual masks
        # do the summing with the input from a file
        check_exist("temp_out.list","w",yes)
        tempout=open("temp_out.list","w")
        for s in crmaskfiles.split(","):
            tempout.write("%s\n" % s)
        tempout.close()
	imsum(input="@temp_out.list",output="crsum.pl",verbose=verbose)
        check_exist("temp_out.list","w",yes)

        # find those pixels that have > bpmthresh flags
	iraf.imcalc("crsum.pl","temp_sum.pl","im1*0",verbose=verbose)
	iraf.imarith("temp_sum.pl","+",bpmthresh,"temp_sum.pl",verbose=verbose)
	iraf.imcalc(input="crsum.pl,temp_sum.pl",output="crsum.pl",
                    equals="if im1 .ge. im2 then 1 else 0",verbose=verbose)
        # combine with BPM
	s1="crsum.pl,"+bpmfile
        s2="crsum2.pl"
        check_exist(s2,"w",yes)
        iraf.imcalc(input=s1,output=s2,
                    equals="if (im1 + im2) .gt. 1 then 1 else 0",
                    verbose=verbose)
        os.rename(bpmfile,bpmfile+".old")
        os.rename(s2,bpmfile)
	hedit(bpmfile,fields="UPDATE",value=s1,add=yes,verify=no,show=verbose)
        if (verbose):
            print "Mask sums written to crsum.pl"
        os.remove("temp_sum.pl")
    if (verbose):
        print "List of output written to %s" % outlist

######################################################################

def irmosaic(ref_pic,picfile,outpfx,dsply,xpaid,starfile,findstars,
             smscale,doauto,docheck,autobox,checksnr,minsnr,makemask,
             compute,align,combine,expand,shift,pixshift,scale,rotate,
             distort,combmode,combtype,reject,lthreshold,hthreshold,
             nlow,nhigh,nkeep,mclip,lsigma,hsigma,rdnoise,gain,snoise,
             sigscale,pclip,grow,subsets,nights,finaltrim,objkey,filtkey,
             nitekey,timekey,expkey,coaddkey,rakey,deckey,rotkey,
             pixscale,rotoffset,rasign,decsign,raindeg,
             clobber=globclob,verbose=globver):

    """ mosaicing"""

    DX={}
    DY={}
    snrs={}
    status={}
    border=-100000
    oldmask=no
    calcsnr=no
    DTOR=0.01745329252 # degrees to radians conversion
    doexpand=expand

    # -------------------------------------------------
    # initial setup

    # get reference file & size
    if (os.path.exists(ref_pic)):
        ref=ref_pic
    else:
        if (os.path.exists(ref_pic+".fits")):
            ref=ref_pic+".fits"
        else:
            ref=ref_pic+".imh"
    if (not (os.path.exists(ref))):
        print "Couldn't open reference image %s" % ref_pic
        sys.exit(1)

    check_exist(picfile,"r")

    # get dimensions of input image
    Xmin=1
    imgets(ref,"i_naxis1")	
    Xmax=float(imgets.value)
    xsize=Xmax
    Ymin=1	
    imgets(ref,"i_naxis2")
    Ymax=float(imgets.value)
    ysize=Ymax
    Xcent=(Xmax+Xmin)/2
    Ycent=(Ymax+Ymin)/2

    if (doexpand):
        Xmax*=2
        Ymax*=2
        # this is the section that we didn't expand yet
        goodsec="[%d:%d,%d:%d]" % (int(Xmax/4),int(3*Xmax/4),int(Ymax/4),int(3*Ymax/4))
        xoffset=Xmax/4
        yoffset=Ymax/4
    else:
        goodsec=""
        xoffset=0
        yoffset=0
                
    # -------------------------------------------------
    # find stars
    if (findstars):

        if (not check_xpa(xpaid)):
            print "XPA ID %s is not valid" % xpaid
            sys.exit(2)
        
        if (starfile != ""):
            check_exist(starfile,"r",clobber)
            starlist=yes
        else:
            starlist=no

        if (not doexpand):
            check_exist(outpfx+ref,"w",clobber)
            iraf.imcopy(ref,outpfx+ref,verbose=verbose)
        else:
            # make image twice original size, for good borders
            if (verbose):
                print "Expanding image %s" % ref
            expand_img(ref,outpfx+ref,scale=2,border=border,clobber=clobber,verbose=no)

        # -------------------------------------------------
        # get reference stars from ref_pic   
        todisp=ref
        ref=outpfx+ref
        check_exist("sm_"+ref,"w",yes)
        # smooth if required
        if (smscale > 0):
            if (doexpand):
                # need to make smoothed version of original image
                check_exist("sm__temp.fits","w",yes)
                iraf.gauss(todisp,"sm__temp.fits",sigma=smscale,ratio=1.0,theta=0,nsigma=3)
                todisp="_temp.fits"
            else:
                todisp="sm_" + todisp

            if (verbose):
                print "Smoothing image with sigma=%d" % smscale
            iraf.gauss(ref,"sm_"+ref,sigma=smscale,ratio=1.0,theta=0,nsigma=3)
        
        # should we display the first image
        cwd=os.getcwd()
        s="xpaset -p %s frame first" % xpaid
        os.system(s)
        if (dsply):
            print "Displaying image"
            if (smscale > 0):
                s="xpaset -p %s file %s/sm_%s" % (xpaid,cwd,todisp)
                os.system(s)                
            else:
                s="xpaset -p %s file %s/%s" % (xpaid,cwd,todisp)
                os.system(s)
            s="xpaset -p %s scale zscale" % xpaid
            os.system(s)
        
	# also write SNR's to file
        check_exist(ref+".snrs","w",clobber)
        if (smscale > 0):            
            centerfile="sm_" + ref
        else:
            centerfile=ref
        # positions are in un-expanded image
        refStar=[]
        snrs={}
        DX={}
        DY={}
        refStar=[]
        if (starlist):
            # get star positions from file
            istars=getlines(starfile)
            for i in range(0,len(istars)):
                s=istars[i]
                [xi,yi]=s.split()
                xi=float(xi)
                yi=float(yi)
                refStar+=[[xi+xoffset,yi+yoffset]]
        else:
            # find stars from regions            
            print "Place green regions around reference stars"
            rsp=raw_input(" (press enter when ready, q to quit) ")
            if (rsp.upper().count("Q") > 0):
                sys.exit(2)
            regs=getregions(xpaid,"green")
            if (verbose):
                print "Read in %d regions" % len(regs)
            check_exist("_temp2.stars","w",yes)
            ofile2=open("_temp2.stars","w")
            # apply offset & write to file
            for i in range(0,len(regs)):
                regs[i][0]+=xoffset
                regs[i][1]+=yoffset
                ofile2.write("%f %f\n" % (regs[i][0],regs[i][1]))
            ofile2.close()
            # center on stars
            iraf.center.interactive=no
            iraf.datapars.scale=1.0
            iraf.datapars.sigma=INDEF
            iraf.datapars.datamin=INDEF
            iraf.datapars.datamax=INDEF
            iraf.datapars.exposure=expkey
            iraf.datapars.airmass="AIRMASS"
            iraf.datapars.filter=filtkey
            iraf.datapars.obstime=nitekey
            iraf.centerpars.calgorithm="centroid"
            iraf.centerpars.cbox=5
            cent=iraf.center(image=centerfile,coords="_temp2.stars",
                             output="",verify=no,update=yes,verbose=yes,
                             Stdout=1,Stderr=1)
            os.remove("_temp2.stars")
            goodstars=0
            for line in cent:                
                if (line.find("ok") >= 0):
                    # centering is OK
                    s=line.split()
                    xi=float(s[3])
                    yi=float(s[4])
                    refStar+=[[xi,yi]]
                    goodstars+=1
                                
        print "Found %d good stars" % len(refStar)
        if (len(refStar) < 1):
            print "Not enough stars: quiting"
            sys.exit(2)
        # print stars, display regions, compute SNRs
        for i in range(0,len(refStar)):
            xi=refStar[i][0]-xoffset
            yi=refStar[i][1]-yoffset
            refStar[i]=[xi,yi]
            if (0):
                s="echo \"image; circle(%f,%f,10) # color=red\" | xpaset %s regions" % (xi,yi,xpaid)
                os.system(s)

            if (calcsnr):
                snr=round(find_snr(ref,xi+xoffset,yi+yoffset,"",verbose=no),1)
            else:
                snr=100
            if (i==0):                
		snrs[ref_pic]=snr
            if (verbose):
                print "Reference star %02d: x=%7.2f, y=%7.2f, snr=%4.2g" % (i,xi,yi,snr)
            hedit(ref,fields="SNR",value=snr,add=yes,update=yes,veri=no,show=no)
        s="xpaset -p %s regions color green" % xpaid
        os.system(s)
        
        # get position info from header
        inf=pyfits.open(ref,"update")
        hdr=inf[0].header
        raref=hdr.get(rakey,"null")
        decref=hdr.get(deckey,"null")
        if (rotkey.upper() != "NONE"):
            rotref=float(hdr.get(rotkey,0))
        else:
            rotref=0
        # offset, specified in arcsec
        raoff=float(hdr.get("RAOFF",0))
        decoff=float(hdr.get("DECOFF",0))

        rotref+=rotoffset
        # convert to decimal
        if (str(raref).find(":") >= 0):
            [raref,decref]=sex2dec(raref,decref)
        # apply offsets
        if (raindeg):
            raref-=raoff*math.cos(decref*DTOR)/3600
        else:
            raref-=raoff*math.cos(decref*DTOR)/3600/15
        decref-=decoff/3600

	# convert to radians
       	rotref*=DTOR

        if (combmode == "exp"):
            # make a new keyword, = Tint*Coadd
            exptime=hdr.get(expkey,1)
            coadds=hdr.get(coaddkey,1)
            exptime*=coadds
            hdr.update("EXPEFF",exptime,comment="Effective exposure")
            
        inf.close()
        check_exist("sm_" + ref,"w",yes)
        status[ref_pic]="ref"
        DX[ref_pic]=0
        DY[ref_pic]=0


        # -------------------------------------------------
        # now get data from the rest of the images
        print " "

	# return Max/Min to original values	
        if (doexpand):
	    Xmax=Xmax/2
	    Ymax=Ymax/2

        sign=-1
        last=""
        files=getlines(picfile)
        nimage=0
        Stars={}

        # Overall pixel drift
        DDX=0
        DDY=0
        
        for fileline in (files):
            s=fileline.split()
            filename=s[0]
            # get offset, if specified
            raoff=0
            decoff=0
            if (len(s) > 1):
                raoff=float(s[1])
            if (len(s) > 2):
                decoff=float(s[2])
            # make sure we don't analyze the reference image again, 
	    # or a blank line
            if (re.search("\S+",filename)):
                check_exist(filename,"r")
            else:
                continue
            if (filename == ref or filename == outpfx+ref or outpfx+filename == ref):
                continue

            nimage+=1
#            if (nimage > 1):
#                continue
            if (verbose):
                print " "
                print "Analyzing image # %d: %s" % (nimage,filename)

            todisp=filename
            if (smscale > 0):
                check_exist("sm_" + filename,"w",yes)
                iraf.gauss(filename,"sm_"+filename,sigma=smscale,ratio=1.0,theta=0,nsigma=3)
                todisp="sm_" + todisp

            # display the image to frame 2
            s="xpaset -p %s frame 2" % xpaid
            os.system(s)
            s="xpaset -p %s file %s/%s" % (xpaid,cwd,todisp)
            os.system(s)                
            s="xpaset -p %s scale zscale" % xpaid
            os.system(s)
            if (combmode == "exp"):
                # make a new keyword, = Tint*Coadd
                inf=pyfits.open(filename,"update")
                hdr=inf[0].header
                exptime=hdr.get(expkey,1)
                coadds=hdr.get(coaddkey,1)
                exptime*=coadds
                hdr.update("EXPEFF",exptime,comment="Effective exposure")
                inf.close()

            if (doauto):
                # try to get shift from header
                inf=pyfits.open(filename)
                hdr=inf[0].header
                ra=hdr.get(rakey,"null")
                dec=hdr.get(deckey,"null")
                if (rotkey.upper() != "NONE"):
                    rot=float(hdr.get(rotkey,0))
                else:
                    rot=0

                # offset, specified in arcsec
                raoff+=float(hdr.get("RAOFF",0))
                decoff+=float(hdr.get("DECOFF",0))
                if abs(raoff)>0 or abs(decoff)>0:
                    print "RAOFF: %3.1f, DECOFF: %3.1f" % (raoff, decoff)

                # convert to decimal
                if (str(ra).find(":") >= 0):
                    [ra,dec]=sex2dec(ra,dec)
                # apply offsets
                if (raindeg):
                    ra-=raoff/(3600*math.cos(dec*DTOR))
                else:
                    ra-=raoff/(15*3600*math.cos(dec*DTOR))
                dec-=decoff/3600
                rot+=rotoffset
                rot*=3.141592654/180

                # figure out offsets, in degrees of arc
                ddec=decref-dec
                dra=raref-ra
                dra*=abs(math.cos((DTOR)*(dec+decref)/2))
                if (not raindeg):
                    dra*=15
                if (rasign != 0):
                    dra*=rasign
                if (decsign != 0):
                    ddec*=decsign
                # scale in arcsec/pixel
		# convert to pixels
		dra*=3600/pixscale
		ddec*=3600/pixscale
                
                # Debugging
                #print "RA-ref: %10.6f, DEC-ref: %11.6f" % (raref, decref)
                #print "RA:     %10.6f, DEC:     %11.6f" % (ra,dec)
                #print "DRA: %.1f, DDEC: %.1f" % (dra, ddec)
                #print "DDX: %.1f,  DDY: %.1f" % (DDX,DDY)

                dxpix=-(dra*math.cos(rot)-ddec*math.sin(rot)) + DDX
                dypix=-(-dra*math.sin(rot)-ddec*math.cos(rot)) + DDY

                # predict where stars should be
       		nstarfound=0
                starstate=range(0,len(refStar))
                xstar=range(0,len(refStar))
                ystar=range(0,len(refStar))
                Star=[]
                snrs[filename]=0
                # more refined shifts
                dx=0
                dy=0
                for i in range(0,len(refStar)):
                    Star.append([-1,-1])
                    if (verbose):
                        print "Trying to find star #%d" % i
                    sign=1
                    # transform ra/dec -> x/y
		    # assumes N is up, E to left
		    rot=rotref
                    xstar[i]=refStar[i][0] + dxpix + dx
                    ystar[i]=refStar[i][1] + dypix + dy
                    if (verbose):
                        print "Offset %5.1f, %5.1f pixels" % \
                              (dxpix+dx,dypix+dy)
                    # see if star is on image
                    offimage=1
                    if (xstar[i] > (Xmax-5) or (xstar[i] < 5)):
                        # state=-1 -> off image
			starstate[i]=-1
                        if (verbose):
                            print "Star #%d is off the image" % i
                        continue  
                    if (ystar[i] > (Ymax-5) or (ystar[i] < 5)):
                        # state=-1 -> off image
			starstate[i]=-1
                        if (verbose):
                            print "Star #%d is off the image" % i
                        continue  

                    # if it got here, must be on image
                    offimage=0
                    if (verbose):
                        print "Star #%d should be at %5.1f,%5.1f" % (i,xstar[i],ystar[i])
                    if (nstarfound == 0):
                        # make a box
                        xl=int(xstar[i])-autobox			
                        xh=int(xstar[i])+autobox
                        yl=int(ystar[i])-autobox
                        yh=int(ystar[i])+autobox
                        # check edges
                        if (xl > Xmax or yl > Ymax):
                            continue
                        if (xh < 1 or yh < 1):
                            continue
                        if (xl < 1):
                            xl=1
                        if (yl < 1):
                            yl=1
                        if (xh > Xmax):
                            xh=Xmax
                        if (yh > Ymax):
                            yh=Ymax
                        boxsec="[%d:%d,%d:%d]" % (xl,xh,yl,yh)
                        s="echo \"image; box(%f,%f,%f,%f) # color=green\" | xpaset %s regions" % (xstar[i],ystar[i],2*autobox,2*autobox,xpaid)
                        os.system(s)

                        # use maximum pixel in box as basis for centering
                        if (smscale > 0):
                            s=iraf.minmax("sm_" + filename + boxsec,update=no,force=yes,verbose=yes,Stdout=1)
                        else:
                            s=iraf.minmax(filename + boxsec,update=no,force=yes,verbose=yes,Stdout=1)                        
                        maxpix=iraf.minmax.maxpix
                        maxpix=re.sub("\[","",maxpix)
                        maxpix=re.sub("\]","",maxpix)
                        [xmax,ymax]=maxpix.split(",")
                        xmax=int(xmax)+xl-1
                        ymax=int(ymax)+yl-1
                    else:
                        xmax=xstar[i]
                        ymax=ystar[i]

                    if (verbose):
                        print "Centering..."
                    check_exist("_temp.star2","w",yes)
                    ofile2=open("_temp2.stars","w")
                    ofile2.write("%d %d" % (xmax,ymax))
                    ofile2.close()
                    iraf.centerpars.calgorithm="centroid"
		    iraf.centerpars.cbox=2*autobox
		    iraf.centerpars.cthreshold=1.0
		    iraf.centerpars.maxshift=autobox
                    centerfile=filename
                    if (smscale > 0):
                        centerfile="sm_" + filename
                    cent=iraf.center(image=centerfile,coords="_temp2.stars",output="",verify=no,update=yes,verbose=yes,Stdout=1,Stderr=1)
                    for line in cent:
                        if (line.find("Warning") >= 0):
                            continue
                        s=line.split()
                    if (s[len(s)-1] == "ok"):
                        xi=float(s[3])
                        yi=float(s[4])
                    else:
                        xi=-1
                        yi=-1
                        continue

                    Star[i]=[xi,yi]
                    os.remove("_temp2.stars")
                    # mark a circle
                    s="echo \"image; circle(%f,%f,%f) # color=red\" | xpaset %s regions" % (xi,yi,autobox/2,xpaid)
                    os.system(s)

                    # This is the first detected star in the image
                    if (nstarfound == 0):
                        # get snr
                        if (calcsnr):
                            snrs[filename]=find_snr(filename,xi,yi,"",verbose=no)
                            hedit(filename,fields="SNR",value=snrs[filename],
                                  add=yes,update=yes,veri=no,show=no)
                        else:
                            snrs[filename]=10
                        # more refined shifts
                        dx=xi-xstar[i]
                        dy=yi-ystar[i]
                    nstarfound+=1
                    
            Stars[filename]=Star
            if (verbose):
                print "Found %d / %d stars for image #%d" % (nstarfound,len(refStar),nimage)
            status[filename]="yes"

            if (abs(rot-rotref) > 2*3.14159/180):
                # too big angle change - exclude
                # exclude from combination
                print "Angle difference too large - excluding"
                status[filename]="no"
            if (nstarfound < 1):
                # too few stars
                status[filename]="no"                    
                
            if (docheck and status[filename] != "no"):
                # check if things worked out
                s=raw_input("Are identifications correct (y/n/c/f/q) [y]? ").upper()
                if (not re.search("\S+",s)):
                    s="Y"
                if (s.find("N") >= 0):
                    # something went wrong
                    status[filename]="no"
                if (s.find("Q") >= 0):
                    # quit
                    sys.exit(2)
                if (s.find("C") >= 0):
                    # stop checking
                    status[filename]="yes"
                    docheck=no
                if (s.find("F") >= 0):
                    # fix manually
                    status[filename]="no"
                    print "Make the region to correct blue, place a cyan region around the correct position"
                    rsp=raw_input(" (press enter when ready, q to quit) ")
                    if (rsp.upper().count("Q") > 0):
                        sys.exit(2)
                    regstocor=getregions(xpaid,"blue")
                    regscor=getregions(xpaid,"cyan")
                    dxcor=regstocor[0][0]-regscor[0][0]
                    dycor=regstocor[0][1]-regscor[0][1]
                    if (verbose):
                        print "Correction offset: %5.1f %5.1f pixels" % (dxcor,dycor)
                    dracor=pixscale*sign*(dxcor*math.cos(rot)-dycor*math.sin(rot))
                    ddeccor=-pixscale*sign*(-dxcor*math.sin(rot)-dycor*math.cos(rot))
                    print "Correction offset: %5.1f %5.1f arcsec" % (dracor,ddeccor)
                    dracor+=raoff
                    ddeccor+=decoff
                    print "Net correction offset: %5.1f %5.1f arcsec" % (dracor,ddeccor)
                    rsp=raw_input(" (press enter to continue, q to quit)")
                    if (rsp.upper().count("Q") > 0):
                        sys.exit(2)

            xshift=0
            yshift=0
            for i in range(0,len(refStar)):
                if (Star[i][0] >= 0 and Star[i][1] >= 0):
                    xshift+=Star[i][0]-refStar[i][0]
                    yshift+=Star[i][1]-refStar[i][1]
            
            if (nstarfound > 0):
                xshift/=nstarfound
                yshift/=nstarfound            

            DX[filename]=xshift
            DY[filename]=yshift

            print "Mean shift: %6.2f, %6.2f pixels" % (xshift,yshift)

            # Write star locations to geomap input file
            ifile=filename + ".inp"
            check_exist(ifile,"w",yes)
            IFILE=open(ifile,"w")
            for i in range(0,len(refStar)):
                if (Star[i][0] >= 0 and Star[i][1] >= 0):
                    s="%f %f %f %f\n" % \
                       (refStar[i][0]+xoffset,refStar[i][1]+yoffset,
                        Star[i][0],Star[i][1])
                    IFILE.write(s)
            IFILE.close()

            # Try to track the systemic drift...
            if status[filename]=="yes":
                DDX+=(xshift-dxpix)
                DDY+=(yshift-dypix)

            if (verbose):
                print "Done with %s" % filename

        # Write listing of files, status, shifts
        shiftfile=picfile+'.shifts'
        check_exist(shiftfile,"w",yes)
        sfout=open(shiftfile,"w")
        for fileline in (files):
            s=fileline.split()
            filename=s[0]
            sfout.write("%s %4s %.2f %.2f %.2f\n" % \
                        (filename,status[filename],snrs[filename],
                         DX[filename],DY[filename]))
        sfout.close()

        # put these back
        if (doexpand):
            Xmax*=2
            Ymax*=2
        
        if (verbose):
            print "Done with list"

    if (verbose):
        print "Done finding stars"            

    # -------------------------------------------------

    # Read shifts from the file if necessary
    if (len(status.keys())==0):
        shiftfile=picfile+'.shifts'
        if (verbose):
            print "Reading in shifts from %s . . ." % shiftfile
        shiftlines=getlines(shiftfile)
        for line in (shiftlines):
            s=line.split()
            filename=s[0]
            status[filename]=s[1]
            snrs[filename]=float(s[2])
            DX[filename]=float(s[3])
            DY[filename]=float(s[4])

    # Organize the list of filenames
    keys=status.keys()
    keys.sort()
        
    # -------------------------------------------------
    # make sure they're good (again)

    if (checksnr):
        if (verbose):
	    print "Checking SNR . . ."
	if (verbose):
	    print "Deleting files with SNR < %f" % minsnr

        print "#\tfile\t\tSNR\t\tTint\t\tCoadd\t\tExpeff\tDX\tDY\n"
        files=[]
        for filename in (keys):
            if (status[filename] != "ref"):                    
                snr=snrs[filename]
            else:
                snr=snrs[0]
            # make sure SNR is OK
            if (snr < minsnr):
                status[filename]="no"

            if (status[filename] == "yes" or status[filename] == "ref"):
                files.append(filename)
                # it's good                
                inf=pyfits.open(filename)
                hdr=inf[0].header
                exptime=hdr.get(expkey,1)
                coadds=hdr.get(coaddkey,1)
                expeff=hdr.get("EXPEFF",1)
                
                s="%d %s\t%f\t%f\t%d\t%f\t%f\t%f" % (len(files),filename,snr,exptime,coadds,expeff,DX[filename],DY[filename])
                print s

        s=raw_input("Exclude any (y/n) [n]? ").upper()
        if (s.find("Y") >= 0):
            # gotta remove some
	    j=1
            while (j != 0):
                j=raw_input("Exclude which number (0 to quit) [0]?")
                if (re.search("\S+",j)):
                    j=int(j)
                else:
                    j=0
                if (j > len(files)):
                    print "%d is greater than the number of files" % j
                else:
                    status[files[j-1]]="no"
                    print "File %d (%s) deleted" % (j,files[j-1])

    # -------------------------------------------------
    # now compute transform	

    if (compute):
        iraf.unlearn(iraf.gregister)
        iraf.unlearn(iraf.geomap)

        # get reference image
        for filename in (keys):
            if (status[filename] == "ref"):
                ref=filename

        # Check for existence of "shifted" ref image
        if not os.path.exists(outpfx+ref):
            if doexpand:
                if verbose:
                    print "Expanding reference image"
                expand_img(ref,outpfx+ref,scale=2,border=border,
                           clobber=clobber,verbose=no)
            else:
                if verbose:
                    print "Copying over reference image"
                iraf.imcopy(ref,outpfx+ref,verbose=verbose)

        # get dimensions of shifted/expanded ref image
        xmin=1
        imgets(outpfx+ref,"i_naxis1")
        xmax=float(imgets.value)
        ymin=1
        imgets(outpfx+ref,"i_naxis2")
        ymax=float(imgets.value)
        
        if (doexpand):
            XMIN=xmax/4
            XMAX=3*xmax/4
            YMIN=ymax/4
            YMAX=3*ymax/4

        # Figure out type of transform
        if (distort):
            tform_type="general"
        elif (shift and pixshift):
            tform_type="pixshift"
        elif (rotate and scale and shift):
            tform_type="rxyscale"
        elif (rotate and shift and not scale):
            tform_type="rotate"
        elif (shift and scale and not rotate):
            tform_type="xyscale"
        elif (shift and not scale and not rotate):
            tform_type="shift"
        else:
            tform_type="shift"

        # Process the rest of the images
        nfiles=0
        keys=status.keys()
        keys.sort()
        for filename in (keys):
            if (status[filename] != "yes"):
                continue
            nfiles+=1
            s="\nAnalyzing image #%d: %s (%f,%f)" % \
               (nfiles,filename,DX[filename],DY[filename])
            print s

	    # output filename
	    ofile=filename + ".tform"
            check_exist(ofile,"w",yes)
            outimage=outpfx + filename
            check_exist(outimage,"w",clobber)

            # Pixshift option, very simple
            if (tform_type=="pixshift"):
                # Create output image
                iscale=2
                if not doexpand:
                    iscale=1
                expand_img(filename,outimage,scale=iscale,border=border,
                           nodata=1,verbose=verbose)
                # Expected destination region
                jxmin=int(round(1+xoffset-DX[filename]))
                jxmax=jxmin+xsize-1
                jymin=int(round(1+yoffset-DY[filename]))
                jymax=jymin+ysize-1
                # Correct for image limits
                (jxmin,djxmin)=check_limit(jxmin,1,lower=yes)
                (jxmax,djxmax)=check_limit(jxmax,Xmax,upper=yes)
                (jymin,djymin)=check_limit(jymin,1,lower=yes)
                (jymax,djymax)=check_limit(jymax,Ymax,upper=yes)
                osec="[%d:%d,%d:%d]" % (jxmin,jxmax,jymin,jymax)
                # Adjust input section, if necessary
                if (djxmin+djxmax+djymin+djymax > 0):
                    if (djxmin > 0):
                        isec="[%d:%d," % (1+djxmin,xsize)
                    elif (djxmax > 0):
                        isec="[%d:%d," % (1,xsize-djxmax)
                    else:
                        isec="[%d:%d," % (1,xsize)
                    if (djymin > 0):
                        isec+="%d:%d]" % (1+djymin,ysize)
                    elif (djymax > 0):
                        isec+="%d:%d]" % (1,ysize-djymax)
                    else:
                        isec+="%d:%d]" % (1,ysize)
                else:
                    isec=""
                iraf.imcopy(filename+isec,outimage+osec,verbose=verbose)
                continue

            # check for existence of star coords file
            ifile=filename + ".inp"
            check_exist(ifile,"r")

            # run geomap
	    iraf.geomap(ifile,ofile,Xmin,Xmax,Ymin,Ymax,
                        fitgeometry=tform_type,interact=no)
            
            # need to get shift
            lines=getlines(ofile)
            for line in lines:
                if (line.find("xshift") >= 0):
                    s=line.split()
                    DX[filename]=float(s[1])+xoffset
                    xshift=DX[filename]
                if (line.find("yshift") >= 0):
                    s=line.split()
                    DY[filename]=float(s[1])+yoffset
                    yshift=DY[filename]
            if (verbose):
                print "Final shifts: %f,%f" % (xshift,yshift)

            if (doexpand):
                if (Xmax/4-xshift < XMIN):
                    XMIN=Xmax/4-xshift
                if (Ymax/4-yshift < YMIN):
                    YMIN=Ymax/4-yshift
                if (3*Xmax/4-xshift > XMAX):
                    XMAX=3*Xmax/4-xshift
                if (3*Ymax/4-yshift > YMAX):
                    YMAX=3*Ymax/4-yshift

            # now do gregister
            if (align):
                check_exist(outimage,"w",clobber)
		iraf.gregister.boundary="constant"
		iraf.gregister.constant=border
		iraf.gregister.geometr="geometric"
		iraf.gregister.xmin="INDEF"   
		iraf.gregister.xmax="INDEF"   
		iraf.gregister.ymin="INDEF"   
		iraf.gregister.ymax="INDEF"   
		iraf.gregister.xscale="INDEF" 
		iraf.gregister.yscale="INDEF"
		iraf.gregister.ncols=Xmax
		iraf.gregister.nlines=Ymax
		iraf.gregister.xsample=1
		iraf.gregister.ysample=1
		iraf.gregister.interpolant="poly3"
                #iraf.gregister.interpolant="linear"
		iraf.gregister.fluxcon=yes
                iraf.gregister(filename,outimage,ofile,transforms=filename + ".inp")

    		# trim border
                if (verbose):
		    print "Trimming border"		
                # get dimensions of  image
		xmin2=1-(xshift-xoffset)
                if ((xmin2-2) <= 1):
		    xmin2=3
		imgets(outimage,"i_naxis1")
		xmax2=float(imgets.value)/2-(xshift-xoffset)
                if (xmax2+2 >= float(imgets.value)):
		    xmax2=float(imgets.value)-2
                ymin2=1-(yshift-yoffset)
                if ((ymin2-2) <= 1):
		    ymin2=3
		imgets(outimage,"i_naxis2")
		ymax2=float(imgets.value)/2-(yshift-yoffset)
                if (ymax2+2 >= float(imgets.value)):
		    ymax2=float(imgets.value)-2
                iraf.imreplace.radius=8
		iraf.imreplace.lower=INDEF
		iraf.imreplace.upper=INDEF
                s="%s[%d:%d,%d:%d]" % (outimage,int(xmin2)-2,int(xmax2)+2,int(ymax2)-2,int(ymax2+2))
		iraf.imreplace(s,border)
                s="%s[%d:%d,%d:%d]" % (outimage,int(xmin2)-2,int(xmax2)+2,int(ymin2)-2,int(ymin2+2))
		iraf.imreplace(s,border)
                s="%s[%d:%d,%d:%d]" % (outimage,int(xmin2)-2,int(xmin2)+2,int(ymin2)-2,int(ymax2+2))
		iraf.imreplace(s,border)
                s="%s[%d:%d,%d:%d]" % (outimage,int(xmax2)-2,int(xmax2)+2,
                                       int(ymin2)-2,int(ymax2+2))
		iraf.imreplace(s,border)
            if (verbose):
		print "Done with transform\n"

            # clean up
            os.remove(ofile)

    # -------------------------------------------------
    # now combine images
    if (combine):
	print "Combining\n"

        # go through list of files to use
        keys=status.keys()
        keys.sort()
        filters={}
        nites={}
        types={}
        typefiltnum={}
        typenitenum={}
        filternum=[]
        nitenum=[]
        for filename in (keys):
            if (not( status[filename] == "ref" or status[filename] == "yes")):
                continue
            if (filename.find(outpfx) != 0):
                filename=outpfx + filename
            # get info about each file        
            # get filter
            inf=pyfits.open(filename)
            hdr=inf[0].header
            filter=hdr.get(filtkey,"none")
            # get night
            nite=hdr.get(nitekey,"none")

            if (filters.has_key(filter)):
                filters[filter] += ",%s" % filename
            else:
                filters[filter] = "%s" % filename
                filternum.append(filter)
            nfilts=filternum.index(filter)+1

            if (nites.has_key(nite)):
                nites[nite] += ",%s" % filename
            else:
                nites[nite] = "%s" % filename
                nitenum.append(nite)
            nnites=nitenum.index(nite)+1
            
            type=""
            if (subsets):            
                type+=filter
            if (nights):
                type+=nite
            
            if (types.has_key(type)):
                types[type] += ",%s" % filename
            else:
                types[type] = "%s" % filename
                typenitenum[type]=nnites
                typefiltnum[type]=nfilts
            keys=types.keys()
            keys.sort()
            ntypes=keys.index(type)+1

            s=filename
            if (not DX.has_key(s)):
                s=re.sub("^" + outpfx,"",s)

            xshift=DX[s]
            yshift=DY[s]

            s="File=%s\tfilter=%s (%d)\tnite=%s (%d) shift=(%f,%f)" % (filename,filter,nfilts,nite,nnites,xshift,yshift)
            print s

        print " "
        # set up imcombine
        # set this to reject boundaries from other images
        iraf.imcombine.lthreshold=border+50000

        iraf.imcombine.reject=reject
        iraf.imcombine.hthreshold=hthreshold
        iraf.imcombine.nlow=nlow
        iraf.imcombine.nhigh=nhigh     
        iraf.imcombine.nkeep=nkeep     
        iraf.imcombine.mclip=mclip     
        iraf.imcombine.lsigma=lsigma    
        iraf.imcombine.hsigma=hsigma    
        iraf.imcombine.sigscale=sigscale  
        iraf.imcombine.pclip=pclip     
        iraf.imcombine.grow=grow      
        iraf.imcombine.combine=combtype

        iraf.imcombine.rdnoise=rdnoise   
        iraf.imcombine.gain=gain      
        iraf.imcombine.snoise=snoise    
        iraf.imcombine.maskval=0
        iraf.imcombine.rejmask=""
        iraf.imcombine.bpmasks=""
        iraf.imcombine.sigma=""
        #iraf.imcombine.zero="median"
        iraf.imcombine.zero="none"
        iraf.imcombine.statsec=""
        # don't mask here
        iraf.imcombine.masktype="none"
        
        keys=types.keys()
        keys.sort()
        for type in keys:
            files=types[type].split(",")
                
            # construct final name
            inf=pyfits.open(files[0])
            object=inf[0].header.get(objkey,"none")
            final_name=object
            if (subsets):
                final_name+="_s%d" % typefiltnum[type]
            if (nights):
                final_name+="_n%d" % typenitenum[type]
        
            if (os.path.exists(final_name + ".fits") and not clobber):
                k=1
                while (os.path.exists("%s_%d.fits" % (final_name,k))):
                    k+=1
                final_name="%s_%d" % (final_name,k)
            final_name+=".fits"                 
            check_exist(final_name,"w",clobber)
            check_exist(final_name + "_exp.pl","w",yes)
            if (verbose):
                print "Writing to %s" % final_name
                
            if (reject=="sigclip"):
                # output sigma file if necessary (for iteration)
                check_exist(final_name + "_sigma.fits","w",yes)
                imcombine.sigma=final_name+"_sigma"

            # determine parameters to imcombine based on type of combination  
            check_exist(".temp.list","w",yes)
            tempout=open(".temp.list","w")
            s=re.sub(",","\n",types[type])
            tempout.write(s)
            tempout.write("\n")
            tempout.close()
            if (combmode=="exp"):
                imcombine("@.temp.list",output=final_name,expmasks=final_name+"_exp",scale="exposure",weight="exposure",expname="EXPEFF",project=no)
            elif (combmode=="sclwts"):
                sclwts("@.temp.list",stars=ref + ".star",rootname=ref)
                imcombine (types[type],output=final_name,expmasks=final_name+"_exp",scale="@"+ref+".scl",weight="@"+ref+".wts")
            elif (combmode=="snr"):
                imcombine ("@.temp.list",output=final_name,expmasks=final_name+"_exp",scale="!SNR",weight="!SNR")
            else:
                print "Be sure STATSEC is correctly set"
                imcombine ("@.temp.list",output=final_name,expmasks=final_name+"_exp",scale="median",weight="median")

            os.remove(".temp.list")
            ntoproc=len(types[type].split(","))
            # convert exposure map to correct sign
            # (as is, it is the # of pixels rejected)
            if (os.path.exists(final_name+"_exp.pl")):
                s0="_temp_exp.pl"
                check_exist(s0,"w",yes)
                iraf.imarith(ntoproc,"-",final_name+"_exp.pl",s0)
                os.remove(final_name+"_exp.pl")
                iraf.rename(s0,final_name+"_exp.pl")

            if (makemask):            
                iraf.unlearn(iraf.imreplace)
                if (verbose):
                    print "Image created, now computing mask"
                check_exist("objmask_" + final_name + ".pl","w")
                check_exist("_f" + final_name + ".fits","w")
                check_exist("_m" + final_name + ".fits","w")
                
                # first we make the mask
                iraf.imreplace.radius=0
                iraf.makemask.prefix="objmask_"
                iraf.makemask.headlist=""
                iraf.makemask.subsample=1
                iraf.makemask.filtsize=15
                iraf.makemask.nsmooth=3
                iraf.makemask.threshtype="nsigma"
                iraf.makemask.nsigthresh=2.0
                iraf.makemask.constthresh=0.0
                iraf.makemask.ngrow=1
                #iraf.makemask.ngrow=20
                iraf.makemask.statsec=goodsec
                iraf.makemask.checklimits=yes
                iraf.makemask.zmin=-32768
                iraf.makemask.zmax=32767
                iraf.makemask.verbose=verbose
                iraf.makemask(inlist=final_name)
                
                [root,ext]=os.path.splitext(final_name)
                print "Updating head of %s: %s -> %s" % (final_name,"OBJMASK","objmask_" +root+ ".pl")
                update_head(final_name,"OBJMASK","objmask_" +root+ ".pl","Name of object mask")
                maskname="objmask_" +root+ ".pl"
                
                # now we read back all the images that were used 
                # to create it, and shift back to their positions
                if (verbose):
                    print "\nShift object masks back to original positions"

                files=types[type].split(",")
                for file in files:
                    s=file
                    if (not status.has_key(s)):
                        s=re.sub("^" + outpfx,"",s)
                    if (status[s] == "ref"):
                        xshift=-xoffset
                        yshift=-yoffset
                    if (status[s] == "yes"):
                        if (shift and pixshift):
                            xshift=round(DX[s])-xoffset
                            yshift=round(DY[s])-yoffset
                        else:
                            xshift=DX[s]-xoffset
                            yshift=DY[s]-yoffset
                    s=re.sub("^" + outpfx,"",s)
                    [root,ext]=os.path.splitext(s)
                    iraf.imshift.shifts=""
                    iraf.imshift.interp="poly3"
                    iraf.imshift.boundary="constant"
                    iraf.imshift.constant=border
                    check_exist("objmask_" + root + ".pl","w",yes)
                    if (verbose):
                        print "Shifting %s by %f,%f for objmask_%s.pl" % (maskname,xshift,yshift,root)

                    # shift mask
                    check_exist("_objmask_" + root + ".pl","w",yes)
                    iraf.imshift(maskname,"_objmask_"+root+".pl",xshift=xshift,yshift=yshift)

                    # trim
                    xmin=1
                    imgets(root,"i_naxis1")
                    xmax=float(imgets.value)
                    ymin=1
                    imgets(root,"i_naxis2")
                    ymax=float(imgets.value)
                    s="[1:%d,1:%d]" % (int(xmax),int(ymax))
                    iraf.imcopy("_objmask_"+root+s,"objmask_"+root+".pl")
                    update_head(root+ext,"OBJMASK","objmask_"+root+".pl","Name of object mask")
                    
                    check_exist("_objmask_" + root + ".pl","w",yes)
                    

            # trim
            if (doexpand and finaltrim):
                if (verbose):
                    print "Trimming final image"
                imgets(final_name,"i_naxis1")
                if (XMIN < 1):
                    XMIN=1
                if (XMAX > int(imgets.value)):
                    XMAX=int(imgets.value)
                imgets(final_name,"i_naxis2")
                if (YMIN < 1):
                    YMIN=1
                if (YMAX > int(imgets.value)):
                    YMAX=int(imgets.value)

                s="[%d:%d,%d:%d]" % (XMIN,XMAX,YMIN,YMAX)
                iraf.imcopy(final_name+s,final_name)
                iraf.imcopy(final_name+"_exp.pl"+s,final_name+"_exp.pl")


        
######################################################################

def load_defaults(instrument, datadir, verbose=globver):
    """ load instrument defaults
    defaults are stored in directory datadir
    """

    def_file=os.path.join(datadir,instrument + ".defaults")
    check_exist(def_file,"r")
    if (verbose):
        print "Loading defaults for %s from file %s" % (instrument,def_file)
    execfile(def_file)

######################################################################
def xfer_masks(inlist,outpfx,origkey,key1,key2,key3,clobber=globclob,verbose=globver):
    """Takes images in inlist
    multiplies masks together, specified by values of key[1-3] in headers
    transfers final mask to file specified by origkey"""

    check_exist(inlist,"r")
    files=getlines(inlist)
    for image in files:
        if (len(image)==0 or image.isspace()):
            continue
        s=image.split()
        image=s[0]
        check_exist(image,"r")
        inf=pyfits.open(image)
        hdr=inf[0].header
        mask1=hdr.get(key1,"null")
        mask2=hdr.get(key2,"null")
        mask3=hdr.get(key3,"null")
        orig=hdr.get(origkey,"null")
        inf.close()

        if (verbose):
            print "Will combine masks %s, %s, %s from file %s for file %s" % (mask1,mask2,mask3,image,orig)
        
        tmask=iraf.mktemp("xfer") + ".pl"
        tmask2=iraf.mktemp("xfer") + ".pl"
        if (mask1 != "null"):
            check_exist(mask1,"r")
        if (mask2 != "null"):
            check_exist(mask2,"r")
        if (mask3 != "null"):
            check_exist(mask3,"r")
        if (mask1 != "null"):
            check_exist(mask1,"r")
            if (mask2 != "null"):
                iraf.imarith(mask1,"+",mask2,tmask,verbose=verbose)
            else:
                iraf.imcopy(mask1,tmask,verbose=verbose)
            if (mask3 != "null"):
                iraf.imarith(tmask,"+",mask3,tmask2,verbose=verbose)
            else:
                iraf.imcopy(tmask,tmask2,verbose=verbose)
        else:
            if (mask2 != "null"):
                if (mask3 != "null"):
                    iraf.imarith(mask2,"+",mask3,tmask2,verbose=verbose)
                else:
                    iraf.imcopy(mask2,tmask2,verbose=verbose)
            else:
                if (mask3 != "null"):
                    iraf.imcopy(mask3,tmask2,verbose=verbose)
                else:
                    print "No masks exist for image: %s" % image
                    sys.exit(1)
        [final,junk]=os.path.splitext(image)
        final=outpfx + final + ".pl"
        check_exist(final,"w",clobber)
        iraf.imcalc(input=tmask2,output=final,equals="if im1 .ge. 1 then 1 else 0",verbose=verbose)
        # delete temporary files
        check_exist(tmask,"w",yes)
        check_exist(tmask2,"w",yes)
        if (orig != "null"):
            update_head(orig,"BPM",final)

######################################################################
def subsky(inlist,outpfx,skykey,stat,lower,upper,clobber=globclob,verbose=globver):
    """ subtract the sky from the images in inlist
        sky is estimated by a constant value
        """

    check_exist(inlist,"r")
    files=getlines(inlist)
    for image in files:
        if (len(image)==0 or image.isspace()):
            continue
        check_exist(image,"r")
        iraf.iterstat(image=image,nsigrej=5,maxiter=10,prin=no,lower=lower,upper=upper,verbose=no)
        if (stat == "mean"):
            skyval=float(iraf.iterstat.mean)
        if (stat == "mode"):
            skyval=float(iraf.iterstat.valmode)
        if (stat == "median"):
            skyval=float(iraf.iterstat.median)
        if (verbose):
            print "Sky for %s is %f" % (image,skyval)
        newimage=outpfx + image
        check_exist(newimage,"w",clobber)
        iraf.imarith(image,"-",skyval,newimage,verbose=verbose)
        try:
            if (len(skykey) > 0 or not (skykey.isspace())):
                iraf.hedit(newimage,fields=skykey,value=skyval,add=yes,verify=no,show=no)
        except:
            x=1
    
######################################################################
def makemask(inlist, prefix, headlist, subsample, filtsize, nsmooth, threshtype,
             nsigthresh,constthresh,ngrow, statsec, checklimits, zmin,zmax,clobber=globclob,verbose=globver):
    """Create object masks for input images.
    
    Image is sky subtracted, either using a constant value, or optionally by
    median filtering the image on a specified spatial scale and subtracting filtered
    image.   If the desired filtering scale is large, the user may want to first
    subsample the image to a smaller size before filtering in order to speed up the
    calculation -- the resulting sky frame is then block replicated back up to the
    original image size before subtracting.
    
    After sky subtraction, the a threshold is applied to the image after optional
    boxcar smoothing.  The threshold may be specified in terms of a number of sky sigma
    above the median sky level as measured using iterative sigma rejection, or as a 
    specified constant number of ADU above the sky level.  The resulting thresholded image
    is turned into a mask with "sky" set to 0 and "objects" set to 1.  The user
    may "grow" additional rings of pixels around masked regions to increase their area.
    
    Finally, the resulting mask may optionally be recorded in the input image header
    using the keyword BPM.  The user may in fact create a mask from one input image and
    add it into the header of another by specifying different lists for inlist and headlist.
    
    Mark Dickinson -- 1992,1993.  
    Latest revision 18 Aug 1993.
    DLK 2002-12-17"""

    # Check validity of boxcar smoothing scale.
    if (nsmooth > 0):
        if (2*int(nsmooth/2) == nsmooth):
            # if nsmooth is even...
            print "parameter nsmooth should be odd."
            return
    # Notify user of effective median filter scale.
    if (subsample > 1 and filtsize > 0):
        realsize = subsample * filtsize
        if (verbose):
            print "After subsampling image, effective median filter scale will be %d" % realsize

    # Size for boxcar smoothing in growing stage.
    nbox = (2*ngrow + 1)
    narea = nbox * nbox

    # Expand input file list.
    infiles=glob.glob(inlist)

    # If parameter headlist != "", we will update image headers with bad pixel list.
    if (headlist != ""):
        headfiles=glob.glob(headlist)
    else:
        headfiles=""

    # Loop through input files.
    for img in infiles:
        # Strip extension off filename if present.
        [img,ext]=os.path.splitext(img)
        if (verbose):
            print "   Working on image %s" % img
        if (verbose):
            print "      Subtracting local sky"

        # If filtsize > 0, median filter image to produce local sky,
        # then subtract that from the image.
        check_exist("_m" + img,"w",yes)
        check_exist("_m" + img + ".fits","w",yes)
        if (filtsize > 0):
            # If subsample > 1, block average input image to smaller size.
            if (subsample > 1):
                if (verbose):
                    print "         Block averaging"
                check_exist("_blk_" + img,"w",yes)
                iraf.blkavg(img,"_blk_" + img, subsample, subsample, option="average")
                workimg="_blk_" + img
            else:
                workimg=img
            check_exist("_f" + workimg,"w",yes)
            check_exist("_f" + workimg + ".fits","w",yes)

            # First, check limits of data range.  Wildly large (positive or negative) data values
            # will screw up fmedian calculation unless checklimits = yes and zmin and zmax are set
            # to appropriate values, e.g. zmin=-32768, zmax=32767.   Note that the old fmedian
            # bug which occurred if zmin=hmin and zmax=hmax has been fixed in IRAF version 2.10.1.
            if (verbose):
                print "        Median filtering sky."
            if (checklimits):
                iraf.minmax(workimg,force=yes,update=yes,verbose=no)
                if (verbose):
                    print "     Data minimum = %f maximum = %f" % (iraf.minmax.minval,iraf.minmax.maxval)
            if (checklimits and (iraf.minmax.minval < zmin or iraf.minmax.maxval > zmax)):
                if (iraf.minmax.minval < zmin):
                    dmin=zmin
                else:
                    dmin=iraf.minmax.minval
                if (iraf.minmax.maxval > zmax):
                    dmax = zmax
                else:                    
                    dmax = iraf.minmax.maxval
                if (verbose):
                    print "     Truncating data range %f to %f" % (dmin,dmax)
                iraf.fmedian(workimg,"_f" + workimg,xw=filtsize,yw=filtsize,
                             boundary="nearest",hmin=-32768,hmax=32767,
                             zmin=dmin,zmax=dmax,Stdout=1)
            else:
                iraf.fmedian(workimg,"_f"+workimg,xw=filtsize,yw=filtsize,boundary="nearest",
                        hmin=-32768,hmax=32767,zmin=INDEF,zmax=INDEF,Stdout=1)

            # If we have block averaged, block replicate median filtered image back to original size.
            if (subsample > 1):
                if (verbose):
                    print  "        Block reproducing."
                check_exist("_f" + img,"w",yes)
                iraf.blkrep ("_f"+workimg,"_f"+img,subsample,subsample)
                iraf.imgets(img,"i_naxis1")
                ix = int(imgets.value)
                iraf.imgets(img,"i_naxis2")
                iy = int(imgets.value)
                cutsec="[1:%d,1:%d]" % (ix,iy)
                iraf.imcopy ("_f"+img+cutsec,"_f"+img,ver=no)
            iraf.imarith(img,"-","_f"+img,"_m"+img)
        else:
            # ...or, just copy the image to a working mask frame.
            iraf.imcopy (img,"_m"+img)

        # Calculate image statistics to determine median sky level and RMS noise.
        if (verbose):
            print "     Computing sky statistics."
        iraf.iterstat.lower="INDEF"
        iraf.iterstat.upper="INDEF"
        iraf.iterstat("_m"+img+statsec,verbose=no,prin=no,Stdout=1)
        if (nsmooth > 0):
            if (verbose):
                print "     Smoothing image before thresholding."
            iraf.boxcar("_m"+img,"_m"+img,nsmooth,nsmooth)

        # Calculate threshold.
        if (threshtype == "constant"):
            thresh = iraf.iterstat.median + constthresh
        else:
            thresh = iraf.iterstat.median + nsigthresh * iraf.iterstat.sigma

        print "     Thresholding image at level %f" % thresh

        # Apply threshold to image, setting "objects" to 1 and "sky" to 0.
        # The order of the imreplace statments must be different if threshold is greater
        # than or less than 1.
        if (thresh >= 0.):
            iraf.imreplace ("_m"+img+ext,0.,upper=thresh,lower=INDEF)
            iraf.imreplace ("_m"+img+ext,1.,lower=thresh,upper=INDEF)
        elif (thresh < 0.):
            iraf.imreplace ("_m"+img+ext,1.,lower=thresh,upper=INDEF)
            iraf.imreplace ("_m"+img+ext,0.,upper=thresh,lower=INDEF)

        # Copy mask to .pl file.
        if (verbose):
            print "     Saving mask as %s%s.pl" % (prefix,img)
        check_exist(prefix + img + ".pl","w",clobber)
        iraf.imcopy("_m"+img,prefix+img+".pl",ver=no)

        # If desired, grow rings around masked objects.
        if (ngrow > 0):
            if (verbose):
                print "     Growing rings around masked objects."
            iraf.imarith(prefix+img+".pl","*",narea,prefix+img+".pl")
            iraf.boxcar(prefix+img+".pl",prefix+img+".pl",nbox,nbox)
            iraf.imreplace(prefix+img+".pl",1,lower=1,upper=INDEF)

        # Record mask name into BPM keyword in image header of files specified 
        # by inlist2.
        if (headlist != ""):
            for headfile in headfiles:
                if (verbose):
                    print "     Recording pixel mask into header of image %s" % headfile
                iraf.hedit(headfile,"BPM",prefix+img+".pl",add=yes,ver=no,update=yes)

        # Clean up.
        if (filtsize > 0):
            iraf.imdelete("_f"+img+ext,ver=no)
            if (subsample > 1 ):
                iraf.imdelete ("_blk_"+img+ext,ver=no)
                iraf.imdelete ("_f"+workimg+ext,ver=no)

        iraf.imdelete("_m"+img+ext,ver=no)

######################################################################
def do_irredux(inlist,darklist,makebpm,make_bpm,bpmsource,bpm,makedark,make_dark,darksubtr,dark_subtr,
               makeflat,make_flat,flatdiv,flat_divide,fixpix,irfixpix,crzap,ircrzap,makemask,mosaic,
               irmosaic,pass1,pass2,defaults,instrument,datadir,clobber=yes,verbose=yes):
    """ Master task to control irredux
        set parameters/defaults
        calls all other tasks
        """

    # defaults
    objkey="OBJECT" 
    filtkey="FILTER"
    nitekey="DATE-OBS"
    timekey="MJD-OBS"
    expkey="TINT"   
    coaddkey="FRMCOADD"
    rakey="RA"      
    deckey="DEC"    
    rotkey="ROTPOSN"
    pixscale=0.15

    # Check for file list inputs
    check_exist(inlist,"r")
    check_exist(darklist,"r")

    # read in key definition
    keydef=os.path.join(datadir, instrument + ".keys")
    check_exist(keydef,"r")
    if (verbose):
        print "Reading key definitions from: %s" % keydef
    lines=getlines(keydef)
    for line in lines:            
        [s1,s2]=line.split("=")
        # evaluate expression
        s="%s=\"%s\"" % (s1,s2)
        exec(s)
        if (s1.find("pixscale") >= 0):
            iraf.irmosaic.pixscale=float(s2)

    # read in defaults
    if (defaults):
        load_defaults(instrument,datadir,verbose)

    if (pass1):
        if (verbose):
            print "First pass..."
        if (makebpm):
            iraf.make_bpm.bpmkey="BPM0"
            iraf.make_bpm.expkey=expkey
            iraf.make_bpm.coaddkey=coaddkey
            iraf.make_bpm.clobber=clobber
            iraf.make_bpm.verbose=verbose
            iraf.make_bpm.modlist=inlist
            if (bpmsource == "dark"):
                iraf.make_bpm(darklist,inlist + ".pl")
            elif (bpmsource == "data"):
                iraf.make_bpm(inlist,inlist + ".pl")
        else:
            # if not making a BPM, but bpm is defined,
            # update the header with that
            if (bpm != ""):
                check_exist(bpm,"r")
                if (verbose):
                    print "Adding BPM0=%s to files in %s" % (bpm,inlist)
                lines=getlines(inlist)
                for line in lines:
                    if (re.search("\S+",line)):
                        update_head(line,"BPM0",bpm,"Original BPM")

        if (makedark):
	    iraf.make_dark.expkey=expkey
	    iraf.make_dark.coaddkey=coaddkey
	    iraf.make_dark.reject="minmax"
	    iraf.make_dark.verbose=verbose
	    iraf.make_dark.clobber=clobber
	    iraf.make_dark.verbose=verbose
	    iraf.make_dark.clobber=clobber
	    iraf.make_dark(objfile=inlist,darkfile=darklist)	

        if (darksubtr):
            check_exist(inlist + ".dark","r")
	    dark_subtr.verbose=verbose
	    dark_subtr.clobber=clobber
	    dark_subtr(inlist=inlist+".dark")
            
        if (makeflat):
            check_exist(inlist+".dark.sub","r")
	    iraf.make_flat.filtkey=filtkey
	    iraf.make_flat.expkey=expkey
	    iraf.make_flat.coaddkey=coaddkey
	    # this will copy BPM0 -> BPM, and use that for imcombine
	    iraf.make_flat.bpmkey="BPM0"
	    iraf.make_flat.clobber=clobber
	    iraf.make_flat.verbose=verbose
	    iraf.make_flat.zero="none"
            iraf.make_flat.prefix="flat"
	    iraf.make_flat(inlist=inlist+".dark.sub")

        if (flatdiv):
            check_exist(inlist+".dark.sub.flat","r")        
            iraf.flat_div.verbose=verbose
	    iraf.flat_div.clobber=clobber
            iraf.flat_div.outpfx="f"
	    iraf.flat_div(inlist=inlist+".dark.sub.flat")
            
        if (fixpix != "no"):
            check_exist(inlist+".dark.sub.flat.div","r")
            if (fixpix == "yes"):
                iraf.irfixpix.verbose=verbose
                iraf.irfixpix.clobber=clobber
                # fixes those pixels in BPM, which should have been set by
                # make_bpm (-> BPM0)
                # then copied by make_flat (->BPM)
                iraf.irfixpix(inlist=inlist+".dark.sub.flat.div")
            else:
                # skip BP fixing
                # just copy file to .fix
                check_exist(inlist+".dark.sub.flat.div.fix","w",yes)
                if (verbose):
                    print "Skipping BP fixing"
                    print "Copying %s to %s" % (inlist+".dark.sub.flat.div",inlist+".dark.sub.flat.div.fix")
                s="/bin/cp %s %s" % (inlist+".dark.sub.flat.div",inlist+".dark.sub.flat.div.fix")
                os.system(s)

        if (crzap != "no"):
            check_exist(inlist+".dark.sub.flat.div.fix","r")
            if (crzap == "yes"):
                iraf.ircrzap.verbose=verbose
                iraf.ircrzap.clobber=clobber
                # set this so it can update the static BPM if needed
                iraf.ircrzap.bpmkey="BPM0"	
                # backup the original BPM
                if (os.path.exists(inlist+".pl")):
                    if (verbose):
                        print "Backup up BPM file %s.pl" % inlist
                    check_exist(inlist + ".back.pl","w",yes)
                    iraf.copy(inlist+".pl",inlist+".back.pl")

                iraf.ircrzap(inlist=inlist+".dark.sub.flat.div.fix")
            else:
                # just copy file to .zap
                check_exist(inlist+".dark.sub.flat.div.fix.zap","w",yes)
                if (verbose):
                    print "Skipping CR zapping"
                    print "Copying %s to %s" % (inlist+".dark.sub.flat.div.fix",inlist+".dark.sub.flat.div.fix.zap")
                s="/bin/cp %s %s" % (inlist+".dark.sub.flat.div.fix",inlist+".dark.sub.flat.div.fix.zap")
                os.system(s)

        if (mosaic):
            check_exist(inlist + ".dark.sub.flat.div.fix.zap","r")		
            lines=getlines(inlist + ".dark.sub.flat.div.fix.zap")
            ref=lines[0].split()[0]

	    iraf.irmosaic.verbose=verbose
	    iraf.irmosaic.clobber=clobber
	    # this will output a thing into OBJMASK
	    iraf.irmosaic.makemask=makemask
	    iraf.irmosaic(ref_pic=ref,picfile=inlist+".dark.sub.flat.div.fix.zap")
            
	if (makemask):
	    # if we've just made a mask and we have one from before, 
	    # we have to combine
	    # if not, we still have to transfer the object mask to the original image
            if (verbose):
		print "\nFixing and combining object masks"
            check_exist(inlist + ".dark.sub.flat.div.fix.zap","r")		
            iraf.xfer_masks.outpfx="mask_"
            iraf.xfer_masks.origkey="NAME1"
            iraf.xfer_masks.key1="BPM0"
            iraf.xfer_masks.key2="CRMASK"
            iraf.xfer_masks.key3="OBJMASK"
            iraf.xfer_masks.clobber=clobber
            iraf.xfer_masks.verbose=verbose
            iraf.xfer_masks(inlist + ".dark.sub.flat.div.fix.zap")
            
    # -----------------------------------
    # end 1st pass
    if (pass2):
        if (verbose):
	    print "Second pass . . ."
        if (makeflat):
	    # for second pass, use object masking
	    # take the masks made from irmosaic
	    # also, make a mask from the combined image
	    # de-shift the 2nd mask to the input image locations
	    # combine it, the individual object masks, and the cr/bad-pix masks
	    # for a final mask
	    # then use that in what follows
            check_exist(inlist + ".dark.sub","r")
	    iraf.make_flat.filtkey=filtkey
	    iraf.make_flat.expkey=expkey
	    iraf.make_flat.coaddkey=coaddkey
	    iraf.make_flat.bpmkey="BPM"
	    iraf.make_flat.clobber=clobber
	    iraf.make_flat.verbose=verbose

	    iraf.make_flat.zero="none"
            check_exist(inlist + ".dark.sub2","w",clobber)
	    iraf.copy(inlist+".dark.sub",inlist+".dark.sub2")
	    iraf.make_flat(inlist=inlist+".dark.sub2",prefix="mflat")

        if (flatdiv):
            check_exist(inlist + ".dark.sub2.flat","r")
	    iraf.flat_div.verbose=verbose
	    iraf.flat_div.clobber=clobber
	    iraf.flat_div(inlist=inlist+".dark.sub2.flat",outpfx="f2")

        if (fixpix != "no"):            
            check_exist(inlist + ".dark.sub2.flat.div","r")
            if (fixpix == "yes"):
                # copy new mask -> BPM1, BPM0 -> BPM
                lines=getlines(inlist + ".dark.sub2.flat.div")
                for line in lines:
                    inf=pyfits.open(line,"update")
                    inf[0].header.update("BPM1",inf[0].header.get("BPM",""))
                    inf[0].header.update("BPM",inf[0].header.get("BPM0",""))                    
                    inf.close()
	       
                iraf.irfixpix.verbose=verbose
                iraf.irfixpix.clobber=clobber
                # fixes those pixels in BPM, which should have been set by
                # make_bpm (-> BPM0)
                # then copied by above routine (->BPM)
                iraf.irfixpix(inlist=inlist+".dark.sub2.flat.div")
            else:
                # skip BP fixing
                # just copy file to .fix
                check_exist(inlist+".dark.sub2.flat.div.fix","w",yes)
                if (verbose):
                    print "Skipping BP fixing"
                    print "Copying %s to %s" % (inlist+".dark.sub2.flat.div",inlist+".dark.sub2.flat.div.fix")
                s="/bin/cp %s %s" % (inlist+".dark.sub2.flat.div",inlist+".dark.sub2.flat.div.fix")
                os.system(s)

        if (crzap != "no"):
            check_exist(inlist+".dark.sub2.flat.div.fix","r")
            if (crzap == "yes"):
                iraf.ircrzap.verbose=verbose
                iraf.ircrzap.clobber=clobber
                # set this so it can update the static BPM if needed
                iraf.ircrzap.bpmkey="BPM0"	
                # backup the original BPM
                if (os.path.exists(inlist+".pl")):
                    if (verbose):
                        print "Backup up BPM file %s.pl" % inlist
                    check_exist(inlist + ".back2.pl","w",yes)
                    iraf.copy(inlist+".pl",inlist+".back2.pl")

                iraf.ircrzap(inlist=inlist+".dark.sub2.flat.div.fix")
            else:
                # just copy file to .zap
                check_exist(inlist+".dark.sub2.flat.div.fix.zap","w",yes)
                if (verbose):
                    print "Skipping CR zapping"
                    print "Copying %s to %s" % (inlist+".dark.sub2.flat.div.fix",inlist+".dark.sub2.flat.div.fix.zap")
                s="/bin/cp %s %s" % (inlist+".dark.sub2.flat.div.fix",inlist+".dark.sub2.flat.div.fix.zap")
                os.system(s)

        if (mosaic):
            check_exist(inlist + ".dark.sub2.flat.div.fix.zap","r")		
            lines=getlines(inlist + ".dark.sub2.flat.div.fix.zap")
            ref=lines[0].split()[0]

	    iraf.irmosaic.verbose=verbose
	    iraf.irmosaic.clobber=clobber
	    # don't output a thing into OBJMASK
	    iraf.irmosaic.makemask=no
	    iraf.irmosaic(ref_pic=ref,picfile=inlist+".dark.sub2.flat.div.fix.zap")

######################################################################
def do_irpanic(inlist,bpm,addbpm,bpmdivide,
               makesky,make_flat,skysub,dark_subtr,makemask,mosaic,
               irmosaic,pass1,pass2,defaults,instrument,datadir,clobber=yes,verbose=yes):
    """ Master task to control irredux
        set parameters/defaults
        calls all other tasks
        """

    # defaults
    objkey="OBJECT" 
    filtkey="FILTER"
    nitekey="DATE-OBS"
    timekey="MJD-OBS"
    expkey="TINT"   
    coaddkey="FRMCOADD"
    rakey="RA"      
    deckey="DEC"    
    rotkey="ROTPOSN"
    pixscale=0.15

    check_exist(inlist,"r")

    # read in key definition
    keydef=os.path.join(datadir, instrument + ".keys")
    check_exist(keydef,"r")
    if (verbose):
        print "Reading key definitions from: %s" % keydef
    lines=getlines(keydef)
    for line in lines:            
        [s1,s2]=line.split("=")
        # evaluate expression
        s="%s=\"%s\"" % (s1,s2)
        exec(s)
        if (s1.find("pixscale") >= 0):
            iraf.irmosaic.pixscale=float(s2)

    # read in defaults
    if (defaults):
        load_defaults(instrument,datadir,verbose)

    if (pass1):
        if (verbose):
            print "First pass..."

        if (addbpm):
            # if not making a BPM, but bpm is defined,
            # update the header with that
            if (bpm != ""):
                check_exist(bpm,"r")
                if (verbose):
                    print "Adding BPM0=%s to files in %s" % (bpm,inlist)
                lines=getlines(inlist)
                for line in lines:
                    if (re.search("\S+",line)):
                        update_head(line,"BPM0",bpm,"Original BPM")

        if (bpmdivide):
            # divide by inverse of BPM
            check_exist(bpm,"r")
            check_exist("bpm_inv.pl","w")
            iraf.imcalc(bpm,"bpm_inv.pl","if (im1 == 1) then 0 else 1")
            check_exist(inlist,"r")
            check_exist(inlist+".flat","w")
            outfile=open(inlist+".flat","w")            
            lines=getlines(inlist)
            for line in lines:
                if (re.search("\S+",line)):
                    check_exist("f"+line,"w",clobber)
                    iraf.imarith(line,"/","bpm_inv.pl","f"+line,divzero=0)
                    outfile.write("f"+line+"\n")
            outfile.close()
                        
        if (makesky):
            check_exist(inlist+".flat","r")
	    iraf.make_flat.filtkey=filtkey
	    iraf.make_flat.expkey=expkey
	    iraf.make_flat.coaddkey=coaddkey
	    # this will copy BPM0 -> BPM, and use that for imcombine
	    iraf.make_flat.bpmkey="BPM0"
	    iraf.make_flat.clobber=clobber
	    iraf.make_flat.verbose=verbose
	    iraf.make_flat.zero="none"
            iraf.make_flat.prefix="sky"
            iraf.make_flat.donorm=no
	    iraf.make_flat(inlist=inlist+".flat")
            check_exist(inlist+"flat.sky","w")
            os.rename(inlist+".flat.flat",inlist+".flat.sky")

        if (skysub):
            check_exist(inlist+".flat.sky","r")        
            iraf.dark_subtr.outpfx="s"
	    iraf.dark_subtr(inlist=inlist+".flat.sky")
                    
        if (mosaic):
            check_exist(inlist + ".flat.sky.sub","r")		
            lines=getlines(inlist + ".flat.sky.sub")
            ref=lines[0].split()[0]

	    iraf.irmosaic.verbose=verbose
	    iraf.irmosaic.clobber=clobber
	    # this will output a thing into OBJMASK
	    iraf.irmosaic.makemask=makemask
	    iraf.irmosaic(ref_pic=ref,picfile=inlist+".flat.sky.sub")
            
	if (makemask):
	    # if we've just made a mask and we have one from before, 
	    # we have to combine
	    # if not, we still have to transfer the object mask to the original image
            if (verbose):
		print "\nFixing and combining object masks"
            check_exist(inlist + ".flat.sky.sub","r")		
            iraf.xfer_masks.outpfx="mask_"
            iraf.xfer_masks.origkey="NAME0"
            iraf.xfer_masks.key1="BPM0"
            iraf.xfer_masks.key2="CRMASK"
            iraf.xfer_masks.key3="OBJMASK"
            iraf.xfer_masks.clobber=clobber
            iraf.xfer_masks.verbose=verbose
            iraf.xfer_masks(inlist + ".flat.sky.sub")
            
    # -----------------------------------
    # end 1st pass
    if (pass2):
        if (verbose):
	    print "Second pass . . ."
        if (makesky):
	    # for second pass, use object masking
	    # take the masks made from irmosaic
	    # also, make a mask from the combined image
	    # de-shift the 2nd mask to the input image locations
	    # combine it, the individual object masks, and the cr/bad-pix masks
	    # for a final mask
	    # then use that in what follows
            check_exist(inlist + ".flat","r")

            iraf.make_flat.filtkey=filtkey
	    iraf.make_flat.expkey=expkey
	    iraf.make_flat.coaddkey=coaddkey
	    # this will copy BPM0 -> BPM, and use that for imcombine
	    iraf.make_flat.bpmkey="BPM"
	    iraf.make_flat.clobber=clobber
	    iraf.make_flat.verbose=verbose
	    iraf.make_flat.zero="none"
            iraf.make_flat.prefix="msky"
            iraf.make_flat.donorm=no

            check_exist(inlist + ".flat2","w",clobber)
	    iraf.copy(inlist+".flat",inlist+".flat2")
            
	    iraf.make_flat(inlist=inlist+".flat2")
            check_exist(inlist+"flat2.sky","w")
            os.rename(inlist+".flat2.flat",inlist+".flat2.sky")

        if (skysub):
            check_exist(inlist+".flat2.sky","r")        
            iraf.dark_subtr.outpfx="s2"
	    iraf.dark_subtr(inlist=inlist+".flat2.sky")

        if (mosaic):
            check_exist(inlist + ".flat2.sky.sub","r")		
            lines=getlines(inlist + ".flat2.sky.sub")
            ref=lines[0].split()[0]

	    iraf.irmosaic.verbose=verbose
	    iraf.irmosaic.clobber=clobber
	    # don't output a thing into OBJMASK
	    iraf.irmosaic.makemask=no
	    iraf.irmosaic(ref_pic=ref,picfile=inlist+".flat2.sky.sub")

######################################################################


######################################################################
# Utility routines                    
######################################################################

def check_limit(xval,xlim,lower=no,upper=no):

    """ check_limit(xval,xlim,dxval,lower=no,upper=no)
    determine if xval violates the upper/lower limit set
    by xlim.  if so, sets dxval to the difference."""

    dxval=0

    if (not lower and not upper):
        print "Must set lower or upper bound in check_limits"
        sys.exit(1)

    if (lower and xval<xlim):
        dxval=xlim-xval
        xval=xlim

    if (upper and xval>xlim):
        dxval=xval-xlim
        xval=xlim

    return (xval,dxval)

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
    for i in range(0,len(lines)):
        if (re.search("\S+",lines[i])):
            lines2.append(lines[i].rstrip())
    return lines2

######################################################################
def expand_img(input,output,scale=2,border=0,
               nodata=0,clobber=globclob,verbose=globver):

    """ takes input image, puts it into central portion of output image, 
    where size(output)=scale*size(input)
    fills remaining portion with value=border"""

    if (not (os.path.exists(input) or os.path.exists(input+".fits"))):
        print "Couldn't open input image %s" % input
        sys.exit(1)
    if (os.path.exists(input+".fits")):
        input+=".fits"
    check_exist(output,"w",clobber)

    scale=int(scale)
    if (scale<1):
        print "Bad scaling in expand_img, assuming scale=2"
        scale=2
        
    fimg=pyfits.open(input)
    hdr=fimg[0].header
    D=fimg[0].data
    (iny,inx)=D.shape

    if (verbose):
        print "Image %s is %d x %d" % (input,inx,iny)
        print "Output image %s will be %d x %d" % (output,scale*inx,scale*iny)

    outx=scale*inx
    outy=scale*iny
    #DX=0*numarray.resize(D,(outy,outx))+border
    DX=0*numpy.resize(D,(outy,outx))+border
    xstart=outx/2-inx/2
    xend=xstart+inx
    ystart=outy/2-iny/2
    yend=ystart+iny
    # print "Xstart=%d Xend=%d Ystart=%d Yend=%d" % (xstart,xend,ystart,yend)
    if (not nodata):
        DX[ystart:yend,xstart:xend]=D

    fimg[0].data=DX
    fimg[0].header.update('EXPND','Expanded by %d from %s' % (scale,input))
    del fimg[0].header['BPM']
    fimg.writeto(output)

######################################################################
def getregions(xpaid,goodcolor=""):
    """regions=getregions(xpaid,goodcolor)
       gets the regions from the XPA point indicated by xpaid
       currently gets only circles or boxes
       if (goodcolor) is specified, will only include those with color=goodcolor
       returns [[x0,y0],[x1,y1],...]
       """

    shape1="circle"
    shape2="box"

    s="xpaget %s regions" % xpaid
    inf=os.popen(s,"r")
    lines=inf.readlines()
    inf.close()
    regions=[]
    defcolor="green"
    for i in range(0,len(lines)):
        lines[i]=lines[i].rstrip()
        if (lines[i].find("global") >= 0):
            j=lines[i].find("color=")
            k=lines[i].find(" ",j)
            defcolor=lines[i][j+6:k]
        if (lines[i].count(shape1) > 0 or lines[i].count(shape2) > 0):
            i1=lines[i].find("(")+1
            i2=lines[i].find(")")
            s=lines[i][i1:i2]
            if (lines[i].count(shape1) > 0):
                [x,y,r]=s.split(",")
            else:
                [x,y,r1,r2,angle]=s.split(",")
            x=float(x)
            y=float(y)
            color=defcolor
            j=lines[i].find("color=")
            if (j >= 0):
                k=lines[i].find(" ",j)
                if (not k > j+6):
                    k=len(lines[i])
                color=lines[i][j+6:k]
            if ((not re.search("\S+",goodcolor)) or (color == goodcolor)):
                regions.append([x,y])
    return regions

######################################################################
def clean_files(root):
    """ remove all temporary files associated with processing
    """

    # get the dark-subtracted files:
    del_files("d" + root + "*.fits")
    # flat-fielded files:
    del_files("fd" + root + "*.fits")
    # bp-fixed files:
    del_files("bfd" + root + "*.fits")
    # cr-fixed files:
    del_files("c*fd" + root + "*.fits")
    # shifted files:
    del_files("S*fd" + root + "*.fits")    
    # pass 2:
    # flat-fielded files:
    del_files("f2d" + root + "*.fits")
    # bp-fixed files:
    del_files("bf2d" + root + "*.fits")
    # cr-fixed files:
    del_files("c*f2d" + root + "*.fits")
    # shifted files:
    del_files("S*f2d" + root + "*.fits")    
    # data files:
    del_files("*fd" + root + "*.fits.tform")    
    # darks
    del_files("dark*.fits")
    # flats
    del_files("flat*.fits*")
    del_files("mflat*.fits*")
    # masks
    del_files("mask*" + root + "*.pl")
    del_files("crmask*" + root + "*.pl")
    del_files("objmask*" + root + "*.pl")
    # others
    del_files("_temp.fits")
    del_files("sm__temp.fits")
    
######################################################################
def del_files(s):
    """ delete the files matching the string in s
    """

    F=glob.glob(s)
    if (len(F) > 0):
        r=raw_input("Delete %s [Y/n]? " % (s)).upper()
        if (r != "N" and r != "NO"):
            for f in F:
                os.remove(f)

######################################################################
def sex2dec(ra,dec):
    """[RA,DEC]=sex2dec(ra,dec)
    converts sexadecimal to decimal
    ra/dec are hh:mm:ss/dd:mm:ss
    RA/DEC are in hours/degrees
    """

    [h,m,s]=ra.split(":")
    RA=float(h)+(float(m)/60.0)+(float(s)/3600)
    [d,m,s]=dec.split(":")
    s=re.sub(",","",s)
    d=abs(float(d))
    DEC=d+(float(m)/60)+(float(s)/3600)
    if (dec.find("-") >= 0):
        DEC=DEC*-1
    return [RA,DEC]
    
######################################################################

def find_snr(image,x,y,statsec,clobber=globclob,verbose=globver):
    """finds noise=sigma of bkg through iterstat
    signal = result of gaussian fit"""

    check_exist(image,"r")
    # get dimesions
    xmin=1
    ymin=1
    imgets(image,"i_naxis1")
    xmax=int(imgets.value)
    imgets(image,"i_naxis2")
    ymax=int(imgets.value)
    if (x < xmin or x > xmax):
        print "X parameter out of range"
        return -1
    if (y < ymin or y > ymax):
        print "Y parameter out of range"
        return -1

    # now find signal
    check_exist("_temp.coord","w",yes)
    OUT=open("_temp.coord","w")
    OUT.write("%f %f\n" % (x,y))
    OUT.close()
    
    check_exist(image + ".temppsf","w",yes)

    iraf.fitpsf.box=10
    iraf.fitpsf.function="radgauss"
    iraf.fitpsf.interac=no
    iraf.fitpsf.verify=no
    iraf.fitpsf.verbose=verbose
    iraf.fitpsf.coords="_temp.coord"
    iraf.fitpsf.output=image+".temppsf"
    iraf.fitpsf.maxiter=50
    iraf.fitpsf.nrej=0
    iraf.fitpsf.krej=3
    iraf.fitpsf.mkbox=no
    x=iraf.fitpsf(image=image,box=10,Stdout=1,Stderr=1)
    os.remove("_temp.coord")

    # read back results
    lines=getlines(image + ".temppsf")
    i=0
    while (re.search("^\#",lines[i])):
        i+=1
    i+=1
    s=lines[i].split()
    signal=float(s[3])
    if (signal <= 0):
        signal=1
    s=lines[i+1].split()
    noise=math.sqrt(float(s[3])**2+float(s[4])**2)
    if (noise <= 0):
        noise=1
    if (verbose):
        print "Signal: %f, Noise: %f, SNR: %f" % (signal,noise,signal/noise)

    os.remove(image + ".temppsf")
    snr=signal/noise
    return snr

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
def check_xpa(xpaid):
    """ checks to see if the XPA ID is valid
    queries for about response
    """

    s="xpaget %s about >& /dev/null" % xpaid
    result=not os.system(s)

    return result

######################################################################

# location of parameter files
_parfile=pardir + "make_bpm.par"
t=iraf.IrafTaskFactory(taskname="make_bpm",value=_parfile,function=make_bpm)
_parfile=pardir + "make_dark.par"
t=iraf.IrafTaskFactory(taskname="make_dark",value=_parfile,function=make_dark)
_parfile=pardir + "dark_subtr.par"
t=iraf.IrafTaskFactory(taskname="dark_subtr",value=_parfile,function=dark_subtr)
_parfile=pardir + "make_flat.par"
t=iraf.IrafTaskFactory(taskname="make_flat",value=_parfile,function=make_flat)
_parfile=pardir + "flat_divide.par"
t=iraf.IrafTaskFactory(taskname="flat_divide",value=_parfile,function=flat_divide)
_parfile=pardir + "irfixpix.par"
t=iraf.IrafTaskFactory(taskname="irfixpix",value=_parfile,function=irfixpix)
_parfile=pardir + "ircrzap.par"
t=iraf.IrafTaskFactory(taskname="ircrzap",value=_parfile,function=ircrzap)
_parfile=pardir + "irmosaic.par"
t=iraf.IrafTaskFactory(taskname="irmosaic",value=_parfile,function=irmosaic)
_parfile=pardir + "expand_img.par"
t=iraf.IrafTaskFactory(taskname="expand_img",value=_parfile,function=expand_img)
_parfile=pardir + "load_defaults.par"
t=iraf.IrafTaskFactory(taskname="load_defaults",value=_parfile,function=load_defaults)
_parfile=pardir + "xfer_masks.par"
t=iraf.IrafTaskFactory(taskname="xfer_masks",value=_parfile,function=xfer_masks)
_parfile=pardir + "subsky.par"
t=iraf.IrafTaskFactory(taskname="subsky",value=_parfile,function=subsky)
_parfile=pardir + "clean_files.par"
t=iraf.IrafTaskFactory(taskname="clean_files",value=_parfile,function=clean_files)
_parfile=pardir + "makemask.par"
t=iraf.IrafTaskFactory(taskname="makemask",value=_parfile,function=makemask)
_parfile=pardir + "do_irredux.par"
t=iraf.IrafTaskFactory(taskname="do_irredux",value=_parfile,function=do_irredux)
_parfile=pardir + "do_irpanic.par"
t=iraf.IrafTaskFactory(taskname="do_irpanic",value=_parfile,function=do_irpanic)


