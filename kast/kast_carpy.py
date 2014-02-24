
# Global python packages
import pyraf
from pyraf import iraf
import os, shutil, math, subprocess
import pyfits

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
KGRATINGS = {"600/4310": {"dispersion": 1.02, "sreg": "[251:1350,*]",
            "slits": [5, 60, 210, 275, 420, 490, 630, 700, 840, 910, 1050]},
            "300/7500": {"dispersion": 4.60, "sreg": "[291:1390,*]",
            "slits": [5, 60, 210, 275, 420, 490, 630, 700, 840, 910, 1050]}}
KARCLIST = "/Users/cenko/iraf/linelists/deimos_linelist.dat"
KSKYLIST = "/Users/cenko/iraf/linelists/skylines.txt"
CCDRNOISE = 2.69 # CCD read noise (e-)
CCDGAIN = 1.28 # CCD e- per DN
BBIASSEC = "[2061:2140,1:4096]" # Overscan section on blue chip
BFLAT = "Flat-B.fits" # Median combined flat field for blue chip
BSHARP = "Sharp-B.fits" # Similar to flat field, but with sharp edges
RBIASSEC = "[2061:2140,1:4096]" # Overscan section on blue chip
RFLAT = "Flat-R.fits" # Median combined flat field for red chip
RSHARP = "Sharp-R.fits" # Similar to flat field, but with sharp edges

#######################

def blue_preproc(image):

    '''Remove overscan and trim image from blue detector.'''

    grisname = get_head(image, "GRISM_N", extn=0)    
    iraf.ccdproc(image, overscan=yes, trim=yes, biassec=BBIASSEC,
                 trimsec=KGRATINGS[grisname]["sreg"])
    
    return

########################

def red_preproc(image):

    '''Remove overscan and trim image from blue detector.'''

    gratname = get_head(image, "GRATING_N", extn=0)    
    iraf.ccdproc(image, overscan=yes, trim=yes, biassec=BBIASSEC,
                 trimsec=KGRATINGS[gratname]["sreg"])
    
    return

########################

def deimos_redcalib(flats, arc):

    '''Calculate y- and x- distortions, create flat-field, calculate wavelength
    solution, for a given grating and central wavelength.'''

    # Pre-process flats
    for image in flats:
        deimos_preproc(image)

    # Create median combination of flats
    iraf.flatcombine(",".join(["rtd%s_B.fits" % x[6:10] for x in flats]),
                     output=BFLAT, combine="median", reject="avsigclip")
    iraf.flatcombine(",".join(["rtd%s_R.fits" % x[6:10] for x in flats]),
                     output=RFLAT, combine="median", reject="avsigclip")

    # Pre-process arc
    graname = get_head(arc, "GRATENAM", extn=0)
    deimos_preproc(arc)
    barc = "rtd%s_B.fits" % arc[6:10]; rarc = "rtd%s_R.fits" % arc[6:10]

    # Generate flats with sharpened edges
    os.system('efits %s "VTKLaplace(i1,numret=1)" %s' % (BFLAT, BSHARP))
    os.system('efits %s "VTKLaplace(i1,numret=1)" %s' % (RFLAT, RSHARP))

    # Calculate y-distortion from sharpened frame
    os.system('getrect -ydist %s -nx 12 -dy 11 -x 1 -y 1 -peakwin 3 -dump' % BSHARP)
    os.system('getrect -ydist %s -nx 12 -dy 11 -x 1 -y 1 -peakwin 3 -dump' % RSHARP)

    # Copy y rectification to flats
    os.system('copyrect -ydist %s %s' % (BSHARP, BFLAT))
    os.system('copyrect -ydist %s %s' % (RSHARP, RFLAT))
    
    # Find slits
    os.system('findslits -ydist %s -th 0.2 -fwhm 2' % BFLAT)
    os.system('findslits -ydist %s -th 0.2 -fwhm 2' % RFLAT)

    # Correct Slit locations
    update_head(BFLAT, ["NSLITS", "CSECT1A", "CSECT1B", "CSECT2A", "CSECT2B",
                        "CSECT3A", "CSECT3B", "CSECT4A", "CSECT4B", "CSECT5A",
                        "CSECT5B"], DGRATINGS[graname]["bslits"])
    update_head(RFLAT, ["NSLITS", "CSECT1A", "CSECT1B", "CSECT2A", "CSECT2B",
                        "CSECT3A", "CSECT3B", "CSECT4A", "CSECT4B", "CSECT5A",
                        "CSECT5B"], DGRATINGS[graname]["rslits"])

    # Create flat fields
    os.system("makeflat -ydist %s" % BFLAT)
    os.system("makeflat -ydist %s" % RFLAT)

    # Apply distortion, slits and flat-fields to arcs
    os.system("copyrect -ydist %s %s" % (BSHARP, barc))
    os.system("copyrect -ydist %s %s" % (RSHARP, rarc))
    os.system("copyslit %s %s -fft -ydist" % (BFLAT, barc))
    os.system("copyslit %s %s -fft -ydist" % (RFLAT, rarc))
    #os.system("flat1d %sslt.fits %s -fft -ydist" % (BFLAT.split(".")[0], barc))
    #os.system("flat1d %sslt.fits %s -fft -ydist" % (RFLAT.split(".")[0], rarc))
    os.system("flat2d %sflt.fits %s" % (BFLAT.split(".")[0], barc))
    os.system("flat2d %sflt.fits %s" % (RFLAT.split(".")[0], rarc))

    # Create new files for x-distortion
    os.system("flat2d %sblz.fits %sf.fits" % (BFLAT.split(".")[0], barc.split(".")[0]))
    os.system("flat2d %sblz.fits %sf.fits" % (RFLAT.split(".")[0], rarc.split(".")[0]))
    os.system("efits %sblz.fits %sff.fits 'i2*greater(i1,0.1)' xdist_B.fits" %
              (BFLAT.split(".")[0], barc.split(".")[0]))
    os.system("efits %sblz.fits %sff.fits 'i2*greater(i1,0.1)' xdist_R.fits" %
              (RFLAT.split(".")[0], rarc.split(".")[0]))

    # Calculate x distortion
    os.system("getrect -xdist -ydist xdist_B.fits -nx 16 -dy 3 -b 3 -x 2 -y 2 -dump -normwt 1e6")
    os.system("getrect -xdist -ydist xdist_R.fits -nx 16 -dy 3 -b 3 -x 2 -y 2 -dump -normwt 1e6")

    # Apply x distortion and calculate wavelength solution
    os.system("copyrect -xdist xdist_B.fits %sf.fits" % barc.split(".")[0])
    os.system("copyrect -xdist xdist_R.fits %sf.fits" % rarc.split(".")[0])
    [grating, cwlen] = get_head(arc, ["GRATENAM", "G3TLTWAV"], extn=0)
    d0 = DGRATINGS[grating]["dispersion"]
    blen = get_head(barc, "NAXIS1"); rlen = get_head(rarc, "NAXIS1")    
    os.system("waverect -xdist -ydist %sf.fits -o 4 -w1 %.2f -w2 %.2f -d0 %.2f -ll %s -arcamp -fwhm 8.5" % (barc.split(".")[0], cwlen - blen * d0, cwlen, d0, DARCLIST))
    os.system("waverect -xdist -ydist %sf.fits -o 4 -w1 %.2f -w2 %.2f -d0 %.2f -ll %s -arcamp -fwhm 8.5" % (rarc.split(".")[0], cwlen, cwlen + rlen * d0, d0, DARCLIST))

    return

########################

def deimos_pipe(science, flats, arc, docalib=yes, doscience=yes):

    '''Take a series of DEIMOS LVM spectra from a given night, process to
    wavelength-calibrated, sky-subtracted, cosmic-ray cleaned 2D images.'''

    if docalib:
        deimos_docalib(flats, arc)

    if doscience:
        deimos_doscience(science, arc)

    return

##########################

def deimos_doscience(science, arc):

    '''Using previously created calibration files, process DEIMOS LVM
    object (and standard) spectra into clean 2D images.'''
    
    barc = "rtd%s_B.fits" % arc[6:10]; rarc = "rtd%s_R.fits" % arc[6:10]
    
    # Pre-process all science images
    scidict = {}
    for image in science:
        [obj, exptime, airmass, skysub, crrej] = get_head(image, ["OBJECT", "EXPTIME",
                                                          "AIRMASS", "SKYSUB", "CRREJ"],
                                                          extn=0)
        if not (skysub == 1):
            skysub = 0
        if not (crrej == 1):
            crrej = 0
        obj = obj.replace(" ", "_")
        if not scidict.has_key(obj):
            scidict[obj] = [[image, exptime, airmass, skysub]]
        else:
            scidict[obj].append([image, exptime, airmass, skysub])
        deimos_preproc(image)

    # Now process all science frame
    for obj, dets in scidict.iteritems():

        for i in range(len(dets)):

            [image, exptime, airmass, skysub] = dets[i]
            
            bsci = "rtd%s_B.fits" % image[6:10]; rsci = "rtd%s_R.fits" % image[6:10]
        
            # Apply y distortion
            os.system("copyrect -ydist %s %s" % (BSHARP, bsci))
            os.system("copyrect -ydist %s %s" % (RSHARP, rsci))
        
            # Copy slits
            os.system("copyslit %s %s -fft -ydist" % (BFLAT, bsci))
            os.system("copyslit %s %s -fft -ydist" % (RFLAT, rsci))
        
            # Apply flats
            #os.system("flat1d %sslt.fits %s -fft -ydist" % (BFLAT.split(".")[0], bsci))
            #os.system("flat1d %sslt.fits %s -fft -ydist" % (RFLAT.split(".")[0], rsci))
            os.system("flat2d %sflt.fits %s" % (BFLAT.split(".")[0], bsci))
            os.system("flat2d %sflt.fits %s" % (RFLAT.split(".")[0], rsci))
            #os.system("flat2d %sblz.fits %sf.fits" % (BFLAT.split(".")[0],
                                                      #bsci.split(".")[0]))
            #os.system("flat2d %sblz.fits %sf.fits" % (RFLAT.split(".")[0],
                                                      #rsci.split(".")[0]))

            # Apply x distortion
            os.system("copyrect -xdist xdist_B.fits %sf.fits" % bsci.split(".")[0])
            os.system("copyrect -xdist xdist_R.fits %sf.fits" % rsci.split(".")[0])

            # Apply wavelength solution
            os.system("copyrect -wdist %sf.fits %sf.fits" % (barc.split(".")[0],
                                                              bsci.split(".")[0]))
            os.system("copyrect -wdist %sf.fits %sf.fits" % (rarc.split(".")[0],
                                                              rsci.split(".")[0]))

            # Tweak based on sky lines (if possible)
            os.system("waverect -wdist -ydist %sf.fits -ll %s -zero 1 -skyamp" % (bsci.split(".")[0], DSKYLIST))
            os.system("waverect -wdist -ydist %sf.fits -ll %s -zero 1 -skyamp" % (rsci.split(".")[0], DSKYLIST))

            # If desired, subtract sky
            os.system("skyrect -ydist %sf.fits -skyint 512 -kx 3 -ky 1 -b 5 -dkx 1.00 -skycut 5.0 -check 15 -pct 40 -g %.2f -rn %.2f" % (bsci.split(".")[0], CCDGAIN, CCDRNOISE))
            os.system("skyrect -ydist %sf.fits -skyint 512 -kx 3 -ky 1 -b 5 -dkx 1.00 -skycut 5.0 -check 15 -pct 40 -g %.2f -rn %.2f" % (rsci.split(".")[0], CCDGAIN, CCDRNOISE))
            if skysub:
                os.system("efits %sf.fits %sfm.fits 'i1-i2' %sfo.fits" %
                          (bsci.split(".")[0], bsci.split(".")[0], bsci.split(".")[0]))
                os.system("efits %sf.fits %sfm.fits 'i1-i2' %sfo.fits" %
                          (rsci.split(".")[0], rsci.split(".")[0], rsci.split(".")[0]))
            else:
                shutil.copy("%sf.fits" % bsci.split(".")[0],
                            "%sfo.fits" % bsci.split(".")[0])
                shutil.copy("%sf.fits" % rsci.split(".")[0],
                            "%sfo.fits" % rsci.split(".")[0])

            # Clean cosmic rays
            if crrej:
                os.system("pyrclean %sfo.fits -th 20.0" % bsci.split(".")[0])
                os.system("pyrclean %sfo.fits -th 25.0" % rsci.split(".")[0])
            else:
                shutil.copy("%sfo.fits" % bsci.split(".")[0],
                            "%sfoc.fits" % bsci.split(".")[0])
                shutil.copy("%sfo.fits" % rsci.split(".")[0],
                            "%sfoc.fits" % rsci.split(".")[0])

            # Unrectify
            os.system("unrect -wdist -xdist -ydist %sfoc.fits" % bsci.split(".")[0])
            os.system("unrect -wdist -xdist -ydist %sfoc.fits" % rsci.split(".")[0])

            # Rename and add / remove keywords
            shutil.copy("%sfocr.fits" % bsci.split(".")[0],
                        "%s_%02i_B.fits" % (obj, i+1))
            update_head("%s_%02i_B.fits" % (obj, i+1), ["EXPTIME", "AIRMASS",
                                                        "SKYSUB", "ONAME", "DISPAXIS"],
                                                       [exptime, airmass, skysub,
                                                        image, 1])
            delete_head("%s_%02i_B.fits" % (obj, i+1), ["LTM1_2", "LTM2_1",
                                                           "CD1_2", "CD2_1"])
            shutil.copy("%sfocr.fits" % rsci.split(".")[0],
                        "%s_%02i_R.fits" % (obj, i+1))
            update_head("%s_%02i_R.fits" % (obj, i+1), ["EXPTIME", "AIRMASS",
                                                        "SKYSUB", "ONAME", "DISPAXIS"],
                                                       [exptime, airmass, skysub,
                                                        image, 1])
            delete_head("%s_%02i_R.fits" % (obj, i+1), ["LTM1_2", "LTM2_1",
                                                           "CD1_2", "CD2_1"])
            
            # Delete raw image
            os.remove(image)


#################################################

def deimos_extract(science, standard, dostandard=yes):

    '''Extract and flux calibrate DEIMOS spectra (assuming they have been
    processed into 2D images using deimos_pipe above).'''

    if dostandard:
        deimos_standard(standard)

    for source in science:

        images = iraffiles("%s_??_B.fits" % source)
        joinstr = ""

        for image in images:

            bimage = image.split(".")[0]; rimage = bimage[:-1] + "R"
            
            # Extract 1D spectra
            if get_head("%s.fits" % bimage, "SKYSUB"):
                iraf.apall(bimage, output="", inter=yes, find=yes, recenter=yes,
                           resize=yes, edit=yes, trace=yes, fittrace=yes, extract=yes,
                           extras=yes, review=yes, background="none",
                           reference="%s_01_B" % standard)
                iraf.apall(rimage, output="", inter=yes, find=yes, recenter=yes,
                           resize=yes, edit=yes, trace=yes, fittrace=yes, extract=yes,
                           extras=yes, review=yes, background="none",
                           reference="%s_01_R" % standard)
            else:
                iraf.apall(bimage, output="", inter=yes, find=yes, recenter=yes,
                           resize=yes, edit=yes, trace=yes, fittrace=yes, extract=yes,
                           extras=yes, review=yes, background="fit",
                           reference="%s_01_B" % standard)
                iraf.apall(rimage, output="", inter=yes, find=yes, recenter=yes,
                           resize=yes, edit=yes, trace=yes, fittrace=yes, extract=yes,
                           extras=yes, review=yes, background="fit",
                           reference="%s_01_R" % standard)

            # Normalize by the standard continuua
            iraf.sarith("%s.ms" % bimage, "/", "s%s_01_B.ms.fits" % standard,
                        "s%s.ms" % bimage)
            iraf.sarith("%s.ms" % rimage, "/", "s%s_01_R.ms.fits" % standard,
                        "s%s.ms" % rimage)

            # Telluric correct
            iraf.telluric("s%s.ms" % bimage, "ts%s.ms" % bimage, "telluric.B.fits",
                          tweakrms=yes, interactive=yes, sample='6850:6950,7575:7700')
            iraf.telluric("s%s.ms" % rimage, "ts%s.ms" % rimage, "telluric.R.fits",
                          tweakrms=yes, interactive=yes, sample='6850:6950,7575:7700')
            
            # Flux calibration
            iraf.calibrate("ts%s.ms" % bimage, "fts%s.ms" % bimage, extinct=yes,
                           flux=yes, extinction="home$extinct/maunakeaextinct.dat",
                           observatory="Keck", sensitivity="%s.B.sens" % standard,
                           airmass='', exptime='')
            iraf.calibrate("ts%s.ms" % rimage, "fts%s.ms" % rimage, extinct=yes,
                           flux=yes, extinction="home$extinct/maunakeaextinct.dat",
                           observatory="Keck", sensitivity="%s.R.sens" % standard,
                           airmass='', exptime='')

            joinstr += "fts%s.ms,fts%s.ms," % (bimage, rimage)
           
        # Combine
        iraf.scombine(joinstr[:-1], "%s.ms.fits" % source, combine="average",
                      reject="avsigclip", w1=INDEF, w2=INDEF, dw=INDEF, nw=INDEF,
                      scale="none", zero="none", weight="none", lsigma=3.0,
                      hsigma=3.0, gain=CCDGAIN, rdnoise=CCDRNOISE)

        # Plot to enable final tweaks
        iraf.splot("%s.ms.fits" % source)

    return

#################################################

def deimos_standard(standard):

    '''Extract standard star and calculate sensitivit function.'''
    
    bstd = "%s_01_B" % standard; rstd = "%s_01_R" % standard

    # Extract standard
    if get_head("%s.fits" % bstd, "SKYSUB"):
        iraf.apall(bstd, output="", inter=yes, find=yes, recenter=yes, resize=yes,
                   edit=yes, trace=yes, fittrace=yes, extract=yes, extras=yes,
                   review=yes, background="none")
        iraf.apall(rstd, output="", inter=yes, find=yes, recenter=yes, resize=yes,
                   edit=yes, trace=yes, fittrace=yes, extract=yes, extras=yes,
                   review=yes, background="none")
    else:
        iraf.apall(bstd, output="", inter=yes, find=yes, recenter=yes, resize=yes,
                   edit=yes, trace=yes, fittrace=yes, extract=yes, extras=yes,
                   review=yes, background="fit")
        iraf.apall(rstd, output="", inter=yes, find=yes, recenter=yes, resize=yes,
                   edit=yes, trace=yes, fittrace=yes, extract=yes, extras=yes,
                   review=yes, background="fit")

    # Fit the continuum
    iraf.continuum("%s.ms" % bstd, output="s%s.ms" % bstd, lines=1, bands=1,
                   type="fit", sample="*", naverage=1, function="chebyshev",
                   order=15, low_reject=3.0, high_reject=5.0, niterate=10, grow=1.0,
                   interac=yes)
    iraf.continuum("%s.ms" % rstd, output="s%s.ms" % rstd, lines=1, bands=1,
                   type="fit", sample="*", naverage=1, function="chebyshev",
                   order=15, low_reject=3.0, high_reject=5.0, niterate=10, grow=1.0,
                   interac=yes)
    
    # Divide through by smoothed standard
    iraf.sarith("%s.ms" % bstd, "/", "s%s.ms" % bstd, "ds%s.ms" % bstd)
    iraf.sarith("%s.ms" % rstd, "/", "s%s.ms" % rstd, "ds%s.ms" % rstd)

    # Create and apply telluric correction
    iraf.splot("%s.ms" % bstd) # Remove telluric, save as a%s.ms
    iraf.splot("%s.ms" % rstd) # Remove telluric, save as a%s.ms
    iraf.sarith("%s.ms" % bstd, "/", "a%s.ms" % bstd, "telluric.B.fits")
    iraf.sarith("%s.ms" % rstd, "/", "a%s.ms" % rstd, "telluric.R.fits")
    iraf.imreplace("telluric.B.fits", 1.0, lower=0.0, upper=0.0)
    iraf.imreplace("telluric.R.fits", 1.0, lower=0.0, upper=0.0)
    iraf.telluric("ds%s.ms" % bstd, "tds%s.ms" % bstd, "telluric.B.fits",
                  xcorr=no, tweakrms=no, interactive=no, sample='6850:6950,7575:7700')
    iraf.telluric("ds%s.ms" % rstd, "tds%s.ms" % rstd, "telluric.R.fits",
                  xcorr=no, tweakrms=no, interactive=no, sample='6850:6950,7575:7700')

    # Define bandpasses for standard star calculation
    iraf.standard("tds%s.ms" % bstd, "%s.B.std" % standard,
                  extinction='home$extinct/maunakeaextinct.dat',
                  caldir='home$standards/', observatory='Keck', interac=yes,
                  star_name=standard, airmass='', exptime='')
    iraf.standard("tds%s.ms" % rstd, "%s.R.std" % standard,
                  extinction='home$extinct/maunakeaextinct.dat',
                  caldir='home$standards/', observatory='Keck', interac=yes,
                  star_name=standard, airmass='', exptime='')

    # Determine sensitivity function
    iraf.sensfunc('%s.B.std' % standard, '%s.B.sens' % standard,
                  extinction='home$extinct/maunakeaextinct.dat',
                  newextinction='extinct.dat', observatory='Keck',
                  function='legendre', order=4, interactive=yes)
    iraf.sensfunc('%s.R.std' % standard, '%s.R.sens' % standard,
                  extinction='home$extinct/maunakeaextinct.dat',
                  newextinction='extinct.dat', observatory='Keck',
                  function='legendre', order=4, interactive=yes)

    return

##########################











