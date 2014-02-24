#! /usr/bin/env python

from iqutils import *
import iqpkg
from mx.DateTime import *
import subprocess
import numpy as np

from pyraf import iraf
iraf.digiphot()
iraf.apphot()

######################

EXPTIMEKEY = "EXPTIME"
OBSDATEKEY = "UTC-OBS"
SIGMA = 2.0

######################

def coaddbyday(inlis):

    # Create a list of dates on which exposures were obtained
    tobs = []
    for image in inlis:
        tob = DateTimeFrom(get_head(image, OBSDATEKEY))
        if not tob.date in tobs:
            tobs.append(tob.date)

    # For each date on which we got exposures
    for odate in tobs:

        # Find all the images from that date
        ims = []; tims = []
        for image in inlis:
            tob = DateTimeFrom(get_head(image, OBSDATEKEY))
            if tob.date == odate:
                ims.append(image)
                tims.append(DateTimeFrom(get_head(image, OBSDATEKEY)).mjd)

        # Swarp cmd
        cmd = "swarp " + " ".join(ims)
        cmd += " -IMAGEOUT_NAME %s.fits " % odate.replace("-", "")
        cmd += "-WEIGHTOUT_NAME %s.weight.fits " % odate.replace("-", "")
        cmd += "-WEIGHT_TYPE MAP_WEIGHT -WEIGHT_SUFFIX .weight.fits"
        subprocess.call(cmd, shell=True)

        # Get midpoint of exposure
        ntims = np.array(tims)
        tmed = np.mean(ntims)
        mtmed = DateTimeFromMJD(tmed).strftime("%Y-%m-%dT%H:%M:%S")
        update_head("%s.fits" % odate.replace("-", ""), "UTC-OBS", mtmed) 
        
    return

###########################

def p48coaddphot(inlis, refstars, ot, filter, outfile):

    outf = open(outfile, "w")
    
    for image in inlis:

        root = image.split(".")[0]
        
        # Get the seeing in the image
        iqpkg.iqobjs(image, SIGMA, get_head(image, "SATURATE"), skyval="0.0")
        seepix = get_head(image, "SEEPIX")

        # Calculate zeropoint
        refstars.wcs2pix(image)
        refstars.set_mag(filter)
        xyfile = open("%s.xy" % root, "w")
        for star in refstars:
            xyfile.write("%10.3f%10.3f\n" % (star.xval, star.yval))
        xyfile.close()
        iraf.phot(image, coords="%s.xy" % root, output="%s.mag" % root,
                  aperture=1.2*float(seepix), interac=no)
        stars = Starlist("%s.mag" % root)
        zp, zpu = stars.zeropt(refstars,method="mean",rejout=0)

        # Measure source
        ot.wcs2pix(image)
        coofile = open("%s.coo" % root, "w")
        coofile.write("%10.3f%10.3f\n" % (ot[0].xval, ot[0].yval))
        coofile.close()
        iraf.phot(image, coords="%s.coo" % root, output="%s.ot" % root,
                  aperture=1.2*float(seepix), calgorithm="none", interac=no)
        stars = Starlist("%s.ot" % root)

        # Just return 99 for non-detection
        if len(stars)==1:
            smag = stars[0].mag
            smagu = stars[0].magu
            mag = smag + zp
        else:
            mag = 99.0
            smagu = 99.0
        mjdobs = DateTimeFrom(get_head(image, OBSDATEKEY)).mjd
        exptime = get_head(image, EXPTIMEKEY)
        outf.write("%15.3f%10.1f%10.3f%10.3f%10.3f\n" % (mjdobs, exptime, mag, smagu, zpu))

    outf.close()

    
