#! /usr/bin/env python

import psycopg2 		# Talk to the database
import pyfits           	# Get header information and trim images
import numpy as np		# Numerical python

import sys, os, shutil, datetime, time, math
import urllib, urllib2
from subprocess import PIPE, Popen

########## Global variables

[IPTFDB, IPTFUSER, IPTFHOST, IPTFPASS] = \
 ["iptf", "iptf_admin", "scidb2.nersc.gov", "ip33d$kyy"]

ROOTDIR = "/global/scratch/sd/cenko/ppqlc"

NPIX = 300; NSUBPIX = 50
BPMDIR = "/project/projectdirs/deepsky/iptf/calib/bias/super"

SEXDIR = "/project/projectdirs/deepsky/iptf/workroom/sbc"
SEXCMD = "/project/projectdirs/deepsky/iptf/scripts/carver/usr/local/bin/sex"
SCAMPDIR = "/project/projectdirs/deepsky/iptf/workroom/sbc"
SCAMPCMD = "/project/projectdirs/deepsky/iptf/scripts/carver/usr/local/bin/scamp"
SWARPDIR = "/project/projectdirs/deepsky/iptf/workroom/sbc"
SWARPCMD = "/project/projectdirs/deepsky/iptf/scripts/carver/usr/local/bin/swarp"
PSFEXDIR = "/project/projectdirs/deepsky/iptf/workroom/sbc"

RSCALE = 2.5; RSSSCALE = 6.0

P48SATURATE = {0: 60000.0, 1: 60000.0, 2: 60000.0, 3: 60000.0,
               4: 60000.0, 5: 60000.0, 6: 60000.0, 7: 60000.0,
               8: 60000.0, 9: 60000.0, 10: 60000.0, 11: 60000.0}

[SUBNSX, SUBNSY, SUBKO, SUBBGO, SUBNORM, SUBCONV, NOISEFILL] = \
 [3, 3, 0, 0, "i", "t", 70000.0]

STARMATCH = 1.5

PTFMARSHALPHOTURL = "http://ptf.caltech.edu/cgi-bin/ptf/transient/edit_phot.cgi"
PTFMARSHALPROGID = 47
PTFMARSHALINSTID = 5
PTFMARSHALUSER = "brad"
PTFMARSHALPASSWD = "iptf374"

########## Routines

def ppqlc(ptfname=None, ra=None, dec=None, dosub=True, filts=['g', 'R'],
          updateMarshal=True):

    if (ptfname==None) and (ra==None or dec==None):
        raise SystemExit("Please specify either a valid PTF name (i.e. 11kly) and/or valid coordinates (i.e. RA and Dec in degrees)")

    # If no PTF name, call the target "None" for now
    if (ptfname==None):
        ptfname="None"

    # Get RA and Dec if not entered by user
    if (ra==None) and (dec==None):
        [ra, dec] = get_coords(ptfname)

    # Create the directory where all this will be stored
    # Structure is: 
    # ROOTDIR/<filt>/refs - reference frames
    # ROOTDIR/<filt>/new - new images
    # ROOTDIR/<filt>/subs - Subtractions
    create_dirs(ptfname, filts)

    for filt in filts:

        # If doing subtractions, find the reference images and stars
        if dosub:
            refpairs = get_p48ref(ra, dec, filt, ptfname)            
        else:
            refpairs = None

        # Pull over all the new images
        newims = get_p48new(ra, dec, filt, ptfname, refpairs)

        # Run the subtractions
        if dosub:
            ppqlc_p48sub(newims, refpairs, ra, dec, ptfname, filt=filt)
        else:
            for i in range(len(newims)):
                image = newims[i]
                shutil.copy("%s/%s/%s/new/new.%s.%i.%i.%05i.fits" % \
                         (ROOTDIR, ptfname, filt, filt, image[2], image[3], i),
                            "%s/%s/%s/subs/new.%s.%i.%i.%05i.sub.fits " % \
                         (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))

        # PSF photometry
        ppqlc_p48phot("%s/%s/%s.%s.ppqlc.dat" % (ROOTDIR, ptfname, ptfname, 
                      filt), newims, ra, dec, ptfname, dosub=dosub, filt=filt)

        if updateMarshal:
            ppqlc_updateMarshal("%s/%s/%s.%s.ppqlc.dat" % (ROOTDIR, ptfname, 
                                ptfname, filt), ptfname, filt=filt)

    return

##########

def ppqlc_updateMarshal(dfile, ptfname, filt="R"):

    """Update the Marshal with Photometry Results"""

    # Open data file
    data = np.loadtxt(dfile)
    if len(data.shape)==1:
        data = [data]

    # Get source ID from database
    db = iptfdb()
    cur = db.cursor()

    cur.execute("SELECT marshal_source_id FROM transients WHERE " + \
                "ptfname = '%s'" % ptfname)
    reply = cur.fetchall()
    sourceid = reply[0][0]

    # Website info
    url = PTFMARSHALPHOTURL 

    # Loop over individual rows
    for row in data:

        # Quantities that don't change
        infodict = {'sourceid': sourceid,
                    'commit': 'yes',
                    'programid': PTFMARSHALPROGID, 
                    'instrumentid': PTFMARSHALINSTID,
                    'issub': 'yes',
                    'reducedby': 'SBC'}

        # Details of given observation
        infodict['jdobs'] = row[0] + 2400000.5
        infodict['filter'] = filt 
        if row[8] == 0:
            infodict['refsys'] = 'SDSS-DR9'
        else:
            infodict['refsys'] = 'IPAC-PTF'

        # If everything failed
        if row[2] == 99.0:
            continue

        # Upper limit
        if row[3] == 99.0:
            infodict['mag'] = 99.0
            infodict['emag'] = 99.0
            infodict['limmag'] = row[2]

        # Otherwise detection
        else:
            infodict['mag'] = row[2]
            infodict['emag'] = row[3]
            infodict['limmag'] = 99.0

        photdata = urllib.urlencode(infodict)
        configure_http()
        flob = urllib2.urlopen(url, photdata)
        s = flob.read()

    return

##########

def configure_http(url=PTFMARSHALPHOTURL,user=PTFMARSHALUSER,
                   passw=PTFMARSHALPASSWD):

     """Configure HTTP to use user/password combo for appropriate URLs"""

     x = urllib2.HTTPPasswordMgrWithDefaultRealm()
     x.add_password(None, url, user, passw)
     auth = urllib2.HTTPBasicAuthHandler(x)
     opener = urllib2.build_opener(auth)
     urllib2.install_opener(opener)

     return

##########

def ppqlc_p48phot(outname, images, ra, dec, ptfname="None", filt="R", 
                  dosub=True):

    """PSF Photometry.  For now will use psfex.  Hope ultimately cjpphot."""

    # Open outfile
    outf = open(outname, "w")

    for i in range(len(images)):

        image = images[i]

        # If no subtracted image, bail
        if not os.path.exists("%s/%s/%s/subs/new.%s.%i.%i.%05i.sub.fits" % \
                       (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)):
            continue

        # Need to update GAIN keyword to float value
        x = pyfits.open("%s/%s/%s/subs/new.%s.%i.%i.%05i.sub.fits" % \
                        (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
        x[0].header['GAIN'] = float(x[0].header['GAIN'])
        x.writeto("%s/%s/%s/subs/new.%s.%i.%i.%05i.sub.fits" % \
                   (ROOTDIR, ptfname, filt, filt, image[2], image[3], i),
                   clobber=True)
        x.close()
        y = pyfits.open("%s/%s/%s/new/new.%s.%i.%i.%05i.fits" % \
                        (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
        y[0].header['GAIN'] = float(y[0].header['GAIN'])
        y.writeto("%s/%s/%s/new/new.%s.%i.%i.%05i.fits" % \
                   (ROOTDIR, ptfname, filt, filt, image[2], image[3], i),
                   clobber=True)
        y.close()

        # Run SExtrator with aperture photometry and vignettes 
        scmd = "%s -c %s/prepsfex.sex " % (SEXCMD, PSFEXDIR) + \
               "-CATALOG_NAME %s/%s/%s/subs/new.%s.%i.%i.%05i.pcat " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
               "-ASSOC_NAME %s/%s/%s/subs/new.%s.%i.%i.%05i.stamps " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
               "%s/%s/%s/new/new.%s.%i.%i.%05i.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) 
        junk = Popen(scmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        # Create PSF and model and extract photometry.  Try first with
        pcmd = "psfex -c %s/psfex.conf -PSF_DIR %s/%s/%s/subs " % \
                (PSFEXDIR, ROOTDIR, ptfname, filt) + \
                "%s/%s/%s/subs/new.%s.%i.%i.%05i.pcat " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)  
        junk = Popen(pcmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        scmd = "%s -c %s/psfex.sex -PSF_NMAX 1 " % (SEXCMD, PSFEXDIR) + \
               "-PSF_NAME %s/%s/%s/subs/new.%s.%i.%i.%05i.psf " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
               "-CATALOG_NAME %s/%s/%s/subs/new.%s.%i.%i.%05i.cat " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
               "%s/%s/%s/new/new.%s.%i.%i.%05i.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) 
   #    "-WEIGHT_IMAGE %s/%s/%s/subs/new.%s.%i.%i.%05i.sub.noise.fits " % \
   #                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
   #           "-WEIGHT_TYPE MAP_RMS -WEIGHT_THRESH %f " % (NOISEFILL - 1) 
        junk = Popen(scmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        scmd = "%s -c %s/psfex.sex -PSF_NMAX 1 " % (SEXCMD, PSFEXDIR) + \
               "-PSF_NAME %s/%s/%s/subs/new.%s.%i.%i.%05i.psf " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
               "-CATALOG_NAME %s/%s/%s/subs/new.%s.%i.%i.%05i.grb " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
               "%s/%s/%s/subs/new.%s.%i.%i.%05i.sub.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
   #    "-WEIGHT_IMAGE %s/%s/%s/subs/new.%s.%i.%i.%05i.sub.noise.fits " % \
   #                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
   #           "-WEIGHT_TYPE MAP_RMS -WEIGHT_THRESH %f " % (NOISEFILL - 1)
        junk = Popen(scmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        # If that failed, increase threshold for sextractor 
        if (not 
         os.path.exists("%s/%s/%s/subs/new.%s.%i.%i.%05i.grb" % \
                       (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))) \
         or os.stat("%s/%s/%s/subs/new.%s.%i.%i.%05i.grb" % \
         (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)).st_size == 0:
            
            scmd = "%s -c %s/psfex.sex " % (SEXCMD, PSFEXDIR) + \
                   "-PSF_NAME %s/%s/%s/subs/new.%s.%i.%i.%05i.psf " % \
                    (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
                "-CATALOG_NAME %s/%s/%s/subs/new.%s.%i.%i.%05i.grb " % \
                    (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
                   "-DETECT_THRESH 3.0 -ANALYSIS_THRESH 3.0 -PSF_NMAX 1 " + \
                   "%s/%s/%s/subs/new.%s.%i.%i.%05i.sub.fits " % \
                    (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) 
    #        "-WEIGHT_IMAGE %s/%s/%s/subs/new.%s.%i.%i.%05i.sub.noise.fits " % \
    #               (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
    #          "-WEIGHT_TYPE MAP_RMS -WEIGHT_THRESH %f " % (NOISEFILL - 1) 
            junk = Popen(scmd, shell=True, stdout=PIPE).stdout
            junk.readlines()

        # Photometric Calibration
        if not os.path.exists("%s/%s/%s/subs/new.%s.%i.%i.%05i.cat" %
                    (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)):
            stars = np.array([])
        else:
            stars = np.loadtxt("%s/%s/%s/subs/new.%s.%i.%i.%05i.cat" % 
                         (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
        refstars = np.loadtxt("%s/%s/%s/subs/new.%s.%i.%i.%05i.stamps" % 
                         (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
        sdss = refstars[0][3]
        
        # For each star, look for nearest match (in pixel space) closer than
        # STARMATCH pixels.  Make sure star is not flagged by SExtractor,
        # and if all looks OK, add the difference to the zpt list
        zpt = []
        for star in stars:
            min = [100.0, 99.0]
            for refstar in refstars:
                dist = np.sqrt(np.power(star[1]-refstar[0],2) + \
                               np.power(star[2]-refstar[1],2))
                if dist < STARMATCH and dist < min[0]:
                    min = [dist, refstar[2]]
            if min[0] != 100.0 and star[5] != 99.0 and star[7] == 0:
                zpt.append(star[5] - min[1])

        # Turn zpt into an array and do some sigma clipping
        zpta = np.array(zpt)
        if len(zpta)>2:
            szpta = sigma_clip(zpta, 2.0, 3)
        else:
            szpta = zpta
        
        # Find the coordinates of the transient/variable source
        wcmd = "sky2xy %s/%s/%s/subs/new.%s.%i.%i.%05i.sub.fits %f %f" % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i, ra, dec)
        coords = Popen(wcmd, shell=True, stdout=PIPE).stdout.readlines()
        [xc, yc] = [float(coords[0].split()[4]), float(coords[0].split()[5])]

        # Grab catalog from subtracted image
        if not os.path.exists("%s/%s/%s/subs/new.%s.%i.%i.%05i.grb" %
                    (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)):
            gstars = np.array([])
        else:
            statinfo = os.stat("%s/%s/%s/subs/new.%s.%i.%i.%05i.grb" %
                         (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
            if statinfo.st_size == 0:
                gstars = np.array([])
            else:
                gstars = np.loadtxt("%s/%s/%s/subs/new.%s.%i.%i.%05i.grb" %
                         (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
                if (len(gstars.shape) == 1) and (len(gstars) != 0):
                    gstars = np.array([gstars])

        # Look for a nearest match in the SExtractor catalog
        ot = [100.0, []]
        for star in gstars:
            dist = np.sqrt(np.power(star[1]-xc,2)+np.power(star[2]-yc,2))
            if dist < STARMATCH and dist < ot[0]:
                ot = [dist, star]

        # Determine magnitude (or limit)
        # If no comparison stars, just return 99
        if len(szpta)==0:
            mag = 99.0
            dmag = 99.0
            ispsf = 0
        # If match and PSF photometry, use this
        elif (ot[0] != 100.0) and (ot[1][3] != 0) and (ot[1][6] < 0.3):
            mag = ot[1][5] - szpta.mean()
            dmag = ot[1][6] # * 1.5
            ispsf = 1
        # Otherwise try aperture photometry
        elif (ot[0] != 100.0) and (ot[1][8] != 0) and (ot[1][11] < 0.3):
            mag = ot[1][10] - szpta.mean()
            dmag = ot[1][11] # * 1.5
            ispsf = 0
        # If no match in star list or no good detection, poor man's upper limit
        else:
            x = pyfits.open("%s/%s/%s/subs/new.%s.%i.%i.%05i.sub.fits" % \
             (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
            cts = np.sum(x[0].data[int(yc)-5:int(yc)+5,int(xc)-5:int(xc)+5])
            sky = np.median(
             x[0].data[int(yc)+15:int(yc)+25,int(xc)+15:int(xc)+25])
            skyerr = np.std(
             x[0].data[int(yc)+15:int(yc)+25,int(xc)+15:int(xc)+25])
            flux = np.maximum(cts - 100.0 * sky, 0.0) 
            df = np.sqrt(100.0 * np.power(skyerr, 2) + flux) 
            if df==0:
                mag = 99.0
                dmag = 99.0
                ispsf = 0
            else:
                mag = -2.5 * np.log10(3 * df + flux) + 27.5 - szpta.mean()
                dmag = 99.0
                ispsf = 0

        # Write results
        if math.isnan(szpta.std()):
            zp = 99.0
        else:
            zp = szpta.std()
 
        outf.write("%15.5f%8.2f%10.3f%10.3f%10.3f%5i%5i%5i%5i\n" % \
                   (x[0].header["OBSMJD"], x[0].header["EXPTIME"], mag,
                    dmag, zp, len(szpta), len(zpta), ispsf, sdss))

    outf.close()
        
    return

##########

def ppqlc_p48sub(images, refpairs, ra, dec, ptfname="None", filt="R"):

    """Run the subtractions"""
 
    # Generate reference catalogs
    for refpair in refpairs:

        scmd = "%s -c %s/scamp.sex " % (SEXCMD, SEXDIR)
        scmd += "-CATALOG_NAME %s/%s/%s/refs/ref.%s.%i.%i.cat " % \
                 (ROOTDIR, ptfname, filt, filt, refpair[0], refpair[1]) 
        scmd += "%s/%s/%s/refs/ref.%s.%i.%i.fits" % \
                 (ROOTDIR, ptfname, filt, filt, refpair[0], refpair[1])
        junk = Popen(scmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

    # Loop over images
    for i in range(len(images)):

        image = images[i]

        # Run sextractor on new images
        scmd = "%s -c %s/scamp.sex " % (SEXCMD, SEXDIR)
        scmd += "-CATALOG_NAME %s/%s/%s/new/new.%s.%i.%i.%05i.cat " % \
                 (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        scmd += "%s/%s/%s/new/new.%s.%i.%i.%05i.fits" % \
                 (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        junk = Popen(scmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        # Align with scamp
        scmd = "scamp -c %s/scamp.conf.cat -VERBOSE_TYPE QUIET " % SCAMPDIR
        scmd += "-ASTREFCAT_NAME %s/%s/%s/refs/ref.%s.%i.%i.cat " % \
                 (ROOTDIR, ptfname, filt, filt, image[2], image[3])
        scmd += "%s/%s/%s/new/new.%s.%i.%i.%05i.cat" % \
                 (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        junk = Popen(scmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        # Configure for swarp
        icmd = "imhead %s/%s/%s/new/new.%s.%i.%i.%05i.fits > " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        icmd += "%s/%s/%s/new/temp.head" % (ROOTDIR, ptfname, filt)
        junk = Popen(icmd, shell=True, stdout=PIPE).stdout
        junk.readlines()
        gcmd = "grep NAXIS %s/%s/%s/new/temp.head > " % \
                (ROOTDIR, ptfname, filt)
        gcmd += "%s/%s/%s/new/ref.%s.%i.%i.%05i.head" % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        junk = Popen(gcmd, shell=True, stdout=PIPE).stdout
        junk.readlines()
        lcmd = "less %s/%s/%s/new/new.%s.%i.%i.%05i.head >> " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        lcmd += "%s/%s/%s/new/ref.%s.%i.%i.%05i.head" % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        junk = Popen(lcmd, shell=True, stdout=PIPE).stdout
        junk.readlines()
        shutil.copy("%s/%s/%s/new/ref.%s.%i.%i.%05i.head" % \
                    (ROOTDIR, ptfname, filt, filt, image[2], image[3], i),
                    "%s/%s/%s/new/ref.%s.%i.%i.%05i.noise.head" % \
                    (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
        shutil.copy("%s/%s/%s/new/ref.%s.%i.%i.%05i.head" % \
                    (ROOTDIR, ptfname, filt, filt, image[2], image[3], i),
                    "%s/%s/%s/new/ref.%s.%i.%i.%05i.uncert.head" % \
                    (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
        os.remove("%s/%s/%s/new/temp.head" % (ROOTDIR, ptfname, filt))

        # Run swarp (including noise and uncertainty images)
        scmd = "swarp -c %s/default.swarp -VERBOSE_TYPE QUIET " % SWARPDIR 
        scmd += "%s/%s/%s/refs/ref.%s.%i.%i.fits -SUBTRACT_BACK N " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3])
        scmd += "-IMAGEOUT_NAME %s/%s/%s/new/ref.%s.%i.%i.%05i.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        scmd += "-WEIGHT_TYPE NONE"
        junk = Popen(scmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        scmd = "swarp -c %s/default.swarp -VERBOSE_TYPE QUIET " % SWARPDIR 
        scmd += "%s/%s/%s/refs/ref.%s.%i.%i.noise.fits -SUBTRACT_BACK N " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3])
        scmd += "-IMAGEOUT_NAME %s/%s/%s/new/ref.%s.%i.%i.%05i.noise.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        scmd += "-WEIGHT_TYPE NONE"
        junk = Popen(scmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        scmd = "swarp -c %s/default.swarp -VERBOSE_TYPE QUIET " % SWARPDIR 
        scmd += "%s/%s/%s/refs/ref.%s.%i.%i.uncert.fits -SUBTRACT_BACK N " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3])
        scmd += "-IMAGEOUT_NAME %s/%s/%s/new/ref.%s.%i.%i.%05i.uncert.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        scmd += "-WEIGHT_TYPE NONE"
        junk = Popen(scmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        # Get info for subtraction
        rhdulist = pyfits.open("%s/%s/%s/refs/ref.%s.%i.%i.fits" % \
                    (ROOTDIR, ptfname, filt, filt, image[2], image[3]))
        nhdulist = pyfits.open("%s/%s" % (image[0], image[1]))

        seepix = max(rhdulist[0].header["SEEING"],
                     nhdulist[0].header["SEEING"])
        r = RSCALE * seepix; rss = RSSSCALE * seepix
    
        refskybkg = rhdulist[0].header["MEDSKY"]
        refskysig = rhdulist[0].header["SKYSIG"]
        newskybkg = nhdulist[0].header["MEDSKY"]
        newskysig = nhdulist[0].header["SKYSIG"]
        il = newskybkg - 10.0 * newskysig; tl = refskybkg - 10.0 * refskysig

        r2hdulist = pyfits.open("%s/%s/%s/new/ref.%s.%i.%i.%05i.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
        n2hdulist = pyfits.open("%s/%s/%s/new/new.%s.%i.%i.%05i.fits" % \
                 (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))

        tu = 0.95 * r2hdulist[0].header["SATURATE"]
        iu = 0.95 * P48SATURATE[image[3]]

        ig = n2hdulist[0].header["GAIN"]
        ir = n2hdulist[0].header["READNOI"]
        tg = ig; tr = ir ### Need to fix this! 

        lmtmag = n2hdulist[0].header["LMT_MG"]

        # Get reference stars for PSF - first open a connection to db
        db = iptfdb()
        cur = db.cursor()
        sdss = 1

        # Try SDSS first
        f = filt.lower()
        cur.execute("SELECT ra, dec, mag FROM sdss_%s WHERE " % f + \
                    "q3c_radial_query(ra, dec, %f, %f, %f) " % \
                     (ra, dec, 2 * NPIX / 3600.0) +
                    "AND NOT q3c_radial_query(ra, dec, %f, %f, %f) " % \
                     (ra, dec, 1.5 * NSUBPIX / 3600.0) +
                    "AND mag < %f " % (lmtmag - 2.0) +
                    "ORDER BY mag ASC ")
        reply = cur.fetchall()
        
        # If no matches, use IPAC database instead
        if len(reply)==0:
            sdss = 0 
            cur.execute("SELECT ra, dec, mag FROM ipac_%s WHERE " % f + \
                        "q3c_radial_query(ra, dec, %f, %f, %f) " % \
                         (ra, dec, 2 * NPIX / 3600.0) + \
                        "AND NOT q3c_radial_query(ra, dec, %f, %f, %f) " % \
                         (ra, dec, 1.5 * NSUBPIX / 3600.0) +
                        "AND ptffield = %i AND ccdid = %i " % \
                         (image[2], image[3]) +
                        " AND mag < %f " % (lmtmag - 2.0) +
                        "ORDER BY mag ASC ")
            reply = cur.fetchall()

        # Convert RA/Dec to pixel coordinates, and write star list
        outf = open("%s/%s/%s/subs/new.%s.%i.%i.%05i.stamps" % \
                 (ROOTDIR, ptfname, filt, filt, image[2], image[3], i), "w")
        for rep in reply:
            wcmd = "sky2xy %s/%s/%s/new/new.%s.%i.%i.%05i.fits %f %f" % \
             (ROOTDIR, ptfname, filt, filt, image[2], image[3], i, rep[0], 
              rep[1])
            coords = Popen(wcmd, shell=True, stdout=PIPE).stdout.readlines()
            [xpix, ypix] = [float(coords[0].split()[4]), 
                            float(coords[0].split()[5])]
            if (np.abs(xpix) < (4 * NPIX)) and (np.abs(ypix) < (4 * NPIX)):
                outf.write("%10.3f%10.3f%10.3f%5i\n" % \
                            (xpix, ypix, rep[2], sdss))
        outf.close()

        # Make sure there are some reference stars to use
        x = np.loadtxt("%s/%s/%s/subs/new.%s.%i.%i.%05i.stamps" % \
                 (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
        nstars = 0
        for s in x:
            if s[0] > 0 and s[0] < 2 * NPIX and s[1] > 0 and s[1] < 2 * NPIX:
                nstars += 1
        if nstars > 3:
            usestars = 1
        else:
            usestars = 0

        # Subtract
        hcmd = "hotpants -inim %s/%s/%s/new/new.%s.%i.%i.%05i.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
                "-tmplim %s/%s/%s/new/ref.%s.%i.%i.%05i.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
                "-outim %s/%s/%s/subs/new.%s.%i.%i.%05i.sub.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
                " -ko %i -bgo %i -tu %f -iu %f -tl %f -il %f -r %f -rss %f " % \
                (SUBKO, SUBBGO, tu, iu, tl, il, r, rss) + \
                " -tni %s/%s/%s/new/ref.%s.%i.%i.%05i.noise.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
                " -imi %s/%s/%s/new/new.%s.%i.%i.%05i.mask.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
                " -nsx %i -nsy %i -tg %f -tr %f -ig %f -ir %f " % \
                (SUBNSX, SUBNSY, float(tg), float(tr), float(ig), float(ir))+ \
                " -savexy %s/%s/%s/subs/new.%s.%i.%i.%05i.xy " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
                " -omi %s/%s/%s/subs/new.%s.%i.%i.%05i.sub.mask.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
                " -oni %s/%s/%s/subs/new.%s.%i.%i.%05i.sub.noise.fits " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i) + \
                " -n %s -c %s -fin %f " % (SUBNORM, SUBCONV, NOISEFILL)
        if usestars:
            hcmd += "-ssf %s/%s/%s/subs/new.%s.%i.%i.%05i.stamps -afssc 0 " % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        junk = Popen(hcmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        # Copy cutout to new image
#       wcmd = "sky2xy %s/%s/%s/new/new.%s.%i.%i.%05i.fits %f %f" % \
#        (ROOTDIR, ptfname, filt, filt, image[2], image[3], i, ra, dec) 
#       coords = Popen(wcmd, shell=True, stdout=PIPE).stdout.readlines()
#       [xc, yc] = [float(coords[0].split()[4]), float(coords[0].split()[5])]
#       nimg = pyfits.open("%s/%s/%s/new/new.%s.%i.%i.%05i.fits" % \
#               (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
#       simg = pyfits.open("%s/%s/%s/subs/new.%s.%i.%i.%05i.sub.fits" % \
#               (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
#       xl = np.max([1, xc - NSUBPIX]); xu = np.min([2 * NPIX, xc + NSUBPIX])
#       yl = np.max([1, yc - NSUBPIX]); yu = np.min([2 * NPIX, yc + NSUBPIX])
#       nimg[0].data[yl:yu,xl:xu] = simg[0].data[yl:yu,xl:xu] + \
#                                   np.median(nimg[0].data[yl:yu,xl:xu]) 
#       nimg.writeto("%s/%s/%s/subs/new.%s.%i.%i.%05i.sphot.fits" % \
#                     (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))

    return
        
##########

def get_p48new(ra, dec, filt="R", ptfname="None", refpairs=None):

    """Retrieve all the new images for a given RA/Dec/filter."""

    # Need to talkto the database for this one
    db = iptfdb()
    cur = db.cursor()
 
    # Grab all the new images available
    cur.execute("SELECT base_dir, filename, ptffield, ccdid from IMAGE " + \
                "WHERE q3c_poly_query(%f, %f, " % (ra, dec) + \
                "ARRAY[ra_ll, dec_ll, ra_lr, dec_lr, ra_ur, dec_ur, " + \
                "ra_ul, dec_ul]) AND filter LIKE '%s' ORDER BY jd" % filt)
    reply = cur.fetchall()
    if len(reply)==0:
        print "No new images at this location: %f, %f" % (ra, dec)
        return []

    # Remove ones where not on FIELD/CCDID of interest
    if refpairs != None:
        nreply = []
        for rep in reply:
            if [rep[2], rep[3]] in refpairs:
                nreply.append(rep)
    else:
        nreply = reply

    # Loop over images
    for i in range(len(nreply)):

        image = nreply[i] 

        # Find the location of transient (in pixels)
        wcmd = "sky2xy %s/%s %f %f" % (image[0], image[1], ra, dec)
        coords = Popen(wcmd, shell=True, stdout=PIPE).stdout.readlines()
        [xc, yc] = [float(coords[0].split()[4]), float(coords[0].split()[5])]
        
        # Find the total number of pixels in each axis
        hdulist = pyfits.open("%s/%s" % (image[0], image[1]))
        [nax1, nax2] = [hdulist[0].header["NAXIS1"],
                        hdulist[0].header["NAXIS2"]]

        # Calculate edges for trimming
        if (xc < 2 * NPIX):
            xl = 1; xu = 2 * NPIX
        elif (xc > (nax1 - 2 * NPIX)):
            xl = nax1 - 2 * NPIX + 1; xu = nax1
        else:
            xl = xc - NPIX + 1; xu = xc + NPIX
        if (yc < 2 * NPIX):
            yl = 1; yu = 2 * NPIX
        elif (yc > (nax2 - 2 * NPIX)):
            yl = nax2 - 2 * NPIX + 1; yu = nax2
        else:
            yl = yc - NPIX + 1; yu = yc + NPIX

        # Trim the image  
        tcmd = "getfits %s/%s %.2f-%.2f %.2f-%.2f -o %s/%s/%s/new/new.%s.%i.%i.%05i.fits" % (image[0], image[1], xl, xu, yl, yu, ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        junk = Popen(tcmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        # Get appropriate bad pixel mask
        t = datetime.datetime.strptime(hdulist[0].header["UTC-OBS"],
                                       "%Y-%m-%dT%H:%M:%S.%f")
        tcmd = "getfits %s/%i%02i%ix/mask_C%02i.fits %.2f-%.2f %.2f-%.2f " % \
            (BPMDIR, t.year, t.month, t.day / 10, image[3], xl, xu, yl, yu) + \
               "-o %s/%s/%s/new/new.%s.%i.%i.%05i.bpm.fits" % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i)
        junk = Popen(tcmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

        # Convert this to a mask for hotpants (i.e., swap 1s and 0s)
        x = pyfits.open("%s/%s/%s/new/new.%s.%i.%i.%05i.bpm.fits" % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))
        x[0].data = np.where(x[0].data == 1, 0.0, 1.0)
        x.writeto("%s/%s/%s/new/new.%s.%i.%i.%05i.mask.fits" % \
                (ROOTDIR, ptfname, filt, filt, image[2], image[3], i))

    return nreply

##########

def get_p48ref(ra, dec, filt="R", ptfname="None"):

    """Retrieve all the reference frames for a given RA/Dec/filter."""

    # Need to talk to the database for this one
    db = iptfdb()
    cur = db.cursor() 
    
    # Find the ccdid and ptffield of all the  processed images where 
    # the ra,dec fall on the image
    cur.execute("SELECT DISTINCT ptffield, ccdid FROM image where " + \
                "q3c_poly_query(%f, %f, " % (ra, dec) + \
                " ARRAY[ra_ll, dec_ll, ra_lr, dec_lr, ra_ur, dec_ur, " + \
                "ra_ul, dec_ul])")
    reply = cur.fetchall()

    # Figure out with ptfield/ccdid pairs we need
    if (len(reply)==0):
        print "No reference images in filter %s" % filt
        cur.close()
        db.close()
        return []
    refpairs = []
    for rep in reply:
        refpairs.append([int(rep[0]), int(rep[1])])
    
    # Grab the appropriate images and trim them
    for refpair in refpairs:

        # Get the info from the database
        cur.execute("SELECT base_dir, filename FROM ref " + \
                    "WHERE ccdid = %i AND ptffield = %i AND filter = '%s'" %
                    (refpair[1], refpair[0], filt))
        reply = cur.fetchall()
        if len(reply)==0:
            print "No reference images for field %i, ccdid %i" % (refpair[0], refpair[1])
            continue

        # Find the location of transient (in pixels)
        wcmd = "sky2xy %s/%s %f %f" % (reply[0][0], reply[0][1], ra, dec)
        coords = Popen(wcmd, shell=True, stdout=PIPE).stdout.readlines()
        [xc, yc] = [float(coords[0].split()[4]), float(coords[0].split()[5])]

        # Find the total number of pixels in each axis
        hdulist = pyfits.open("%s/%s" % (reply[0][0], reply[0][1]))
        [nax1, nax2] = [hdulist[0].header["NAXIS1"],
                        hdulist[0].header["NAXIS2"]]
        
        # Calculate edges for trimming
        if (xc < 4 * NPIX):
            xl = 1; xu = 4 * NPIX
        elif (xc > (nax1 - 4 * NPIX)):
            xl = nax1 - 4 * NPIX + 1; xu = nax1
        else:
            xl = xc - 2 * NPIX + 1; xu = xc + 2 * NPIX
        if (yc < 4 * NPIX):
            yl = 1; yu = 4 * NPIX
        elif (yc > (nax2 - 4 * NPIX)):
            yl = nax2 - 4 * NPIX + 1; yu = nax2
        else:
            yl = yc - 2 * NPIX + 1; yu = yc + 2 * NPIX

        # Trim the images 
        tcmd = "getfits %s/%s %.2f-%.2f %.2f-%.2f -o %s/%s/%s/refs/ref.%s.%i.%i.fits" % (reply[0][0], reply[0][1], xl, xu, yl, yu, ROOTDIR, ptfname, filt, filt, refpair[0], refpair[1])
        junk = Popen(tcmd, shell=True, stdout=PIPE).stdout
        junk.readlines()
  
        # Trim noise and uncertainty, too
        tcmd = "getfits %s/%s.noise.fits %.2f-%.2f %.2f-%.2f -o %s/%s/%s/refs/ref.%s.%i.%i.noise.fits" % (reply[0][0], reply[0][1].split(".")[0], xl, xu, yl, yu, ROOTDIR, ptfname, filt, filt, refpair[0], refpair[1])       
        junk = Popen(tcmd, shell=True, stdout=PIPE).stdout
        junk.readlines()
        tcmd = "getfits %s/%s %.2f-%.2f %.2f-%.2f -o %s/%s/%s/refs/ref.%s.%i.%i.uncert.fits" % (reply[0][0], reply[0][1].replace("refimg", "uncert") , xl, xu, yl, yu, ROOTDIR, ptfname, filt, filt, refpair[0], refpair[1])
        junk = Popen(tcmd, shell=True, stdout=PIPE).stdout
        junk.readlines()

    cur.close()
    db.close()
    return refpairs

##########

def create_dirs(ptfname, filts):

    "Create the directory structure where all the data will be stored."""

    # Make sure rootdir is there
    if not os.path.exists(ROOTDIR):
        os.mkdir(ROOTDIR)
    
    # First move the old directory (if there is one)
    if os.path.exists("%s/%s" % (ROOTDIR, ptfname)):
       #ndirs = subprocess.check_output("ls %s | grep %s | wc" % 
       #                                 (ROOTDIR, ptfname), shell=True)
       #shutil.move("%s/%s" % (ROOTDIR, ptfname), "%s/%s.%i" %
       #            (ROOTDIR, ptfname, int(ndirs.split()[0])))
       shutil.rmtree("%s/%s" % (ROOTDIR, ptfname))
    os.mkdir("%s/%s" % (ROOTDIR, ptfname))

    # Then create all the new ones
    for filt in filts:
        os.mkdir("%s/%s/%s" % (ROOTDIR, ptfname, filt))
        os.mkdir("%s/%s/%s/refs" % (ROOTDIR, ptfname, filt))
        os.mkdir("%s/%s/%s/subs" % (ROOTDIR, ptfname, filt))
        os.mkdir("%s/%s/%s/new" % (ROOTDIR, ptfname, filt))

    return

##########

def get_coords(ptfname):

    """For a given PTF name, return the RA and Dec in decimal degrees."""

    # Need a cursor here
    db = iptfdb()
    cur = db.cursor()

    # Query here to get ra and dec based on name
    cur.execute("select ra, dec from transients where ptfname = '%s'" % ptfname)
    reply = cur.fetchall()
    if len(reply)==0 or len(reply[0])!=2:
        raise SystemExit("Please enter a valid PTF Name!: %s" % ptfname)
    radeg = float(reply[0][0])
    dcdeg = float(reply[0][1])
    
    # Close and sign off
    cur.close()
    db.close()

    return [radeg, dcdeg]
    
##########

def iptfdb(dbname=IPTFDB, user=IPTFUSER, host=IPTFHOST, passw=IPTFPASS):

    """Connect to the iPTF database."""

    db = psycopg2.connect("dbname= %s user=%s host=%s password=%s" % \
          (dbname, user, host, passw))

    return db

##########

def sigma_clip(data, clipsig=2.0, maxiter=3):

    """Iterative sigma clipping."""

    nloop = 0
    while nloop < maxiter:

        good = np.where(np.abs(data - np.mean(data)) < clipsig * np.std(data))
        data = data[good]
        nloop += 1

    return data

##########

if __name__ == "__main__":

    ppqlc(ptfname=sys.argv[1])

