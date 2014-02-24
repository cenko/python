# Global packages

# Local packages
from iqutils import *
import iqpkg

def id_ot(inlis, refstars, tmass, ra, dec, radius, output="Cands.reg", 
          sigma=3.0, satval=60000.0, wtroot="/Users/scenko/CALIB/p48/mask_C", 
          pixtol=5.0):

     fchips = {}     # Stores field / chip dictionary 
     center = astrocoords(ra, dec)     # Center of error region
     newcands = []     # New candidate OTs
     
     for image in inlis:
         
         # Identify all field and chip combinations
         [field, chip] = get_head(image, ["PTFFIELD", "CCDID"])
         if not fchips.has_key(field):
             fchips[field] = [chip]
         else:
             try:
                 ind = fchips[field].index(chip)
             except ValueError:
                 fchips[field].append(chip)

         # Also grab objects
         iqpkg.iqobjs(image, sigma, satval, skyval="0.0", wtimage="%s%02i.fits" %
                (wtroot, int(chip)), wtcut=0.1, fwhm=1.5, pix=1.0)

     for (field, chips) in fchips.iteritems():

         for chip in chips:

             nonmatches = []
             for image in inlis:

                 # Only look at images for the given field/chip pair
                 [nfield, nchip] = get_head(image, ["PTFFIELD", "CCDID"])
                 if (nfield == field) and (nchip == chip):

                     # Grab all the objects from this image that don't appear in
                     # reference star list
                     stars = Starlist(get_head(image, "STARFILE"))
                     refstars.wcs2pix(image)
                     new = stars.nomatch(refstars, tol=pixtol)

                     # Make sure the new objects fall inside error circle
                     new.pix2wcs(image)
                     new_incircle = []
                     for star in new:
                         temp = astrocoords(star.radeg(), star.dcdeg())
                         dist = temp.diff(center, degree=1)
                         if sqrt(pow(dist[0],2) + pow(dist[1],2)) < radius:
                             new_incircle.append(star)

                     # If first image from field/chip pair, add everything to
                     # nonmatches.  Otherwise, only add objects if they have shown
                     # up already
                     if len(nonmatches) == 0:
                         nonmatches = new_incircle
                     else:
                         s1 = Starlist(stars=new_incircle)
                         sr = Starlist(stars=nonmatches)
                         s1.wcs2pix(image)
                         sr.wcs2pix(image)
                         temp, junk = sr.match(s1, maxnum=10000, useflags=no)
                         nonmatches = temp.stars

             # When done with a field/chip pair, add candidates to list
             for star in nonmatches:
                 newcands.append(star)

     # Final comparison with 2MASS 
     temp = Starlist(stars=newcands)
     temp.wcs2pix(image)
     tmass.wcs2pix(image)
     result = temp.nomatch(tmass, tol=pixtol)
     result.write(output)
     
                     

