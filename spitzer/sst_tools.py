import pyfits, os
import numpy as np
from pyraf import iraf

from iqutils import *

EDGETOL = 100.0
ROOTDIR = "/Users/scenko/OBS/spitzer"

def do_all(field, epochs=[2], copyref=True):

    cwd = os.getcwd()

    # Create directories
    if not os.path.exists("%s/refs/%s" % (ROOTDIR, field)):
        os.mkdir("%s/refs/%s" % (ROOTDIR, field))
        os.mkdir("%s/refs/%s/raw" % (ROOTDIR, field))
        os.mkdir("%s/refs/%s/conv" % (ROOTDIR, field))
    if not os.path.exists("%s/%s" % (ROOTDIR, field)):
        os.mkdir("%s/%s" % (ROOTDIR, field))

    if copyref:
        # Copy reference images to relevant location
        shutil.copy("%s/%s_epoch1.fits" % (field, field),
                    "%s/refs/%s/raw/ref.fits" % (ROOTDIR, field))
        shutil.copy("%s/%s_cov_epoch1.fits" % (field, field),
                    "%s/refs/%s/raw/ref.cov.fits" % (ROOTDIR, field))
        shutil.copy("%s/%s_unc_epoch1.fits" % (field, field),
                    "%s/refs/%s/raw/ref.unc.fits" % (ROOTDIR, field))
        shutil.copy("%s/stamps.reg" % field, 
                    "%s/refs/%s/raw/" % (ROOTDIR, field))
        shutil.copy("%s/2mass.reg" % field,
                    "%s/refs/%s/raw/" % (ROOTDIR, field))
#       shutil.copy("%s/%s_epoch1_highres.fits" % (field, field),
#                   "%s/refs/%s/conv/ref.fits" % (ROOTDIR, field))
#       shutil.copy("%s/%s_cov_epoch1_highres.fits" % (field, field),
#                   "%s/refs/%s/conv/ref.cov.fits" % (ROOTDIR, field))
#       shutil.copy("%s/%s_unc_epoch1_highres.fits" % (field, field),
#                   "%s/refs/%s/conv/ref.unc.fits" % (ROOTDIR, field))
        #shutil.copy("%s/stamps.conv.reg" % field,
        #           "%s/refs/%s/conv/" % (ROOTDIR, field))
#       shutil.copy("%s/2mass.reg" % field,
#                   "%s/refs/%s/conv/" % (ROOTDIR, field))

    for ep in epochs:

        dobs = get_head("%s/%s_epoch%i.fits" % (field, field, ep), "DATE_OBS")
        dobs = dobs[:10].replace("-", "")
        os.mkdir("%s/%s/%s" % (ROOTDIR, field, dobs))
        os.mkdir("%s/%s/%s/raw" % (ROOTDIR, field, dobs))
        os.mkdir("%s/%s/%s/conv" % (ROOTDIR, field, dobs))

        # Copy new files
        shutil.copy("%s/%s_epoch%i.fits" % (field, field, ep),
                    "%s/%s/%s/raw/new.fits" % (ROOTDIR, field, dobs))
        shutil.copy("%s/%s_cov_epoch%i.fits" % (field, field, ep),
                    "%s/%s/%s/raw/new.cov.fits" % (ROOTDIR, field, dobs))
        shutil.copy("%s/%s_unc_epoch%i.fits" % (field, field, ep),
                    "%s/%s/%s/raw/new.unc.fits" % (ROOTDIR, field, dobs))
#       shutil.copy("%s/%s_epoch%i_highres.fits" % (field, field, ep),
#                   "%s/%s/%s/conv/new.fits" % (ROOTDIR, field, dobs))
#       shutil.copy("%s/%s_cov_epoch%i_highres.fits" % (field, field, ep),
#                   "%s/%s/%s/conv/new.cov.fits" % (ROOTDIR, field, dobs))
#       shutil.copy("%s/%s_unc_epoch%i_highres.fits" % (field, field, ep),
#                   "%s/%s/%s/conv/new.unc.fits" % (ROOTDIR, field, dobs))

        # Copy old files
        os.system("cp %s/refs/%s/raw/* %s/%s/%s/raw/" % (ROOTDIR, field, 
                   ROOTDIR, field, dobs))
        os.system("cp %s/refs/%s/conv/* %s/%s/%s/conv/" % (ROOTDIR, field,
                   ROOTDIR, field, dobs))

        # Run subs and html
        os.chdir("%s/%s/%s/raw" % (ROOTDIR, field, dobs))
        spitzer_sub("new.fits", "new.unc.fits", "new.cov.fits",
                    "ref.fits", "ref.unc.fits", "ref.cov.fits",
                    "sub.fits", stamps="stamps.reg", tmass="2mass.reg") 
        sub2html("cands.html", Starlist("cands.reg"), "new.fits", 
                 "ref.fits", "sub.fits")
        sub2html("cands_all.html", Starlist("cands_all.reg"), "new.fits",
                 "ref.fits", "sub.fits")
        os.chdir("../conv")
        #spitzer_sub("new.fits", "new.unc.fits", "new.cov.fits",
        #           "ref.fits", "ref.unc.fits", "ref.cov.fits",
        #           "sub.fits", stamps="stamps.conv.reg", tmass="2mass.reg") 
        #sub2html("cands.html", Starlist("cands.reg"), "new.fits", 
        #        "ref.fits", "sub.fits")
        #sub2html("cands_all.html", Starlist("cands_all.reg"), "new.fits",
        #        "ref.fits", "sub.fits")
        os.chdir(cwd)

    return

##########

def spitzer_sub(new, newunc, newcov, ref, refunc, refcov, out,
                stamps=None, tmass=None):

    # Remove nan from mosaics
    x = pyfits.open(new)
    y = np.nan_to_num(x[0].data)
    x[0].data = y
    new2 = "n%s" % new
    x.writeto(new2)

    # Find stars in new image
    os.system("$REDUCTION/runsex.pl %s 5.0 -weight %s" % (new2, newcov))

    # Find stars in reference
    os.system("$REDUCTION/runsex.pl %s 5.0 -weight %s" % (ref, refcov))

    # Create file for geotran
    stars = Starlist("%s.stars" % new2)
    refstars = Starlist("%s.stars" % ref)
    refstars.pix2wcs(ref)
    refstars.wcs2pix(new2)
    a,b = stars.match(refstars, maxnum=1000)

    # Create file for geotran
    b.pix2wcs(new2)
    b.wcs2pix(ref)
    refroot = ref.split(".")[0]
    outf = open("%s.match" % refroot, "w")
    for i in range(len(a)):
        outf.write("%10.3f%10.3f%10.3f%10.3f\n" % (a[i].xval, a[i].yval,
                   b[i].xval, b[i].yval))
    outf.close()

    # Geomap and geotran
    [naxis1, naxis2] = get_head(new, ["NAXIS1", "NAXIS2"])
    iraf.geomap("%s.match" % refroot, "%s.db" % refroot, 1, naxis1, 1, naxis2,
                fitgeometry="rotate", interactive=no)
    iraf.geotran(ref, "t%s" % ref, "%s.db" % refroot, "%s.match" % refroot)
    iraf.geotran(refunc, "t%s" % refunc, "%s.db" % refroot, 
                 "%s.match" % refroot)

    # Get stars for PSF matching
    if stamps != None:

        psfstars = Starlist(stamps)
        psfstars.pix2wcs(ref)
        psfstars.wcs2pix(new2)
        outf = open("stamps.lis", "w")
        for star in psfstars:
            outf.write("%10.3f%10.3f\n" % (star.xval, star.yval))
        outf.close()

    # Appropriate parameters
    iraf.iterstat(ref)
    update_head(ref, ["MEDSKY", "SKYSIG"], [iraf.iterstat.median,
                                            iraf.iterstat.sigma])
    [refskybkg, refskysig] = get_head(ref, ["MEDSKY", "SKYSIG"])
    tl = refskybkg - 10 * refskysig; tu = 30000.0
    iraf.iterstat(new)
    update_head(new, ["MEDSKY", "SKYSIG"], [iraf.iterstat.median, 
                                            iraf.iterstat.sigma])
    [newskybkg, newskysig] = get_head(new, ["MEDSKY", "SKYSIG"])
    il = newskybkg - 10 * newskysig; iu = 30000.0
    
    # Run hotpants
    hpcmd = "hotpants -inim %s -tmplim t%s -outim %s -tni t%s -ini %s -nsx 3 -nsy 3 -savexy %s.xy -ko 0 -bgo 0 -oni u%s -n t -tl %.2f -tu %.2f -il %.2f -iu %.2f -r 7.5 -rss 18.0" % (new2, ref, out, refunc, newunc, new2, out, tl, tu, il, iu)
    if stamps != None:
        hpcmd += " -ssf stamps.lis -afssc 0"
    os.system(hpcmd)
    #iraf.imarith(new2, "-", "t%s" % ref, out)

    # Create appropriate weight image
    #iraf.imexpr("sqrt(a**2 + b**2)", "u%s" % out, a=newunc, b=refunc)
    iraf.imarith(1, "/", "u%s" % out, "temp1.fits")
    iraf.imarith("temp1.fits", "/", "u%s" % out, "w%s" % out)
    os.remove("temp1.fits")

    # Pick up candidates
    os.system("$REDUCTION/runsex.pl %s 2.0 -weight w%s" % (out, out))
    stars = Starlist("%s.stars" % out)
    cands = []
    [nax1, nax2] = get_head(ref, ["NAXIS1", "NAXIS2"])
    stars.pix2wcs(new)
    stars.wcs2pix(ref)
    refu = pyfits.open(refunc)
    for star in stars:
        if star.xval > EDGETOL and star.xval < (nax1 - EDGETOL) and star.yval > EDGETOL and star.yval < (nax2 - EDGETOL) and refu[0].data[star.yval,star.xval] != 0 and star.fwhmw > 0.5 and star.fwhmw < 10.0:
            cands.append(star)
    scands = Starlist(stars=cands)

    # Filter out bright stars
    if tmass==None:
        scands.write("cands.reg")
    else:
        stmass = Starlist(tmass)
        stmass.wcs2pix(ref)
        s2cands = scands.nomatch(stmass)
        s2cands.write("cands.reg")
    
    # More comprehensive list
    os.system("$REDUCTION/runsex.pl %s 1.5" % out)
    stars = Starlist("%s.stars" % out)
    cands = []
    [nax1, nax2] = get_head(ref, ["NAXIS1", "NAXIS2"])
    stars.pix2wcs(new)
    stars.wcs2pix(ref)
    refu = pyfits.open(refunc)
    for star in stars:
        if star.xval > EDGETOL and star.xval < (nax1 - EDGETOL) and star.yval > EDGETOL and star.yval < (nax2 - EDGETOL) and refu[0].data[star.yval,star.xval] != 0 and star.fwhmw > 0.5 and star.fwhmw < 10.0:
            cands.append(star)
    scands = Starlist(stars=cands)
    scands.write("cands_all.reg")

    return

def sub2html(outname, cands, newimg, refimg, subimg, dx=50, dy=50):

    outf = open(outname, "w")

    # Write main observation info
    [obj, dobs, exptime] = get_head(newimg, ["OBJECT", "DATE_OBS", "EXPTIME"])
    outf.write("<html><head><title>Spitzer Observations of %s</title></head>\n" % obj)
    outf.write("<body bgcolor=#ffffff text=#000000>\n")
    outf.write("<h2>Spitzer Subtraction Results for %s</h2>\n" % obj)
    outf.write("<h5>Observation Date: %s</h5>\n" % dobs)
    outf.write("<h5>Exposure Time: %s</h5>\n" % exptime)

    # Start table
    outf.write("<p><p><p>\n")
    outf.write("<table border=0>\n")
    outf.write("<tr><td>Cand. Name</td><td>New Image</td><td>Reference Image</td>")
    outf.write("<td>Subtracted Image</td><td>Coordinates (J2000.0)</td>\n")
    outf.write("</tr>\n")

    cands.wcs2pix(subimg)
    
    # Loop over candidates
    for i in range(len(cands)):

        cand = cands[i]
        [cra, cdec] = cand.coo.deg()

        # Create new jpg's for display
        [xval, yval] = [cand.xval, cand.yval]
        [nax1, nax2] = get_head(subimg, ["NAXIS1", "NAXIS2"])
        trimreg = "[%i:%i,%i:%i]" % (max(xval-dx,1), min(xval+dx,nax1),
                                     max(yval-dy,1), min(yval+dy,nax2))
        iraf.imcopy("%s%s" % (subimg, trimreg), "stemp.fits")
        clist = Starlist(stars=[cand])
        clist.write("ctemp.reg")
        #os.system("ds9 -zscale stemp.fits -align yes -width 100 -height 100 -colorbar no -saveimage sub%03i.jpg -exit" % i)

        [[xval, yval]] = imwcs2pix(newimg, cra, cdec)
        [nax1, nax2] = get_head(newimg, ["NAXIS1", "NAXIS2"])
        news = pyfits.open(newimg)
        x = news[0].data[max(yval-dy,1):min(yval+dy,nax2),max(xval-dx,1):min(xval+dx,nax1)]
        news[0].data = x
        news.writeto("ntemp.fits")
        #os.system("ds9 -zscale ntemp.fits -align yes -width 100 -height 100 -colorbar no -saveimage new%03i.jpg -exit" % i)

        [[xval, yval]] = imwcs2pix(refimg, cra, cdec)
        [nax1, nax2] = get_head(refimg, ["NAXIS1", "NAXIS2"])
        refs = pyfits.open(refimg)
        x = refs[0].data[max(yval-dy,1):min(yval+dy,nax2),max(xval-dx,1):min(xval+dx,nax1)]
        refs[0].data = x
        refs.writeto("rtemp.fits")
        #os.system("ds9 -zscale rtemp.fits -align yes -width 100 -height 100 -colorbar no -saveimage ref%03i.jpg -exit" % i)

        os.system("rm ?temp.fits")
        os.system("rm ctemp.reg")

        outf.write("<tr>\n")
        outf.write("<td>%s-%03i</td>\n" % (newimg[:-5], i))
        outf.write('<td><img src="new%03i.jpg" width="%i"></td>\n' % (i, 4*dx))
        outf.write('<td><img src="ref%03i.jpg" width="%i"></td>\n' % (i, 4*dx))
        outf.write('<td><img src="sub%03i.jpg" width="%i"></td>\n' % (i, 4*dx))
        outf.write('<td>RA: %.5f, Dec: %.5f</td>\n' % (cra, cdec))

    outf.write("</table>\n")
    outf.write("</body></html>\n")
    outf.close()











    
