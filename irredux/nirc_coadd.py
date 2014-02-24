from iqutils import *
import math
import time

def nirc_coadd(inlis, outfile):

    refimage = inlis[0]
    [raref,dcref,rotref]=get_head(refimage,["RA","DEC","ROTPOSN"])
    rotref = (rotref - 180.52) * math.pi / 180.0
    iraf.imcopy(refimage, "s%s" % refimage[6:])
    iraf.imcopy("nirc_coaddbpm.pl", "s%s.pl" % refimage.split('.')[0][6:])
    iraf.hedit("s%s" % refimage[6:], "BPM", "s%s.pl" %
               refimage.split(".")[0][6:], add=no, addonly=no, verify=no, 
               show=no, update=yes)
    coadd,exptime=get_head("s%s" % refimage[6:], ["FRMCOADD", "TINT"])
    iraf.hedit("s%s" % refimage[6:], "EXPTIME", coadd*exptime, 
               add=no, addonly=no, verify=no, show=no, update=yes)

    for i in range(1, len(inlis)):
        image = inlis[i]
        [ra,dc,rot] = get_head(image, ["RA","DEC","ROTPOSN"])
        rot = (rot - 180.52) * math.pi / 180.0
        ddec = dcref - dc; dra = raref - ra
        dra *= abs(math.cos((math.pi / 180.0)*(dc + dcref)/2.0))
        dra *= 3600 / 0.15; ddec *= 3600 / 0.15
        dxpix=(dra*math.cos(rot)-ddec*math.sin(rot))
        dypix=(-dra*math.sin(rot)-ddec*math.cos(rot))
        iraf.imshift(image, "s%s" % image[6:], dxpix, dypix, 
                     boundary_type="constant", constant=-100001.0)
        iraf.imexpr("a < -2000.0 ? 1.0 : 0.0", "s%s.pl" %
                    image.split(".")[0][6:], a="s%s" % image[6:])
        iraf.hedit("s%s" % image[6:], "BPM", "s%s.pl" % 
                   image.split(".")[0][6:], add=no, addonly=no, verify=no, 
                   show=no, update=yes)
        coadd,exptime=get_head("s%s" % image[6:], ["FRMCOADD", "TINT"])
        iraf.hedit("s%s" % image[6:], "EXPTIME", coadd*exptime, add=no,
                   addonly=no, verify=no, show=no, update=yes)

    #time.sleep(30.0)

    infiles = ""
    for image in inlis:
        infiles += "s%s," % image[6:]
    iraf.imcombine(infiles, outfile, combine="median", reject="avsigclip",
                   masktype="goodvalue", maskvalue=0, scale="exposure", 
                   zero="none", weight="exposure", expname="'EXPTIME'", 
                   rdnoise=100.0, gain=6.0) 
                   
    

