#!/usr/bin/env python

# D. Fox, June 2004
# Adapted from "cda.py" by A. Ptak, Aug 2001
# Based on code developed for XAssist, see http://xassist.pha.jhu.edu
# added USNOB/NOMAD and other switches,  A.J. Pickles, Jan 2005

# $Id: getusnob.py,v 1.1 2004/06/11 22:15:55 derekfox Exp derekfox $

import sys, string, os, traceback
import getopt, re, time, glob

from types import *
from httplib import HTTP

# Defaults / Globals
def_orefcat="asc_ref.cat"
def_catusno="usnob"
#def_catusno= "nomad"
def_magclr="R2"
def_bright="5.0"
def_faint="19.0"
def_epoch="2000.0"
def_silent="off"
def_radius="20"
def_email_base="astro.caltech.edu"
MAXTRY = 15
# USNO web site
domain = "www.nofs.navy.mil"
#ajp: USNO updates the script occasionally. 
#     Check by running search on web-page, then view source info.
script = "/cgi-bin/tfch3tI.cgi"


##############################

class dummy_ui:
    def error_mesg(self, txt):
        print txt

    def writeln(self, txt):
        print txt

##################################################
    
def download_web_page(domain, url):
    try:
        h = HTTP(domain)
        h.putrequest("GET", url)
        h.putheader('Accept', 'text/html')
        h.putheader('Accept', 'text/plain')
        h.endheaders()
    except:
        return(None)

    try:
        errcode, errmsg, headers = h.getreply()
    except:
        sys.stderr.write("Error in receiving response from " + domain + \
                         '\n')
        return None
    if errcode != 200:
        sys.stderr.write("Error in receiving response from " + domain + \
                         '\n')
        return None
    results = h.getfile().read()
    return(results)

##############################

def dec2sxg(coord,ra=0):

    try:
        coox=float(coord)
    except:
        sys.stderr.write("Bad coordinate string in dec2sxg\n")
        sys.exit(4)

    if ra:
        rah=int(coox/15.0)
        ram=int(60*(coox/15.0-rah))
        ras=60*(60*(coox/15.0-rah)-ram)
        return ("%02d" % rah,"%02d" % ram,"%07.4f" % ras)

    dcsgn='+'
    if re.search("^-",coord):
        dcsgn='-'

    absdec=abs(coox)
    dcd=int(absdec)
    dcm=int(60*(absdec-abs(dcd)))
    dcs=60*(60*(absdec-abs(dcd))-dcm)

    return (dcsgn+"%02d" % dcd,"%02d" % dcm,"%06.3f" % dcs)

##################################################

def main():

    # Parse Command line
    try:
        opts, args = getopt.getopt( sys.argv[1:], 
                     "hdqc:o:m:b:f:e:s:u:r:",
              ["help","debug","qonly","catusno=","orefcat=","magclr=","bright=","faint=",
	       "epoch=","silent=","user=","retry=","regionf=","mail="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    # Defaults
    debug = 0
    qonly = 0

    # User-controllable values
    orefcat=def_orefcat
    catusno=def_catusno
    magclr=def_magclr
    bright=def_bright
    faint =def_faint
    regfile=""
    epoch=def_epoch
    silent=def_silent
    radius=def_radius
    user=os.environ['USER']
    email=user+'@'+def_email_base
    retry=MAXTRY

    # Process details
    (xdir,xname)=os.path.split(sys.argv[0])
    pid=os.getpid()

    # Options parsing
    for opt, val in opts:
        # Usage info only
        if opt in ("-h", "--help"):
            usage()
            sys.exit(1)
        # Debug mode
        elif opt in ("-d", "--debug"):
            debug=1
        # Echo the web query only
        elif opt in ("-q", "--qonly"):
            qonly=1
            print "Only querying, no data files will be downloaded"
        # Name of input USNO catalog
        elif opt in ("-c", "--catusno"):
            catusno=val
        # Name of output reference catalog
        elif opt in ("-o", "--orefcat"):
            orefcat=val
        # Name of USNO/NOMAD colormag
        elif opt in ("-m", "--magclr"):
            magclr=val
        # Bright mag limit
        elif opt in ("-b", "--bright"):
            bright=val
        # Faint mag limit
        elif opt in ("-f", "--faint"):
            faint=val
        # Epoch of coordinates
        elif opt in ("-e", "--epoch"):
            epoch=val
        # Silent Output
        elif opt in ("-s", "--silent"):
            silent=val
        # Username
        elif opt in ("-u", "--user"):
            user=val
            email=user+'@'+def_email_base
        # Number of retry attempts
        elif opt in ("-r", "--retry"):
            retry=val
        # Region file
        elif opt in ("-g --regfile"):
            regfile=val
        # User's email address
        elif opt in ("--mail"):
            email=val
        else:
            sys.stderr.write("Unmatched option %s\n" % opt)

    # Try to guess a good region file name
    if len(regfile)==0:
        if orefcat.find('.cat')>0:
            regfile=orefcat.replace('.cat','.reg')

    # Parse arguments
    if len(args) < 2:
        usage()
        sys.exit(2)
    elif len(args)==2:
        (rain,dcin)=args[0:2]
    elif len(args)==3:
        (rain,dcin,radius)=args[0:3]
    else:
        sys.stderr.write("Ignoring extra arguments\n")
        (rain,dcin,radius)=args[0:3]

    # Parse coordinates (RA)
    re1=re.search("(\d+)[:\s](\d+)[:\s]([\d\.]+)",rain)
    re2=re.search("([\d\.]+)",rain)
    if re1:
        rah="%02d" % int(re1.group(1))
        ram="%02d" % int(re1.group(2))
        ras="%07.4f" % float(re1.group(3))
    elif re2:
        (rah,ram,ras)=dec2sxg(rain,ra=1)
    else:
        sys.exit("Failed to parse RA string %s\n" % rain)
    raout='+'.join([rah,ram,ras])
    ranice=':'.join([rah,ram,ras])

    # Parse coordinates (Dec)
    re1=re.search("([\+\-]?\d+)[:\s](\d+)[:\s]([\d\.]+)",dcin)
    re2=re.search("([\+\-]?[\d\.]+)",dcin)
    if re1:
        dsgn='+'
        if re.search("^-",dcin):
            dsgn='-'
        dcd=dsgn+("%02d" % abs(int(re1.group(1))))
        dcm="%02d" % int(re1.group(2))
        dcs="%06.3f" % float(re1.group(3))
    elif re2:
        (dcd,dcm,dcs)=dec2sxg(dcin,ra=0)
    else:
        sys.exit("Failed to parse Dec string %s\n" % dcin)
    dcout='+'.join([dcd,dcm,dcs])
    dcnice=':'.join([dcd,dcm,dcs])
        
    # First request: (colbits=cb_flg&) removed below, as hex output harder to read back in
    url = ("%s?ra=%s&dec=%s&equinox=J2000&epoch=%s&cextract=rect&" + \
           "rawid=%s&decwid=%s&wunits=Minutes&cat=%s&surims=None&" + \
           "getcat=yes&colbits=cb_id&colbits=cb_ra&slf=ddd.ddd/dd.ddd&" + \
           "colbits=cb_sigra&colbits=cb_mura&colbits=cb_smural&colbits=cb_fitpts&" + \
           "colbits=cb_mag&clr=%s&skey=mag&bri=%s&fai=%s&" + \
           "gzf=No&cftype=ASCII") % \
          (script,raout,dcout,epoch,radius,radius,catusno,magclr,bright,faint)

    if silent in ("off"):
        print "Sending %s request for %s, %s (%s')" % (catusno,ranice,dcnice,radius)

    if qonly:
        print "Full URL of initial query:"
        print url
        sys.exit(0)
    
    results = download_web_page(domain, url)

    if debug:
        test1="%s-%d-1.html" % (xname,pid)
        fil=open(test1, 'w')
        fil.write(results)

    if type(results)==NoneType or len(results)<100:
        sys.exit("Something wrong, results too short\n")

    # Successful first request:  Get link for "progress" web page
    urltrk=""
    inrange=0
    lines=string.split(results,'\n')
    for line in lines:
        re1=re.search("refresh.+(http:.+\.html)",line)
        if re1:
            urltrk=re1.group(1)
            break

    if len(urltrk)>1:
        urlcat=urltrk.replace("fch.html",catusno)
        if silent in ("off"):
           print "Tracking results through %s" % urltrk
           print "Catalog should appear as %s" % urlcat
    else:
        sys.exit("Trouble identifying tracking URL\n")

    re2=re.search("http://([^/]+)(/.+)$",urlcat)
    if re2:
        dom2=re2.group(1)
        url2=re2.group(2)
    else:
        sys.exit("Problem parsing target URL\n")

    # Wait a bit
    time.sleep(10)

    # Start trying to retrieve the catalog
    itry=0
    while 1:
        results2 = download_web_page(dom2, url2)

        if results2:
            lines=string.split(results2,'\n')
            if len(lines)>=26:
                break

        itry +=1
        if itry>retry:
            print "Giving up"
            return

        time.sleep(10)

    # Add comment flags to header
    for i in range(len(lines)):
        if re.search('^\s*\d\d\d\d',lines[i]):
            break
        lines[i]='#'+lines[i]
    nhead=i

    # Write the catalog
    fcat=open(orefcat, 'w')
    for line in lines:
        # Remove leading spaces
        while line.startswith(' '):
            line=line[1:]
        fcat.write(line+'\n')
    fcat.close()

    # Write the region file (if requested)
    if len(regfile)>0:
        freg=open(regfile,'w')
        freg.write("# Region file format: DS9 version 3.0\n")
        freg.write("global color=blue\n")
        for line in lines[nhead:]:
            els=line.split()
            if len(els)<16:
                continue
            rstar=5.0*(20.0/float(els[16]))**2
            regline=("fk5;circle(%s,%s,%.2f\") # text={%s} " + \
                     "B1MAG={%s} R1MAG={%s} B2MAG={%s} R2MAG={%s} " + \
                     "I2MAG={%s}\n") % (els[1],els[2],rstar,els[0],
                                   els[12],els[13],els[14],els[15],els[16])
            freg.write(regline)
        freg.close()

    # Successful completion
    if silent in ("off"):
	    print "Done."

    sys.exit(0)

##################################################    

def usage():

    (xdir,xname)=os.path.split(sys.argv[0])

    print "Usage: %s [-c catusno] [-o orefcat] [-m magclr] [-b brtmag] [-f fntmag] [-e epoch] [-s flag] [-r retry] [-g regfile] <ra> <dec> [size]" % xname
    print "    <ra> and <dec> are sexagesimal hours/deg or decimal deg/deg"
    print "    <size> is length of a side in arcmin (%s)" % def_radius
    print "    -c <usnob/nomad> gives alternate input catalog (%s)" % def_catusno
    print "    -o <orefcat> gives alternate name for ouput reference catalog (%s)" % def_orefcat
    print "    -m <R2/K> gives alternate USNO or NOMAD sort magnitude (%s)" % def_magclr
    print "    -b <val> sets bright search magnitude (%s)" % def_bright
    print "    -f <val> sets  faint search magnitude (%s)" % def_faint
    print "    -e <epoch> gives the desired input coordinate epoch (%s)" % def_epoch
    print "    -s <off/on> sets the silent output flag (%s)" % def_silent
    print "    -r <val> sets the number of retries (%s)" % MAXTRY
    print "    -g <regfile> gives alternate name for ouput region file (%s)" % def_orefcat.replace('.cat','.reg')

##############################

if __name__ == "__main__":

    main()

