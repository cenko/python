#!/usr/bin/env python

# D. Fox, Jan 2004
# Adapted from "cda.py" by A. Ptak, Aug 2001
# Based on code developed for XAssist, see http://xassist.pha.jhu.edu

# $Id: getxdss.py,v 1.2 2004/03/23 02:41:24 derekfox Exp derekfox $

import sys, string, os, traceback
import getopt, re, time, glob

from httplib import HTTP
from ftplib import FTP

# Defaults / Globals
def_source="XDSS"
def_epoch="2000"
def_radius="20"
def_email="derekfox@astro.caltech.edu"

##############################

class dummy_ui:
    def error_mesg(self, txt):
        print txt

    def writeln(self, txt):
        print txt
        

class xa_ftp:
    def __init__(self, xassist=None, host=None):
        self.xassist = xassist
        if not xassist:
            self.ui = dummy_ui()
        else:
            self.ui = xassist.get_ui()
        self.host = host
        self.connected = 0
        self.ftp = None
        if host != None:
            self.open()

    def close(self):
        self.ftp = None
        self.connected = 0

    def get(self, dir, filename, localfile=None):
        self.ui.writeln("Trying to download %s/%s..." % (dir, filename))
        if not self.connected:
            if not self.open():
                return(1)
        if not localfile:
            localfile = filename
        try:
            fil = open(localfile, "wb")
        except:
            self.ui.error_mesg("Could not open " + localfile + " for writing")
            return(1)

        try:
            self.ftp.retrbinary("RETR " + dir + '/' + filename,
                                fil.write, 8192)
            self.ui.writeln("Downloaded " + filename)
        except:
            self.ui.error_mesg("Error downloading " + filename)
            #traceback.print_exc()
            #self.ui.writeln("Error was " + sys.exc_info()[1])
            fil.close()
            return(1)
        return(0)
    
    def glob(self, dir, wildcards):
        if not(self.connected):
            if self.open():
                return(1)
        lines = []
        maxlines = 1000
        for wildcard in wildcards:
            cmd = "LIST " + dir + "/" + wildcard
            print cmd
            conn = self.ftp.transfercmd(cmd)
            fp = conn.makefile('rb')
            i = 0
            while i < 1000:
                line = fp.readline()
                if not line: break
                #if line[-1:] == '\n': line = line[:-1]
                #if line[-1:] == '\r': line = line[:-1]
                line = string.strip(line)
                vals = string.split(line)
                file = os.path.basename(vals[-1])
                print "file = " + file
                lines.append(file)
                i = i + 1
            fp.close();
            conn.close();
            self.reopen()
        return(lines)

    def mget(self, dir, wildcards):
        flist = self.glob(dir, wildcards)
        if flist == 1:
            self.ui.writeln("No files found")
            return(1)
        #print flist
        status = 0
        for file in flist:
            if self.get(dir, file):
                status = 1
        return(status)
            
    def open(self):
        self.ui.writeln("Opening connection to " + self.host)
        try:
            self.ftp = FTP(self.host)
        except:
            self.ui.error_mesg("Initialization of FTP object failed")
            return(1)
	self.ftp.set_pasv(1) # Add 9/17/01, A.P.
        try:
            self.ftp.login()
        except:
            self.ui.error_mesg("Log into " + self.host + " failed")
            return(1)
        self.connected = 1
        return(0)
    
    def reopen(self):
        self.close()
        self.open()

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
    print "Placed request, downloading web page..."
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

def dec2sxg(coord,ra=1):

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
    dcm=int(60*(absdec-dcd))
    dcs=60*(60*(absdec-dcd)-dcm)

    return (dsgn+"%02d" % dcd,"%02d" % dcm,"%06.3f" % dcs)

##################################################

def main():

    # CADC www site, it may work to change this to a mirror URL
    domain = "cadcwww.dao.nrc.ca"
    script = "/cadcbin/getdss/1"

    # Parse Command line
    try:
        opts, args = getopt.getopt( sys.argv[1:], 
                     "hdqn:e:u:",
                    ["help","debug","qonly","name=","epoch=","user="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    # Defaults
    debug = 0
    qonly = 0
    next = 1
    source=def_source
    epoch=def_epoch
    radius=def_radius
    email=def_email
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
        # Name of source
        elif opt in ("-n", "--name"):
            source=val
        # Epoch of coordinates
        elif opt in ("-e", "--epoch"):
            epoch=val
        # User's email address
        elif opt in ("-u", "--user"):
            email=val
        else:
            sys.stderr.write("Unmatched option %s\n" % opt)

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

    # Parse coordinates
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
    raout='%20'.join([rah,ram,ras])
    ranice=':'.join([rah,ram,ras])

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
        (dcd,dcm,dcs)=dec2sxg(dcin)
    else:
        sys.exit("Failed to parse Dec string %s\n" % dcin)
    dcout='%20'.join([dcd,dcm,dcs])
    dcnice=':'.join([dcd,dcm,dcs])
        
    # First request
    url = ("%s?ra=%s&dec=%s&epoch=%s&radius2=%s&dss_select=XDSS&" + \
           "output_mode=FITS&colour_selection=ALL&email_user=%s&" + \
           ".cgifields=colour_selection&.cgifields=dss_select") % \
           (script,raout,dcout,epoch,radius,email)

    print "Sending request for %s, %s" % (ranice,dcnice)

    if qonly:
        print "Full URL of initial query:"
        print url
        sys.exit(0)
    
    results = download_web_page(domain, url)

    if debug:
        test1="%s-%d-1.html" % (xname,pid)
        fil=open(test1, 'w')
        fil.write(results)

    if len(results)<100:
        sys.exit("Something wrong, results too short\n")

    # Successful first request:  Get link for "progress" web page
    url2=""
    inrange=0
    lines=string.split(results,'\n')
    for line in lines:
        if inrange:
            re1=re.search("(http:[^<]+)</A>",line)
            if re1:
                url2=re1.group(1)
                break
        elif string.find(line,"progress")>0:
            inrange=1

    if len(url2)>1:
        print "Tracking results through %s" % url2
    else:
        sys.exit("Trouble finding url2\n")

    re2=re.search("http://([^/]+)(/.+)$",url2)
    if re2:
        dom2=re2.group(1)
        url2=re2.group(2)
    else:
        sys.exit("Problem parsing target URL\n")

    # Wait for successful processing of query
    time.sleep(10)
    
    while 1:
        results2 = download_web_page(dom2, url2)

        if debug:
            if results2:
                test2="%s-%d-2.html" % (xname,pid)
                fil=open(test2,'w')
                fil.write(results2)
                fil.close()
            
        if results2:
            if string.find(results2,"Process terminated with success")>0:
                break
            if string.find(results2,"Process terminated")>0:
                sys.stderr.write("Request has failed:\n")
                sys.stderr.write(results2)
                sys.exit(5)

        time.sleep(10)

    if len(results2)<100:
        sys.exit("Something wrong, results2 too short\n")

    # Successful second request: Parse for "data" link
    lines=string.split(results2,'\n')
    for line in lines:
        if string.find(line,"compressed tar file")>0:
            re1=re.search("(ftp://.+\.gz)\">",line)
            if re1:
                ftpurl=re1.group(1)
                break
            else:
                sys.exit("Trouble finding ftpurl\n")
                
    print "Downloading files as %s" % ftpurl

    re2=re.search("ftp://([^/]+)(/.+/)([^/]+)$",ftpurl)
    if re2:
        ftphost=re2.group(1)
        ftpdir=re2.group(2)
        ftpfile=re2.group(3)
    else:
        sys.exit("Problem parsing target ftp URL\n")

    ftp = xa_ftp(host=ftphost)
    for i in range(1, 25):
        print "Download attempt %d" % i
        status = ftp.get(ftpdir, ftpfile)
        if status:
            print "Download failed"
            # Wait for 30 seconds
            print "Sleeping for 30 seconds..."
            time.sleep(30)
            ftp.reopen()
        else:
            break

    # Uncompress/Untar file
    os.system("tar -xzf %s" % ftpfile)

    # Rename file(s)
    files=glob.glob("coord_xdss*.fits")
    if len(files)>0:
        os.remove(ftpfile)
        for file in files:
            re1=re.search("xdss_(.)_",file)
            band=re1.group(1)
            os.rename(file,"%s_%s.fits" % (source,band))

    # Successful completion
    print "Done."
    sys.exit(0)

##################################################    

def usage():

    (xdir,xname)=os.path.split(sys.argv[0])

    print "Usage: %s [-n name] [-e epoch] [-u email] <ra> <dec> [size]" % xname
    print "    <ra> and <dec> are sexagesimal hours/deg or decimal deg/deg"
    print "    <size> is length of a side in arcmin (%s)" % def_radius
    print "    -n name gives the source name for output files (%s)" % def_source
    print "    -e epoch gives the epoch of the coordinates (%s)" % def_epoch
    print "    -u email gives the user email address for problems (%s)" % def_email

##############################

if __name__ == "__main__":

    main()

