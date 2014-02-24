#! /usr/bin/env python

import psycopg2
import sys, datetime, os
from urllib2 import urlopen

# Grab CandID from argv
if len(sys.argv) != 2:
    print "Wrong number of arguments!  Format should be:"
    print "pymp_ptf_sub.py <subtraction ID>"
    sys.exist(1)
else:
    subid = sys.argv[1]

# Get RA, Dec, and UT from database
conn = psycopg2.connect(database="subptf", user="subptf", 
        password="p33d$kyy", host="scidb2.nersc.gov", port=5432)
cur = conn.cursor()

cur.execute("select s.ra_c, s.dec_c, p.utc_obs, s.height, s.width from subtraction as s, proc_image as p where s.id = %s and s.proc_image_id = p.id;" % subid)
reply = cur.fetchall()
[radeg, dcdeg, tobs, height, width] = reply[0]
cur.close()
#print reply

# Convert to appropriate formats
rah=int(radeg/15.0)
ram=int(60*(radeg/15.0-rah))
ras=60*(60*(radeg/15.0-rah)-ram)

dcsgn="+"
if cmp(dcdeg,0) < 0:
    dcsgn="-"
absdec = abs(dcdeg)
dcd=int(absdec)
dcm=int(60*(absdec-abs(dcd)))
dcs=60*(60*(absdec-abs(dcd))-dcm)

rasxg = "%02i:%02i:%07.4f" % (rah,ram,ras)
dcsxg = dcsgn+"%i:%02i:%06.3f" % (dcd,dcm,dcs)

radius = max(float(height) / 2.0, float(width) / 2.0) * 3600.0 * 1.10

# PyMPChecker will host the regions file at a url determined by the input 
# parameters. Convert the input parameters to what PyMPChecker uses in the url.
raf = rasxg.replace(":", "%3A")
dcf = dcsxg.replace(":", "%3A")
tobss = "%04i-%02i-%02i %02i:%02i:%02i" % (tobs.year, tobs.month, tobs.day,
         tobs.hour, tobs.minute, int(tobs.second))
tobsf = str(tobss).replace("-", "%2F").replace(":", "%3A").replace(" ", "+")
input_url = ("http://dotastro.org/PyMPC/PyMPC/?in_1=" + raf + "&in_2="
    + dcf + "&in_3=" + tobsf + "&in_4=" + str(int(radius)))
#print input_url

# Submit the parameters to PyMPChecker in the form of the url. It will run the 
# calculations (takes a few seconds). Then the client retrieves the html, 
# searches for the link to the regions file, and uses wget to pull down the 
# file.
response = urlopen(input_url)
html = response.read()
regions_line_start = html.find("/mpc_data/ds9_regions/")
regions_line_end = html.find("\">DS9 Region File<")
ds9_regions_url = ("http://dotastro.org" + 
    html[regions_line_start:regions_line_end])
regions_file = "%s.reg" % subid
os.system("wget -q -O " + regions_file + " " + ds9_regions_url)

# Check resulting region file
fout = open("%s.reg" % subid)
lines = fout.readlines()
if len(lines)==2:
    # No asteroid matches
    print 0
else:
    for i in range(2,len(lines)):
        print lines[i].strip()

os.remove("%s.reg" % subid)
