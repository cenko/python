#!/usr/bin/env python

"""
Python tools for Palomar Transient Factory.

SBC 2013/02/06
"""

# "Standard" Python packages
import pgdb, datetime, urllib2, os, re, shutil
from lxml import etree
from mx.DateTime import *
from types import *

try:
    import MySQLdb
except ImportError:
    print "Could not import MySQLdb: Can't communicate with PTEL database"

try:
    import sqlite3
except ImportError:
    print "Could not import sqlite3: Can't update IAU/CRTS databases"

# Custom Python packages
import oscar_defs

try:
    from iqutils import *
except:
    print "Could not import iqutils: Some functionality may be lost!"

################################
# Global Variable Definitions
################################

# PTF Database Properties
[PTFdbname, PTFdbuser, PTFdbpass, PTFdbip] = \
 ["ptfcands", "tcp", "classify", "131.215.193.69"]

# P60 Target list default properties
[P60PTFoffset, P60PTFroi, P60PTFproperty] = \
 [[180.0, -180.0], [1025, 1, 1024, 2048], {"NOFOCUS":1}]

P60Programs = {0: ["2009A-I0113", "'PTF: Follow-Up'", "'Kulkarni'", 1.00, "white"],
               1: ["2009B-I0001", "'PTF: P60 Transient Vetting'", "'Kulkarni'", 
                   1.015, "red"],
               2: ["2009B-I0002", "'PTF: Transients in the Local Universe'", 
                   "'Kasliwal'", 1.1, "orange"],
               4: ["2009B-I0004", "'PTF: Spectral Vetting'", "'Kulkarni'",
                   0.00, "white"],
               6: ["2009B-I0006", "'PTF: Luminous Supernovae'", "'Quimby'",
                   1.02, "yellow"],
               7: ["2009B-I0007", "'PTF: No Transient Left Behind'", "'Nugent'",
                   0.00, "white"],
               8: ["2009B-I0008", "'PTF: II-P Cosmology'", "'Poznanski'",
                   0.00, "white"],
               9: ["2009B-I0009", "'PTF: Type Ia SNe in the Low-z Universe'",
                   "'Sullivan'", 0.00, "white"],
               10: ["2009B-I0010", "'PTF: Blindly Selected SNe'", "'Quimby'",
                    0.00, "white"],
               11: ["2009B-I0011", "'PTF: Cataclysmic Variables'", "'Shara'",
                    0.00, "white"],
               12: ["2009B-I0012", "'PTF: Core-collapse SNe'", "'Gal-Yam'",
                    1.02, "cyan"],
               13: ["2009B-I0013", "'PTF: Evolution and Dispersion in the Spectral Properties of Ia SNe'",
                    "'Ellis'", 0.00, "white"],
	       14: ["2009B-I0014", "'PTF: Examining evolution in the properties of Type Ia SNe'",
                    "'Sullivan'", 0.00, "white"],
               15: ["2009B-I0015", "'PTF: LCOGT Light Curves'", "'Howell'",
                    0.00, "white"],
               16: ["2009B-I0016", "'PTF: Evolution of the spectral properties of SN Ia'", 
		    "'Sullivan'", 0.00, "white"],
	       17: ["2009B-I0017", "'PTF: Tidal Disruption Candidate'", 
		    "'Bloom'", 0.00, "white"],
	       18: ["2009B-I0018", "'PTF: GALEX Seredipitous Observations'",
	            "'Ofek'", 1.02, "green"],
	       19: ["2009B-I0019", "'PTF: Follow-Up of Swift-UVOT Targets'",
	            "'Ofek'", 0.00, "white"],
	       20: ["2009B-I0020", "'PTF: Low-Resolution Spectral Vetting and Follow-Up'",
		    "'Quimby'", 0.00, "white"],
	       21: ["2009B-I0021", "'LSNe Follow-Up'", "'Quimby'", 0.00, "white"], 
               22: ["2009B-I0022", "'NIR Photometry'", "'Cenko'", 0.00, "white"]
               }

P60Priority = {5: 80000, 4: 60000, 3: 55000, 2:50000, 1: 40000, 0: 30000}
P60Weights = {"Age": 1.0, "Magnitude": 0.01} 
P60RefWeights = {"airmass": [1.7,0.25], "moondeg": [120.0,0.1]}
P60StdWeights = {"airmass": [2.5,0.25], "moondeg": [150.0,0.1]}
P60VetWeights = {"airmass": [2.5,0.25], "moondeg": [170.0,0.1]}
P60Filters = {"g": "gpr", "r": "rpr", "i": "ipr", "z": "zpr",
              "b": "B", "v": "V"}
P60_OSS = "oscar_auto.py -v -r -o -f 4.0 -s 10000.0"
#P60_OSS = "/usr/local/bin/python /home/p60/oss/oscar_auto.py -v -r -o -f 4.0 -s 10000.0"

# PAIRITEL Filters
PTELFilters = {"J": "J", "H": "H", "K": "Ks"}

# PAIRITEL Default parameters for PTF targets
[PTELPTFprojid, PTELPTFpriority, PTELPTFexptime] = \
 ["'PTF'", 90.0, 1800.0]

# PAIRITEL Default parameters for all targets
[PTELseebad, PTELtranbad, PTELairbad, PTELdith, PTELmno, PTELrasz, PTELdecsz,
 PTELmaxoff] = [5.0, 0.2, 3.0, 1, 1.0, 10.0, 10.0, 50.0]

# PAIRITEL Scheduler database info
[PTELdbip, PTELdb, PTELdbuser, PTELdbpass] = \
 ["192.33.141.15", "Mt_Hopkins", "obs", "ilove2mass"] 

# PTF Follow-up Marshal Website
#[MarshalURL, MarshalBaseURL, MarshalUser, MarshalPass] = \
#  ["http://navtara.caltech.edu/cgi-bin/ptf/list_phot2.cgi",
#  "http://navtara.caltech.edu/cgi-bin/ptf/", "ptf", "discover"]
[MarshalURL, MarshalBaseURL, MarshalUser, MarshalPass] = \
  ["http://ptf.caltech.edu/cgi-bin/ptf/transient/list_phot2.cgi",
  "http://ptf.caltech.edu/cgi-bin/ptf/transient", "p60", "iptf060"]

# Color / Priority Mapping (for HTML)
#ColorMap = {5: "red", 4: "orange", 3: "yellow", 2: "green", 1: "cyan"}

# IAU SNe Table
[IAUSNURL, IAUSNtable] = \
   ["http://www.cbat.eps.harvard.edu/lists/RecentSupernovae.html",
    "/home/cenko/python/PTF/IAUSNe.sql"]

# CRTS OT Table
[CRTSURL, CRTStable] = \
    ["http://nesssi.cacr.caltech.edu/catalina/Allns.html",
     "/home/cenko/python/PTF/CRTS.sql"]

#########################################################################

class MarshalRequest:

    """Object to handle all the details of a marshal observation request"""

    ##########################

    def __init__(self, element=None):

        # Define the entire set of variables
        self.name = "PTF"
        self.ra = None
        self.dec = None
        self.filter = None
        self.mag_guess = None
        self.snr = None
        self.mag_lim = None
        self.nvisit = 1
        self.spacing = DateTimeDelta(1/24.0)
        self.group = None
        self.classification = "Unknown"
        self.disc_date = utc()
        self.programs = {}
        self.comment = None
        self.ids = None

        # Parse the etree element if requested
        if (element != None):
            self.parse(element)

    #######################

    def parse(self, element):

        """Parse a row from the list_phot.html page (in etree format) and 
        fill in appropriate request variables."""

        # Request ID
        self.id = element[0][0].attrib["value"]

        # Name
        self.name = "PTF%s" % element[1].text

        # RA / Dec
        self.ra = float(element[2].text)
        self.dec = float(element[3].text)

        # Filter
        self.filter = element[4].text

        # Magnitude Guess
        if element[5].text=="N/A":
            self.mag_guess = 19.0
        else:
            self.mag_guess = float(element[5].text.strip(">"))

        # SNR
        self.snr = float(element[6].text)

        # Limiting Magnitude
        self.mag_lim = float(element[7].text)

        # Number of visits
        self.nvisit = int(element[8].text)

        # Visit spacing
        self.spacing = DateTimeDeltaFrom(element[9].text) 

        # Group
        self.group = int(element[10].text)

        # Classification
        self.classification = element[11].text

        # Discovery Date
        self.disc_date = DateTimeFrom(element[12].text)

        # Program / Priority
        #if (element[13].text != None):
            #programs = element[13].text.split(","); programs.pop()
            #priorities = element[14].text.split(","); priorities.pop()
            #while len(programs):
                #self.programs[int(programs.pop())] = int(priorities.pop())
        self.program = int(element[13].text)
        self.priority = abs(int(element[14].text)) # Sets -1 -> 1

        # Comment
        self.comment = element[15].text

#########################################################################

def update_tlist(ptffile='ptf/p60-ptf.lis', outfile='targets/p60-targets.lis',
                 url=MarshalURL):

    """Get new PTF targets, insert them into P60 list (after removing old
       one), and respond to marshal with list of checked out targets."""

    configure_http()
    targets=fetch_marshal_targets(filters=P60Filters)
    new_p60_targets(targets, outfile=ptffile)
    targlis=oscar_defs.oscar_target_list(catalog=outfile)
    i=0
    while i < len(targlis):
        for key,item in P60Programs.iteritems():
            if targlis[i].keyword['P60PRID'] == item[0]:
                targlis.__delitem__(i)
                i-=1
        i+=1

    targlis.write_targets(outfile)
    os.system("less %s >> %s" % (ptffile, outfile))
    #p60_checkout(outfile, targets)

    return

#########################################################################

def p60_checkout(p60file, targets, mjd=None, standards="targets/p60-standard.lis",
                 tmpfile="ptf/oscar_tmp.lis", url=MarshalURL):

    """Simulate list of P60 observations for the night, and check-out those
       targets likely to be observed from marshal."""

    # If MJD not specified, round to the nearest one
    if mjd==None:
        mjd = round(utc().mjd)

    # Simulate P60 target list for that night
    cmd = "%s -d -t %i %s %s > %s" % (P60_OSS, mjd, p60file, standards, tmpfile)
    scmd = os.popen(cmd, 'r', -1)
    scmd.readlines()

    # Parse simulated target list
    otargs = []
    tmp = open(tmpfile, 'r')
    lines = tmp.readlines()
    i = 0
    for i in range(len(lines)):
        line = lines[i]
        re1 = re.search("Asked to observe '(PTF\d\d[a-z]+)'", line)
        if re1:
            name = re1.group(1)
            line2 = lines[i+1]
            re2 = re.search("Exposure:\s+(\d)x\s+(\w+)\s+for\s+(\w+)", line2)
            otargs.append([name, int(re2.group(1)), re2.group(2), 
                           float(re2.group(3))])
            i+=3
        else:
            i+=1

    # Match lists and retrieve IDs
    ids = []
    for [name, number, filter, exptime] in otargs:
        for target in targets:
            if (name == target.name)  and (filter == P60Filters[target.filter]):
                ids.append(target.id)

    # Alert marshal to check out targets
    data = "instrumentid=6&commit=yes"
    for id in ids:
        data += "&ids=%s" % id
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req)
    ret = response.read()

    return
    
#########################################################################

def ptel_do_all(url=MarshalURL, filters=PTELFilters, projid=PTELPTFprojid, 
                priority=PTELPTFpriority, exptime=PTELPTFexptime):

    """Get new PTF targets, insert them into the PAIRITEL scheduler
       database."""

    configure_http()
    targets = fetch_marshal_targets(url=url, filters=filters)
    update_ptel_queue(targets, projid=projid, priority=priority, 
                      exptime=exptime)
    add_lowptel_targets()

    return

#########################################################################

def new_p60_targets(targets, outfile="p60-ptf.lis"):

    """Take a target list formatted as elementtree elements and  
       convert to P60 format."""

    # Need to set current time
    tnow = utc()

    # Calculate P60 Priorities
    targetdict = {}
    for target in targets:
        if P60Programs.has_key(target.program):
            target.p60priority = max(P60Priority[target.priority] * \
                                     P60Programs[target.program][3] - \
                                     P60Weights["Age"] * \
                                     (tnow - target.disc_date).days - \
                                     P60Weights["Magnitude"] * \
                                     target.mag_guess, 0)
        else:
            target.p60priority = 0.0

        # Create Target Dictionary
        if targetdict.has_key(target.name):
            targetdict[target.name].append(target)
        else:
            targetdict[target.name] = [target]

    # New list to save targets in
    p60targets = []

    # Loop through sources
    for name, exposures in targetdict.iteritems():

        tpriority = -1000000.0
        mexposure = None
        oexposures = {}
        reference = False

        for exposure in exposures:

            if exposure.p60priority > tpriority:
                tpriority = exposure.p60priority
                mexposure = exposure
            if exposure.mag_lim > 22.0:
                reference = True

            if not oexposures.has_key(exposure.filter):
                oexposures[exposure.filter] = exposure
            elif oexposures.has_key(exposure.filter) and exposure.priority > oexposures[exposure.filter].priority:
                oexposures[exposure.filter] = exposure

        p60exposures = []
        for null, oexposure in oexposures.iteritems():
            
            if (oexposure.mag_lim > 22.0):
                p60exposures.append(oscar_defs.oscar_exposure(P60Filters[oexposure.filter], 180.0, 3))
            elif oexposure.mag_guess < 18:
                p60exposures.append(oscar_defs.oscar_exposure(P60Filters[oexposure.filter], 60.0, 1))
            elif oexposure.mag_guess < 20:
                p60exposures.append(oscar_defs.oscar_exposure(P60Filters[oexposure.filter], 120.0, 1))
            else:
                p60exposures.append(oscar_defs.oscar_exposure(P60Filters[oexposure.filter], 180.0, 1))


        # Create oscar_target
        if reference:
            weights = P60RefWeights
        elif (mexposure.program==1) and (mexposure.priority==5):
            weights = P60VetWeights
        else:
            weights = P60StdWeights

        p60target = oscar_defs.oscar_target(mexposure.name, mexposure.ra, mexposure.dec, 
                                            p60exposures, priority=tpriority, 
                                            timing={"monitor":['%.2f' %
                                            (12.0 / mexposure.spacing.hours)]},
                                            repeat=mexposure.nvisit, offset=P60PTFoffset, 
                                            roi=P60PTFroi, keyword={"P60PRID":
                                            P60Programs[mexposure.program][0], "P60PRNM":
                                            P60Programs[mexposure.program][1], "P60PRPI":
                                            P60Programs[mexposure.program][2], "OBJTYPE":
                                            "Transient"}, property=P60PTFproperty,
                                            weight=weights)

        # Add to target list
        p60targets.append(p60target)

    # Write out new target list (for now)
    p60lis=oscar_defs.oscar_target_list()
    p60lis.targets=p60targets
    p60lis.write_targets(outfile)

    return

##########################################################################

def configure_http(url=MarshalBaseURL,user=MarshalUser,passw=MarshalPass):

     """Configure HTTP to use user/password combo for appropriate URLs"""

     x = urllib2.HTTPPasswordMgrWithDefaultRealm()
     x.add_password(None, url, user, passw)
     auth = urllib2.HTTPBasicAuthHandler(x)
     opener = urllib2.build_opener(auth)
     urllib2.install_opener(opener)

     return

##########################################################################

def fetch_marshal_targets(url=MarshalURL, filters=None):

     """Query PTF Marshal for targets requiring photometric follow-up.
     Will only grab targets for which request filter is a key in the 
     specified filters dictionary."""

     # Grab HTML from Marshal
     flob = urllib2.urlopen(url)
     s = flob.read()
     html = etree.HTML(s)

     # Convert to useful format (etree Element)
     targets = etree.Element("Targets")
     trows = html.findall('.//tr')
     for row in trows:
         if (len(row)==16) and (len(row[0])>0):
             targets.append(row)

     # Add targets of interest
     newtargets=[]
     for target in targets:

         # Only want targets with appropriate filters
         if target[4].text in filters:
             try:
                 newtarget = MarshalRequest(element=target)
                 newtargets.append(newtarget)
             except:
                 print "Problem with target %s" % target[1].text

     newtargets.sort()
     return newtargets

#########################################################################

def add_lowptel_targets(url=MarshalURL, filters=P60Filters, mag_cut=18):

    """Search through all targets requiring P60 photometry for bright (< 18 mag)
    sources.  Add these to the PTEL queue with low priority."""

    newtargs = []
    alltargs = fetch_marshal_targets(url=url, filters=filters)

    tnow = utc()
    
    for target in alltargs:

        if target.mag_guess < mag_cut:

            target.priority = 1
            newtargs.append(target)

    update_ptel_queue(newtargs, priority=10.0)

    return

#########################################################################

def update_ptel_queue(targets, projid=PTELPTFprojid, priority=PTELPTFpriority,
                      exptime=PTELPTFexptime):

    """Take a list of MarshalRequests and insert them into the PAIRITEL
    queue.  For now, if source is already in the queue, will just ignore.
    If source is not in the queue, will request 20 observations spaced
    by 3 days.  In the future the details of the requests should be 
    handled by the marshal and all individual requests will be added 
    to the queue."""

    # Keep track of what we've scheduled already
    done_targets = []

    # Need to talk to PAIRITEL database
    try:
        db=MySQLdb.connect(db=PTELdb,user=PTELdbuser,host=PTELdbip,
                           passwd=PTELdbpass)
    except:
        print "Error connecting to PTEL database at Mt. Hopkins"
        return
    cursor=db.cursor()

    # Loop through targets
    for target in targets:

        # If already seen source before, move on
	if "%s_Marshal" % target.name in done_targets:
            continue

        # Check and see if already in the database
        cursor.execute("SELECT ObjID FROM object WHERE Name = '%s_Marshal'" %
                       target.name)
        data = cursor.fetchall()
        
        # If already in the database, skip for now
        if len(data)>0:
            done_targets.append("%s_Marshal" % target.name)
            continue

        # Project ID 
        if target.priority == 5:
            projid = "'PTF1'" 
        elif target.priority == 4:
            projid = "'PTF2'"
        elif target.priority == 3:
            projid = "'PTF3'"
        elif target.priority == 2:
            projid = "'PTF4'"
        else:
            projid = "'PTF5'"

        # Otherwise, first add object
        objid = add_ptel_object(cursor, "%s_Marshal" % target.name, target.ra, 
                                target.dec, projid=projid, priority=priority, 
                                active="'yes'")

        # Add first observation of the object
        add_ptel_obs(cursor, "'"+objid.replace("'","")+".1'", objid, projid,
                     "%s_Marshal" % target.name, exptime=exptime, 
                     sequence="'yes'", seq_after=0, seq_after_con=48.0, 
                     active="'yes'")

        # Add another 19, with appropriate sequence constraints
        for i in range(1,20):
            add_ptel_obs(cursor, "'"+objid.replace("'","")+".%i'" % (i+1),
                         objid, projid, "%s_Marshal" % target.name, 
                         exptime=exptime, sequence="'yes'",
                         seq_after="'"+objid.replace("'","")+".%i'" % i,
                         seq_after_con=48.0, active="'yes'")

        # Update list of completed targets 
        done_targets.append("%s_Marshal" % target.name)

    # Close database
    db.close()
    
    return
    
#####################################################################

def add_ptel_obs(cursor, obsid, objid, projid, name, exptime=PTELPTFexptime,
                 sequence="'no'",seq_after=0,seq_after_con=0,active="'yes'"):

    """Insert a new observation in the PAIRITEL obs table.  No return."""

    cursor.execute("INSERT INTO obs(ObsID,ObjID,ProjID,IsActive,Name,Exptime_req,is_in_sequence,seq_after_id,seq_after_constraint) values(%s,%s,%s,%s,'%s',%s,%s,%s,%s)" % (obsid,objid,projid,active,name,exptime,sequence,seq_after,seq_after_con))

    return       

######################################################################

def add_ptel_object(cursor, name, ra, dec, projid=PTELPTFprojid, 
                    priority=PTELPTFpriority, active="'yes'"):

    """Insert a new source into the PAIRITEL OBJECT table.  Returns
    ID of the new object."""
   
    # Need to find highest object number in the given project
    cursor.execute("SELECT ObjID FROM object WHERE ProjID LIKE %s" % projid)
    data = cursor.fetchall()

    # If there's nothing there, start with 1
    if len(data) == 0:
        objid = "'"+projid.replace("'","")+'.1'+"'"

    # Otherwise need to find the highest number ID 
    else:
        ids = []
        for line in data:
            ids.append(int(line[0].split(".")[1]))
        objid = "'"+projid.replace("'","")+".%i'" % (max(ids)+1)

    # Insert the new object
    cursor.execute("INSERT INTO object(ObjID,ProjID,Name,ra,decl,user_priority,IsActive) values(%s,%s,'%s',%s,%s,%s,%s)" % (objid, projid, name, ra, dec, priority, active))    

    return objid

#######################################################################

def p60sked2html(file=None, url=MarshalURL, filters=P60Filters, 
                 ptffile="p60-ptf.lis", oldfile="p60-targets.lis",
                 newfile="p60-targets.new", mjd=None):

    """Update P60 target list.  Simulate nightly observations.  Output
    nicely formatted html for viewing."""

    configure_http()
    targets = fetch_marshal_targets(url=url, filters=filters)
    new_p60_targets(targets, outfile=ptffile)
    targlis=oscar_defs.oscar_target_list(catalog=oldfile)
    i=0
    while i < len(targlis):
        for key,item in P60Programs.iteritems():
            if targlis[i].keyword['P60PRID'] == item[0]:
                targlis.__delitem__(i)
                i-=1
        i+=1

    targlis.write_targets(newfile)
    os.system("less %s >> %s" % (ptffile, newfile))
    obstargets = p60_simulate(newfile, "p60-targets.sim", mjd=mjd)
    p60_write_html(targlis, targets, obstargets, file=file, mjd=mjd)
    return

#######################################################################

def p60_write_html(oldtargets, newtargets, obstargets, file=None, mjd=None):

    """Write (in html format) the simulated P60 target list for a 
    given night."""

    # If MJD not specified, round to the nearest one
    if mjd==None:
        mjd = round(utc().mjd)
    tsked = DateTimeFromMJD(mjd)

    # Dictionary for total times
    expdict={}
    
    # Start with the header
    if file==None:
        print """
        <HTML>
        <HEAD>
        <LINK type="text/css" rel="stylesheet" href="style.cgi">

        <TITLE>P60 Simulated Schedule</TITLE>

        </HEAD>
        <BODY>

        """
    else:
        ofile=open(file, 'w')
        ofile.write('<HTML>\n<HEAD>\n<LINK type="text/css" rel="stylesheet" href="style.cgi">\n\n<TITLE>P60 Simulated Schedule</TITLE>\n\n</HEAD>\n<BODY>\n\n')

    # Table Properties
    if file==None:
        print """
        <center>

        <h2>P60 Simulated Schedule for %s</h2>
        <table class="schedule">
        <tr>
          <th>Target Name</th>
          <th>P60 Program</th>
          <th>Filter(s)</th>
          <th>Exposure Time</th>
          <th>Priority</th>
        </tr>

        """ % tsked.strftime("%d %B %Y")
    else:
        ofile.write('<center>\n\n<h2>P60 Simulated Schedule for %s</h2>\n<table class="schedule">\n<tr>\n  <th>Target Name</th>\n  <th>P60 Program</th>\n  <th>Filter(s)</th>\n  <th>Exposure Time</th>\n  <th>Priority</th>\n</tr>\n\n' % tsked.strftime("%d %B %Y"))

    # Rows
    for obstarget in obstargets:

        program=None
        priority=0

        # Get program information
        for newtarget in newtargets:

            if obstarget[0]==newtarget.name:
                if newtarget.priority > priority:
                    program=newtarget.program
                    priority=newtarget.priority
                elif newtarget.priority == priority:
                    if P60Programs[newtarget.program][3] > P60Programs[program][3]:
                        program=newtarget.program
                        priority=newtarget.priority


        # If non-PTF target
        if program==None:

            if file==None:
                print """
                  <tr>
                    <td align=center>%s</td>
                    <td align=center>Caltech Internal</td>
                """ % obstarget[0]

            else:
                ofile.write('  <tr>\n    <td align=center>%s</td>\n    <td align=center>Caltech Internal</td>\n' % obstarget[0])

        # For PTF Targets
        else:

            if file==None:
                print """
                  <tr bgcolor=%s>
                    <td align=center><a href="http://navtara.caltech.edu/cgi-bin/ptf/view_source.cgi?name=%s">%s</a></td>
                    <td align=center>%s</td>
                """ % (P60Programs[program][4], obstarget[0][3:], obstarget[0], P60Programs[program][1].strip("'"))

            else:
                ofile.write('  <tr bgcolor=%s>\n    <td align=center><a href="http://navtara.caltech.edu/cgi-bin/ptf/view_source.cgi?name=%s">%s</a></td>\n    <td align=center>%s</td>\n' % (P60Programs[program][4], obstarget[0][3:], obstarget[0], P60Programs[program][1].strip("'")))

        # Observation details
        filters = obstarget[1][0][1]
        exptime = int(obstarget[1][0][0]) * float(obstarget[1][0][2])
        for j in range(1,len(obstarget[1])):
            if not (obstarget[1][j][1]==obstarget[1][j-1][1]):
                filters+=", %s" % obstarget[1][j][1]
                exptime+= int(obstarget[1][j][0]) * float(obstarget[1][j][2])

        if file==None:
            print """
                <td align=center>%s</td>
                <td align=center>%.1f</td>
                <td align=center>%i</td>
               </tr>

            """ % (filters, exptime, priority)

        else:
            ofile.write('    <td align=center>%s</td>\n    <td align=center>%.1f</td>\n    <td align=center>%i</td>\n  </tr>\n\n' % (filters, exptime, priority))

        if not expdict.has_key(program):
            expdict[program]=exptime
        else:
            expdict[program]+=exptime

    # Final allocations
    totexp=0.0
    for key, value in expdict.iteritems():
        totexp+=value

    allstr = "<p>\nEstimated nightly usage: "
    for key, value in expdict.iteritems():
        if key==None:
            program="Other"
        else:
            program=P60Programs[key][1]
        allstr += "%s: %.2f hr (%.1f%%); " % (program.strip("'"), (value / 3600.0), (100.0 * value / totexp))

    if file==None:
        print allstr
    else:
        ofile.write(allstr)

    # Footer
    if file==None:
        print """
        </table>
        </center>
        </body>
        </html>
        """

    else:
        ofile.write("</table>\n</center>\n</body>\n</html>\n")
        ofile.close()

    return
 
#######################################################################

def p60_simulate(targetfile, simfile, standards="p60-standard.lis", mjd=None):

    """Simulate nightly P60 observations.   Return an ordered list of 
    targets observed.  ***Note: Because simulation will overwrite target 
    list, we need to be careful with filenames."""

    # If MJD not specified, round to the nearest one
    if mjd==None:
        mjd = round(utc().mjd)

    # Rename target list
    shutil.copy(targetfile, "junk.lis")

    # Simulate P60 target list for that night
    cmd = "%s -d -t %i junk.lis %s > %s" % (P60_OSS, mjd, standards, simfile)
    scmd = os.popen(cmd, 'r', -1)
    scmd.readlines()

    # Delete temporary files
    os.remove("junk.lis")
    if os.path.exists("junk.lis.1"):
        os.remove("junk.lis.1")

    # Parse simulated target list
    obstargs = []
    tmp = open(simfile, 'r')
    lines = tmp.readlines()
    tmp.close()
    i = 0
    for i in range(len(lines)):
        line = lines[i]
        re1 = re.search("Asked to observe '(\S+)'", line)
        if re1:
            name = re1.group(1)
            line2 = lines[i+1]
            re2 = re.search("Exposure:\s+(\d)x\s+(\w+)\s+for\s+(\w+)", line2)
            if (len(obstargs)>0) and (name == obstargs[len(obstargs)-1][0]):
                obstargs[len(obstargs)-1][1].append([int(re2.group(1)),
                                                     re2.group(2),
                                                     float(re2.group(3))])
            else:
                obstargs.append([name, [[int(re2.group(1)), re2.group(2),
                                        float(re2.group(3))]]])

            i+=3
        else:
            i+=1

    return obstargs

#######################################################################

def create_iausn_table(filename=IAUSNtable):

    conn = sqlite3.connect(filename)
    c = conn.cursor()

    c.execute('''CREATE TABLE IAUSN (name text, ra real, dec real,
                                     host text, type text, date date,
                                     reference int)''')

    conn.commit()
    c.close()

    return

#######################################################################

def create_crts_table(filename=CRTStable):

    conn = sqlite3.connect(filename)
    c = conn.cursor()

    c.execute('''CREATE TABLE CRTS (id text, ra real, dec real, date date,
                                    mag real, otclass text, url text)''')

    conn.commit()
    c.close()

    return

#######################################################################

def update_crts_table(filename=CRTStable, url=CRTSURL):

    # Connect to database
    conn = sqlite3.connect(filename)
    c = conn.cursor()

    # Grab HTML from CRTS
    flob = urllib2.urlopen(url)
    s = flob.read()
    html = etree.HTML(s)

    # Grab all the lines in the table
    rows = html.findall('.//tr')

    # Parse all the entries
    for i in range(1, len(rows)):

        row = rows[i]
        name = row[0][0].text.strip()
        
        # Check to see if entry already exists
        c.execute("SELECT COUNT(*) FROM CRTS WHERE id = ?", (name,))
        row2 = c.fetchall()

        if not row2[0][0]:

            s1 = "(http://[a-zA-Z./]*\d{8}/\d+.html)"
            r1 = re.search(s1, row[0][0].attrib['onclick'])
            if not r1:
                url = "Unknown"
            else:
                url = r1.group(1)
            ra = row[1].text.strip()
            dec = row[2].text.strip()
            ddate = "%s-%s-%s" % (row[3].text.strip()[0:4],
                                  row[3].text.strip()[4:6],
                                  row[3].text.strip()[6:8])
            mag = row[4].text.strip()
            if (len(row[12])>0):
                otclass = row[12][0].text
            else:
                otclass = row[12].text.split()[0]
                
            print "%s\t%10.3f%10.3f%12s%6.2f%10s\t%s" % (name, float(ra),
                                                         float(dec), ddate,
                                                         float(mag), otclass,
                                                         url)

            # Insert into table
            c.execute("INSERT INTO CRTS VALUES (?, ?, ?, ?, ?, ?, ?)",
                      (name, ra, dec, ddate, mag, otclass, url)) 

    # Commit changes and close cursor
    conn.commit()
    c.close()
            
    return

#######################################################################

def update_iausn_table(filename=IAUSNtable, url=IAUSNURL):

    # Connect to database
    conn = sqlite3.connect(filename)
    c = conn.cursor()

    # Grab HTML from IAU SNe URL
    flob = urllib2.urlopen(url)
    s = flob.read()
    html = etree.HTML(s)

    # Run through all the lines in the HTML table
    pre = html.find('.//pre')
    i=0

    while (i<len(pre)):

        # SN Name
        name = pre[i].text

        # Check to see if already in database
        c.execute("SELECT COUNT(*) FROM IAUSN WHERE name = ?", (name,))
        row = c.fetchall()

        if not row[0][0]:

            # Parse first line to get host and discovery date
            s1 = "^\s+(\S+\s+\S*)\s+(\d{4}\s\d{2}\s\d{2})\s+(\d{1,2}\s+\d{1,2}\.\d\s+[+-]\d{1,2}\s+\d{1,2})\s+\S*\s*\S*\s*(\d{1,2}\.?\d?)"
            r1 = re.match(s1, pre[i].tail)
            if r1:
                host = r1.group(1).strip()
                ddate = r1.group(2).replace(" ", "-")
            else:
                host = "Unknown"
                ddate = "1900-01-01"

            # Parse second line for coordinates and reference
            i+=1
            s2 = "^\s+(\d{1,2}\s+\d{1,2}\s+\d{1,2}\.?\d?\d?)\s+([+-]\s*\d{1,2}\s+\d{1,2}\s+\d{1,2}\.?\d?)"
            r2 = re.match(s2, pre[i].tail)
            if r2:
                coo = astrocoords(r2.group(1).replace(" ", ":"),
                                  r2.group(2).replace(" ", ":"))
                ra = coo.radeg()
                dec = coo.dcdeg()
            else:
                ra = -99.0
                dec = -99.0
            ref = pre[i][0].tail.strip()
            
            # Parse third line for type
            i+=1
            s3 = "^\s+(\S+)\s+(\d{4}\w+)"
            r3 = re.match(s3, pre[i].tail)
            if r3:
                sntype = r3.group(1)
            else:
                sntype = "Unknown"

            # Print results
            print "%i\t%s\t%10.3f%10.3f%15s%5s%15s%5s" % (i, name, ra, dec, host, sntype, ddate, ref)
            
            # Insert into table
            c.execute("INSERT INTO IAUSN VALUES (?, ?, ?, ?, ?, ?, ?)",
                      (name, ra, dec, host, sntype, ddate, ref)) 

            # Increment i once more
            i+=1
            
        # If match, increment i
        else:
            i+=3

    # Commit changes and close cursor
    conn.commit()
    c.close()

    return
    
#######################################################################

def p60_realtime(ptffile='ptf/p60-ptf.lis', outfile="oscar_new_target.lis",
                 oldfile="targets/p60-targets.lis", url=MarshalURL):

    configure_http()
    targets=fetch_marshal_targets(filters=P60Filters)
    new_p60_targets(targets, outfile=ptffile)
    newts=oscar_defs.oscar_target_list(catalog=ptffile)
    oldts=oscar_defs.oscar_target_list(catalog=oldfile)
    otargs=[]
    for newt in newts:
        done=0
        for oldt in oldts:
            if (newt.name==oldt.name) and (newt.keyword["P60PRID"]==oldt.keyword["P60PRID"]):
                done=1
                continue
        if not done:
            otargs.append(newt)

    otargets=oscar_defs.oscar_target_list()
    otargets.targets=otargs
    otargets.write_targets('ptf/p60-rt.lis')
    if (os.path.getsize('ptf/p60-rt.lis')>0):
        os.system('less ptf/p60-rt.lis >> oscar_new_target.lis')
        os.system('less ptf/p60-rt.lis >> targets/p60-targets.lis')

    return

#######################################################################
# OLD ROUTINES - Not used since institution of marshal
#######################################################################

###################################################

def add_pteltarget(name,ra,dec,projid=PTELPTFprojid,priority=PTELPTFpriority,
                   objid=None,active="'yes'",obsid=None,exptime=PTELPTFexptime,
                   seebad=PTELseebad,tranbad=PTELtranbad,airbad=PTELairbad,
                   dith=PTELdith,mno=PTELmno,rasz=PTELrasz,decsz=PTELdecsz,
                   maxoff=PTELmaxoff):

    """Add object and observation entry into database for PTEL target"""

    # Connect to Mt Hopkins database
    try:
        db=MySQLdb.connect(db=PTELdb,user=PTELdbuser,host=PTELdbhost,
                           passwd=PTELdbpass)
    except:
        print "Error connecting to PTEL database at Mt. Hopkins"
        return
    cursor=db.cursor()

    # Check to see if object already exists
    cursor.execute("SELECT ObjID FROM object WHERE Name LIKE %s" % name) 
    dat=cursor.fetchall()
    if (len(dat)>1):
        objid=str((dat[0])[0]).strip()

    # Otherwise need to create new object
    else:
        if (objid==None):
            cursor.execute("SELECT ObjID FROM object WHERE ProjID LIKE %s" 
                           % projid)
            dat=cursor.fetchall()
            if len(dat) == 0:
                objid="'"+projid.replace("'","")+'1'+"'"
            else:
                ids=[]
                for i in range(len(dat)):
                    temp=(((dat[i])[0]).split("."))[1]
                    if int(temp)<1000:
                        ids.append(int(temp))
                objid="'"+projid.replace("'","")+"."+str(max(ids)+1)+"'"
        cursor.execute("INSERT INTO object(ObjID,ProjID,Name,ra,decl,user_priority,IsActive) values(%s,%s,%s,%s,%s,%s,%s)" % (objid,projid,name,ra,dec,priority,active)) 

    # Update obsid if necessary
    if (obsid==None):
        cursor.execute("SELECT ObsID FROM obs WHERE ObjID LIKE %s" % objid)
        dat=cursor.fetchall()
        if len(dat)==0:
            obsid=objid.rstrip("'")+".1'"
        else:
            ids=[]
            for i in range(len(dat)):
                temp=(((dat[i])[0]).split("."))[2]
                if int(temp)<1000:
                    ids.append(int(temp))
            obsid=objid.rstrip("'")+"."+str(max(ids)+1)+"'"

    # Insert obs row into table
    cursor.execute("INSERT INTO obs(ObsID,ObjID,ProjID,ExpTime_req,worst_seeing,worst_trans,worst_air,Name,dither_pattern,frame_overlap,ra_sz,dec_sz,step_size_arcsec) values(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)" % (obsid,objid,projid,exptime,seebad,tranbad,airbad,name,dith,mno,rasz,decsz,maxoff))

    # Close database
    db.close()

#######################################################################

def update_ptf_p60(outfile="p60-targets.lis", 
                   filters=['gpr','rpr','ipr','zpr']):

    """Query PTF photometry database and update P60 target list based 
    on results of monitoring campaigns."""

    # Begin with empty new target list
    p60targets=[]

    # Connect to database and establish cursor
    ptfdb=pgdb.connect(database=PTFdbname,user=PTFdbuser,password=PTFdbpass,
                       host=PTFdbip)
    ptfc=ptfdb.cursor()
 
    # Want list of all targets with P60 photometry
    query="SELECT DISTINCT ON (shortname) shortname FROM photometry WHERE telescope='P60'"
    ptfc.execute(query)
    rows=ptfc.fetchall()
    if (len(rows)==0):
        print "No PTF/P60 targets in the photometry database!"
        ptfdb.close()
        return

    # Loop over targets
    for row in rows:

        name=row[0]
        exposures=[]
        nodetect=True
        
        # Make sure target is 'Surely Transient'
        query="SELECT type2 FROM saved_cands WHERE shortname = '%s'" % name
        ptfc.execute(query)
        row=ptfc.fetchone()
        if (row==None) or (row[0]!="SurelyTransient"):
            continue

        # Loop over filters
        for filter in filters:

            # Below detection threshold?
            if ptf_p60_nodetect(name, filter):

                # Add reference if not obtained yet
                if not ptf_p60_hasref(name,filter):
                    exposures.append(oscar_defs.oscar_exposure(filter,180.0,5))

            # Faint?
            elif ptf_p60_faint(name, filter):
                exposures.append(oscar_defs.oscar_exposure(filter,180.0,1))
                nodetect=False    

            # Otherwise keep observing
            else:
                exposures.append(oscar_defs.oscar_exposure(filter,120.0,1))
                nodetect=False

        # If exposures requested, append target
        if not exposures==[]:
        
            # Need RA and Dec first
            query="SELECT ra,dec FROM saved_cands WHERE shortname='%s'" % name
            ptfc.execute(query)
            row=ptfc.fetchone()
 
            # Priority decreased for reference images
            if nodetect:
                priority=1.05*P60PTFpriority
                keyword=P60PTFrefkeyword
            else:
                priority=1.10*P60PTFpriority
                keyword=P60PTFfollowkeyword

            # Create and append target
            target=oscar_defs.oscar_target('PTF%s' % name, row[0], row[1],
                                           exposures, priority=priority,
                                           offset=P60PTFoffset,
                                           roi=P60PTFroi,
                                           keyword=keyword,
                                           property=P60PTFproperty)
            p60targets.append(target)

    # Close database
    ptfdb.close()

    # Write out target list
    if (len(p60targets)==0):
        print "No P60/PTF Targets require updating"
        return
    else:
        tlis=oscar_defs.oscar_target_list()
        tlis.targets=p60targets
        tlis.write_targets(outfile)

############################################################################

def ptf_p60_nodetect(name, filter, tlim=14.0, minobs=3, snrlim=3.0, 
                     maglim=23.0, magerrlim=0.50):

    """Determine if a PTF/P60 target has gone undetected for a specified
    period of time in a given filter.  Requires at least minobs non-detections
    in that time period (otherwise just unobserved)."""

    # Connect to database and establish cursor
    ptfdb=pgdb.connect(database=PTFdbname,user=PTFdbuser,password=PTFdbpass,
                       host=PTFdbip)
    ptfc=ptfdb.cursor()

    # Determine current mjd and cutoff window
    tnow=gmt()
    tcut=tnow.mjd-tlim

    # Grab all photometry within desired time window
    query="SELECT mag, magerr, snr, isupperlimit FROM photometry WHERE shortname='%s' and filter='%s' and mjd>%f" % (name, filter, tcut)
    ptfc.execute(query)
    rows=ptfc.fetchall()
    
    # Close database
    ptfdb.close()

    # Require at least minobs observations
    if (len(rows)<minobs):
        print "Not enough observations of PTF%s in filter %s" % (name, filter)
        return False

    # Loop over rows
    for row in rows:
        
        # If upper limit, continue
        if (row[3] == 'True'):
            continue

        # Otherwise, check to see if we consider this a detection
        if (row[0] < maglim) and (row[1] < magerrlim) and (row[2] > snrlim):
            return False

    # If we make it through all observations, then not detected
    return True

#############################################################################

def ptf_p60_hasref(name, filter, mintime=1, minulim=20.0):

    """Determine if a given PTF/P60 target has reference images in a given
    filter.  Need 5 observations within mintime days of each other,
    all with limiting magnitudes greater than minulim."""

    # Connect to database and establish cursor
    ptfdb=pgdb.connect(database=PTFdbname,user=PTFdbuser,password=PTFdbpass,
                       host=PTFdbip)
    ptfc=ptfdb.cursor()

    # Query to look for minobs observations within mintime, all with minulum
    query = "SELECT p1.shortname FROM photometry AS p1, photometry AS p2, photometry AS p3, photometry AS p4, photometry AS p5 WHERE p1.shortname = '%s' AND p1.shortname=p2.shortname AND p1.shortname=p3.shortname AND p1.shortname=p4.shortname AND p1.shortname=p5.shortname AND p1.filter='%s' AND p1.filter=p2.filter AND p1.filter=p3.filter AND p1.filter=p4.filter AND p1.filter=p5.filter AND ABS(p1.mjd-p2.mjd) < %f AND ABS(p1.mjd-p3.mjd) < %f AND ABS(p1.mjd-p4.mjd) < %f AND ABS(p1.mjd-p5.mjd) < %f AND p1.isupperlimit='True' AND p2.isupperlimit='True' AND p3.isupperlimit='True' AND p4.isupperlimit='True' AND p5.isupperlimit='True' AND p1.limmag > %f AND p2.limmag > %f AND p3.limmag > %f AND p4.limmag > %f and p5.limmag > %f" % (name, filter, mintime, mintime, mintime, mintime, minulim, minulim, minulim, minulim, minulim)
    ptfc.execute(query)
    rows=ptfc.fetchall()

    if (len(rows)==0):
        return False
    else:
        return True

#############################################################################

def ptf_p60_faint(name, filter, tmin=2.0, maglim=20.0):

    """Determine if a given PTF/P60 target is considered 'faint'.  Grab
    all observations within the last tmin days.  If magnitude is < maglim,
    or if no detection and limiting magnitude is < maglim, target is
    considered faint."""

    # Connect to database and establish cursor
    ptfdb=pgdb.connect(database=PTFdbname,user=PTFdbuser,password=PTFdbpass,
                       host=PTFdbip)
    ptfc=ptfdb.cursor()

    # Determine current mjd and cutoff window
    tnow=gmt()
    tcut=tnow.mjd-tmin

    # Query to get all photometry within last tmin days
    query = "SELECT mag, isupperlimit, limmag FROM photometry WHERE shortname='%s' AND filter='%s' AND mjd>%f" % (name, filter, tcut)
    ptfc.execute(query)
    rows=ptfc.fetchall()

    # Close database
    ptfdb.close()

    # Loop over rows
    for row in rows:

        # If non-detection and upper limit deep enough
        if (row[0]==-99.0) and (row[2]>maglim):
            return True

        # Otherwise if detection
        if (row[0]>maglim) and (row[2]>maglim):
            return True

    # If we make it all the way through, return False
    return False

##########################################################################     
     
def get_new_ptf_4p60(outfile):

     """Program to query PTF candidate database, and return objects
        lacking P60 multi-color observations."""

     # Connect to database and establish cursor
     ptfdb=pgdb.connect(database=PTFdbname,user=PTFdbuser,password=PTFdbpass,
                        host=PTFdbip)
     ptfc=ptfdb.cursor()

     # Query for objects that are surely transients but lack P60 obs
     p60query="SELECT DISTINCT ON (shortname) shortname, ra, dec, mag, rundate FROM (SELECT * FROM saved_cands AS s1 INNER JOIN (SELECT DISTINCT ON (s2.shortname) s2.shortname FROM saved_cands AS s2 EXCEPT (SELECT DISTINCT ON (p1.shortname) p1.shortname FROM photometry AS p1, photometry AS p2, photometry AS p3, photometry AS p4 WHERE p1.shortname = p2.shortname AND p1.shortname=p3.shortname AND p1.shortname = p4.shortname AND UPPER(p1.telescope) = 'P60' AND UPPER(p2.telescope) = 'P60' and UPPER(p3.telescope) = 'P60' and UPPER(p4.telescope) = 'P60' AND UPPER(p1.filter) LIKE 'GPR' AND UPPER(p2.filter) LIKE 'RPR' AND UPPER(p3.filter) LIKE 'IPR' AND UPPER(p4.filter) LIKE 'ZPR' AND ABS(p1.mjd - p2.mjd) < 1.0 AND ABS(p1.mjd - p3.mjd) < 1.0 AND ABS(p1.mjd - p4.mjd) < 1.0)) AS Temp USING (shortname) ORDER BY s1.rundate DESC) AS Temp2 WHERE shortname IS NOT NULL AND shortname <> 'None' AND shortname <> '' and type2='SurelyTransient'"
     ptfc.execute(p60query)
     rows=ptfc.fetchall()

     # If no rows, no objects require immediate P60 observations
     if (len(rows)==0):
          print "No new PTF objects require P60 observations!"
          ptfdb.close()
          return

     # Otherwise, loop over rows and create targets
     p60targets=[]
     tnow=datetime.date.today()

     for row in rows:

          # If last observation > 10 days, use long exposure
          if ((tnow-datetime.date(int(row[4][:4]),int(row[4][4:6]),
                                  int(row[4][6:8]))).days > 10):
               target=oscar_defs.oscar_target('PTF%s' % row[0], row[1], row[2],
                       [oscar_defs.oscar_exposure('gpr',180.0,1),
                        oscar_defs.oscar_exposure('rpr',180.0,1),
                        oscar_defs.oscar_exposure('ipr',180.0,1),
                        oscar_defs.oscar_exposure('zpr',180.0,1)],
                       priority=1.15*P60PTFpriority, offset=P60PTFoffset,
                       roi=P60PTFroi, keyword=P60PTFnewkeyword,
                       property=P60PTFproperty)

          # If last mag > 20, make exposures 120 s
          elif (row[3] > 20.0):
               target=oscar_defs.oscar_target('PTF%s' % row[0], row[1], row[2],
                       [oscar_defs.oscar_exposure('gpr',120.0,1),
                        oscar_defs.oscar_exposure('rpr',120.0,1),
                        oscar_defs.oscar_exposure('ipr',120.0,1),
                        oscar_defs.oscar_exposure('zpr',120.0,1)],
                       priority=1.15*P60PTFpriority, offset=P60PTFoffset,
                       roi=P60PTFroi, keyword=P60PTFnewkeyword,
                       property=P60PTFproperty)

          # If mag > 18, make exposures 90 s
          elif (row[3] > 18.0):
               target=oscar_defs.oscar_target('PTF%s' % row[0], row[1], row[2],
                       [oscar_defs.oscar_exposure('gpr',90.0,1),
                        oscar_defs.oscar_exposure('rpr',90.0,1),
                        oscar_defs.oscar_exposure('ipr',90.0,1),
                        oscar_defs.oscar_exposure('zpr',90.0,1)],
                       priority=1.15*P60PTFpriority, offset=P60PTFoffset,
                       roi=P60PTFroi, keyword=P60PTFnewkeyword,
                       property=P60PTFproperty)

          # Else make exposures 60 s
          else:
               target=oscar_defs.oscar_target('PTF%s' % row[0], row[1], row[2],
                       [oscar_defs.oscar_exposure('gpr',60.0,1),
                        oscar_defs.oscar_exposure('rpr',60.0,1),
                        oscar_defs.oscar_exposure('ipr',60.0,1),
                        oscar_defs.oscar_exposure('zpr',60.0,1)],
                       priority=1.15*P60PTFpriority, offset=P60PTFoffset,
                       roi=P60PTFroi, keyword=P60PTFnewkeyword,
                       property=P60PTFproperty)

          p60targets.append(target)

     # Now create and write the list
     tlis=oscar_defs.oscar_target_list()
     tlis.targets=p60targets
     tlis.write_targets(outfile)

     # Close database connection
     ptfdb.close()
     return

############################################

def get_new_ptf_4ptel():

     """Program to query PTF candidate database, and return objects
        lacking PAIRITEL multi-color observations."""

     # Connect to database and establish cursor
     ptfdb=pgdb.connect(database=PTFdbname,user=PTFdbuser,password=PTFdbpass)
     ptfc=ptfdb.cursor()

     # Query for objects that are surely transients but lack PAIRITEL obs
     ptelquery="SELECT DISTINCT ON (shortname) shortname, ra, dec, mag, rundate FROM (SELECT * FROM saved_cands AS s1 INNER JOIN (SELECT DISTINCT ON (s2.shortname) s2.shortname FROM saved_cands AS s2 EXCEPT (SELECT DISTINCT ON (p1.shortname) p1.shortname FROM photometry AS p1, photometry AS p2, photometry AS p3 WHERE p1.shortname = p2.shortname AND p1.shortname=p3.shortname AND UPPER(p1.telescope) = 'PAIRITEL' AND UPPER(p2.telescope) = 'PAIRITEL' and UPPER(p3.telescope) = 'PAIRITEL' AND UPPER(p1.filter) LIKE 'J%' AND UPPER(p2.filter) LIKE 'H%' AND UPPER(p3.filter) LIKE 'K%' AND ABS(p1.mjd - p2.mjd) < 1.0 AND ABS(p1.mjd - p3.mjd) < 1.0 AS Temp USING (shortname) ORDER BY s1.rundate DESC) AS Temp2 WHERE shortname IS NOT NULL AND shortname <> 'None' AND shortname <> ''"
     ptfc.execute(ptelquery)
     rows=ptfc.fetchall()

     # If no rows, no objects require immediate PAIRITEL observations
     if (len(rows)==0):
          print "No new PTF objects require P60 observations!"
          ptfdb.close()
          return
     
     tnow=datetime.date.today()
     
     for row in rows:

          # Only worthwhile to look at things brighter than 18 mag
          if (row[3] > 18):
               print "PTF%s is too faint for PAIRITEL observations" % row[1]
               continue

          # If very old, put in for long exposure and increased priority
          if ((tnow-datetime.date(int(row[4][:4]),int(row[4][4:6]),
                                  int(row[4][6:8]))).days > 10):
               add_pteltarget("PTF%s" % row[0], row[1], row[2],
                              projid=PTELPTFprojid, priority=2*PTELPTFpriority,
                              objid=None, active="'yes'", obsid=None,
                              exptime=2*PTELPTFexptime)

          # Otherwise put in for standard PTF request
          else:
               add_pteltarget("PTF%s" % row[0], row[1], row[2],
                              projid=PTELPTFprojid, priority=PTELPTFpriority,
                              objid=None, active="'yes'", obsid=None,
                              exptime=PTELPTFexptime)

     # Close db
     ptfdb.close()
     return

#########################################################################

def fetch_marshal_targets2(outfile='p60-ptf.lis', url=MarshalURL):

     """Script to query Follow-Up Marshal for new targets to observe.
        Return formatted as ElementTree."""

     # Grab HTML from Marshal
     flob = urllib2.urlopen(url)
     s = flob.read()
     html = etree.HTML(s)
     # Convert to useful format (etree Element)
     targets = etree.Element("Targets")
     trows = html.findall('.//tr')
     for row in trows:
          if (len(row)==16) and (len(row[0])>0):
               targets.append(row)
     # Need today's date
     tnow = now()

     # Convert to oscar_target format
     otargs=[]
     i = 0
     while i < len(targets):
          
          target = targets[i]
          tname = "PTF%s" % target[1].text
          tra = float(target[2].text)
          tdec = float(target[3].text)
          tfilt = P60filts[target[4].text]
          if target[5].text=="N/A":
               tmag=19.0
          elif re.search(">", target[5].text):
               tmag=99.0
          else:
               tmag = float(target[5].text)
          try:
               tlim = float(target[7].text)
          except ValueError:
               tlim = 20.0
          tclass = target[11].text
          tdisc = DateTimeFrom(target[12].text)
          if type(target[13].text)==NoneType:
               tprog = ['0']
          else:
               tprog = target[13].text.split(',')
               tprog.pop()
          if type(target[14].text)==NoneType:
               tpri = ['0']
          else:
               tpri = target[14].text.split(',')
               tpri.pop()
    
          # Exposure time will depend on magnitude
          if (tlim > 22.0):
               texp = [oscar_defs.oscar_exposure(tfilt, 180.0, 3)]
          elif tmag < 18:
               texp = [oscar_defs.oscar_exposure(tfilt, 60.0, 1)]
          elif tmag < 20:
               texp = [oscar_defs.oscar_exposure(tfilt, 120.0, 1)]
          else:
               texp = [oscar_defs.oscar_exposure(tfilt, 180.0, 1)]

          # Get other exposure requests for same target
          while (i+1<len(targets)) and (targets[i+1][1].text == tname[3:]):
               if targets[i+1][5].text=="N/A":
                    t2mag=19.0
               elif (re.search(">", targets[i+1][5].text)):
                    t2mag=99.0
               else:
                    t2mag = float(targets[i+1][5].text)
               t2filt = P60filts[targets[i+1][4].text]
               if (float(targets[i+1][7].text) > 22.0):
                    texp.append(oscar_defs.oscar_exposure(t2filt, 180.0, 3))
               elif t2mag < 18:
                    texp.append(oscar_defs.oscar_exposure(t2filt, 60.0, 1))
               elif t2mag < 20:
                    texp.append(oscar_defs.oscar_exposure(t2filt, 120.0, 1))
               else:
                    texp.append(oscar_defs.oscar_exposure(t2filt, 180.0, 1))
               i+=1

          # Set priority
          try:
              bestpri=0; bestprog=0
              for j in range(len(tprog)):
                  (prog,pri)=(float(tprog[j]), float(tpri[j]))
                  if P60Programs.has_key(prog):
                      tpriority = max(P60Priority[pri] * P60Programs[prog][3] \
                                  - P60Weights["Age"] * (tnow-tdisc).days - \
                                  P60Weights["Magnitude"] * tmag, 0)
                  else:
                      tpriority = 0
                  if tpriority > bestpri:
                      bestpri = tpriority
                      bestprog = prog
          except:
              bestpri = 0.0 
              bestprog = 0

          # Create oscar_target
          otarg = oscar_defs.oscar_target(tname, tra, tdec, texp, 
                                          priority=bestpri,
                                          offset=P60PTFoffset,
                                          roi=P60PTFroi, 
                                          keyword={"P60PRID": 
                                                   P60Programs[bestprog][0],
                                                   "P60PRNM":
                                                   P60Programs[bestprog][1],
                                                   "P60PRPI":
                                                   P60Programs[bestprog][2],
                                                   "OBJTYPE":
                                                   "Transient"},
                                          property=P60PTFproperty)
          
          # Add to target list
          otargs.append(otarg)

          # Increment i
          i+=1

     # Write out new target list (for now)
     olis=oscar_defs.oscar_target_list()
     olis.targets=otargs
     olis.write_targets(outfile)

##########################################################################
