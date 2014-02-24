#!/usr/bin/env python

from oscar_session import *

######################################################################

def main():

    # Parse Command line
    try:
        opts, args = getopt.getopt( sys.argv[1:], 
                     "hdvorb:f:s:t:z:",
                    ["help","dummy","verbose","observe","rawranks",
                     "boost=","fwait=","swait=","time=","zero="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    ####################

    # Defaults
    verbose=0
    dummy=0
    skipstart=0
    rawranks=0
    boostnm=[]
    fwait=2.0
    stdwait=6.0
    time=0.0
    zeronm=[]

    # Options parsing
    for opt, val in opts:
        # Usage info only
        if opt in ("-h", "--help"):
            # help request
            usage()
            sys.exit(1)
        elif opt in ("-d","--dummy"):
            # dummy operation / no socket connection
            dummy=1
        elif opt in ("-v","--verbose"):
            # verbose operation
            verbose += 1
        elif opt in ("-o","--observe"):
            # commence with observations right away
            skipstart=1
        elif opt in ("-r","--rawranks"):
            # prioritize targets by raw priorities only
            rawranks=1
        elif opt in ("-b","--boost"):
            # boost the priority of specified target
            boostnm.extend(val.upper().split(','))
        elif opt in ("-f","--fwait"):
            # set the focus-waiting time
            try:
                fwait=float(val)
            except:
                print "Error converting fwait setting '%s' to float" % val
        elif opt in ("-s","--swait"):
            # set the standard-field waiting time
            try:
                stdwait=float(val)
            except:
                print "Error converting swait setting '%s' to float" % val
        elif opt in ("-t","--time"):
            # set the "current time" as MJD
            try:
                time=1.0*float(val)-MJDEPHEP
            except:
                print "Error converting time setting '%s' to float" % val
        elif opt in ("-z","--zero"):
            # zero the priority of specified target
            zeronm.extend(val.upper().split(','))
        else:
            sys.stderr.write("Unmatched option %s\n" % opt)

    # Command-line arguments
    if len(args)<2:
        usage()
        sys.exit(1)
    [catfile,stdfile]=args[0:2]

    ####################

    # Define the oscar_session object
    I=oscar_session(dummy=dummy,verbose=verbose,time=time,
                    skipstart=skipstart,amdomes=1,
                    allnight=1,endnight=1,writetgts=1)

    # Waiting time between standard field observations, in hours
    I.stdwait=stdwait

    # But we need to focus pretty regularly
    I.fwait=fwait

    # Number of stellar images in each focus run
    I.saonseq=9

    # Read in target catalog
    I.read_targets(catfile)

    # Set "rawranks" property as requested
    I.targets.rawranks=rawranks

    # Zero target priority for "zero" targets
    if len(zeronm)>0:
        for tgt in I.targets:
            if tgt.name.upper() in zeronm:
                tgt.priority=0

    # Boost target priority for "boost" targets
    if len(boostnm)>0:

        # Find the "highest-priority" target
        hipri=0.0
        for tgt in I.targets:
            if tgt.priority>hipri:
                hipri=float(tgt.priority)

        # Set the "boost" targets to 10x hipri
        for tgt in I.targets:
            if tgt.name.upper() in boostnm:
                tgt.priority=10*hipri
    
    # Read in standards catalog
    I.read_standards(stdfile)

    # Say hello to OCS
    I.connect()

    # Update ephemerides
    I.sunstat()

    # Do a full night
    I.full_night()

    # Quit
    I.done()
    sys.exit(0)
    
######################################################################

def usage():

    (xdir,xname)=os.path.split(sys.argv[0])
    print "Usage:  %s [-h] [-d -t <MJDSTART>] [-m] [-o] [-r] [-b <boosttgt>] [-z <zerotgt>] <catalog> <stdcatalog>" % xname

######################################################################

if __name__=="__main__":

    main()

