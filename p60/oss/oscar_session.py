
from oscar_defs import *

######################################################################
# 
# To-Do and Thoughts
#
#  1 - Write interesting record files
# x2 - Various levels of verbosity
# 
# Simple stuff:
# 
# x1 - Request biases in the morning if they weren't done in the
#      afternoon
# x2 - Adjust priority of targets during the night
# 
# Note:  See oscar_defs.py for global variables
#        oscar_session.py,v.old archives older versions
# 
# $Id: oscar_session.py,v 1.1 2007/06/20 18:19:57 derekfox Exp $
# 
######################################################################

class oscar_session:

    def __init__(self,dummy=0,skipstart=0,
                 verbose=0,time=0.0,nbias=5,
                 ndome=9,fstep=0.05,nfocus=10,nfine=5,
                 saonseq=9,amdomes=1,allnight=1,endnight=1,
                 writetgts=1):

        # Initialization arguments
        self.nosock=dummy
        self.verbose=verbose

        # Some constants
        self.alphabet=list(string.uppercase)

        # Establish basic properties
        self.ocssock=0
        self.ocshost="198.202.125.198"
        self.ocsport=1326
        self.lastillegal=0
        self.lastbadfilt=0
        self.lastocsretval=''
        self.lasttgtname=''
        self.obstarg={}
        self.obsfilt={}

        self.catfile=''
        self.targets=[]
        self.standards=[]
        self.filters=[]
        self.badfilt=[]
        self.afternoonfilt=[]
        self.morningfilt=[]
        self.scistat={}
        self.scitime={}

        # Default defocus setting
        self.defocus=defocus_default

        # Define some reserved keywords
        self.binroikey='BINROI'
        self.mscsizekey='MSCSIZE'
        self.mscposnkey='MSCPOSN'
        self.reserved_keys=[self.binroikey,self.mscsizekey,
                            self.mscposnkey]

        # Full-chip, single-pixel
        self.fullbinroi=[[1,1],[1,1,2048,2048]]
        self.binroi=[self.fullbinroi]
        self.binroicode={binroi2str(self.fullbinroi):'A'}

        # Additional defaults
        for binroi in binroi_defaults:
            if binroi in self.binroi:
                # no double-booking
                continue
            self.binroi.append(binroi)
            self.binroicode[binroi2str(binroi)]= \
                 self.alphabet[len(self.binroicode)]
            if len(self.binroicode)>=25:
                print "Can't have more than 26 default %s settings" % \
                      self.binroikey
                break

        # Status
        self.didcals=0
        self.didstart=0
        self.ocsgoodnight=0

        # Queue-observation properties
        self.skipstart=skipstart
        self.amdomes=amdomes
        self.allnight=allnight
        self.endnight=endnight
        self.writetgts=writetgts

        # Waiting time between dome status checks
        self.twait=60

        # Sun limit (degrees, negative for below horizon)
        self.sunlim=-12.0

        # Afternoon calibrations
        self.nbias=nbias
        self.ndome=ndome

        # New targets input
        self.newtgtfile=NEWTGTFILE
        self.tmptgtfile=TMPTGTFILE

        # New priority input
        self.newprifile=NEWPRIFILE
        self.tmpprifile=TMPPRIFILE

        # Science status information
        self.scistatfile=SCISTATFILE

        # Preferred telescope stow position
        self.stow_az=0.0      # Azimuth
        self.stow_elev=60.0   # Elevation
        self.stow_da=220.0    # Dome Angle

        # Focus loop and default focus exposure
        self.fstep=fstep
        self.nfocus=nfocus
        self.nfine=nfine
        self.ffilt='R'
        self.ftime=10.0
        self.saoftime=4.0
        self.saonseq=saonseq

        # Standard field exposures
        self.standard_exp=[oscar_exposure('R',30.0,1)]

        # Preferred interval for focus / standard star observations
        self.fwait=2.0        # Hours
        self.lastfocmjd=0.0   # MJD
        self.stdwait=6.0      # Hours
        self.lastnstd=0       # Index of last standard field
        self.laststdmjd=0.0   # MJD

        # Whether to use SAOFOCUS routine to focus on each field
        self.usesao=1
        # Whether to use a focus loop to focus on each field
        self.usefloop=0

        # Initialize Observer object
        self.observer=palomar

        # Establish the time_now setting
        self.tnow=time_now(fake=self.nosock,time=time)
        self.observer.date=self.now()
        self.set_twilight()

        # Initialize random number sequence
        random.seed(self.now())

    ##############################

    def connect(self):

        if not self.nosock:
            if self.verbose:
                print "Waiting to avoid bind() errors..."
            time.sleep(30)
            
        if self.verbose:
            print "Establishing communication with ocssock"

        if not self.nosock:
            self.ocssock=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
            self.ocssock.connect((self.ocshost,self.ocsport))
            self.ocssock.send('HI')
            ocsretval=self.ocssock.recv(2048)
            if not re.search('READY',ocsretval):
                print "Got strange return value "+ocsretval
                sys.exit(1)

        if self.verbose:
            print "Successful communication with ocssock"

    ##############################

    def disconnect(self):

        if self.verbose:
            print "Attempting to disconnect from OCS"

        ok=self.ocssend("DISCONNECT\n",wait=1,output=0)
        if not ok:
            print "Disconnect failed"
        else:
            print "Disconnected from OCS"

        if self.verbose:
            print "Please wait 60 seconds before attempting to reconnect"

    ##############################

    def done(self):

        if self.verbose:
            print "Sending DONE to OCS"

        if not self.nosock:
            ok=self.ocssend("DONE\n",wait=1,output=0)
        else:
            ok=1

        if not ok:
            print "DONE failed"
        else:
            print "Okay, DONE"

        if self.verbose:
            print "Please wait 60 seconds before attempting to reconnect"

    ##############################

    def read_targets(self,file):

        self.catfile=file
        self.targets=oscar_target_list(file,tnow=self.tnow,
                           badfilt=self.badfilt,sunlim=self.sunlim,
                           observer=self.observer)
        self.update_filters()
        self.update_binroi()
        self.update_standard_exp()
        
        if self.verbose:
            print "Read in %d targets from %s" % \
                  (len(self.targets),file)

    ##############################

    def read_standards(self,file):
        
        self.standards=oscar_target_list(file,tnow=self.tnow,
                              badfilt=self.badfilt,sunlim=self.sunlim,
                              observer=self.observer)
        
        if self.verbose:
            print "Read in %d standard targets from %s" % \
                  (len(self.standards),file)

    ##############################

    def write_targets(self,file=""):

        if len(file)==0:
            file=self.catfile
            
        # Back up the old target list
        catbackup=backupfile(file)
        shutil.copy(file,catbackup)
        # Write the new target list
        self.targets.write_targets(file)

    ##############################

    def update_filters(self):

        self.filters=self.targets.list_of_filters()

    ##############################

    def update_binroi(self):

        for t in self.targets:
            tbin=copy.deepcopy(t.binning)
            troi=copy.deepcopy(t.roi)
            if [tbin,troi] not in self.binroi:
                self.binroi.append([tbin,troi])
                if len(self.binroicode)<26:
                    self.binroicode[binroi2str([tbin,troi])]= \
                        self.alphabet[len(self.binroicode)]
                else:
                    self.binroicode[binroi2str([tbin,troi])]= \
                        self.alphabet[int(len(self.binroicode)/26)-1]+ \
                        self.alphabet[len(self.binroicode) % 26]
        
    ##############################

    def update_standard_exp(self):

        stdexp=[]
        for filt in self.filters:
            if standard_exposures.has_key(filt):
                stime=standard_exposures[filt]
                stdexp.append(oscar_exposure(filt,stime,1))

        self.standard_exp=stdexp

    ##############################

    def sunalt(self):
        [sunaz,sunalt]=self.tnow.sunpos(self.observer)
        return sunalt

    def sunstat(self):

        sunalt=self.sunalt()

        print "Sun is currently at altitude %.1f degrees" % sunalt

        if sunalt>self.sunlim:
            print "Currently cannot observe"

    ##############################

    def utdate(self):
        return self.tnow.date()

    def now(self):
        return self.tnow.now()

    def mjdnow(self):
        return self.tnow.mjd()

    def lstnow(self):
        return self.tnow.lst(self.observer)

    def lstrange(self):
        return self.tnow.lstrange(self.observer,sunalt=0.0)

    def lstmid(self):
        return self.tnow.lstmid(self.observer,sunalt=0.0)

    def set_twilight(self):
        self.twilight=self.tnow.twilight(self.observer,sunalt=self.sunlim)
        return self.twilight

    ##############################
    ##############################

    def afternoon(self,filters=[],nbias=0,ndome=0,wait=0,output=0):

        if len(filters)>0:
            dfilt=filters
        else:
            dfilt=self.filters

        # We require a request of >=1 bias
        if nbias<1:
            nbias=self.nbias

        # We require a request of >=1 exposure in >=1 filter
        if ndome<1:
            ndome=self.ndome

        if not self.nosock:

            ocsmsg = "AFTERNOON\n"

            # Biases:  Specify one for each binning/ROI combo
            for [bin,roi] in self.binroi:
                ocsmsg+="BIAS %d %d %d %d %d %d %d\n" % \
                         (bin[0],bin[1],roi[0],roi[1],roi[2],roi[3],
                          nbias)
                binroicode=self.binroicode[binroi2str([bin,roi])]
                ocsmsg+="KEYWORD %s %s\n" % (self.binroikey,binroicode)

            # Domeflats:  Request by filters
            for filt in dfilt:
                if dome_exposures.has_key(filt):
                    exptime=dome_exposures[filt]
                    ocsmsg+="%s %d %.1f\n" % (filt,ndome,exptime)
                    ocsmsg+="KEYWORD %s A\n" % self.binroikey
                else:
                    print "Bad filter request for domeflats: '%s'" % filt

            # And send!
            self.ocssend(ocsmsg,wait=wait,output=output)

            # Keep track of filters observed
            self.afternoonfilt=dfilt

            if self.verbose and not wait:
                print "Spawned afternoon calibrations request"

        else:
            # Advance clock for calibrations
            self.didcals=1
            tcalibrations=ndome*(1+len(dfilt))*TREADOUT
            self.tnow += tcalibrations/SECSPERDAY

            if self.verbose:
                print "Advanced clock for calibrations"

    ##############################

    def morning(self,filters=[],nbias=0,ndome=0,exclude=[],
                wait=0,output=0):

        if len(filters)>0:
            dfilt=filters
        else:
            dfilt=self.filters

        # Excluding specified filters
        if len(exclude)>0:
            for xfilt in exclude:
                if xfilt in dfilt:
                    dfilt.remove(xfilt)

        # We require a request of >=1 exposure in >=1 filter
        if ndome<1:
            ndome=self.ndome

        if not self.nosock:

            ocsmsg = "MORNING\n"

            if nbias>0:
                # Biases:  Specify one for each binning/ROI combo
                for [bin,roi] in self.binroi:
                    ocsmsg+="BIAS %d %d %d %d %d %d %d\n" % \
                             (bin[0],bin[1],roi[0],roi[1],roi[2],roi[3],
                              nbias)
                    binroicode=self.binroicode[binroi2str([bin,roi])]
                    ocsmsg+="KEYWORD %s %s\n" % (self.binroikey,binroicode)

            # Domeflats:  Request filters
            for filt in dfilt:
                if dome_exposures.has_key(filt):
                    exptime=dome_exposures[filt]
                    ocsmsg+="%s %d %.1f\n" % (filt,ndome,exptime)
                    ocsmsg+="KEYWORD %s A\n" % self.binroikey
                else:
                    print "Bad filter request for domeflats: '%s'" % filt

            # And send!
            self.ocssend(ocsmsg,wait=wait,output=output)

            # Keep track of filters observed
            self.morningfilt=dfilt

            if self.verbose and not wait:
                print "Spawned morning calibrations request"

        else:
            # Advance clock for calibrations
            self.didcals=1
            tcalibrations=ndome*len(dfilt)*TREADOUT
            self.tnow += tcalibrations/SECSPERDAY

            if self.verbose:
                print "Advanced clock for calibrations"

    ##############################

    def wait_for_night(self,sunlim=None):

        # Parse argument
        if type(sunlim)==NoneType:
            sunlim=self.sunlim

        # Wait for nightfall
        while self.sunalt()>=sunlim:
            if self.nosock:
                self.tnow += 120.0/SECSPERDAY
            else:
                time.sleep(120)

    ##############################

    def wait_for_targets(self):

        # Wait until we get a nonzero target
        # (but not past sunrise!)
        ntarg=-1
        while self.sunalt()<self.sunlim:
            if self.nosock:
                self.tnow += 120.0/SECSPERDAY
            else:
                time.sleep(120)
            ntarg=self.targets.rank_targets()
            if ntarg>=0:
                break

        return ntarg

    ##############################

    def check_ready(self,output=0):

        # "check_ready" messages were saturating the output for
        # simulations, so we needed to tone them down
        if self.nosock:
            return 1

        return self.ocssend("CHECKREADY\n",wait=1,output=output)

    ##############################

    def move_and_focus(self,tgt,fmed=0,nfocus=-1,fstep=-1.0,
                       filter="",exptime=-10.0,wait=0,output=0):

        # Focus run inputs or defaults
        if fmed<=0:
            fcenter="CURRENT"
        else:
            fcenter="%.3f" % fmed

        if nfocus<0:
            nfocus=self.nfine

        if fstep<0:
            fstep=self.fstep

        if len(filter)==0:
            filter=self.ffilt

        if exptime<0:
            exptime=self.ftime

        # Basic error check
        if not hasattr(tgt,'coords'):
            print "move_and_focus must be passed the target itself"
            return

        [radeg,decdeg]=tgt.coodeg()
        eqnx=tgt.equinox

        # Construct the message for the OCS

        ocsmsg  = "DOFOCUSLOOP\n"
        ocsmsg += "%s\n" % tgt.name
        # Send "0.0 0.0" for offset
        ocsmsg += "%.4f %.4f %.3f 0.0 0.0\n" % \
                  (radeg/15.0,decdeg,eqnx)

        ocsmsg += "%s %s %d %.3f %.3f\n" % \
                  (filter,fcenter,nfocus,fstep,exptime)
        ocsmsg += "KEYWORD %s A\n" % self.binroikey
        
        ok=self.ocssend(ocsmsg,wait=wait,output=output)

        self.lastfocmjd=self.mjdnow()
        return ok

    ##############################

    def move_and_saofocus(self,tgt,fmed=0,nfocus=-1,fstep=-1.0,
                          filter="",exptime=-1.0,wait=0,output=0):

        # Focus run inputs or defaults
        if fmed<=0:
            fcenter="CURRENT"
        else:
            fcenter="%.3f" % fmed

        if nfocus<0:
            nfocus=self.saonseq

        if fstep<0:
            fstep=self.fstep

        if len(filter)==0:
            filter=self.ffilt

        if exptime<0:
            exptime=self.saoftime

        # Basic error check
        if not hasattr(tgt,'coords'):
            print "move_and_saofocus must be passed the target itself"
            return

        [radeg,decdeg]=tgt.coodeg()
        eqnx=tgt.equinox

        # Construct the message for the OCS

        ocsmsg  = "DOSAOFOCUS\n"
        ocsmsg += "%s\n" % tgt.name
        # Send "0.0 0.0" for offset
        ocsmsg += "%.4f %.4f %.3f 0.0 0.0\n" % \
                  (radeg/15.0,decdeg,eqnx)

        ocsmsg += "%s %s %d %.3f %.3f\n" % \
                  (filter,fcenter,nfocus,fstep,exptime)
        ocsmsg += "KEYWORD %s A\n" % self.binroikey
        
        ok=self.ocssend(ocsmsg,wait=wait,output=output)

        # Check for a redo request
        if re.search('REDO',self.lastocsretval):
            re1=re.search('REDO\s+([\d\.]+)',self.lastocsretval)
            newfoc=re1.group(1)

            ocsmsg  = "DOSAOFOCUS\n"
            ocsmsg += "%s\n" % tgt.name
            # Send "0.0 0.0" for offset
            ocsmsg += "%.4f %.4f %.3f 0.0 0.0\n" % \
                      (radeg/15.0,decdeg,eqnx)
            ocsmsg += "%s %s %d %.3f %.3f\n" % \
                      (filter,newfoc,nfocus,fstep,exptime)
            ocsmsg += "KEYWORD %s A\n" % self.binroikey

            ok=self.ocssend(ocsmsg,wait=wait,output=output)

        self.lastfocmjd=self.mjdnow()
        return ok

    ##############################

    def start_of_night(self,nstd=0,ffilt='',ftime=10.0):

        ok=1

        if nstd==0:
            nstd=self.choose_standard()

        self.lastnstd=nstd
        topstd=self.standards[nstd]
        self.lastfocmjd=self.mjdnow()
        self.laststdmjd=self.mjdnow()

        if len(ffilt)>0:
            self.ffilt=ffilt
            self.ftime=ftime

        if not self.nosock:

            # OCS Startnight Command
            ocsmsg = "STARTNIGHT\n"
            # Target Name
            ocsmsg += topstd.name+"\n"
            # RA DEC
            radec=topstd.coords.deg()
            eqnx=topstd.coords.equinox()
            # Currently NOT sending "0 0" for Proper Motions
            ocsmsg += "%.6f %.6f %.2f\n" % \
                      (radec[0]/15,radec[1],eqnx)
            # Filter
            ocsmsg += self.ffilt+"\n"
            # Exptime
            ocsmsg += "%.2f\n" % self.ftime
            # Number of focus steps
            ocsmsg += "%d\n" % self.nfocus
            # Size of focus steps
            ocsmsg += "%.3f\n" % self.fstep
            # Subsequent standard-field exposures
            for stdexp in self.standard_exp:
                ocsmsg += "%s %.2f\n" % (stdexp.filter,stdexp.exptime)
            # And send!
            ok=self.ocssend(ocsmsg,wait=0,output=0)

            if self.verbose:
                print "Spawned start of night request"

        else:

            # Advance clock for start of night
            self.didstart=1
            tstartofnight=self.nfocus*TREADOUT
            self.tnow += tstartofnight/SECSPERDAY

            if self.verbose:
                print "Advanced clock for start of night request"

        return ok

    ##############################

    def choose_target(self,boostactive=3.0):

        [iexp,jra,kdec]=[0,0,0]

        # Factor by which to multiply scores of partially-complete targets
        if boostactive<1:
            boostactive=1.0

        # Loop until we get a useful target
        ntarg=-1
        while ntarg<0:

            # Get the highest-ranked target
            ntarg=self.targets.rank_targets(boostactive=boostactive)

            # Check for out-of-targets condition
            if ntarg<0:
                sys.stderr.write("Out of targets to observe\n")
                break

            # Best target
            toptarg=self.targets[ntarg]

            # Choose specific exposure
            [iexp,jra,kdec]=toptarg.choose_exposure(badfilt=self.badfilt)

            # Move to next target if we get no exposure
            if iexp<0:
                toptarg.check_done()
                ntarg=-1

        return [ntarg,iexp,jra,kdec]

    ##############################

    def choose_standard(self):

        nstd=self.standards.rank_targets()

        if nstd<0:
            # no good standards?
            sys.stderr.write("Couldn't find standard target to observe\n")

        return nstd

    ##############################

    def observe_target(self,ntargin,iexp=-1,jra=-1,kdec=-1,wait=1):

        if len(ntargin)==4:
            [ntarg,iexp,jra,kdec]=ntargin
        else:
            ntarg=ntargin

        toptarg=self.targets[ntarg]

        if iexp<0 or jra<0 or kdec<0:
            [iexp,jra,kdec]=toptarg.choose_exposure(badfilt=self.badfilt)

        self.submit_target(toptarg,iexp,jra,kdec,wait=0,output=0)
        topexp=toptarg.exposures[iexp]

        if self.verbose:
            print "Asked to observe '%s'" % toptarg.name,
            if toptarg.nr>1 or toptarg.nd>1:
                print " at (%d,%d)" % (jra,kdec),
            print ""
            print "Exposure:  "+str(topexp),
            if topexp.imult>1:
                print "/ #%d" % (1+topexp.obsi[jra,kdec]),
            print ""

        if wait:

            # Wait for OCS to complete last command
            retval=self.wait()

            if retval:

                # Successful observation
                self.targets[ntarg].update(iexp,jra,kdec)

                # Record-keeping
                self.obstarg[ntarg]=1
                self.obsfilt[topexp.filter]=1

                if self.verbose:
                    print "Successful observation of %s" % toptarg.name

            elif self.lastillegal:

                # Update target status to prevent rerequest
                self.targets[ntarg].property['ILLEGAL']=1
                # Clear status flag
                self.lastillegal=0

                if self.verbose:
                    print "Target %s request illegal" % toptarg.name

            elif self.lastbadfilt:

                # Add to list of bad filters
                badfilt=topexp.filter
                self.badfilt.append(badfilt)
                # Remove from tonight's filter list,
                # thus preventing domeflat requests
                if badfilt in self.filters:
                    self.filters.remove(badfilt)
                    self.update_standard_exp()
                # Clear status flag
                self.lastbadfilt=0

                if self.verbose:
                    print "Filter '%s' not currently in filter wheel" % \
                          badfilt

        else:

            if self.verbose:
                print "Spawned observation request for %s" % toptarg.name

        if self.nosock:

            # Advance clock for the observation
            ttarget=toptarg.exposures[iexp].exptime+TREADOUT
            self.tnow += ttarget/SECSPERDAY

        # Update our records
        self.lasttgtname=copy.deepcopy(toptarg.name)

        return

    ##############################

    def submit_target(self,targin,iexp,jra,kdec,wait=0,output=0):

        # Control behavior: Whether to apply the offset or not
        applyoffset=0

        # Exposure in question
        xpr=targin.exposures[iexp]

        # Switch targets if exposure belongs to a subtarget
        if type(xpr.target)!=StringType:
            # Switch to subtarget
            target=xpr.target
        else:
            # Retain target that was given to us
            target=targin

        # Mosaic properties for this exposure
        obsi=xpr.obsi
        #[jmax,kmax]=obsi.getshape()
        [jmax,kmax]=obsi.shape

        # Begin constructing message for OCS
        ocsmsg = "SCIENCE\n"
        # Target Name
        ocsmsg += target.name+"\n"

        ###################
        ## Exposure Info ##
        ###################

        [rasxg,dcsxg]=target.coosxg()
        coords=astrocoords(rasxg,dcsxg,equinox=target.equinox,
                           system=target.system)

        # Adjust target coords for mosaics
        jx=jra-(jmax-1)/2
        ky=kdec-(kmax-1)/2
        coords.shift([jx*0.5*FOV,ky*0.5*FOV],arcsec=1)

        # Get a copy of the default offset
        offset=copy.deepcopy(target.offset)

        # Dither corrections
        dither=[0.0,0.0]
        if len(xpr.edither)>0:
            # Adjust target coords for exposure-specific dither
            npos=len(xpr.edither) # number of dither positions
            ipos=obsi[jra,kdec] % npos
            dithers=xpr.edither[ipos].split(',')
            dither=[float(dithers[0]),float(dithers[1])]
        elif len(target.dither)==2:
            # Random dither in a box (default behavior)
            dither[0]=(random.random()-0.5)*target.dither[0]
            dither[1]=(random.random()-0.5)*target.dither[1]
        # Adjust offset for the dither
        offset[0] += dither[0]
        offset[1] += dither[1]

        # Apply the offset to the coordinates
        if applyoffset:
            coords.shift(offset,arcsec=1)

        # Extract final coordinates for exposure
        radec=coords.deg()
        eqnx=coords.equinox
        # Syntax:  <RA hours> <Dec deg> <Eqnx> <RA offset> <Dec offset>
        #   (No sending of proper motions)
        if applyoffset:
            ocsmsg += "%.6f %.6f %.2f 0.0 0.0\n" % \
                          (radec[0]/15,radec[1],eqnx)
        else:
            ocsmsg += "%.6f %.6f %.2f %.2f %.2f\n" % \
                          (radec[0]/15,radec[1],eqnx,offset[0],offset[1])

        # Filter
        ocsmsg += xpr.filter+"\n"
        # Exptime
        ocsmsg += "%.2f\n" % xpr.exptime

        # Defocus?
        if target.property.has_key('DEFOCUS'):
            if target.property['DEFOCUS']:
                pass

        #####################
        ## Header Keywords ##
        #####################

        # Binning
        tbin,troi=target.binning,target.roi
        ocsmsg += "%d %d\n" % (tbin[0],tbin[1])

        # ROI
        ocsmsg += "%d %d %d %d\n" % (troi[0],troi[1],troi[2],troi[3])

        # BINROI keyword
        binroicode=self.binroicode[binroi2str([tbin,troi])]
        ocsmsg += "KEYWORD %s %s\n" % (self.binroikey,binroicode)

        # Mosaic setting
        ocsmsg += "KEYWORD %s %d,%d\n" % (self.mscsizekey,jmax,kmax)
        ocsmsg += "KEYWORD %s %d,%d\n" % (self.mscposnkey,jra,kdec)

        # Generic header keywords facility
        for keywd in target.keyword.keys():
            if keywd in self.reserved_keys:
                continue
            ocsmsg += "KEYWORD %s %s\n" % (keywd,target.keyword[keywd])

        # And send!
        self.ocssend(ocsmsg,wait=wait,output=output)

    ##############################

    def observe_standard(self,nstd=0,fmed=0,nfocus=-1,fstep=-1.0,
                         ffilt="",ftime=-10.0,retry=0,wait=1,output=0):

        # Return value
        ok=1

        # Parse arguments/keywords
        if nstd==0:
            nstd=self.choose_standard()

        topstd=self.standards[nstd]

        if topstd.score()<=0:
            sys.stderr.write("Zero-score standard field requested\n")
            return

        if nfocus<0:
            nfocus=self.nfine

        if fstep<0:
            fstep=self.fstep

        if len(ffilt)==0:
            ffilt=self.ffilt

        if ftime<0:
            ftime=self.ftime

        # Submit for observation
        if not self.nosock:

            if self.verbose:
                print "Taking standard filter images."

            # OCS "Standard" Command
            ocsmsg = "STANDARD\n"
            # Target Name
            ocsmsg += topstd.name+"\n"
            # RA DEC
            radec=topstd.coords.deg()
            eqnx=topstd.coords.equinox
            # Currently sending "0.0 0.0" for offset
            ocsmsg += "%.6f %.6f %.2f 0.0 0.0\n" % \
                      (radec[0]/15,radec[1],eqnx)
            # Subsequent standard-field exposures
            for stdexp in self.standard_exp:
                ocsmsg += "%s %.2f\n" % (stdexp.filter,stdexp.exptime)
            # BINROI
            ocsmsg += "KEYWORD %s A\n" % self.binroikey
            
            # And send!
            ok=self.ocssend(ocsmsg,wait=wait,output=output)

            if self.verbose:
                if wait:
                    print "Done"
                else:
                    print "Spawned STANDARD command"
            
        else:

            # Advance clock for the observations
            tstd=(ftime+TREADOUT)*nfocus

            for exp in self.standard_exp:
                tstd += exp.exptime+TREADOUT

            self.tnow += tstd/SECSPERDAY

            if self.verbose:
                print "Advanced clock for focus loop + " + \
                      "standards observation of %s" % topstd.name

        # Update our records
        self.lastnstd=nstd
        self.laststdmjd=self.mjdnow()
        self.lasttgtname=copy.deepcopy(topstd.name)

        return ok

    ##############################

    def check_new_targets(self):
        
        # Check for new targets
        if not os.path.exists(self.newtgtfile):
            return

        os.rename(self.newtgtfile,self.tmptgtfile)
        if self.verbose:
            print "Reading new targets file %s" % self.newtgtfile

        newtgt=oscar_target_list(self.tmptgtfile,tnow=self.tnow,
                     badfilt=self.badfilt,sunlim=self.sunlim,
                     observer=self.observer)
        if self.verbose:
            print "Found %d new target(s)" % len(newtgt)
        for newt in newtgt:
            self.targets.append(newt)
        os.remove(self.tmptgtfile)

        # Update filter list, etc.
        self.update_filters()
        self.update_binroi()
        self.update_standard_exp()

    ##############################

    def check_new_priority(self):

        # Check for new targets
        if not os.path.exists(self.newprifile):
            return
        
        os.rename(self.newprifile,self.tmpprifile)
        if self.verbose:
            print "Reading new priorities from %s" % self.newprifile

        prilines=getlines(self.tmpprifile)

        for line in prilines:
            # Comment line
            if re.search('^#',line):
                continue
            # Blank line
            if not re.search('\S+',line):
                continue
            els=line.split()
            # Only one word in line?
            if len(els)<2:
                continue
            renm=els[0].upper()
            try:
                newpri=float(els[1])
            except:
                print "Failed to convert priority '%s' to float" % els[1]
                continue
            # Reset priorities for matching targets
            for tgt in self.targets:
                tgtname=tgt.name.upper()
                if tgtname.find(renm)>=0:
                    tgt.priority=newpri
                    if self.verbose:
                        print "Reset target '%s' to priority %.1f" % \
                              (tgt.name,newpri)

        # Done with changes; delete file
        os.remove(self.tmpprifile)

        # Re-sort target list
        self.targets.sort_targets()

    ##############################

    def check_status(self):
        
        # Update "Science Status" information
        if not os.path.exists(self.scistatfile):
            return

        # Status keys to use for seeing, etc.
        seekey='SEEING_R'
        extkey='EXTINCT_R'
        skykey='SKYBKG_R'

        statlines=getlines(self.scistatfile)
        for line in statlines:
            reln=re.search("^(\d+)\s+(\S+)\s+(.+)",line)
            time,keywd,keyval=reln.groups()
            self.scistat[keywd]=keyval
            self.scitime[keywd]=time

        if self.scistat.has_key(seekey):
            oscar_values['seeing']=self.scistat[seekey]
        if self.scistat.has_key(extkey):
            oscar_values['extinction']=self.scistat[extkey]
        if self.scistat.has_key(skykey):
            oscar_values['sky']=self.scistat[skykey]

    ##############################

    def queue_observe(self,amdomes=None,allnight=None,endnight=None,
                      writetgts=None):

        # Parse inputs
        if type(amdomes)!=NoneType:
            self.amdomes=amdomes
        if type(allnight)!=NoneType:
            self.allnight=allnight
        if type(endnight)!=NoneType:
            self.endnight=endnight
        if type(writetgts)!=NoneType:
            self.writetgts=writetgts

        # Error checking - Sun location
        if self.sunalt()>=self.sunlim:
            if self.allnight:
                self.wait_for_night()
            else:
                print "Sun is too high to observe, suggest wait_for_night"
                return

        # Error checking - Target list
        [ntarg,iexp,jra,kdec]=self.choose_target()
        if ntarg<0:
            print "Can't find a good target to observe right now"
            if self.allnight:
                print "Waiting for targets to come along..."
                self.wait_for_targets()
            else:
                print "Please augment target list and retry"
                return

        try:

            while self.sunalt()<self.sunlim:

                # Check facility status
                if not self.check_ready():
                    if self.verbose:
                        print self.lastocsretval
                    # We continue after waiting for these reasons:
                    #   (1) Failed check_ready takes 10 minutes, which
                    #       is pretty long already
                    #   (2) Give the routine a chance to quit at end
                    #       of night
                    time.sleep(self.twait)

                # Time for a standard field?
                if self.stdwait>0 and \
                       24*(self.mjdnow()-self.laststdmjd)>self.stdwait:

                    if self.laststdmjd>0 and \
                           self.standards[self.lastnstd].score()>0:
                        # Use the same field as before
                        nstd=self.lastnstd
                    else:
                        # Time for a new standard field
                        nstd=self.choose_standard()

                    # Observe a standard field
                    self.observe_standard(nstd=nstd,wait=1)

                # Check priorities, targets, status
                self.check_new_priority()
                self.check_new_targets()
                self.check_status()

                # Choose an exposure for observation
                [ntarg,iexp,jra,kdec]=self.choose_target()

                # No targets left!
                if ntarg<0:

                    # If less than 1 hour remaining, go to standard field:
                    if 24*(self.twilight[1]-self.now())<1:
                        print "Ending night on standard field"
                        if self.laststdmjd>0 and \
                               self.standards[self.lastnstd].score()>0:
                            # Use the same field as before
                            nstd=self.lastnstd
                        else:
                            # Time for a new standard field
                            nstd=self.choose_standard()

                        # Observe a standard field
                        self.observe_standard(nstd=nstd,wait=1)

                        # Lather, rinse, repeat
                        continue

                    # Lots of time left -- close dome (temporarily) or quit
                    if self.allnight:
                        self.close_dome(wait=1,output=0)
                        ntarg=self.wait_for_targets()
                        if ntarg>=0:
                            continue
                        else:
                            break
                    else:
                        break

                # Pull out some useful items
                toptgt=self.targets[ntarg]
                topexp=toptgt.exposures[iexp]

                # Time to re-focus?
                focusnow=self.fwait>0 and \
                       24*(self.mjdnow()-self.lastfocmjd)>self.fwait
                   
                # Is it a new field we need to focus on?
                if toptgt.name != self.lasttgtname:
                    focusnow=1
                    if toptgt.property.has_key('NOFOCUS'):
                        if toptgt.property['NOFOCUS']:
                            focusnow=0

                # Focus on field if it's a new field
                if focusnow:
                    # Requested filter
                    ffilt=topexp.filter
                    # Check that this filter is reasonable for focusing
                    if not standard_exposures.has_key(ffilt):
                        # No good - use standard focus filter
                        ffilt=self.ffilt
                    # Move to the target & focus
                    if self.usesao:
                        self.move_and_saofocus(toptgt,filter=ffilt,wait=1)
                    if self.usefloop:
                        self.move_and_focus(toptgt,nfocus=self.nfine,
                                            filter=ffilt,wait=1)

                # Observe the top target
                self.observe_target([ntarg,iexp,jra,kdec],
                                    wait=1)

        finally:

            # Shutdown the facility if appropriate
            if self.endnight or self.sunalt()>self.sunlim:
                self.sunstat()
                # ENDNIGHT closes dome & stows telescope
                self.end_night()

        if not self.endnight and ntarg<0:
            print "Out of targets! (not closing)"
            self.sunstat()
            return

        if self.writetgts:
            # Write the new target list
            self.write_targets()

        # Morning domeflats if requested
        if self.amdomes:
            # Request bias images, too, if we didn't get them already
            if len(self.afternoonfilt)>0:
                nbias=0
            else:
                nbias=self.nbias
            # We do most of the filters now
            self.morning(nbias=nbias,exclude=self.afternoonfilt,
                         wait=1)

        # The End
        
    ##############################

    def full_night(self,skipstart=None,amdomes=None,
                   allnight=None,endnight=None,writetgts=None):

        # Parse inputs
        if type(skipstart)!=NoneType:
            self.skipstart=skipstart
        if type(amdomes)!=NoneType:
            self.amdomes=amdomes
        if type(allnight)!=NoneType:
            self.allnight=allnight
        if type(endnight)!=NoneType:
            self.endnight=endnight
        if type(writetgts)!=NoneType:
            self.writetgts=writetgts

        # Afternoon calibrations only if the sun is up
        if self.sunalt()>0:
            if self.amdomes:
                # Minimal set of domeflats in afternoon
                self.afternoon(filters=['R'],wait=1)
            else:
                # Full set of calibrations in afternoon
                self.afternoon(wait=1)

        # Wait for sunset
        if self.sunalt()>-5:
            # Open the dome
            self.open_dome(wait=1)
            # Wait for sunset
            self.wait_for_night(sunlim=-5)

        # Open the dome (another try -- no harm if already open)
        self.open_dome(wait=1)

        # Wait for astronomical twilight
        if self.sunalt()>=self.sunlim:
            self.wait_for_night()

        # Wait for facility-ready status
        self.sunstat()
        if not self.check_ready():
            if self.verbose:
                print self.lastocsretval
            time.sleep(self.twait)

        # Make sure we are focusing in some fashion
        if not self.usesao and not self.usefloop:
            self.usesao=1

        # Start of night standard observation
        nstd=self.choose_standard()
        topstd=self.standards[nstd]
        if self.skipstart:
            # Pretend like we are observing a standard field
            self.lastnstd=nstd
            self.laststdmjd=self.mjdnow()
            self.lasttgtname=copy.deepcopy(topstd.name)
        else:
            # Focus first
            if self.usesao:
                self.move_and_saofocus(topstd,nfocus=self.saonseq+2,wait=1)
            else:
                self.move_and_focus(topstd,nfocus=self.nfocus,wait=1)
            # Observe the standard field
            self.observe_standard(nstd=nstd,wait=1)

        # Observe
        self.queue_observe()

    ##############################

    def focus_test(self,filters=[],endnight=0,nloop=-1):

        iloop=0

        if len(filters)==0:
            filters=self.ffilt

        try:

            while self.sunalt()<self.sunlim:

                # Check facility status
                if not self.check_ready():
                    if self.verbose:
                        print self.lastocsretval
                    time.sleep(self.twait)
                    continue

                # Select all Score>0 Standards
                goodstd=[]
                for std in self.standards:
                    if std.score()>0:
                        goodstd.append(std)
                        
                if len(goodstd)<1:
                    print "No standard fields visible?"
                    break

                # Pick one at random
                topstd=random.choice(goodstd)
                nstd=self.standards.index(topstd)
                self.laststdmjd=self.mjdnow()

                # Observe it
                if not self.nosock:

                    if self.verbose:
                        print "Moving telescope to standard field"

                    self.move_std(nstd,wait=1,output=0)

                    for filt in filters:

                        if self.verbose:
                            print "Executing focus loop in %s" % filt

                        self.focus_loop(0,self.nfine,self.fstep,filt,
                                        self.ftime,wait=1,output=0) 

                    if self.verbose:
                        print "Done"

                else:

                    # Advance clock for the observations
                    tstd=(self.ftime+TREADOUT)*self.nfine*len(filters)

                    self.tnow += tstd/SECSPERDAY

                    if self.verbose:
                        print "Advanced clock for focus loops " + \
                              "on target %s" % topstd.name

                # Check if we are done
                iloop+=1
                if nloop>0 and iloop>=nloop:
                    break

        finally:

            # Shutdown the facility if appropriate
            if endnight or self.sunalt()>self.sunlim:
                self.sunstat()
                self.end_night()

        if not endnight and len(goodstd)<1:
            print "No good standards available -- not closing"
            self.sunstat()
            return

    ##############################

    def wait(self,output=0):

        retval=0
        self.lastillegal=0
        self.lastocsretval=''

        # No messing around if we are doing simulations
        if self.nosock:
            return 1

        if self.verbose:
            print "Waiting for return message on ocssock..."
        
        ocsretval=self.ocssock.recv(2048)
        self.lastocsretval=ocsretval
        
        if re.search('SUCCESS',ocsretval):
            retval=1
            if self.verbose:
                print "Request completed successfully"
        elif re.search('FAILURE',ocsretval):
            if re.search('NOFILT',ocsretval):
                self.lastbadfilt=1
                if self.verbose:
                    print "Request failed (bad filter)"
            else:
                if self.verbose:
                    print "Request failed"
        elif re.search('ILLEGAL',ocsretval):
            # Set "illegal" flag
            self.lastillegal=1
            if self.verbose:
                print "Request failed (illegal target)"
        elif re.search('REDO',ocsretval):
            # Redo a focus run
            if self.verbose:
                print "OCS requested Redo"
        elif re.search('DISCONNECT',ocsretval):
            # Disconnect request confirmation
            retval=1
            if self.verbose:
                print "OCS is disconnecting"
        elif re.search('GOODNIGHT',ocsretval):
            # OCS is going bye-bye
            self.ocsgoodnight=1
            if self.verbose:
                print "OCS is going away, time to quit"
            if self.writetgts:
                self.write_targets()
            # Just quit, I guess...
            sys.exit(0)
        else:
            print "Strange return from oscsock:"
            if not output:
                print ocsretval

        if output:
            print ocsretval

        return retval

    ##############################

    def end_night(self,wait=1):

        if self.verbose:
            print "Sending end of night command"

        if not self.nosock:

            ocsmsg="ENDNIGHT\n"
            retval=self.ocssend(ocsmsg,wait=1,output=0)

            if retval:
                print "Successful end of night"
            else:
                print "Failed end of night"

    ##############################
    #   Individual OCS routines  #
    ##############################

    def ocssend(self,message,wait=0,output=0):

        ok=1
        if not self.nosock:
            self.ocssock.send(message)
        elif self.verbose>1:
            print "Sending this message to OCS:"
            print message,

        if wait:
            ok=self.wait(output=output)

        return ok

    ##############################

    def stow_tel(self,azimuth=-1,elevation=-1,domeangle=-1,
                 wait=0,output=0):

        if azimuth<0 or elevation<0 or domeangle<0:
            [azimuth,elevation,domeangle] = \
                 [self.stow_az,self.stow_elev,self.stow_da]

        ocsmsg="STOWTEL %.3f %.3f %.3f\n" % \
                (azimuth,elevation,domeangle)

        self.ocssend(ocsmsg,wait=wait,output=output)

    ##############################

    def stow_tel_dflat(self,wait=0,output=0):

        [azimuth,elevation,domeangle] = \
             [0., 90., 90.]

        ocsmsg="STOWTEL %.3f %.3f %.3f\n" % \
                (azimuth,elevation,domeangle)

        self.ocssend(ocsmsg,wait=wait,output=output)

    ##############################

    def lamp_on(self,wait=1,output=0):

        self.ocssend("LAMPON\n",wait=wait,output=output)

    ##############################

    def lamp_off(self,wait=1,output=0):

        print "LAMPOFF command not yet implemented"

    ##############################

    def open_dome(self,wait=0,output=0):

        if self.verbose:
            print "Sending OPENDOME command."

        if not self.nosock:
            self.ocssend("OPENDOME\n",wait=wait,output=output)
            
    ##############################

    def close_dome(self,wait=0,output=0):

        if self.verbose:
            print "Sending CLOSEDOME command."

        if not self.nosock:
            self.ocssend("CLOSEDOME\n",wait=wait,output=output)

    ##############################

    def move_tel(self,ntarg,wait=0,output=0):

        tgt=self.targets[ntarg]
        [radeg,decdeg]=tgt.coords.deg()
        eqnx=tgt.coords.equinox

        ocsmsg="MOVETEL %.4f %.4f %.3f %s\n" % \
                    (radeg/15.0,decdeg,eqnx,tgt.name)

        self.ocssend(ocsmsg,wait=wait,output=output)

    ##############################

    def move_std(self,nstd,wait=0,output=0):

        tgt=self.standards[nstd]
        [radeg,decdeg]=tgt.coords.deg()
        eqnx=tgt.coords.equinox

        ocsmsg="MOVETEL %.4f %.4f %.3f %s\n" % \
                    (radeg/15.0,decdeg,eqnx,tgt.name)

        self.ocssend(ocsmsg,wait=wait,output=output)

    ##############################

    def move_focus(self,focus,wait=1,output=0):

        ocsmsg="MOVEFOCUS %.3f\n" % focus
        self.ocssend(ocsmsg,wait=wait,output=output)

    ##############################

    def adj_focus(self,dfocus,wait=1,output=0):

        ocsmsg="ADJFOCUS %.3f\n" % dfocus
        self.ocssend(ocsmsg,wait=wait,output=output)
        
    ##############################

    def set_focus(self,filter,focus,wait=1,output=0):

        ocsmsg="SETFOCUS %s %.3f\n" % (filter,focus)
        self.ocssend(ocsmsg,wait=wait,output=output)
        
    ##############################

    def get_focus(self,filter):

        focus=0.0
        self.ocssend("GETFOCUS %s\n" % filter)
        ok=self.wait()
        if ok:
            ocsretval=self.lastocsretval
            rexp=re.search('is\s+([\d\.]+)',ocsretval)
            if rexp:
                try:
                    focus=float(rexp.group(1))
                except:
                    print "Failed to parse focus return string as float"
            else:
                print "Could not read focus value from return string"

        return focus
        
    ##############################

    def take_bias(self,nbias,wait=0,output=0):

        ocsmsg="TAKEBIAS %d\n" % nbias
        self.ocssend(ocsmsg,wait=wait,output=output)
        
    ##############################

    def take_dflat(self,ndome,filter,exptime,wait=0,output=0):

        ocsmsg="DOMEFLAT %s %d %.3f\n" % (filter,ndome,exptime)
        self.ocssend(ocsmsg,wait=wait,output=output)
        
    ##############################

    def take_image(self,nimage,exptime,wait=0,output=0):

        ocsmsg="TAKEIMAGES %d %.3f\n" % (nimage,exptime)
        self.ocssend(ocsmsg,wait=wait,output=output)
        
    ##############################

    def move_filter(self,filter,wait=1,output=0):

        ocsmsg="MOVEFILT %s\n" % filter
        self.ocssend(ocsmsg,wait=wait,output=output)

    ##############################

    def focus_loop(self,fmed=0,nfocus=-1,fstep=-1.0,
                   filter="",exptime=-10.0,wait=0,output=0):

        if fmed<=0:
            fcenter="CURRENT"
        else:
            fcenter="%.3f" % fmed

        if nfocus<0:
            nfocus=self.nfocus

        if fstep<0:
            fstep=self.fstep

        if len(filter)==0:
            filter=self.ffilt

        if exptime<0:
            exptime=self.ftime

        ocsmsg="FOCUSLOOP %s %s %d %.3f %.3f\n" % \
                (filter,fcenter,nfocus,fstep,exptime)
        self.ocssend(ocsmsg,wait=wait,output=output)

    ##############################

    def get_loop_focus(self,filter="",wait=1,output=1):

        if len(filter)==0:
            filter=self.ffilt

        ocsmsg="GETLOOPFOCUS %s\n" % filter
        self.ocssend(ocsmsg,wait=wait,output=output)

    ##############################

    def get_offset(self,filter,wait=1,output=1):

        ocsmsg="GETOFFSET %s\n" % filter
        self.ocssend(ocsmsg,wait=wait,output=output)

    ##############################

    def set_offset(self,filter,offset,wait=1,output=1):

        ocsmsg="SETOFFSET %s %.3f\n" % (filter,offset)
        self.ocssend(ocsmsg,wait=wait,output=output)

    ############################
    #   Generic Communication  #
    ############################

    def generic_tcs(self,command,wait=0,output=0):

        ocsmsg="TCS %s\n" % command
        self.ocssend(ocsmsg,wait=wait,output=output)

    ##############################

    def generic_dhe(self,command,wait=0,output=0):

        ocsmsg="DHE %s\n" % command
        self.ocssend(ocsmsg,wait=wait,output=output)

    ##############################

    def generic_flt(self,command,wait=0,output=0):

        ocsmsg="FLT %s\n" % command
        self.ocssend(ocsmsg,wait=wait,output=output)

######################################################################
