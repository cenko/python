#! /usr/bin/env python
#$Id: ocs.py,v 2.6 2007/07/14 19:12:02 cenko Exp $

from ocs_globals import *
import sys, socket, time, os, re, string, shutil, ephem
import telnetlib, getopt, oscar_session, mx.DateTime
import smtplib
from email.MIMEText import MIMEText
from types import *
from math import *
from iqutils import *
from badflats import *
from check_keywords import *
from oscar_defs import *
from time_now import *

############################################################################
# A first effort at an object-oriented ocs
############################################################################

class ocs_ooriented:

     #######################################################################
     # __init__: Initialization routine.  Always run when an instance of 
     # ocs_ooriented is executed.  Sets all the variables that will be used 
     # by the instances of ocs_ooriented, and little else.  Takes four 
     # optional arguments, TEST and PROC, which determine whether to run in
     # a 'test' directory (1 = Yes, 2 = No), whether to start automatic 
     # data processing pipeline (1 = Yes, 2 = No), whether to run in 
     # asynchronous mode, and whether user is conducting daytime tests.
     #######################################################################

     def __init__(self, TEST=0, PROC=1, ASYNC=1, DAY=0, VERBOSE=0, \
                  SKIPSTART=1, AUTO=1):

          ##################################################################
          # Global variables.  These are defined in the file ocs_globals.py,
          # since there should be easy access to change them if necessary.
          ##################################################################

          # OSS Variables
          self.oss_ip = oss_ip		# IP address of OSS host machine
          self.oss_port = oss_port      # Port for OSS connection
          self.oss_x = oss_x			# OSS executable filename

          # Arcview Variables
          self.aview_ip = aview_ip	# IP Address of Arcview host
          self.comreply_port = comreply_port 	# Port for Command/Reply Chan.
          self.async_port = async_port	# Port for Asynchronous channel
          self.comreply_short_tout = comreply_short_tout
          self.comreply_long_tout = comreply_long_tout
          self.async_tout = async_tout  # Default async timeout
          self.longreply = longreply    # Expected initial reply from long com.

          # Filter List 
          self.filterlist = filterlist  # List of filters CURRENTLY in wheel
          self.rfilterlist = rfilterlist # Reverse filter list
          self.home_filter = home_filter	# Filter of home position

          # String variables
          self.done = done
          self.success = success
          self.fail = fail
          self.ok = ok
          self.rawtxt = rawtxt
          self.error = error
	  self.NOT = NOT
          self.null = null
          self.nofilt = nofilt
          self.redo = redo
          self.wetness = wetness
          self.paramrange = paramrange
          self.illegal = illegal
          self.fixed = fixed
          self.goodnight = goodnight

          # Email Notification Variables
          self.email_subject = email_subject
          self.email_from = email_from
          self.email_to = email_to

          # Standard Stow Position
          self.stow_ha = stow_ha
          self.stow_dec = stow_dec
          self.stow_dome = stow_dome
          self.stow_dome_safe = stow_dome_safe
          self.stow_error = stow_error  # Allow stow pointing error (')
          self.dome_error = dome_error	# Allow dome error (deg)
 
          # Executable Files
          self.p60keys_x = p60keys_x    # p60keys executable file
          self.p60keys1_x = p60keys1_x  # p60keys ex. file (1x thru)
          self.p60redux_x = p60redux_x  # Reduction Pipeline Executable
          self.p60redux1_x = p60redux1_x # Reduction ex file (1x thru)
          self.p60transient_x = p60transient_x # Transfer PTF data->transient 
          self.p60fastfoc_x = p60fastfoc_x # Fast focus reduction script
          self.p60copyfoc_x = p60copyfoc_x # Copy focus files to hipfocus 
          self.datatran_x = datatran_x	# Data Transfer executable

          # Misc. Global Variables
          self.rootdir = rootdir	# Root directory for ocs operations
          self.p60manual = p60manual	# Manual commands file
          self.foctable = foctable	# Current Focus table file
          self.baugiroot = baugiroot	# Root directory on baugi
          self.targetlist = targetlist	# Path to target list
          self.standardlist = standardlist # Path to standard list
          self.telescope_error = telescope_error # Allow pointing error (")
          self.focus_error = focus_error # Allow Focus Error (mm)
          self.offset_error = offset_error # Allow Offset Error (")
          self.image_status = image_status # Image Status File
          self.facility_status = facility_status # Facility Status File
          self.defocus = defocus	# Amount to offset for defoc. images

          ##################################################################
          # State Variables (denoted by '_s' suffix).  These are divided into
          # 7 groups: Telescope Status, Filter Status, CCD Status, Target
          # Status, Socket Status, Process Status, and Misc.  Each is 
          # described in detail
          ##################################################################

          ##################################################################
          # 1) Telescope status
          ##################################################################

          # ready_s: Overall Telescope Readiness to observe.  'READY',  
          # 'NOT_READY', 'UNKNOWN' 
          self.ready_s = 'UNKNOWN' 

          # power_s: Telescope Power. 'READY', 'NOT_READY', 'UNKNOWN'
          self.power_s = 'UNKNOWN'

          # oil_pad_s: Oil Pad Status: 'READY', 'NOT_READY', 'UNKNOWN'
          self.oil_pad_s = 'UNKNOWN' 

          # control_s: Remote Control Status.  Will be one of four values:
          # 'CONSOLE', 'AVAILABLE', 'REMOTE', 'UNKNOWN'
          self.control_s = 'UNKNOWN'

          # telinit_s: 1 => Telescope successfully initialized, 0 => Not yet
          self.telinit_s = 0

          # dome_shutters_s: 'OPEN', 'CLOSED', or 'UNKNOWN'
          self.dome_shutters_s = 'UNKNOWN'

          # dome_motion_s: 'MOVING', 'TRACKING', 'STOPPED', 'UNKNOWN', 
          # 'GO_TO'
          self.dome_motion_s = 'UNKNOWN'

          # focus_position_s: Floating point number in mm
          self.focus_position_s = -1.0

          # focus_motion_s: 'STATIONARY', 'MOVING', 'UNKNOWN'
          self.focus_motion_s = 'UNKNOWN'
       
          # telescope_ra_s: Floating point, decimal hours 
          self.telescope_ra_s = -1.0

          # telescope_dec_s: Floating point, decimal degrees
          self.telescope_dec_s = -1.0

          # telescope_equinox_s: Coordinate Equinox
          self.telescope_equinox_s = 2000.0

          # telescope_mflag_s: 0=proper motion, 1=RA spatial rate, 2=RA 
          # angular rate (1 and 2 non-sidereal motion)
          self.telescope_mflag_s = -1.0 

          # telescope_rarate_s: Rate of RA proper motion, floating point
          self.telescope_rarate_s = -1.0

          # telescope_decrate_s: Rate of DEC proper motion, floating point
          self.telescope_decrate_s = -1.0

          # telescope_raoff_s: RA Offset (arcmin):
          self.telescope_raoff_s = -1.0
          
          # telescope_decoff_s: DEC Offset (arcmin):
          self.telescope_decoff_s = -1.0

          # telescope_motion_s: 'MOVING', 'TRACKING', 'STOPPED', 'UNKNOWN'
          # 'IN_POSITION'
          self.telescope_motion_s = 'UNKNOWN'

          # telescope_hourang_s: Hour angle of telescope 
          self.telescope_hourang_s = -1

          # Current Airmass
          self.telescope_airmass_s = -1.0

          # dome_azimuth_s: If telescope is stowed, angle of dome
          # shutters.  Otherwise -1
          self.dome_azimuth_s = -1

          # weather_s: 'OKAY', 'NOT_OKAY', 'UNKNOWN'
          self.weather_s = 'UNKNOWN' 

          # lamp_s: 'ON', 'OFF', 'UNKNOWN'
          self.lamp_s = 'UNKNOWN'

          # remote_close_s: 'OKAY', 'NOT_OKAY', 'UNKNOWN' 
          self.remote_close_s = 'UNKNOWN' 

          # sunlight_s: 'OKAY', 'NOT_OKAY', 'UNKNOWN' 
          self.sunlight_s = 'UNKNOWN' 

          # mirror_temp_s: Temperature of Primary Mirror
          self.mirror_temp_s = -1

          # floor_temp_s = Temperature of Dome Floor
          self.floor_temp_s = -1

          # bot_tube_temp = Temperature at bottom of telescope tube
          self.bot_tube_temp_s = -1

          # mid_tube_temp_s = Temperature at middle of telescope tube
          self.mid_tube_temp_s = -1

          # top_bute_temp_s = Temperature at top of telescope tube
          self.top_tube_temp_s = -1

          # top_air_temp_s = Air above tube
          self.top_air_temp_s = -1

          # primary_cell_temp_s 
          self.primary_cell_temp_s = -1

          # secondary_cell_temp_s
          self.secondary_cell_temp_s = -1

	  # Outside air temperature
	  self.outside_air_temp_s = -1

          # Pointing coordinates
          self.point_ra_s = self.point_dec_s = self.point_equinox_s = -1.0

          ##################################################################
          # 2) Filter Status
          ##################################################################
          
          # filter_s: Current Filter (or 'UNKNOWN')
          self.filter_s = 'UNKNOWN'

          # filter_motion_s: 'MOVING', 'STATIONARY', 'UNKNOWN'
          self.filter_motion_s = 'UNKNOWN'

          ##################################################################
          # 3) CCD Status
          ##################################################################
 
          # ccd_s: 'IDLE', 'EXPOSE', 'READOUT', 'UNKNOWN'
          self.ccd_s = 'UNKNOWN'
      
          # ccd_temp_s: CCD Temperature as Float
          self.ccd_temp_s = -1.0

          # dhe_temp_s: Electronics Temperature as float
          self.dhe_temp_s = -1.0

          # neck_temp_s: Neck of can temp as float
          self.neck_temp_s = -1.0

          # can_temp_s: Bottom Dewar Can temp as float
          self.can_temp_s = -1.0

          # CCD Idle status
          self.ccd_idle_s = 1

          # CCD Autoclear status
          self.ccd_autoclear_s = 1

	  # CCD Binning
	  self.ccd_bin_s = [-1, -1]

	  # CCD ROI
	  self.ccd_roi_s = [-1, -1, -1, -1]

          ##################################################################
          # 4) Target Status 
          ##################################################################

          # object_s: Object Being Observed
          self.object_s = 'UNKNOWN'
           
          # object_ra_s: Object RA
          self.object_ra_s = -1.0

          # object_dec_s: Object DEC
          self.object_dec_s = 1.0

          # object_equinox_s: Object Equinox
          self.object_equinox_s = 2000.0 

          # object_rarate_s: Object RA Rate
          self.object_rarate_s = -1.0
 
          # object_decrate_s: Object DEC Rate
          self.object_decrate_s = -1.0

          # object_raproper_s: Object RA Proper Motion
          self.object_raproper_s = -1.0

          # object_decproper_s: Object Dec Proper Motion
          self.object_decproper_s = -1.0

          # imgtype_s: Current Image Type: 'BIAS', 'DOMEFLAT', 'PHOTOSTD', 
          # 'FOCUS', 'SAOFOCUS', 'SCIENCE', 'UNKNOWN'
          self.imgtype_s = 'UNKNOWN'

          # p60prid_s: Proposal ID
          self.p60prid_s = 'UNKNOWN'

          # p60prnm_s: Proposal Name
          self.p60prnm_s = 'UNKNOWN'

          # p60prtm_s: Proprietary period in months.  Default = 24
          self.p60prtm_s = 24

          # p60prpi_s: Proprosal PI
          self.p60prpi_s = 'UNKNOWN'

          # moondeg_s: 180 - degrees from moon
          self.moondeg_s = -1.0
 
          # mscsize_s: Size of object mosaic: 'RA steps, DEC steps'.  Each
          # step is 1/2 a FOV
          self.mscsize_s = '1,1'

          # mscposn_s: Position within object mosaic. '0,0' => no mosaic
          self.mscposn_s = '0,0'

          # objtype_s: Object type keyword, from list compiled by AGY
          self.objtype_s = 'UNKNOWN'

          # binroi_s: Binning/ROI setting (used for processing)
          self.binroi_s = 'A'
 
          # defocus_s: Should the images be defocussed? 1 => yes, 0 => no
          self.defocus_s = 0

          # autotrig_s: Is this observation an automated trigger?
          self.autotrig_s = 'F'

          # trigsour_S: If autotrig == 'T', where did the trigger come from?
          self.trigsour_s = 'NONE'

          # generic_keys_s: Generic Keyword facility
          self.generic_keys_s = ''

          #################################################################
          # 5) Socket Status (1 => Connected, 0 => Not Connected)
          #################################################################

          # comreply_s: Command/Reply Channel
          self.comreply_s = 0

          # async_s: Asynchronous Channel
          self.async_s = 0

          # osssock_s: Connection to OSS
          self.osssock_s = 0

          ##################################################################
          # 6) Processes Status (1 => Running, 2 => Not Running)
          #################################################################
 
          # oss_s: Spawned OSS Process
          self.oss_s = 0

          # keywords_s: Process to add keywords and change filename
          self.p60keys_s = 0

          # proc_s: Data Reduction Pipeline
          self.p60redux_s = 0

          # p60transient_s: Transfer PTF data to transient
          self.p60transient_s = 0

          # fastfoc_s: Fast Focus Reduction Process
          self.p60fastfoc_s = 0

          # copyfoc_s: Copy focus files to hipfocus dir
          self.p60copyfoc_s = 0

          # datatran_s: Data Transfer Process
          self.datatran_s = 0

          # arcview_s: Arcview Process
          self.arcview_s = 0

          #################################################################
          # 7) Misc Status Variables
          #################################################################

          # logfile_s: Is the logfile open? 1 => Yes, 0 => No
          self.logfile_s = 0

          # alert_s: Are we observing a GRB alert? 1 => Yes, 0 => No
          self.alert_s = 0 

          #################################################################
          # Non-status Global Variables: Variables needed by various parts
          # of the ocs during execution, divided up into logical groups
          # as well as I could
          #################################################################

          ###############################################################
          # Flag variables.  Set by flag calls to class.  Described above
          ###############################################################

          self.TEST = TEST		# All data to test directory
          self.PROC = PROC		# Automatically Begin data pipeline?
          self.ASYNC = ASYNC            # Run in asynchronous mode?
          self.DAY = DAY	 	# Daytime testing?
          self.VERBOSE = VERBOSE	# Verbose mode?
          self.SKIPSTART = SKIPSTART	# Skip initial standard observation
          self.AUTO = AUTO		# Automated mode (ie spawn oss)?
          
          ################################################################
          # Socket Variables.  All variables will eventually be sockets.  
          # When not active, set to 0
          ################################################################

          self.sock = 0			# Server for osssock
          self.sock2 = 0		# Server for async
          self.osssock = 0		# OSS Socket Connection
 	  self.comreply = 0		# Command/Reply Arcview Channel
	  self.async = 0	        # Async Arcview Channel

          #################################################################
          # Filter Information
          #################################################################
      
          self.focus = ['R', 13.0, 10.0, 1.0] # Current Filter / Focus pair
          self.offset = {}		# Filter offset dictionary

          ################################################################
          # Spawned Processes.  Will contain pointers to processes, 0 for now
          #################################################################

          self.arcview = 0		# Arcview 
          self.oss = 0			# Spawned OSS
          self.p60keys = 0		# Keywords/File rename process
          self.p60redux = 0		# Reduction Pipeline
          self.p60transient = 0		# PTF Data transfer 
          self.p60fastfoc = 0		# Fast Focus pipeline
          self.p60copyfoc = 0		# Copy focus files to hipfocus dir
          self.datatran = 0		# Data Transfer

          #################################################################
          # Misc. Global Variables
          #################################################################

          self.logfile = 0		# Nightly logfile
          self.rawdict = {}		# Dictionary of Transferred Raw Images
          self.procdict = {}		# Dictionary of Processed Files
          self.homedir = ''		# Directory to Write data to
          self.curfocus = 0		# Current [Focus, Filter] Pair
          self.loopnum = 0		# Focus Loop Iteration ID
	  self.stdloop = 0              # Standard Loop Iteration ID
          self.focusnum = 0		# # of images in focus loop
          self.focuscurnum = 0		# # of current image in focus loop
          self.lowfocus = 0		# Lowest focus value in focus loop
          self.focusstep = 0		# Focus Increment in focus loop

          #################################################################
          # Ephemeris Information for Palomar.  Needed to calculate moondeg
          #################################################################

          self.palomar = ephem.Observer()
          self.palomar.long, self.palomar.lat, self.palomar.elev = \
               '-116.863',        '33.3560',        1706.0
          self.palomar.date = ephem.now()
          self.moon = ephem.Moon()
          self.moon.compute(self.palomar)
	  self.tnow = time_now() 

#############################################################################

     #######################################################################
     # comreply_send(self, command, async=1, short=1, wait=1, timeout=-1): 
     # Routine
     # to send 'command' to arcview over the command/reply channel, then deal
     # appropriately with the response.  If short command, final response
     # comes over command/reply channel, and routine will exit.  If long 
     # command, routine will send command over comreply, then wait for 
     # acknowledgement.  If wait = 1, routine will wait for final return 
     # over async channel (with 'timeout' set), otherwise will just return
     # with acknowledgement reply.
     ########################################################################
           
     def comreply_send(self, command, short=1, wait=1, timeout=-1,\
      source='', asynccommand=''):

          if (self.comreply_s == 0):
               print 'comreply_send: Not Connected to ArcVIEW'
               self.log_write('comreply_send: %s\n%s: Not Connected to \
                ArcVIEW' % (command, self.fail))
               return '%s: Not Connected to Arcview' % self.fail

          try:
               if (self.ASYNC == 1) or (short == 1):
                    self.comreply.settimeout(self.comreply_short_tout)
               elif (timeout == -1):
                    self.comreply.settimeout(self.comreply_long_tout)
               else:
                    self.comreply.settimeout(timeout)

               self.comreply.send(command) 
               reply = self.comreply.recv(1024)
               self.log_write('comreply_send: %s\n%s' % (command, reply))

               if (short == 1) or (self.ASYNC == 0) or (wait == 0):
                    return reply 
               elif re.search(self.ok, reply):
                    reply = self.async_recv(source=source,\
                     command=asynccommand,timeout=timeout)
                    return reply 
               else:
                    return '%s: %s' % (self.fail, reply)

          except socket.timeout:
               print 'Timeout from Comreply'
               self.log_write('comreply_send:%s\n%s: TIMEOUT FROM COMREPLY' \
                % (command, self.fail))
               return '%s: TIMEOUT FROM COMREPLY' % self.fail 

     #######################################################################
     # async_recv(self, timeout=-1): Receive reply over asynchronous channel
     # If timeout is specified, will use given timeout, otherwise uses
     # default value of async_tout
     #######################################################################
 
     def async_recv(self, source='', command='', timeout=-1):

          if (self.async_s == 0):
               print 'Not Connected to async Channel'
               self.log_write('async_recv: %s: Not Connected to async channel'\
                % (self.fail)) 
               return '%s: Not Connected to async Channel' % self.fail

          try:
               if (timeout == -1):
                    self.async.settimeout(self.async_tout)
               else:
                    self.async.settimeout(timeout)

               reply = self.async.recv(1024)

               done = 0
               replys = self.parse_async(reply)
               for line in replys:
                    self.log_write('async_recv: %s' % line)
                    self.handle_async(line)
                    if (((source == '') and (command == '')) or \
                     ((re.search('\('+source+'.*\)',line)) and \
                     (re.search(command,line)))) and not \
                     (re.search('_ASYNC',line)):
                         done = 1

               if (done == 1):
                    return self.success
               else:
                    reply = self.async_recv(source=source, \
                                            command=command, timeout=timeout)
                    return reply

          except socket.timeout:
               print 'Timeout from async'
               self.log_write('async_recv: %s: TIMEOUT FROM ASYNC' % self.fail)
               return '%s: TIMEOUT FROM ASYNC' % self.fail

     #####################################################################
     # parse_async(self, reply)
     #####################################################################

     def parse_async(self, message):
     
          try:
               lines = message.splitlines()
               replys = []
               reply = ''
               for line in lines:
                    if re.search('\(TCS:.*\)',line) or \
                     re.search('\(FLT:.*\)',line) or \
                     re.search('\(DHE:.*\)',line) or \
                     re.search('_ASYNC',line): 
                         reply += '%s\n' % line
 	                 replys.append(reply)
                         reply = ''
                    elif (line == ''):
                         pass
                    else:
                         reply += '%s\n' % line
               return replys

          except LookupError:
               print 'Error Parsing Message: %s' % message
               self.log_write('parse_async: Error Parsing Message %s' % \
                message)
               return []          

     #######################################################################
     # handle_async(self, message): Returns null if only received 
     # asynchronous messages.  Otherwise strips off asynchronous messages,
     # deals with them, and returns callback messages
     #######################################################################
     
     def handle_async(self, message):

          # For now, pass on true async messages
          if re.search('_ASYNC', message): 
                    self.log_write('handle_async: %s' % message)

          # Only DHE messages expected are from EXPOSE and INIT
          elif re.search('\(DHE:.*\)', message): 

               if re.search(self.error, message):
                    pass

               elif (re.search('EXPOSE', message)) or \
                (re.search('INIT',message)): 
                    if (re.search(self.done, message)): 
                         self.ccd_s = 'IDLE'
                    else:
                         self.ccd_s = 'UNKNOWN'

               # Otherwise unknown
               else:
                    print 'Unknown Return from Async: %s' % message
                    self.log_write('handle_async: Unknown return: %s' % \
                     message) 
                    self.ccd_s = 'UNKNOWN'

          # Only FLT messages expected are from MOVE and INIT
          elif re.search('\(FLT:.*\)',message): 

               if re.search('INIT',message) or \
                re.search('MOVE', message):
                    if re.search(self.done, message): 
                         reply = self.get_flt_status()
                    else:
                         self.filter_s = 'UNKNOWN'
                         self.filter_motion_s = 'UNKNOWN'

               # Otherwise unknown
               else:
                    print 'Unknown return from async: %s' % message
                    self.log_write('handle_async: Unknown return: %s' % \
                     message) 
                    self.filter_s = 'UNKNOWN'
                    self.filter_motion_s = 'UNKNOWN'

          # Whole slew of TCS messages to deal with
          elif re.search('\(TCS:.*\)', message): 
 
               # Do we have an unknown message?
               unknown = True 

               # telinit
               if re.search('TELINIT', message): 
                    unknown = False 
                    if re.search(self.done, message): 
                         self.telinit_s = 1
                    else:
                         self.telinit_s = 0

               # Opening / Closing Dome 
               if re.search('CLOSE', message) or \
                re.search('OPEN', message): 
                    unknown = False
                    [self.dome_shutters_s, self.dome_motion_s] = \
                     self.get_tcs_status(['Dome_Shutter_Status',\
                     'Dome_Motion_Mode'], '?STATUS')
                    [self.dome_shutters_s, self.dome_motion_s] = \
                     self.set_type([self.dome_shutters_s, \
                     self.dome_motion_s], StringType)

               # Moving the Telescope
               if re.search('GOPOS', message) or re.search('GOREF',message)\
	        or re.search('MOVE', message): 
                    unknown = False
                    # First grab dome information
                    [self.dome_shutters_s, self.dome_motion_s] = \
                     self.get_tcs_status(['Dome_Shutter_Status',\
                     'Dome_Motion_Mode'], '?STATUS')
                    [self.dome_shutters_s, self.dome_motion_s] = \
                     self.set_type([self.dome_shutters_s, self.dome_motion_s],\
                     StringType)
                    # Check to see if the telescope is still settling
                    [self.telescope_motion_s] = self.get_tcs_status( \
                     ['Telescope_Motion_Status'], '?POS')
                    [self.telescope_motion_s] = self.set_type(\
                     [self.telescope_motion_s], StringType)
                    while (self.telescope_motion_s == 'SETTLING'):
                         time.sleep(3.0)
                         [self.telescope_motion_s] = self.get_tcs_status( \
                          ['Telescope_Motion_Status'], '?POS')
                         [self.telescope_motion_s] = self.set_type(\
                          [self.telescope_motion_s], StringType)
                    # When done settling, get all positional info
                    [self.telescope_ra_s, self.telescope_dec_s, \
                     self.telescope_raoff_s, self.telescope_decoff_s, \
                     self.telescope_equinox_s, \
                     self.telescope_motion_s, self.telescope_airmass_s] = \
                     self.get_tcs_status(['Telescope_RA', \
                     'Telescope_Dec', 'Telescope_RA_Offset', \
                     'Telescope_Dec_Offset', 'Telescope_Equinox',\
                     'Telescope_Motion_Status', 'Telescope_Airmass'], '?POS') 
                    # Convert positional info to appropriate type
                    [self.telescope_raoff_s, self.telescope_decoff_s, \
                     self.telescope_airmass_s] = \
                     self.set_type([self.telescope_raoff_s, \
                     self.telescope_decoff_s, self.telescope_airmass_s], \
                     FloatType)
                    [self.dome_shutters_s, self.dome_motion_s, \
                     self.telescope_ra_s, self.telescope_dec_s, \
                     self.telescope_equinox_s, self.telescope_motion_s] = \
                     self.set_type([self.dome_shutters_s, self.dome_motion_s,\
                     self.telescope_ra_s, self.telescope_dec_s, \
                     self.telescope_equinox_s, self.telescope_motion_s], \
                     StringType)
                    self.set_coords()

               # Stow Telescope
               if re.search('STOW', message): 
                    unknown = False
                    # Check to see if the telescope is still settling
                    [self.telescope_motion_s] = self.get_tcs_status( \
                     ['Telescope_Motion_Status'], '?POS')
                    [self.telescope_motion_s] = self.set_type(\
                     [self.telescope_motion_s], StringType)
                    while (self.telescope_motion_s == 'SETTLING'):
                         time.sleep(3.0)
                         [self.telescope_motion_s] = self.get_tcs_status( \
                          ['Telescope_Motion_Status'], '?POS')
                         [self.telescope_motion_s] = self.set_type(\
                          [self.telescope_motion_s], StringType)
                    # When done Settling, get all positional info
                    [self.telescope_hourang_s,\
                     self.telescope_dec_s,\
                     self.dome_azimuth_s] = \
                     self.get_tcs_status(\
                     ['Telescope_HA','Telescope_Dec',\
                     'Dome_Azimuth'], '?POS')
                    [self.telescope_hourang_s,\
                     self.telescope_dec_s] = self.set_type(\
                     [self.telescope_hourang_s,\
                     self.telescope_dec_s], StringType)
                    [self.dome_azimuth_s] = self.set_type(\
                     [self.dome_azimuth_s], FloatType)
                    self.set_coords()

               # Focus Motion
               if re.search('ADJFOCUS',message) or \
                re.search('GOFOCUS',message):
                    unknown = False
                    [self.focus_position_s] = self.get_tcs_status(\
                     ['Focus_Position'], '?POS')
                    [self.focus_position_s] = \
                     self.set_type([self.focus_position_s], FloatType)
                    if re.search(self.done,message):
                         self.focus_motion_s = 'STATIONARY'
                    else:
                         self.focus_motion_s = 'UNKNOWN'

               # Lamp On
               if re.search('LAMPON',message):
                    unknown = False
                    [self.lamp_s] = self.get_tcs_status(\
                     ['Lamp_Status'], '?STATUS')
                    [self.lamp_s] = self.set_type([self.lamp_s], StringType)

               # Otherwise unknown TCS command
               if (unknown == True):
                    print 'Unknown return from async: %s' % message
                    self.log_write('handle_async: Unknown return: %s' % \
                     message) 

          # Unknown command altogether
          else:
               print 'Unknown return from async: %s' % message
               self.log_write('handle_async: Unknown return: %s' % message) 

     #######################################################################
     # osssock_send(self, reply)
     #######################################################################

     def osssock_send(self, reply):

          try: 
               if (re.search(self.nofilt, reply)):
                    self.osssock.send('%s: %s' % (self.fail, self.nofilt))
               elif (re.search(self.redo, reply)):
                    self.osssock.send(reply)
               elif (re.search(self.paramrange, reply)):
                    self.osssock.send(self.illegal)
               elif (re.search(self.fail, reply)):
                    self.osssock.send(self.fail)
               elif (re.search(self.success, reply)):
                    self.osssock.send(self.success)
               else:
                    self.osssock.send(reply)
               print 'Sent %s to Scheduler' % reply
               return

          except (socket.error), e:
               self.log_write('Error Communicating with OSS: %s' % e)
               self.osssock_s = 0
               return

     #######################################################################
     # osssock_recv(self)
     #######################################################################

     def osssock_recv(self):

          try:
               reply = self.osssock.recv(1024)
               return reply
          except (socket.error), e:
               self.log_write('Error Communicating with OSS: %s' % e)
               self.osssock_s = 0
               return

     #######################################################################
     # get_tcs_status(self, strlist, command): 
     # Routine parses through the 
     # return from 'command', looking for strings in the list strlist.
     # If it finds the string, it returns the given value.  If it doesn't,
     # then it returns a -1 for that string.  Returns all values in a 
     # corresponding string list
     #######################################################################

     def get_tcs_status(self, strlist, command):

          retlist = []

          response = self.comreply_send('TCS %s' % command,short=1,wait=1)

          lines = response.splitlines()

          if (len(lines) < 2):
               for i in range(len(strlist)):
                    retlist.append(-1)
               print 'get_tcs_status: Bad Return from ArcVIEW'
               self.log_write('get_tcs_status: Bad Return from ArcVIEW')
               return retlist

          for inp in strlist:
               for j in range(len(lines)):
                    if re.search(inp, lines[j]): 
                         null, status = lines[j].split('=')
                         retlist.append(status)
                         break
                    elif (j == len(lines) - 1):
                         retlist.append(-1)
                         break
                    else:
                         pass
           
          self.log_write('get_tcs_status: %s\n%s' % (strlist, retlist))
          return retlist

     #######################################################################
     # get_flt_status(self):
     #######################################################################

     def get_flt_status(self):

          reply = self.comreply_send('FLT GET position', short=1)
          if not re.search('FILTER1:', reply):
               print 'Bad Return from Filter Status'
               self.log_write('get_flt_status: %s: Bad return from FLT: %s' % \
                (self.fail, reply))
               self.filter_s = 'UNKNOWN'
               self.filter_motion_s = 'UNKNOWN'
               return self.fail
          else:
               try:
                    self.log_write('get_flt_status: %s' % reply)
                    self.filter_s = reply.split()[1].split('(')[0]
                    self.filter_motion_s = 'STATIONARY'
                    return self.success
               except (LookupError,AttributeError,TypeError,ValueError):
                    print 'Bad Return from Filter Status'
                    self.log_write('get_flt_status: %s: Bad Return: %s' % \
                     (self.fail, reply))
                    self.filter_s = 'UNKNOWN'
                    self.filter_motion_s = 'UNKNOWN'
                    return self.fail

     #######################################################################
     # get_temp_status(self, strlist):
     #######################################################################

     def get_temp_status(self, strlist):

          retlist = []

          for name in strlist: 
               reply = self.comreply_send('DHE TP get %s' % name, short=1) 
               if re.search(self.error, reply):
                    retlist.append(-1)
                    self.log_write('get_temp_status: %s: No Temp %s' % \
                     (self.fail, name))
               else:
                    try:
                         temp = reply.split()[0]
                         [temp] = self.set_type([temp], FloatType)
                         retlist.append(temp)
                         self.log_write('get_temp_status: %s = %.2f' % (name, 
                          temp))
                    except (LookupError,ValueError,TypeError,AttributeError):
                         print 'Error Reading Temp %s' % temp
                         retlist.append(-1)
                         self.log_write('get_temp_status: %s: Error Reading\
                          Temp %s: %s' % (self.fail, name, reply))

          return retlist
          
                    
     #########################################################################
     # sanity_check(self, day=0): Basic Check to make sure things are OK for 
     # observing.  First simple networking checks, then verifies disk space
     # on baugi and mimir.  Checks existence of focus table file.
     # Finally, if not day time observations, makes
     # sure telescope and oil pads are powered on.  Returns 'SUCCESS' or
     # 'FAILURE: ' + Error Message
     #########################################################################

     def sanity_check(self):

          print 'Performing Basic Sanity Check ... ',

          check1 = os.system('ping -c 1 espn.com')
          check2 = os.system('ping -c 1 %s' % self.aview_ip)

          if not (check1 == 0) or not (check2 == 0):
               print 'Fatal Networking Error'
               return '%s: Fatal Networking Error' % self.fail

          try:
               check3 = os.popen('df', 'r', -1)
               #usage = check3.readlines()[5].split()[4].split('%')[0]
               usage = check3.readlines()[1].split()[4].split('%')[0]
               if (int(usage) > 80):
                    print 'Not Enough Disk Space on Mimir'
                    return '%s: Not Enough Disk Space on Mimir' % self.fail

               check4 = os.popen('ssh p60@%s \'df\'' % self.aview_ip, 'r', -1) 
               usage = check4.readlines()[3].split()[4].split('%')[0]
               if (int(usage) > 90):
                    print 'Not Enough Disk Space on Baugi'
                    return 'FAILURE: Not Enough Disk Space on Baugi'

          except (LookupError,TypeError,AttributeError,ValueError):
               print 'Bad Reply from Disk Space Check'
               return '%s: Not Enough Disk Space' % self.fail
             

          if not os.path.exists(self.foctable):
               return self.fail + ': Could not find Focus Table'

          if (self.DAY == 0):
               
               [self.power_s, self.oil_pad_s] =  self.get_tcs_status(\
                  ['Telescope_Power_Status','Oil_Pad_Status'], '?STATUS')
               [self.power_s, self.oil_pad_s] = self.set_type(\
                [self.power_s, self.oil_pad_s], StringType)

               if not (self.power_s == 'READY'):
                    print 'No Telescope Power'
                    return '%s: Fatal Telescope Power Error' % self.fail

               elif not (self.oil_pad_s == 'READY'): 
                    print 'No Telescope Oil'
                    return '%s: Fatal Oil Pad Error' % self.fail

          self.log_write('sanity_check: Basic Sanity Check Passed')
          print 'Passed'
          return self.success 

     ########################################################################
     # get_full_status(self): Update all status variables from TCS and CCD
     # that are possible.  To be used when starting up for the night and as a
     # diagnostic to track down errors.  No return.
     # *UNSET STATUS VARIABLES*: telinit, focus_motion, filter,
     # filter_motion, ccd, ccdtemps
     ########################################################################

     def get_full_status(self):

          [self.ready_s, self.power_s, self.oil_pad_s, self.control_s, \
           self.dome_shutters_s, self.dome_motion_s, self.weather_s, \
           self.lamp_s, self.remote_close_s, self.sunlight_s] = \
           self.get_tcs_status(['Telescope_Ready_Status', \
           'Telescope_Power_Status', \
           'Oil_Pad_Status', 'Telescope_Control_Status', \
           'Dome_Shutter_Status', 'Dome_Motion_Mode', 'Weather_Status', \
           'Lamp_Status', 'Remote_Close_Status', 'Sunlight_Status'], \
           '?STATUS')

          [self.focus_position_s, self.telescope_ra_s, self.telescope_dec_s, \
           self.telescope_equinox_s, \
           self.telescope_rarate_s, self.telescope_decrate_s, \
           self.telescope_raoff_s, self.telescope_decoff_s, \
           self.telescope_motion_s, \
           self.telescope_hourang_s, self.telescope_airmass_s,\
           self.dome_azimuth_s] = self.get_tcs_status( \
           ['Focus_Position', 'Telescope_RA', 'Telescope_Dec', \
           'Telescope_Equinox', \
           'Telescope_RA_Rate','Telescope_Dec_Rate','Telescope_RA_Offset',\
           'Telescope_Dec_Offset', 'Telescope_Motion_Status', \
           'Telescope_HA','Telescope_Airmass','Dome_Azimuth'], '?POS')

          [self.mirror_temp_s, self.floor_temp_s, self.bot_tube_temp_s, \
           self.mid_tube_temp_s, self.top_tube_temp_s, self.top_air_temp_s, \
           self.primary_cell_temp_s, self.secondary_cell_temp_s, \
           self.outside_air_temp_s] = \
           self.get_tcs_status(['Mirror_Temp', 'Floor_Temp', 'Bot_Tube_Temp', \
           'Mid_Tube_Temp', 'Top_Tube_Temp', 'Top_Air_Temp', \
           'Primary_Cell_Temp', 'Secondary_Cell_Temp', 'Outside_Air_Temp'], \
           '?WEATHER')

          [self.focus_position_s, self.telescope_rarate_s, \
           self.telescope_decrate_s, self.telescope_raoff_s, \
           self.telescope_decoff_s, \
           self.dome_azimuth_s, \
           self.mirror_temp_s, self.floor_temp_s, self.bot_tube_temp_s, \
           self.mid_tube_temp_s, self.top_tube_temp_s, self.top_air_temp_s, \
           self.primary_cell_temp_s, self.secondary_cell_temp_s,\
           self.outside_air_temp_s, self.telescope_airmass_s] = \
           self.set_type([self.focus_position_s, self.telescope_rarate_s, \
           self.telescope_decrate_s, self.telescope_raoff_s, \
           self.telescope_decoff_s, \
           self.dome_azimuth_s, \
           self.mirror_temp_s, self.floor_temp_s, self.bot_tube_temp_s, \
           self.mid_tube_temp_s, self.top_tube_temp_s, self.top_air_temp_s, \
           self.primary_cell_temp_s, self.secondary_cell_temp_s, \
           self.outside_air_temp_s, self.telescope_airmass_s], FloatType)

          [self.ready_s, self.power_s, self.oil_pad_s, self.control_s, \
           self.dome_shutters_s, self.dome_motion_s, self.weather_s, \
           self.lamp_s, self.remote_close_s, self.sunlight_s, \
           self.telescope_ra_s, self.telescope_dec_s, \
           self.telescope_equinox_s, \
           self.telescope_motion_s,self.telescope_hourang_s] = self.set_type( \
          [self.ready_s, self.power_s, self.oil_pad_s, self.control_s, \
           self.dome_shutters_s, self.dome_motion_s, self.weather_s, \
           self.lamp_s, self.remote_close_s, self.sunlight_s, \
           self.telescope_ra_s, self.telescope_dec_s, \
           self.telescope_equinox_s, \
           self.telescope_motion_s, \
           self.telescope_hourang_s], StringType)
 
          self.set_coords()

          null = self.get_flt_status()

          # Removed to improve efficiency 5 June 2007
#         [self.ccd_temp_s, self.dhe_temp_s, self.neck_temp_s, \
#          self.can_temp_s] = self.get_temp_status(['CCDTEMP', 'DHETEMP', \
#          'NECKTEMP', 'CANTEMP'])

     #######################################################################
     # write_full_status(self)
     #######################################################################

     def write_full_status(self):

          try:
               outfile = open(self.facility_status, 'a+')
               outfile.write('%s\t%s\t\t%s\n' % (time.asctime(), \
                'Telescope_Ready_Status', self.ready_s))
               outfile.write('%s\t%s\t\t%s\n' % (time.asctime(), \
                'Remote_Close_Status', self.remote_close_s))
               outfile.write('%s\t%s\t\t\t%s\n' % (time.asctime(), \
                'Weather_Status', self.weather_s))
               outfile.write('%s\t%s\t%s\n' % (time.asctime(), \
                'Telescope_Control_Status', self.control_s))
               outfile.write('%s\t%s\t\t%s\n' % (time.asctime(), \
                'Dome_Shutter_Status', self.dome_shutters_s))
               outfile.write('%s\t%s\t\t%s\n' % (time.asctime(), \
                'Dome_Motion_Mode', self.dome_motion_s))
               outfile.write('%s\t%s\t\t\t%s\n' % (time.asctime(), \
                'Lamp_Status', self.lamp_s))
               outfile.write('%s\t%s\t\t\t%s\n' % (time.asctime(), \
                'Sunlight_Status', self.sunlight_s))
               outfile.write('%s\t%s\t\t\t%.4f\n' % (time.asctime(), \
                'Telescope_RA', self.telescope_ra_s))
               outfile.write('%s\t%s\t\t\t%.4f\n' % (time.asctime(), \
                'Telescope_Dec', self.telescope_dec_s))
               outfile.write('%s\t%s\t\t%.2f\n' % (time.asctime(), \
                'Telescope_Equinox', self.telescope_equinox_s))
               outfile.write('%s\t%s\t\t%.2f\n' % (time.asctime(), \
                'Telescope_RA_Offset', self.telescope_raoff_s))
               outfile.write('%s\t%s\t\t%.2f\n' % (time.asctime(), \
                'Telescope_Dec_Offset', self.telescope_decoff_s))
               outfile.write('%s\t%s\t\t%s\n' % (time.asctime(), \
                'Telescope_Motion_Status', self.telescope_motion_s))
               outfile.write('%s\t%s\t\t\t%s\n' % (time.asctime(), \
                'Object_Name', self.object_s))
               outfile.write('%s\t%s\t\t\t%.4f\n' % (time.asctime(), \
                'Object_RA', self.object_ra_s))
               outfile.write('%s\t%s\t\t\t%.4f\n' % (time.asctime(), \
                'Object_Dec', self.object_dec_s))
               outfile.write('%s\t%s\t\t\t%.2f\n' % (time.asctime(), \
                'Object_Equinox', self.object_equinox_s))
               outfile.write('%s\t%s\t\t\t%.4f\n' % (time.asctime(), \
                'Telescope_HA', self.telescope_hourang_s))
               outfile.write('%s\t%s\t\t\t%.2f\n' % (time.asctime(), \
                'Dome_Azimuth', self.dome_azimuth_s))
               outfile.write('%s\t%s\t\t\t%.1f\n' % (time.asctime(), \
                'Bot_Tube_Temp', self.bot_tube_temp_s))
               outfile.write('%s\t%s\t\t\t%.2f\n' % (time.asctime(), \
                'CCD_Temp', self.ccd_temp_s))
               outfile.write('%s\t%s\t\t\t%.2f\n' % (time.asctime(), \
                'DHE_Temp', self.dhe_temp_s))
               outfile.write('%s\t%s\t\t\t%.2f\n' % (time.asctime(), \
                'Neck_Temp', self.neck_temp_s))
               outfile.write('%s\t%s\t\t\t%.2f\n' % (time.asctime(), \
                'Can_Temp', self.can_temp_s))
               outfile.write('\n')
               outfile.close()
          except TypeError:
               outfile.write('\n')
               outfile.close()

     #######################################################################
     # check_ready(self)
     #######################################################################

     def check_ready(self):

          [self.ready_s, self.weather_s, self.remote_close_s] = \
           self.get_tcs_status(['Telescope_Ready_Status', 'Weather_Status', \
           'Remote_Close_Status'], '?STATUS')
          [self.ready_s, self.weather_s, self.remote_close_s] = self.set_type(\
           [self.ready_s, self.weather_s, self.remote_close_s], StringType)

          if not (self.weather_s == 'OKAY'):
               print 'Weather Status Not Okay'
               return self.fail
          elif not (self.remote_close_s == 'OKAY'):
               print 'Remote Close Status Not Okay'
               return self.fail
          elif not (self.ready_s == 'READY'):
               reply = self.tel_init(wait=1)
               [self.ready_s] = self.get_tcs_status(\
                ['Telescope_Ready_Status'], '?STATUS')
               if re.search(self.fail, reply) or not (self.ready_s == 'READY'):
                    print 'Telescope not Ready to Observe'
                    return self.fail
               
          print 'Ready to Observe'
          return self.success

     ########################################################################
     # arcview_connect(self): Connects to ARCView 
     # (if not connected already)
     # and attempts to verify connection. No return 
     ########################################################################

     def arcview_connect(self):

          print 'Connecting to ArcVIEW Server',

          try:
               self.comreply = socket.socket(socket.AF_INET, \
                                             socket.SOCK_STREAM)
               self.comreply.settimeout(self.comreply_short_tout)
               self.comreply.connect((self.aview_ip, self.comreply_port))
               self.comreply_s = 1
               self.log_write('arcview_connect: COMREPLY Connected')
          except socket.timeout:
               self.comreply = self.comreply_s = 0
               print 'Timeout Connecting Comreply Channel'
               self.log_write('arcview_connect: %s: Timeout Connecting \
                Comreply Channel' % self.fail) 
               return self.fail

          if (self.ASYNC == 1):
               reply = self.comreply_send('COMSTCP blocking off', short=1) 
          else:
               reply = self.comreply_send('COMSTCP blocking on', short=1) 

          if not re.search(self.done, reply):
               self.comreply = self.comreply_s = 0
               print 'Error Setting Async Status'
               self.log_write('arcview_connect: %s: Error setting async \
                status' % self.fail) 
               return self.fail

          if (self.ASYNC == 1):
               self.sock2 = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
               self.sock2.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
               self.sock2.bind(('', self.async_port))
               self.sock2.listen(5)

               time.sleep(2.0)
               if (self.DAY == 1):
                    reply = self.comreply_send('FLT INIT',short=0,wait=0)
               else:
                    reply = self.take_control()
                    if re.search(self.fail, reply):
                         self.log_write('initialize: %s: Error Taking Control \
                          of Telescope: %s' % (self.fail, reply))
                         return self.fail
                    reply = self.comreply_send('TCS ADJFOCUS 0.01',short=0,\
                     wait=0)
                    time.sleep(5.0)
               self.sock2.settimeout(self.async_tout)

               try:
                    self.async, self.async_address = self.sock2.accept()
                    reply2 = self.async.recv(1024)
                    while re.search('_ASYNC', reply2):
                         reply2 = self.async.recv(1024)
                    if re.search(self.done, reply2):
                         self.async_s = 1
                         self.log_write('arcview_connect: ASYNC Connected')
                    else:
                         self.async = self.async_s = 0
                         print 'Error Connecting Async Channel'
                         self.log_write('arcview_connect: %s: ASYNC NOT \
                          CONNECTED' % self.fail)
                         return self.fail
               except socket.timeout:
                    self.async = self.async_s = 0
                    self.log_write('arcview_connect: %s: ASYNC NOT \
                     CONNECTED' % self.fail)
                    return self.fail

          return self.success
                    
     #########################################################################
     # arcview_spawn(self): If it isn't already running, spawns ArcVIEW
     # process, and verifies succes.  Returns success or failure. 
     ######################################################################## 

     def arcview_spawn(self):

          print 'Spawning ArcVIEW process on baugi',

          test = os.popen('ssh p60@%s \'ps xa | grep ArcVIEW\'' % \
           self.aview_ip, 'r', -1) 
          lines = test.readlines()

          if (len(lines) >= 3):
             self.arcview_s = 1
             self.log_write('arcview_spawn: ArcVIEW already running')
             return self.success
          else:
               self.arcview1 = os.popen('ssh p60@%s start_arcviewVNC_remote' \
                % self.aview_ip, 'r', -1)
               # Below command is new for baugi2
               time.sleep(30)
               self.arcview2 = os.popen('ssh -f p60@%s \
                 \'start_ArcVIEW_Server\'' % self.aview_ip, 'r', -1)
               time.sleep(30)
               
               if (self.comreply_s == 1):
                    self.comreply.close()
                    self.comreply = self.comreply_s = 0
               if (self.async_s == 1):
                    self.async.close()
                    self.async = self.async_s = 0
               if not (self.sock2 == 0):
                    self.sock2.close()
                    self.sock2 = 0

               test = os.popen('ssh p60@%s \'ps xa | grep ArcVIEW\'' % \
                self.aview_ip, 'r', -1) 
               lines = test.readlines()

               if (len(lines) < 3):
                    self.arcview_s = self.arcview = 0
                    print 'Error Spawning ArcVIEW'
                    self.log_write('arcview_spawn: %s Spawning ArcVIEW' % \
                     self.fail)
                    return self.fail
               else:
                    self.arcview_s = 1
                    self.log_write('arcview_spawn: %s: ArcVIEW Spawned' % \
                     self.success)
                    return self.success

     ########################################################################
     # arcview_close(self): Kills ArcVIEW process (if possible). No return 
     ########################################################################      
     def arcview_close(self): 
    
          print 'Closing ArcVIEW process on baugi ...',

          test = os.popen('ssh p60@%s \'ps xa | grep ArcVIEW\'' % \
           self.aview_ip, 'r', -1)
          lines = test.readlines()
                                                                                
          if (len(lines) < 3):
               pass
          else:
               #os.system('ssh p60@%s \'APP SHUTDOWN\'' % self.aview_ip)
               # Below command is new for baugi2
               os.system('ssh p60@%s shutdown_ArcVIEW_Server' % self.aview_ip)
               time.sleep(30)

          self.arcview_s = 0
          if not (self.arcview == 0):
               self.arcview.close()
               self.arcview = 0

          if (self.comreply_s == 1):
               self.comreply.close()
               self.comreply = self.comreply_s = 0
          if (self.async_s == 1):
               self.async.close()
               self.async = self.async_s = 0
          if not (self.sock2 == 0):
               self.sock2.close()
               self.sock2 = 0

          print 'Done'
          self.log_write('arcview_close: Done')
          return

     #########################################################################
     # arcview_disconnect(self): No return
     #########################################################################

     def arcview_disconnect(self):

          print 'Disconnecting from ArcVIEW Server ...', 

          if not (self.comreply == 0):
               self.comreply.shutdown(2)
               self.comreply.close()
               self.comreply = self.comreply_s = 0
          if not (self.async == 0):
               self.async.shutdown(2)
               self.async.close()
               self.async = self.async_s = 0
          if not (self.sock2 == 0):
               self.sock2.shutdown(2)
               self.sock2.close()
               self.sock2 = 0

          print 'Done'
          self.log_write('arcview_disconnect: %s' % self.success)
          return

     ########################################################################
     # take_break(self, stime=900.0): No return
     ########################################################################

     def take_break(self, stime=900.0):
          
          print 'Powering down ArcVIEW and CCD for a rest ...' 

          # Close ArcVIEW
          self.arcview_disconnect()
          self.arcview_close()

          # Telnet to tps
          try:
               tps = telnetlib.Telnet('198.202.125.197')
               tps.write('/Off 1\r')
               tps.write('/X\r')
               tps.close()
               time.sleep(stime)
               tps = telnetlib.Telnet('198.202.125.197')
               tps.write('/On 1\r')
               tps.write('/X\r')
               tps.close()
               self.log_write('take_break: %s' % self.success)
               print 'Done'
               return

          # Simply report failure to logfile if unsuccessful
          except socket.error:
               self.log_write('take_break: %s: Could not connect to TPS' %
                self.fail)
               print 'Failed.'
               return
                
     ########################################################################
     # recycle_power(self): No return
     ########################################################################

     def recycle_power(self, dev='both'):
     
          print 'Recycle Power in Progress ...',

          try:
               tps = telnetlib.Telnet('198.202.125.197')
               if (dev == 'dhe') or (dev == 'both'):
                    tps.write('/BOOT 1\r')
                    time.sleep(20)
                    response = tps.read_eager()
               if (dev == 'flt') or (dev == 'both'):
                    tps.write('/BOOT 2\r')
                    time.sleep(20)
                    response = tps.read_eager()
               tps.write('/X\r')
               time.sleep(5)
               tps.close()
               self.log_write('recycle_power: %s' % self.success)
               print 'Done'
               return

          except socket.error:
               self.log_write('recycle_power: %s: Could not connect to TPS' %\
                self.fail) 
               return

     ######################################################################
     # set_coords(self, command):
     ######################################################################

     def set_coords(self):
          
          try:

               if (self.telescope_equinox_s == 'UNKNOWN') or \
                (self.telescope_ra_s == 'UNKNOWN') or \
                (self.telescope_dec_s == 'UNKNOWN') or \
                (self.telescope_hourang_s == 'UNKNOWN'): 
                    self.telescope_equinox_s = -1.0
                    self.telescope_ra_s = -1.0
                    self.telescope_dec_s = -1.0
                    self.hourang_s = -1.0
                    return

               if (type(self.telescope_equinox_s) == StringType):
                    self.telescope_equinox_s = \
                     float(self.telescope_equinox_s[1:])
               if (type(self.telescope_ra_s) == FloatType):
                    self.telescope_ra_s = self.telescope_ra_s * 15.0
               if (type(self.telescope_hourang_s) == FloatType):
                    self.telescope_hourang_s = \
                     self.telescope_hourang_s * 15.0

               telescope1 = astrocoords(self.telescope_ra_s, \
                self.telescope_dec_s, equinox=self.telescope_equinox_s)
               self.telescope_ra_s = telescope1.radeg() / 15.0
               self.telescope_dec_s = telescope1.dcdeg()

               telescope2 = astrocoords(self.telescope_hourang_s, \
                self.telescope_dec_s)
               self.telescope_hourang_s = telescope2.radeg() / 15.0

               return

          except (TypeError,ValueError,LookupError,AttributeError):
               print 'Error with Coordinates: %s' % command
               self.log_write('set_coords: %s: Error with command %s' % \
                (self.fail, command))
               return

     ########################################################################
     # set_type(self, inlist, type)
     ########################################################################

     def set_type(self, inlist, intype):
          
          retlist = []
          for object in inlist:
               if (type(object) == intype):
                    pass
               elif (intype == FloatType):
                    try:
                         object = float(object)
                    except (ValueError, TypeError):
                         print 'Error Converting %s to float' % object
                         object = -1.0
               elif (intype == StringType):
                    try:
                         object = string(object)
                    except (ValueError, TypeError):
                         print 'Error Converting %s to string' % object
                         object = 'UNKNOWN'
               elif (intype == IntType):
                    try:
                         object = int(object)
                    except (ValueError, TypeError):
                         print 'Error Converting %s to integer' % object      
                         object = -1
               else:
                    print 'Unexpected Type: %s: No Change Made' % object
               retlist.append(object) 

          return retlist
               
     ########################################################################
     # reset_controller(self): Returns 'SUCCESS' or 'FAILURE'
     ########################################################################

     def reset_controller(self):

          print 'Resetting Electronics Controller ... ',

          reply = self.comreply_send('DHE CLOSE', short=1)
          if not re.search(self.done, reply):
               print 'Could Not Close Shutter'
               self.log_write('reset_controller: %s: Could Not Close Shutter' \
                % self.fail)
               return '%s: %s' % (self.fail, reply)
          
          time.sleep(3.0)

          reply = self.comreply_send('DHE INIT', \
           short=0, wait=1, source='DHE', asynccommand='INIT')

          if (self.ASYNC == 0):
               if re.search(self.error,reply):
                    print 'Count not Reinitialize Controller'
                    self.ccd_s == 'UNKNOWN'
                    self.log_write('reset_controller: %s: Could not \
                     Reinitialize Controller: %s' % (self.fail, reply))
                    return '%s: %s' % (self.fail, reply)
               else:
                    self.ccd_s == 'IDLE'
		    self.ccd_idle_s = self.ccd_autoclear_s = 1
		    self.ccd_bin_s = [1, 1]
		    self.ccd_roi_s = [1, 1, 2048, 2048]
                    print 'Done'
                    return self.success
          elif (self.ASYNC == 1):
               if (self.ccd_s == 'IDLE'):
		    self.ccd_idle_s = self.ccd_autoclear_s = 1
		    self.ccd_bin_s = [1, 1]
		    self.ccd_roi_s = [1, 1, 2048, 2048]
                    print 'Done'
                    return self.success
               else:
                    print 'Count not Reinitialize Controller'
		    self.ccd_idle_s = self.ccd_autoclear_s = -1
		    self.ccd_bin_s = [-1, -1]
		    self.ccd_autoclear_s = [-1, -1, -1, -1]
                    self.log_write('reset_controller: %s: Could not \
                     Reinitialize Controller: %s' % (self.fail, reply))
                    return '%s: %s' % (self.fail, reply)

     #######################################################################
     # get_homedir(self): Sets homedir global variable, which specifies
     # current UT year, month, and day.  No return
     #######################################################################

     def get_homedir(self):

          if (self.DAY == 1) or (self.TEST == 1):
               self.homedir = 'test'
          else:
               uttime = mx.DateTime.gmtime()
               date = uttime.date
               homedir = date.split('-')
               self.homedir = homedir[0] + homedir[1] + homedir[2]

          self.log_write('get_homedir: %s' % self.homedir)
          return
              
     ########################################################################
     # create_dirs(self): Creates nightly directories.  No return
     ########################################################################
   
     def create_dirs(self):
  
          if not os.path.exists('%s%s' % (self.rootdir, self.homedir)):
               os.mkdir('%s%s' % (self.rootdir, self.homedir))
          if not os.path.exists('%s%s/raw' % (self.rootdir, self.homedir)):
               os.mkdir('%s%s/raw' % (self.rootdir, self.homedir))
          if not os.path.exists('%s%s/calib' % (self.rootdir, self.homedir)): 
               os.mkdir('%s%s/calib' % (self.rootdir, self.homedir))
          if not os.path.exists('%s%s/proc' % (self.rootdir, self.homedir)):
               os.mkdir('%s%s/proc' % (self.rootdir, self.homedir))
          if not os.path.exists('%s%s/misc' % (self.rootdir, self.homedir)):
               os.mkdir('%s%s/misc' % (self.rootdir, self.homedir))
          if not os.path.exists('%s%s/hipfocus' % (self.rootdir,self.homedir)):
               os.mkdir('%s%s/hipfocus' % (self.rootdir, self.homedir))
          if not os.path.exists('%s%s/flats' % (self.rootdir, self.homedir)):
               os.mkdir('%s%s/flats' % (self.rootdir, self.homedir))

          os.system('chmod -R 775 %s%s' % (self.rootdir, self.homedir))

          os.system('ssh p60@%s \'mkdir %s%s\'' % (self.aview_ip, \
           self.baugiroot, self.homedir))

          self.log_write('create_dirs: %s: Created directories in %s%s and \
           %s' % (self.success, self.rootdir, self.homedir, self.baugiroot))
          return

     ########################################################################
     # read_foctable(self): Reads in focus table (including current focus 
     # value pair and filter offsets) and writes them to the variables
     # 'focus' and 'offset'.  No return
     ########################################################################

     def read_foctable(self):
   
          print 'Reading in Focus Table ... ',

          try:
               foctable = open(self.foctable, 'r', -1)
               null = foctable.readline()

               filter, focus, temp, airmass = foctable.readline().split()
               self.focus = [filter, float(focus), float(temp), float(airmass)]

               lines = foctable.readlines()
               for i in range(len(lines)):
                    filter, offset = lines[i].split()
                    self.offset[filter] = float(offset)

               foctable.close()
               print 'Done'
               self.log_write('read_foctable: %s' % self.success) 
               return

          except (IOError,ValueError,TypeError,LookupError,AttributeError):
               print 'Error Reading Focus Table'
               self.log_write('read_foctable: %s' % self.fail)
               return 

     #######################################################################
     # write_foctable(self)
     #######################################################################

     def write_foctable(self):

          print 'Writing Focus Table ... ',

          try:
               foctable = open(self.foctable, 'w+', -1)
               foctable.write('%s\n' % time.asctime())
               foctable.write('%s\t%.2f\t%.1f\t%.2f\n' % \
                (self.filterlist[self.focus[0]], \
                self.focus[1], self.focus[2], self.focus[3]))
               for filter in self.offset:
                    foctable.write('%s\t%.2f\n' %(filter, self.offset[filter]))
               foctable.close()
               print 'Done'
               self.log_write('write_foctable: %s' % self.success)
               return
          except (IOError,ValueError,TypeError,LookupError,AttributeError):
               print 'Error Writing Focus Table'
               self.log_write('write_foctable: %s: Error writing focus table'\
                % self.fail)
               return

     #######################################################################
     # open_logfile(self): Opens nightly log file.  No return.
     #######################################################################

     def open_logfile(self):

          print 'Opening Logfile ... ',

          os.chdir('%s%s' % (self.rootdir, self.homedir))
          self.logfile = open('%s.log' % self.homedir, 'a+', -1)
          self.logfile_s = 1
          os.chdir(self.rootdir)
          self.log_write('open_logfile: %s' % self.success)
          print 'Done'
          return

     #######################################################################
     # close_logfile(self)
     #######################################################################

     def close_logfile(self):
     
          print 'Closing Logfile ... ',

          if (self.logfile_s == 1):
               self.logfile.close()
               self.logfile_s = 0
          self.log_write('close_logfile: %s' % self.success)
          print 'Done'
          return

     #######################################################################
     # spawn_keys(self): Spawns process to write additional FITS header
     # keywords to image files (and rename them).  No return.
     #######################################################################

     def spawn_keys(self):
       
          print 'Spawning Keyword Process ... ',
          os.chdir('%s%s/raw' % (self.rootdir, self.homedir))
          self.p60keys = os.popen(self.p60keys_x, 'r', -1)
          self.p60keys_s = 1
          os.chdir(self.rootdir)          
          self.log_write('spawn_keys: %s' % self.success)
          print 'Done'
          return

     ######################################################################
     # spawn_redux(self)
     ######################################################################

     def spawn_redux(self):

          print 'Spawing Reduction Pipeline ... ', 

          # First Spawn the Process
          os.chdir('%s%s/proc' % (self.rootdir, self.homedir))
          self.p60redux = os.popen(self.p60redux_x, 'r', -1)
          self.p60redux_s = 1

          self.p60transient = os.popen(self.p60transient_x, 'r', -1)
          self.p60transient_s = 1

          # Create Symbolic Link for OSS
          try:
               os.remove('%s%s' % (self.rootdir, self.image_status))
          except:
               pass
               
          os.symlink('%s%s/proc/%s' % (self.rootdir, self.homedir,\
                                       self.image_status),
                     '%s%s' % (self.rootdir, self.image_status))

          os.chdir(self.rootdir)
          self.log_write('spawn_redux: %s' % self.success)
          print 'Done'
          return

     ######################################################################
     # spawn_fastfoc(self)
     ######################################################################

     def spawn_fastfoc(self):

          print 'Spawing Fast Focus Pipeline ... ',

          # First Spawn the Process
          os.chdir('%s%s/hipfocus' % (self.rootdir, self.homedir))
          self.p60fastfoc = os.popen(self.p60fastfoc_x, 'r', -1)
          self.p60fastfoc_s = 1

          # And the copy process
          os.chdir('%s%s/raw' % (self.rootdir, self.homedir))
          self.p60copyfoc = os.popen(self.p60copyfoc_x, 'r', -1)
          self.p60copyfoc_s = 1

          os.chdir(self.rootdir)
          self.log_write('spawn_fastfoc: %s' % self.success)
          print 'Done'
          return

     ######################################################################
     # close_redux(self)
     ######################################################################

     def close_redux(self):

          print 'Closing Reduction Pipeline ... ',

          if (self.p60redux_s == 1):
               self.p60redux.close()
               self.p60redux_s = 0
          if (self.p60transient_s == 1):
               self.p60transient.close()
               self.p60transient_s = 0
          if (self.p60keys_s == 1):
               self.p60keys.close()
               self.p60keys_s = 0
          if (self.p60fastfoc_s == 1):
               self.p60fastfoc.close()
               self.p60fastfoc_s = 0
          if (self.p60copyfoc_s == 1):
               self.p60copyfoc.close()
               self.p60copyfoc_s = 0

          self.log_write('close_redux: %s' % self.success)
          print 'Done'
          return

     ######################################################################
     # open_dict(self): Opens directory dictionary (used to keep track of
     # which images have already been transferred from baugi).  No return
     ######################################################################
           
     def open_dict(self):
     
          print 'Opening Raw File Dictionary ... ',
          dircommand = 'ssh p60@%s \'ls %s%s/*.fits\'' % (self.aview_ip, \
           self.baugiroot, self.homedir) 
          temp = os.popen(dircommand, 'r', -1)
          lines = temp.readlines()
          for line1 in lines:
               try:
                    self.rawdict[line1.split('\n')[0]] = 1
               except (IndexError,ValueError):
                    self.log_write('open_dict: Bad File from baugi')
          print 'Done'

          print 'Opening Processed File Dictionary ...',
          os.chdir('%s%s/proc' % (self.rootdir, self.homedir))
          dir = os.listdir('.')
          for file in dir:
               if re.search('20\d*p.fits',file) and not \
                re.search('stars', file):
                    self.procdict[file] = 1
          os.chdir(self.rootdir)
          print 'Done'

          self.log_write('open_dict: %s' % self.success)
          return

     ######################################################################
     # set_obs_params(self, obsname, imagenumber=-1, imagestoread=1)
     # Sets observation parameters for a given exposure.  Returns either 
     # success of failure
     ######################################################################

     def set_obs_params(self, obsname, imagenumber=-1, imagestoread=1, \
                        exptime=-1):

          print 'Setting Observation Parameters ... ',

          setcommand = 'DHE SET multipleextensions yes, displayimage =\
                       no, imagestoread %i, rootname %s%s/%s' % (imagestoread,\
                       self.baugiroot, self.homedir, obsname)
          if not (imagenumber == -1):
               setcommand += ', imagenumber %i' % imagenumber
          if not (exptime == -1):
               setcommand += ', exposuretime = %.2f[s]' % exptime

          reply = self.comreply_send(setcommand, short=1) 
          if re.search(self.done, reply):
               self.log_write('set_obs_params: %s\n%s' % (self.success,  
                setcommand))
               print 'Done'
               return self.success
          else:
               print 'Error Setting Parameters'
               self.ccd_s = 'UNKNOWN'
               self.log_write('set_obs_params: %s\n%s' % (self.fail, \
                setcommand))
               return '%s: %s' % (self.fail, reply)

     ######################################################################
     # flt_init(self): Initialize filter wheel.  Returns success or failure
     ######################################################################

     def flt_init(self, wait=0):

          print 'Initializing Filter Wheel ... ',

          self.filter_motion_s = 'MOVING'
          reply = self.comreply_send('FLT INIT', short=0, \
           wait=wait, timeout=120, source='FLT', asynccommand='INIT')

          if (self.ASYNC == 0) and (re.search(self.done,reply)): 
               self.filter_s = self.home_filter
               self.log_write('flt_init: %s' % self.success)
               print 'Done'
               return self.success
          elif (self.ASYNC==1) and (wait == 0) and (re.search(self.ok,reply)):
               print 'Okay'
               return self.success
          elif (self.ASYNC == 1) and (wait == 1):
               if (self.filter_s == self.home_filter) and \
                (self.filter_motion_s == 'STATIONARY'):
                    print 'Done'
                    return self.success
               else:
                    print 'Error Initializing Filter Wheel'
                    self.log_write('flt_init: %s\n%s' % (self.fail, reply))
                    return '%s: %s' % (self.fail, reply)
          else:
               print 'Error Initializing Filter Wheel'
               self.filter_s = 'UNKNOWN'
               self.filter_motion_s = 'UNKNOWN'
               self.log_write('flt_init: %s\n%s' % (self.fail, reply))
               return '%s: %s' % (self.fail, reply)

     ######################################################################
     # oss_start(self): Starts scheduler.  If in automated mode (self.AUTO
     # == 1, will spawn the scheduler in a new process.  Otherwise, waits
     # for user to begin scheduler and then connect to ocs.  Performs 
     # preliminary systems to check to verify communication
     ######################################################################
     
     def oss_start(self, raw=0, fwait=2.0, swait=6.0):

          print 'Starting Scheduler' 

          if (self.AUTO == 1):
               self.oss_spawn(raw=raw, fwait=fwait, swait=swait)
          else:
               print 'Please initiate scheduler to begin observations'

          # Open Socket Connection to OSS
          try:
               self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
               self.sock.bind((self.oss_ip, self.oss_port))
               self.sock.listen(5)

               self.osssock, self.oss_address = self.sock.accept()
               self.osssock_s = 1

               # Test Connection to OSS
               greeting = self.osssock.recv(1024)
               if (greeting == 'HI'):
                    print 'Good Handshake with OSS'
                    self.osssock.send('READY')
                    self.log_write('oss_start: %s' % self.success) 
                    return self.success
               else:
                    print 'Bad Handshake with OSS'
                    self.log_write('oss_start: %s' % self.fail)
                    return self.fail
 
          except socket.error:
               print 'Error Connecting to OSS'
               self.log_write('oss_start: %s: Socket Error' % self.fail)
               return self.fail

     ######################################################################
     # oss_spawn(self): Spawns scheduler (to be used in automated mode).
     # No return.   
     ######################################################################

     def oss_spawn(self, raw=0, fwait=2.0, swait=6.0):

          print 'Spawning Automated Scheduling Process ... ',

          os.chdir(self.rootdir)
          osscommand = self.oss_x
          if (raw == 1):
               osscommand += ' -r'
          if (self.SKIPSTART == 1):
               osscommand += ' -o'
          osscommand += ' -f %.1f -s %.1f' % (fwait, swait)
          osscommand += ' %s' % self.targetlist
          osscommand += ' %s' % self.standardlist
          osscommand += ' >> /dev/null'

          self.oss = os.popen(osscommand, 'r', -1)
          self.oss_s = 1
          self.log_write('oss_spawn: %s' % self.success)
          print 'Done'
          return

     ########################################################################
     # oss_disconnect(self)
     ########################################################################

     def oss_disconnect(self):

          print 'Closing Down Scheduler'

          try:
               if (self.osssock_s == 0):
                    self.log_write('oss_disconnect: OSS Already disconnected')
                    return self.success
               else:
                    self.osssock.send('DISCONNECT')
                    self.osssock.close()
                    self.sock.close()
                    self.osssock = 0
                    self.osssock_s = 0
                    self.sock = 0
                    time.sleep(120)
                    self.log_write('oss_disconnect: %s' % self.success)
                    return self.success
          except socket.error:
               self.log_write('oss_disconnect: %s' % self.fail)
               return self.fail

     #######################################################################
     # oss_close(self)
     #######################################################################

     def oss_close(self):

          try:
               if (self.oss == 0):
                    self.log_write('oss_close: OSS Process already Closed')
                    return self.success
               else:
                    self.osssock.close()
                    self.sock.close()
                    self.oss.close()
                    self.oss = 0
                    self.osssock = 0
                    self.osssock_s = 0
                    self.sock = 0
                    self.log_write('oss_close: %s' % self.success)
                    return self.success
          except socket.error:
               self.log_write('oss_close: %s' % self.fail)
               return self.fail

     ########################################################################
     # log_write(self, str): Prints messages to screen and, if applicable, 
     # writes them to the logfile
     #######################################################################

     def log_write(self, str):
 
          if (self.VERBOSE == 1):
               print str
          if (self.logfile_s == 1):
               self.logfile.write('%s\n' % str)
               self.logfile.write('%s\n\n' % time.asctime(time.gmtime()))

     #######################################################################
     # take_control(self)
     ########################################################################

     def take_control(self):

          print 'Taking Control of the Telescope ... ',

          reply = self.comreply_send('TCS OBTAIN', short=1)

          if (re.search(self.done, reply)): 
	       print 'Done'
	       self.control_s = 'REMOTE'
	       self.log_write('take_control: %s' % self.success)
	       return self.success 
          else:
               self.log_write('take_control: %s' % self.fail)
	       self.control_s = 'UNKNOWN'
               print 'Error Taking Control'
               return '%s: %s' % (self.fail, reply)

     #######################################################################
     # tel_init(self, wait=1)
     #######################################################################

     def tel_init(self, wait=1):

          print 'Initializing Telescope ... ',

          reply = self.comreply_send('TCS TELINIT', short=0, \
           wait=wait, timeout=240, source='TCS', asynccommand='TELINIT')
          
          if (self.ASYNC == 0) and re.search(self.done, reply):
               self.telinit_s = 1
               self.log_write('tel_init: %s' % self.success)
               print 'Done'
               return self.success 
          elif (self.ASYNC==1) and (wait == 0) and (re.search(self.ok,reply)):
               print 'Okay'
               return self.success
          elif (self.ASYNC == 1) and (wait == 1):
               if (self.telinit_s == 1):
                    print 'Done'
                    return self.success
               else:
                    print 'Error Initializing Telescope'
                    return '%s: TCS TELINIT: %s' % (self.fail, reply)
          else:
               print 'Error Initializing Telescope'
               self.telinit_s = 0
               self.log_write('tel_init: %s\n%s' % (self.fail, reply)) 
               return '%s: TCS TELINIT: %s' % (self.fail, reply)

     ######################################################################
     # misc_init(self)
     ######################################################################

     def misc_init(self):

          self.curfocus = 0
          self.loopnum = 0
          self.stdloop = 0
          self.focusnum = 0
          self.focuscurnum = 0

          print 'Copying Files From Previous Morning ... ',

          if os.path.exists('%scalib/BPM.pl' % self.rootdir):
               shutil.copyfile('%scalib/BPM.pl' % self.rootdir, \
                '%s%s/proc/BPM.pl' % (self.rootdir, self.homedir))

          if os.path.exists('%smorning/proc/Bias.fits' % self.rootdir):
               shutil.copyfile('%smorning/proc/Bias.fits' % self.rootdir, \
                '%s%s/proc/Bias.fits' % (self.rootdir, self.homedir))
          elif os.path.exists('%scalib/Bias.fits' % self.rootdir):
               shutil.copyfile('%scalib/Bias.fits'  % self.rootdir, \
                '%s%s/proc/Bias.fits' % (self.rootdir, self.homedir))

          for filter in self.filterlist:
               if os.path.exists('%smorning/proc/Flat-%s.fits' % \
                (self.rootdir, filter)):
                    shutil.copyfile('%smorning/proc/Flat-%s.fits' % \
                     (self.rootdir, filter), '%s%s/proc/Flat-%s.fits' \
                     % (self.rootdir, self.homedir, filter))
               elif os.path.exists('%scalib/Flat-%s.fits' % (self.rootdir, \
                filter)):
                    shutil.copyfile('%scalib/Flat-%s.fits' % (self.rootdir, \
                     filter), '%s%s/proc/Flat-%s.fits' \
	             % (self.rootdir, self.homedir, filter))

          os.system('cp %smorning/raw/20*r.fits %s%s/raw' % (self.rootdir, \
           self.rootdir, self.homedir))

          print 'Done'
          return

     #######################################################################
     # initialize(self)
     #######################################################################

     def initialize(self, aclose=0):

          if (aclose==1):
               self.arcview_disconnect()
               self.arcview_close()
          reply = self.arcview_spawn()
          if re.search(self.fail, reply):
               self.log_write('initialize: Could not spawn ArcVIEW')
               return self.fail
          null = self.arcview_disconnect()
          reply = self.arcview_connect()
          if re.search(self.fail, reply):
               self.log_write('initialize: Could not connect to ArcVIEW')
               return self.fail
          self.ccd_s = 'IDLE'
          #reply = self.reset_controller()
          #if re.search(self.fail, reply):
               #self.log_write('initialize: Could not reset controller')
               #reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
               #if re.search(self.fail, reply2):
                    #print 'Could Not Recover from error'
                    #return self.fail
               #else:
                    #print 'Recovered from Error'

          reply = self.sanity_check()
          if re.search(self.fail, reply):
               self.log_write('initialize: %s: Failed Sanity Check: %s' % \
                (self.fail, reply))
               return self.fail

          self.get_homedir()
          self.create_dirs()
          self.read_foctable()
          self.open_logfile()
          self.open_dict()
          self.spawn_keys()
          self.misc_init()

          obsnum = self.read_obsnum()
          reply = self.set_obs_params('null', imagenumber=obsnum, \
                                       imagestoread=1)
          if re.search(self.fail, reply):
               self.log_write('initialize: Could not set night observational \
                              parameters')
               return self.fail
    
          reply = self.flt_init(wait=1)
          if re.search(self.fail, reply):
               self.log_write('initialize: %s: Error Initializing filter \
                Wheel: %s' % (self.fail, reply)) 
               return self.fail

          self.get_full_status()
          self.write_full_status()

          reply = self.take_control()
          if re.search(self.fail, reply):
               self.log_write('initialize: %s: Error Taking Control of \
                Telescope: %s' % (self.fail, reply)) 
               return self.fail

          if (self.DAY == 0) and (self.TEST == 0) and (self.weather_s=='OKAY')\
           and (self.remote_close_s == 'OKAY') and \
           (not self.ready_s == 'READY'):
               reply = self.tel_init(wait=1)
               if re.search(self.fail, reply):
                    self.log_write('initialize: %s: Could not initialize \
                     Telescope: %s' % (self.fail, reply))
                    return self.fail
          elif not (self.ready_s == 'READY'):
               print 'Warning: Not Going to Initialize Telescope'

          self.log_write('initialize: %s' % self.success)
          return self.success

   
     #########################################################################
     # stow_tel(self, az, elev, dome, wait=0: Returns success or fail
     #######################################################################

     def stow_tel(self, hourang=stow_ha, dec=stow_dec, \
      dome=stow_dome, wait=1):

          print 'Stowing Telescope ... ',

          self.telescope_motion_s = 'MOVING'
          reply = self.comreply_send('TCS STOW %.2f %.2f %.2f' % (hourang, \
           dec, dome), short=0, wait=wait, timeout=240, source='TCS', \
           asynccommand='STOW') 

          if (self.ASYNC == 0) and (re.search(self.done,reply)):
               self.telescope_motion_s = 'IN_POSITION'
               self.telescope_hourang_s = hourang 
               self.telescope_dec_s = dec 
               self.dome_azimuth_s = dome
               self.log_write('stow_tel: %s: %.2f %.2f %.2f' % (self.success, \
                hourang, dec, dome))
               print 'Done'
               return self.success
          elif (self.ASYNC==1) and (wait==0) and (re.search(self.ok,reply)):
               print 'Okay'
               return self.success
          elif (self.ASYNC==1) and (wait==1):
               if (not re.search(self.fail, reply)) and \
                self.check_stow(hourang, dec, dome):
                    print 'Done'
                    return self.success
               else:
                    print 'Error Stowing Telescope'
                    self.log_write('stow_tel: %s: %s' % (self.fail, reply))
                    return '%s: TCS STOW: %s' % (self.fail, reply)
          else:
               print 'Error Stowing Telescope'
               self.telescope_motion_s = 'UNKNOWN'
               self.telescope_hourang_s = -1.0
               self.telescope_dec_s = -1.0
               self.dome_azimuth_s = -1.0
               self.log_write('stow_tel: %s: %s' % (self.fail, reply))
               return '%s: TCS STOW: %s' % (self.fail, reply)

     #####################################################################
     # lamp_on(self): Returns success or fail
     #####################################################################

     def lamp_on(self):

          print 'Turning Lamp On ... ',

          reply = self.comreply_send('TCS LAMPON',short=0,\
           wait=1,timeout=30,source='TCS',asynccommand='LAMPON')

          if (self.ASYNC==0) and (re.search(self.done,reply)):
               self.lamp_s = 'ON'
               print 'Done'
               return self.success
          elif (self.ASYNC == 1):
               if (self.lamp_s == 'ON'): 
                    print 'Done'
                    self.log_write('lamp_on: %s' % self.success)
                    return self.success
               else:
                    print 'Error Turning Lamp On'
                    self.log_write('lamp_on: %s: %s' % (self.fail, reply))
                    return '%s: %s' % (self.fail, reply)
          else:
               print 'Error Turning Lamp On'
               self.lamp_s = 'UNKNOWN'
               self.log_write('lamp_on: %s: %s' % (self.fail, reply))
               return '%s: %s' % (self.fail, reply)

     #####################################################################
     # lamp_off(self): Returns success or fail
     #####################################################################
                                                                                
     def lamp_off(self):
                                                                                
          print 'Turning Lamp Off ... ',

          reply = self.comreply_send('TCS LAMPOFF', short=1,\
           wait=1)

          if re.search(self.done, reply): 
               print 'Done'
	       self.lamp_s = 'OFF'
               self.log_write('lamp_off: %s' % self.success)
               return self.success
          else:
               print 'Error Turning Off Lamp'
	       self.lamp_s = 'UNKNOWN'
               self.log_write('lamp_off: %s: %s' % (self.fail, reply))
               return '%s: %s' % (self.fail, reply) 

     ######################################################################
     # move_filter(self, filter, wait=0): Returns success or fail
     ######################################################################

     def move_filter(self, filter, wait=0):

          print 'Moving Filter ... ',

          self.filter_motion_s = 'MOVING'
          reply = self.comreply_send('FLT MOVE %s' % self.filterlist[filter], \
           short=0, wait=wait, timeout=45, source='FLT') 

          if (self.ASYNC == 0) and (re.search(self.done, reply)):
               self.filter_motion_s = 'STATIONARY'
               self.filter_s = self.filterlist[filter]
               print 'Done'
               self.log_write('move_filter: %s: %s' % (self.success, filter))
               return self.success
          elif (self.ASYNC==1) and (wait==0) and (re.search(self.ok,reply)):
               print 'Okay'
               self.log_write('move_filter: %s: %s' % (self.success, filter))
               return self.success
          elif (self.ASYNC == 1) and (wait == 1):
               if (self.filter_motion_s == 'STATIONARY') and (self.filter_s \
                ==self.filterlist[filter]) and not re.search(self.fail,reply):
                    print 'Done'
                    self.log_write('move_filter: %s: %s' %(self.success,filter))
                    return self.success
               else:
                    print 'Error Moving Filter'
                    self.log_write('move_filter: %s: %s' % (self.fail,reply))
                    return '%s: %s' % (self.fail, reply)
          else:
               print 'Error Moving Filter'
               self.filter_s == 'UNKNOWN'
               self.filter_motion_s == 'UNKNOWN'
               self.log_write('move_filter: %s: %s' % (self.fail, reply))
               return '%s: %s' % (self.fail, reply) 

     ########################################################################
     # take_images(self, number, exptime, name, wait=0): Returns
     # success or failure
     #######################################################################
                    
     def take_images(self, number, exptime, obsname, foc=0, wait=0):

          print 'Taking Images ... '

          if (foc == 0) and (not self.ccd_idle_s == 1):
               reply = self.set_idle('ON')
               if not re.search(self.success,reply):
                    self.log_write('take_images: %s: %s' % (self.fail, reply))
                    return '%s: %s' % (self.fail, reply)
          
          if (foc == 0) and (not self.ccd_autoclear_s == 1):
               reply = self.set_autoclear('ON')
               if not re.search(self.success,reply):
                    self.log_write('take_images: %s: %s' % (self.fail, reply))
                    return '%s: %s' % (self.fail, reply)

          reply = self.set_obs_params(obsname, imagestoread=number, \
                                      exptime=exptime)
          if re.search(self.fail, reply):
               self.log_write('take_images: %s: %s' % (self.fail, reply))
               reply2 = self.handle_error(level=1, errstr='DHE: %s' % reply)
               if re.search(self.fail, reply2):
                    print 'Could Not Recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail

          reply = self.pre_process(number, obsname)
	  if re.search(self.fail, reply):
	       self.ccd_s = 'UNKNOWN'
	       self.log_write('take_images: %s: %s' % (self.fail, reply))
               return '%s: %s' % (self.fail, reply)
	       
          # Check shutter status
          reply = self.direct_cmd('DHE GET autoshutter', short=1)
          self.log_write('Autoshutter Status: %s' % reply)

          tout = number * (exptime + 30.0) + 10.0
          self.ccd_s = 'EXPOSE'
          reply = self.comreply_send('DHE EXPOSE', short=0, wait=wait, \
           timeout=tout, source='DHE', asynccommand='EXPOSE') 

          if (self.ASYNC == 0) and (re.search(self.done,reply)):
               self.ccd_s = 'IDLE'
               print 'Done'
               self.log_write('take_images: %s' % self.success)
               return self.success
          elif (self.ASYNC==1) and (wait==0) and (re.search(self.ok,reply)):
               time.sleep(exptime + 1.0)
               self.ccd_s = 'READOUT'
               print 'Okay'
               self.log_write('take_images: %s: Reading Out' % self.success)
               return self.success
          elif (self.ASYNC==1) and (wait == 1):
               if (self.ccd_s == 'IDLE') and (not re.search(self.fail,reply)): 
                    print 'Done'
                    self.log_write('take_images: %s' % self.success)
                    return self.success
               else:
                    self.log_write('take_images: %s: %s' % (self.fail,reply))
                    reply2 = self.handle_error(level=1,errstr='DHE %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could Not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail
                         
          else:
               self.ccd_s = 'UNKNOWN'
               print 'Error Taking Image'
               self.log_write('take_images: %s: %s' % (self.fail, reply))
               reply2 = self.handle_error(level=1,errstr='DHE %s' % reply)
               if re.search(self.fail,reply2):
                    print 'Could not Recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from error'
                    return self.fail

     ########################################################################
     # open_dome(self, wait=1): Returns success or failure
     ########################################################################

     def open_dome(self, wait=1):

          print 'Opening Dome ... ',

          self.dome_motion_s = 'MOVING'
          reply = self.comreply_send('TCS OPEN', short=0, \
           wait=wait, timeout=300, source='TCS', asynccommand='OPEN')

          if (self.ASYNC == 0) and (re.search(self.done,reply)):
               self.dome_shutters_s = 'OPEN'
               self.dome_motion_s = 'TRACKING'
               print 'Done'
               self.log_write('open_dome: %s' % self.success)
               return self.success
          elif (self.ASYNC==1) and (wait == 0) and (re.search(self.ok,reply)):
               self.log_write('open_dome: %s: Opening Dome' % self.success)
               print 'Okay'
               return self.success
          elif (self.ASYNC == 1) and (wait == 1):
               if (self.dome_shutters_s == 'OPEN') and ((self.dome_motion_s \
                == 'TRACKING') or (self.dome_motion_s == 'GO_TO') or \
                (self.dome_motion_s == 'ANTICIPATE')) and not \
                re.search(self.fail, reply):
                    print 'Done'
                    self.log_write('open_dome: %s' % self.success)
                    return self.success
               else:
                    print 'Error Opening Dome'
                    self.log_write('open_dome: %s: %s' % (self.fail, reply))
                    return '%s: TCS OPEN: %s' % (self.fail, reply)
          else:
               print 'Error Opening Dome'
               self.dome_shutters_s = 'UNKNOWN'
               self.dome_motion_s = 'UNKNOWN'
               self.log_write('open_dome: %s: %s' % (self.fail, reply))
               return '%s: TCS OPEN: %s' % (self.fail, reply)
           
     ########################################################################
     # close_dome(self, wait=1): Returns success or failure
     ########################################################################
                                                                                
     def close_dome(self, wait=1):
                                                                                
          print 'Closing Dome ... ',

          self.dome_motion_s = 'MOVING'
          reply = self.comreply_send('TCS CLOSE', short=0, \
           wait=wait, timeout=300, source='TCS', asynccommand='CLOSE')

          if (self.ASYNC == 0) and (re.search(self.done,reply)):
               self.dome_shutters_s = 'CLOSED'
               self.dome_motion_s = 'OFF'
               print 'Done'
               self.log_write('close_dome: %s' % self.success)
               return self.success
          elif (self.ASYNC==1) and (wait == 0) and (re.search(self.ok,reply)):
               self.log_write('close_dome: %s: Closing Dome' % self.success)
               print 'Okay'
               return self.success
          elif (self.ASYNC == 1) and (wait == 1):
               if (self.dome_shutters_s == 'CLOSED') and not \
                re.search(self.fail, reply):
                    print 'Done'
                    self.log_write('close_dome: %s' % self.success)
                    return self.success
               else:
                    print 'Error Closing Dome'
                    self.log_write('close_dome: %s: %s' % (self.fail, reply))
                    return '%s: TCS CLOSE: %s' % (self.fail, reply) 
          else:
               print 'Error Closing Dome'
               self.dome_shutters_s = 'UNKNOWN'
               self.dome_motion_s = 'UNKNOWN'
               self.log_write('close_dome: %s: %s' % (self.fail, reply))
               return '%s: TCS CLOSE: %s' % (self.fail, reply)

     ########################################################################
     # move_tel(self, name, ra, dec, equinox=2000, mflag=0, rarate=0, 
     # decrate=0)
     #######################################################################

     def move_tel(self, name, ra, dec, equinox=2000.0, mflag=0, rarate=0, \
                  decrate=0, wait=0): 

          print 'Moving Telescope ... ',

          command = 'TCS INPOS \"%s\" %.6f %.6f %.2f' % (name, ra, dec, equinox)
          if not ((rarate == 0) and (decrate == 0)):
               command += ' %i %.2f %.2f' % (mflag, rarate, decrate)
          reply = self.comreply_send(command, short=1)
          if not re.search(self.done, reply):
               print 'Error Inputting Coordinates'
               self.telescope_ra_s = self.telescope_dec_s = \
	        self.telescope_equinox_s = self.telescope_mflag_s = \
		self.telescope_rarate_s = self.telescope_decrate_s = \
		self.point_ra_s = self.point_dec_s = -1.0 
               self.telescope_motion_s = 'UNKNOWN'
               self.log_write('move_tel: %s: TCS INPOS: %s' % (self.fail, \
                reply))
               return '%s: TCS INPOS: %s' % (self.fail, reply)

          self.telescope_motion_s = 'SLEWING'
          reply = self.comreply_send('TCS GOPOS', short=0, \
           wait=wait, timeout=240, source='TCS', asynccommand='GOPOS')

	  # Update where telescope should be pointing
	  self.point_ra_s = ra
	  self.point_dec_s = dec
          self.point_equinox_s = equinox

          if (self.ASYNC == 0) and (re.search(self.done,reply)):
               self.telescope_ra_s = ra
               self.telescope_dec_s = dec
               self.telescope_equinox_s = equinox
               self.telescope_mflag_s = mflag
               self.telescope_rarate_s = rarate
               self.telescope_decrate_s = decrate
               self.telescope_motion_s = 'TRACKING'
               print 'Done'
               self.log_write('move_tel: %s' % self.success)
               return self.success
          elif (self.ASYNC==1) and (wait == 0) and (re.search(self.ok,reply)):
               print 'Okay'
               self.log_write('move_tel: %s: Moving Telescope' % self.success)
               return self.success
          elif (self.ASYNC == 1) and (wait == 1):
               if (not re.search(self.fail, reply)) and self.check_coords():
                    print 'Done'
                    self.log_write('move_tel: %s' % self.success)
                    return self.success
               else:
                    print 'Error Moving Telescope'
                    self.log_write('move_tel: %s: %s' % (self.fail, reply))
                    return '%s: TCS GOPOS: %s' % (self.fail, reply)
          else:
               print 'Error Moving Telescope'
	       self.telescope_ra_s = self.telescope_dec_s = \
	        self.telescope_equinox_s = self.telescope_mflag_s = \
		self.telescope_rarate_s = self.telescope_decrate_s = \
		self.point_ra_s = self.point_dec_s = -1.0
               self.telescope_motion_s = 'UNKNOWN'
               self.log_write('move_tel: %s: TCS GOPOS: %s' % (self.fail, \
                reply))
               return '%s: TCS GOPOS: %s' % (self.fail, reply)

     ########################################################################
     # move_focus(position, wait=0):
     ########################################################################

     def move_focus(self, position, wait=0):

          print 'Moving focus ... ',

          self.focus_motion_s = 'MOVING'
          reply = self.comreply_send('TCS GOFOCUS %.2f' % position, short=0, \
           wait=wait, timeout=45, source='TCS', asynccommand='GOFOCUS') 

          if (self.ASYNC == 0) and (re.search(self.done, reply)):
               self.focus_position_s = position
               self.focus_motion_s = 'STATIONARY'
               print 'Done'
               self.log_write('move_focus: %s' % self.success)
               return self.success
          elif (self.ASYNC==1) and (wait == 0) and (re.search(self.ok,reply)):
               self.log_write('move_focus: %s: Moving Focus' % self.success)
               print 'Okay'
               return self.success
          elif (self.ASYNC==1) and (wait==1):
               if (not self.focus_position_s == -1.0) and \
                (abs(self.focus_position_s - position) < self.focus_error) and \
                (self.focus_motion_s == 'STATIONARY') and \
                (not re.search(self.fail, reply)):
                    print 'Done'
                    self.log_write('move_focus: %s' % self.success)
                    return self.success
               else:
                    print 'Error Moving Focus'
                    self.log_write('move_focus: %s: %s' % (self.fail, reply))
                    return '%s: TCS GOFOCUS: %s' % (self.fail, reply)
          else:
               print 'Error Moving Focus'
               self.focus_position_s = -1.0
               self.focus_motion_s = 'UNKNOWN'
               self.log_write('move_focus: %s: %s' % (self.fail, reply))
               return '%s: TCS GOFOCUS: %s' % (self.fail, reply)
               
     #######################################################################
     # adj_focus(deltafoc, wait=0):
     #######################################################################

     def adj_focus(self, deltafoc, wait=0):

          print 'Adjusting Focus Position ... ',

          self.focus_motion_s = 'MOVING'

          old_position = self.focus_position_s
          reply = self.comreply_send('TCS ADJFOCUS %.2f' % deltafoc, \
           short=0, wait=wait, timeout=45, source='TCS', \
           asynccommand='ADJFOCUS')
          
          if (self.ASYNC == 0) and (re.search(self.done,reply)):
               self.focus_position_s = old_position + deltafoc
               self.focus_motion_s = 'STATIONARY'
               print 'Done'
               self.log_write('adj_focus: %s' % self.success)
               return self.success
          elif (self.ASYNC==1) and (wait == 0) and (re.search(self.ok,reply)):
               print 'Okay'
               self.log_write('adj_focus: %s: Moving Focus' % self.success)
               return self.success
          elif (self.ASYNC==1) and (wait==1):
               if (not self.focus_position_s == -1.0) and \
                (abs(self.focus_position_s - (old_position + deltafoc)) < \
                self.focus_error) and (self.focus_motion_s == 'STATIONARY') \
                and (not re.search(self.fail, reply)):
                    print 'Done'
                    self.log_write('adj_focus: %s' % self.success)
                    return self.success
               else:
                    print 'Error Adjusting Focus'
                    self.log_write('adj_focus: %s: %s' % (self.fail, reply))
                    return '%s: TCS ADJFOCUS: %s' % (self.fail, reply)
          else:
               print 'Error Adjusing Focus'
               self.focus_position_s = -1.0
               self.focus_motion_s = 'UNKNOWN'
               self.log_write('adj_focus: %s: %s' % (self.fail, reply))
               return '%s: TCS ADJFOCUS: %s' % (self.fail, reply)

     ########################################################################
     # check_coords(name, ra, dec, epoch=2000.0, mflag=0, rarate=0, 
     #  decrate=0): Check to see if telescope is in proper
     #  position.  Returns True or False
     #######################################################################

     def check_coords(self):

          try:
	       target = astrocoords(self.point_ra_s*15.0, self.point_dec_s)
               targeteph = ephem.readdb("%s,f|S,%s,%s,15.0,%.2f" % ('None', \
                target.sxg()[0], target.sxg()[1], self.point_equinox_s))
               targeteph.compute(self.palomar)
     
               telescope = astrocoords(self.telescope_ra_s*15.0,\
                self.telescope_dec_s)
               teleph = ephem.readdb("%s,f|S,%s,%s,15.0,%.2f" % ('None', \
                telescope.sxg()[0], telescope.sxg()[1], 2000.0))
               teleph.compute(self.palomar)

               # Now calculate offset
               offset = ephem.separation(targeteph, teleph)
               diff = float(offset) * 180.0 / math.pi * 3600.0

          except (ArithmeticError,LookupError,ValueError,AttributeError):
               print 'Error Verifying Coordinates'
               self.log_write('check_coords: %s: Error Verifying Coordinates' \
                % self.fail)
               return False

          if (self.telescope_motion_s == 'TRACKING') and \
	   (diff < self.telescope_error):
               self.log_write('check_coords: True')
               return True
          else:
               self.log_write('check_coords: False')
               return False

     ########################################################################
     # check_stow(self)
     #######################################################################

     def check_stow(self, ha=stow_ha, dec=stow_dec, \
      dome=stow_dome):

          try:
               target = astrocoords(ha * 15.0, dec)
               object = astrocoords(self.telescope_hourang_s * 15.0, \
                self.telescope_dec_s) 
               offset = target.diff(object,arcmin=1,degree=0)
               diff = sqrt( power(offset[0],2) + power(offset[1],2) )
          except (ArithmeticError,LookupError,ValueError,AttributeError):
               print 'Error Verifying Coordinates'
               self.log_write('check_stow: %s: Error Verifying Coordinates' \
                % self.fail)
               return False

          if ((self.telescope_motion_s == 'IN_POSITION') or \
           (self.telescope_motion_s == 'STOPPED')) and \
           (abs(self.dome_azimuth_s - dome) < self.dome_error) and \
           (diff < self.stow_error):
               self.log_write('check_stow: True')
               return True
          else:
               self.log_write('check_stow: False')
               return False

     ########################################################################
     # get_focus(self, filter)
     ########################################################################

     def get_focus(self, filter, airmass=-1, btubtemp=-1):

          self.read_foctable()
          if not (self.filterlist.has_key(filter)):
               print 'No Filter %s in Wheel' % filter
               self.log_write('get_focus: %s: No Filter %s in Wheel' % \
                (self.fail, filter))
               return -1.0
          
          try:
               if (self.filterlist[filter] == self.focus[0]):
                    focus = self.focus[1]
               else:
                    focus = self.focus[1] - self.offset[self.focus[0]] + \
                     self.offset[self.filterlist[filter]]

               if (btubtemp == -1):
                    btubtemp = self.bot_tube_temp_s
               if (airmass == -1):
                    airmass = self.telescope_airmass_s

               if not (btubtemp == -1.0):
                    tempdiff = btubtemp - self.focus[2]
                    tempoff = -0.037 * tempdiff
                    focus += tempoff
               if not (airmass == -1.0):
                    amdiff = airmass - self.focus[3]
                    # Changed 8/28/07 by SBC w/ new focus data
                    #amoff = -0.033 * amdiff
                    amoff = -0.156 * amdiff
                    focus += amoff

               return focus 
               
          except (ArithmeticError,LookupError,AttributeError,TypeError):
               print 'Error Determing Focus Position'
               self.log_write('get_focus: %s: Error Determing Best Focus \
                Position' % self.fail)
               return -1.0

     #######################################################################
     # read_obsnum(self)
     #######################################################################

     def read_obsnum(self):

          obsnum = 100
          command = os.popen('ssh p60@%s \'ls %s%s/*.fits\'' % \
           (self.aview_ip, self.baugiroot, self.homedir), 'r', -1)
          dir = command.readlines()
          for i in range(len(dir)):
               try:
                    name = dir[i].split('/')[4].split('.')[0]
                    curnum = int(name[len(name)-4:])
                    if (curnum > obsnum):
                         obsnum = curnum
               except (LookupError,AttributeError,TypeError,ValueError):
                    self.log_write('Bad File: %s' % dir[i])

          return (obsnum + 1)
          
     #######################################################################
     # get_obsnum(self)
     #######################################################################

     def get_obsnum(self):

          reply = self.comreply_send('DHE GET imagenumber', short=1) 
          if (re.search(self.fail, reply)):
               return -1
          try:
               obsnum, null = reply.split('\n')
               if re.match('\d{1,4}', obsnum):
                    return int(obsnum)
               else:
                    return -1
          except (LookupError,AttributeError,TypeError,ValueError):
               return -1

     ########################################################################
     # pre_process(self, number)
     #######################################################################

     def pre_process(self, number, obsname):

          self.get_moondeg()

          obsnum = self.get_obsnum()
          if (obsnum == -1):
               self.log_write('pre_process: %s: Error getting obsnum' % \
                self.fail)
               return self.fail

          for i in range(number):
               filename = obsname 
               if (filename[len(filename)-1].isdigit()):
                    filename += '_'
               filename += '%04d' % (obsnum + i)
               self.write_keys(filename)

          return self.success

     #######################################################################
     # get_moondeg(self)
     #######################################################################

     def get_moondeg(self):
 
          if (self.telescope_ra_s == -1.0) or (self.telescope_dec_s == -1.0) \
           or (self.telescope_equinox_s == -1.0):
               self.moondeg = -1.0
               return

          try:
               target = astrocoords(self.telescope_ra_s, self.telescope_dec_s,\
                self.telescope_equinox_s)

               targeteph = ephem.readdb("%s,f|S,%s,%s,15.0,%.2f" % ('None', \
                target.sxg()[0], target.sxg()[1], self.telescope_equinox_s))
               self.palomar.date = ephem.now()
               self.moon.compute(self.palomar)
               targeteph.compute(self.palomar)
               self.moondeg_s = 180.0 - (ephem.separation(self.moon,targeteph)\
                / 0.0174532925199433)
               return
        
          except (ArithmeticError,AttributeError,LookupError,ValueError):
               self.log_write('get_moondeg: %s: Error Calculating moondeg'  %\
                self.fail)
               return

     #######################################################################
     # set_binning(self, xbin=1, ybin=1)
     ######################################################################

     def set_binning(self, xbin=1, ybin=1):

          print 'Setting CCD Binning ... ',

          reply = self.comreply_send('DHE SET binning %d %d' % (xbin, ybin), \
           short=1)
          if (re.search(self.done, reply)):
               print 'Done'
	       self.ccd_bin_s = [xbin, ybin]
               return self.success
          else:
               self.ccd_s = 'UNKNOWN'
	       self.ccd_bin_s = [-1, -1]
               print 'Error Setting Binning'
               return '%s: %s' % (self.fail, reply)

     ########################################################################
     # set_roi(self, x1=1, y1=1, xl=2048, yl=2048)
     #######################################################################

     def set_roi(self, x1=1, y1=1, xl=2048, yl=2048):

          print 'Setting CCD ROI ... ',

          reply = self.comreply_send('DHE SET roi %d %d %d %d' % (x1, y1, \
           xl, yl), short=1)
          if (re.search(self.done, reply)):
               print 'Done'
	       self.ccd_roi_s = [x1, y1, xl, yl]
               return self.success
          else:
               self.ccd_s = 'UNKNOWN'
	       self.ccd_roi_s = [-1, -1, -1, -1]
               print 'Error Setting ROI'
               return '%s: %s' % (self.fail, reply)
     
     ######################################################################
     # set_shutter(self, command)
     ######################################################################

     def set_shutter(self, command):

          print 'Setting CCD Shutter Status ' + command + ' ... ',

          reply = self.comreply_send('DHE DO shutter %s' % command, \
           short=1)
          if (re.search(self.done, reply)):
               print 'Done'
               return self.success
          else:
               self.ccd_s = 'UNKNOWN'
               print 'Error Setting Shutter Status'
               return '%s: %s' % (self.fail, reply)

     ######################################################################
     # set_sampling(self, command)
     ######################################################################

     def set_sampling(self, command):

          print 'Setting CCD Sampling to %s' % command

          reply = self.comreply_send('DHE SET sampling %s' % command, \
           short=1)
          if re.search(self.done, reply):
               return self.success
          else:
               return self.fail

     ######################################################################
     # set_idle(self, command)
     ######################################################################

     def set_idle(self, command):

          print 'Setting CCD Idle Mode ' + command + ' ... ',

          reply = self.comreply_send('DHE SET idle %s' % command, \
           short=1)
          if (re.search(self.done, reply)):
               print 'Done'
	       if re.search('ON', command):
	            self.ccd_idle_s = 1
	       elif re.search('OFF', command):
	            self.ccd_idle_s = 0
               return self.success
          else:
               self.ccd_s = 'UNKNOWN'
               print 'Error Setting Idle Mode'
               self.ccd_idle_s = -1
               return '%s: %s' % (self.fail, reply)

     #####################################################################
     # set_autoclear(self, command)
     ######################################################################

     def set_autoclear(self, command):
      
          print 'Setting CCD Autoclear Mode ' + command + ' ... ',

          reply = self.comreply_send('DHE SET autoclear %s' % command, \
           short=1)
          if (re.search(self.done, reply)):
               print 'Done'
	       if re.search('ON', command):
	            self.ccd_autoclear_s = 1
	       elif re.search('OFF', command):
	            self.ccd_autoclear_s = 0
               return self.success
          else:
               self.ccd_s = 'UNKNOWN'
	       self.ccd_autoclear_s = -1
               print 'Error Setting Autoclear'
               return '%s: %s' % (self.fail,  reply)

     ########################################################################
     # parallel_shift(self)
     ########################################################################

     def parallel_shift(self):

          print 'Performing Parallel Shift ... ',

          reply = self.comreply_send('DHE MEMORY manualcommand tim none ' + \
           'none none none none PSH', short=1)

          if re.search(self.error, reply):
               print 'Error with Parallel Shift'
               self.ccd_s = 'UNKNOWN'
               self.log_write('parallel_shift: ERROR: %s' % reply)
               return '%s: %s' % (self.fail, reply)
          else:
               print 'Done'
               return self.success

     ######################################################################
     # find_hip(self, ra, dec)
     #######################################################################

     def find_hip(self, ra, dec):

          print 'Searching for Nearby Hipparcos Star'

          searchcmd = 'find_Hipparcos %f %f %shipparcos/' % (ra * 15.0, dec,\
           self.rootdir)
          cmd = os.popen(searchcmd, 'r', -1)
          reply = cmd.readlines()

          ####################
	  # Error Checking ??
          ####################

          if (len(reply) == 0):
               return -1
          else:
               reply = reply[0].rstrip().split()
          if not (len(reply) == 3):
               return -1
          else:
               distance = reply[0]
               [ra] = self.set_type([reply[1]], FloatType)
               [dec] = self.set_type([reply[2]], FloatType)
               if (ra == -1.0) or (dec == -1.0):
                    return -1
               else:
                    return [ra / 15.0, dec]

     #######################################################################
     # write_keys(self, filename)
     #######################################################################

     def write_keys(self, filename):

          os.chdir(self.rootdir + self.homedir + '/raw')
          
          try:
               keysfile = open('%s%s/raw/%s.xkeys' % (self.rootdir, \
                self.homedir, filename), 'w+', -1) 
          except IOError:
               self.log_write('write_keys: %s: Error Opening File %s' % \
                (self.fail, filename))

          keysfile.write('IAXIS     2\n')
          keysfile.write('IAXIS1     2048\n')
          keysfile.write('IAXIS2     2048\n')
          keysfile.write('INSTRMNT     CCDcam\n')

          keysfile.write('AUTOTRIG     %s\n' % self.autotrig_s)
          keysfile.write('TRIGSOUR     %s\n' % self.trigsour_s)
          keysfile.write('IMGTYPE     %s\n' % self.imgtype_s)
          keysfile.write('P60PRNM     %s\n' % self.p60prnm_s)
          keysfile.write('P60PRID     %s\n' % self.p60prid_s)
          keysfile.write('P60PRPI     %s\n' % self.p60prpi_s)
          keysfile.write('P60PRTM     %d\n' % self.p60prtm_s)
          keysfile.write('MOONDEG     %.2f\n' % self.moondeg_s)
          keysfile.write('MSCSIZE     %s\n' % self.mscsize_s)
          keysfile.write('MSCPOSN     %s\n' % self.mscposn_s)
          keysfile.write('OBJTYPE     %s\n' % self.objtype_s)
          keysfile.write('BINROI     %s\n' % self.binroi_s)
          keysfile.write('OBJECT     %s\n' % self.object_s)
          keysfile.write('FILTER     %s\n' % self.rfilterlist[self.filter_s])
          keysfile.write('EQUINOX     %s\n' % self.telescope_equinox_s)
          keysfile.write('TARGRA     %.5f\n' % (self.object_ra_s*15.0))
          keysfile.write('TARGDEC     %.5f\n' % self.object_dec_s)

          if not (self.generic_keys_s == ''):
               keysfile.write(self.generic_keys_s)

          if (self.imgtype_s == 'BIAS'):
               keysfile.write('DARK     TRUE\n')
          else:
               keysfile.write('DARK     FALSE\n')

          if (self.imgtype_s == 'FOCUS'):
               keysfile.write('LOOPNUM     %d\n' % self.loopnum)
               keysfile.write('FOCUSSEQ     %dof%d\n' % (self.focuscurnum, \
                self.focusnum))
          elif (self.imgtype_s == 'PHOTOSTD'):
               keysfile.write('STDLOOP     %d\n' % self.stdloop)
          elif (self.imgtype_s == 'SCIENCE'):
               keysfile.write('DEFOCUS     %i\n' % self.defocus_s)
          elif (self.imgtype_s == 'SAOFOCUS'):
               keysfile.write('SAOSHIFT     32\n')
               keysfile.write('RAOFF     90.0\n')
               keysfile.write('DECOFF     -90.0\n')
               keysfile.write('LOOPNUM     %d\n' % self.loopnum)
               keysfile.write('NUMFOC     %d\n' % self.focusnum)
               for i in range(1, self.focusnum+1):
                    keysfile.write('SAOFOC%02d     %.2f\n' % (i, \
                     self.lowfocus + (i-1)*self.focusstep))

          os.chdir(self.rootdir)

          return

     ########################################################################
     # ftp_files(self)
     ########################################################################

     def ftp_files(self):

          if (os.path.exists('%s%s/raw' % (self.rootdir, self.homedir))):
               os.chdir('%s%s/raw' % (self.rootdir, self.homedir))
          else:
               print 'No Such Directory: %s%s/raw' %(self.rootdir,self.homedir)
               self.log_write('ftp_files: %s: No Such Directory' % self.fail)
               return

          command = os.popen('ssh p60@%s \'ls %s%s/*.fits\'' % (self.aview_ip,\
	   self.baugiroot, self.homedir), 'r', -1)
          lines = command.readlines()
          for i in range(len(lines)):
               try:
                    file = lines[i].split('\n')[0]
                    if not self.rawdict.has_key(file):
                         os.system('scp p60@%s:%s .' % (self.aview_ip, file)) 
                         self.rawdict[file] = 1
               except (LookupError,ValueError,AttributeError,TypeError):
                    self.log_write('ftp_files: Bad file from baugi')

          os.chdir(self.rootdir)
          return

     #######################################################################
     # email_notify(self, str)
     #######################################################################

     def email_notify(self, str):

          try:
               msg = MIMEText(str)
               msg['Subject'] = self.email_subject
               msg['From'] = self.email_from
               msg['To'] = self.email_to
               s = smtplib.SMTP()
               s.connect()
               s.sendmail(self.email_from, [self.email_to], msg.as_string())
               s.close()
               return
          except email.Errors:
               self.log_write('email_notify: %s: Error Sending Email' %\
                self.fail)
               return

     #######################################################################
     # direct_cmd(self, command, wait=1, timeout=-1, source='',
     #  asynccommand='')
     #######################################################################

     def direct_cmd(self, command, short=0, wait=1, timeout=-1, \
      source='', asynccommand=''):

          reply = self.comreply_send(command, short=short, \
           wait=wait, timeout=timeout, source=source, asynccommand=asynccommand)

          return reply

     #######################################################################
     # take_bias(self, number, xbin=1, ybin=1, x1=1, y1=1, xl=2048, 
     #  yl=2048, binroi='A')
     #######################################################################

     def take_bias(self, number, xbin=1, ybin=1, x1=1, y1=1, xl=2048, \
      yl=2048, binroi='A'):

          print 'Taking Bias Frames'

          # First set Bias keywords
          self.imgtype_s = 'BIAS'
          self.p60prnm_s = 'CALIBRATION'
          self.p60prid_s = 'CALIBRATION'
          self.p60prtm_s = 0
          self.p60prpi_s = 'CALIBRATION'
          self.mscsize_s = '1,1'
          self.mscposn_s = '0,0'
          self.objtype_s = 'BIAS'
          self.object_s = 'Bias'
          self.binroi_s = binroi

          # Remove Processed Bias Images
          if (binroi == 'A'):
               biasroot = 'Bias'
          else:
               biasroot = 'Bias-%s' % binroi

          if os.path.exists('%s%s/proc/%s.fits' % (self.rootdir, \
           self.homedir, biasroot)):
               os.remove('%s%s/proc/%s.fits' % (self.rootdir, \
                self.homedir, biasroot))
          if os.path.exists('%s%s/proc/%s.lis' % (self.rootdir, \
           self.homedir, biasroot)):
               os.remove('%s%s/proc/%s.lis' % (self.rootdir, \
                self.homedir, biasroot))

          # Set Binning/ROI
          reply = self.set_binning(xbin=xbin, ybin=ybin)
          if (re.search(self.fail, reply)):
               self.log_write('take_bias: ' + reply)
               shutil.copyfile('%s/calib/%s.fits' % (self.rootdir, biasroot), \
                '%s%s/proc/%s.fits' % (self.rootdir, self.homedir, biasroot))
               reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
               if re.search(self.fail,reply2):
                    print 'Could not Recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail

          reply = self.set_roi(x1=x1, y1=y1, xl=xl, yl=yl)
          if (re.search(self.fail, reply)):
               self.log_write('take_bias: ' + reply)
               shutil.copyfile('%s/calib/%s.fits' % (self.rootdir, biasroot), \
                '%s%s/proc/%s' % (self.rootdir, self.homedir, biasroot))
               reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
               if re.search(self.fail,reply2):
                    print 'Could not Recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail

          # Take Images
          reply = self.take_images(number, 0.0, self.object_s, \
           wait=1)

          if not (re.search(self.success, reply)):
               self.log_write('take_bias: %s: %s' % (self.fail, reply))
               shutil.copyfile('%s/calib/%s.fits' % (self.rootdir, biasroot), \
                '%s%s/proc/%s.fits' % (self.rootdir, self.homedir, biasroot))
               return reply
          else:
               print 'Done Taking Bias Frames'
               self.log_write('take_bias: %s' % self.success)
               self.ftp_files()
               return self.success

     #######################################################################
     # dome_flat(self, filter, exptime, number, wait=1)
     #######################################################################

     def dome_flat(self, filter, number, exptime):

          print 'Taking Dome Flats'

          self.imgtype_s = 'DOMEFLAT'
          self.p60prid_s = 'CALIBRATION'
          self.p60prnm_s = 'CALIBRATION'
          self.p60prtm_s = 0
          self.p60prpi_s = 'CALIBRATION'
          self.objtype_s = 'DOMEFLAT'
          self.object_s = 'Flat-' + filter
          self.mscsize_s = '1,1'
          self.mscposn_s = '0,0'

          # Stow Telescope
          if not (self.check_stow()):
               reply = self.stow_tel(wait=1)
               if re.search(self.fail, reply):
                    print 'Error Stowing Telescope'
                    self.log_write('dome_flat: %s: %s' % (self.fail, reply))
                    reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
                    if re.search(self.fail,reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Turn on dome lamp
          if not (self.lamp_s == 'ON'):
               reply = self.lamp_on()
               if re.search(self.fail, reply):
                    print 'Error Turning on Flat Lamp'
                    self.log_write('dome_flat: %s: %s' % (self.fail, reply))
                    reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
                    if re.search(self.fail,reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Move filter
          reply = self.move_filter(filter, wait=1)
          if re.search(self.fail, reply):
               print 'Error Moving Filter'
               self.log_write('dome_flat: %s: %s' % (self.fail, reply))
               reply2 = self.handle_error(level=1,errstr='FLT: %s' % reply)
               if re.search(self.fail,reply2):
                    print 'Could not Recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail

          # Remove Old Flat
          if (os.path.exists('%s%s/proc/%s.fits' % (self.rootdir, \
           self.homedir, self.object_s))):
               os.system('rm -f %s%s/proc/%s*.fits' % (self.rootdir, \
                self.homedir, self.object_s))

          # Take Images
          reply = self.take_images(number, exptime, self.object_s, \
           wait=1)
          if not (re.search(self.success, reply)):
               if os.path.exists('%scalib/%s.fits' % (self.rootdir, \
                self.object_s)):
                    shutil.copyfile('%scalib/%s.fits' % (self.rootdir, \
                     self.object_s), '%s%s/proc/%s.fits' % (self.rootdir, \
                     self.homedir, self.object_s)) 
               self.log_write('dome_flat: %s: %s' % (self.fail, reply))
               return reply
          else:
               print 'Done Taking Dome Flats'
               self.ftp_files()
               self.log_write('dome_flat: %s' % self.success)
               return self.success

     ########################################################################
     # science_obs(name, ra, dec, epoch, raoff=0, decoff=0, filter, exptime
     #######################################################################

     def science_obs(self, name, filter, exptime, ra, dec, equinox=2000.0, \
      raoff=0, decoff=0, xbin=1, ybin=1, x1=1, y1=1, xl=2048, yl=2048, \
      wait=0):

          print 'Executing Science Observation'

          # Get Full Status Update
          self.get_full_status()
          self.write_full_status()
          self.imgtype_s = 'SCIENCE'

          if (not self.remote_close_s == 'OKAY') or \
           (not self.weather_s == 'OKAY') or (not self.ready_s == 'READY'):
               print 'Telescope Not Ready to Observe'
               self.log_write('science_obs: %s: Telescope Not Ready to \
                Observe' % self.fail)
               return self.fail

          # Spawn reduction if necessary
          if (self.PROC == 1) and (self.p60redux_s == 0):
               self.spawn_redux()

          # Open Dome
          if not ((self.dome_shutters_s == 'OPEN') and ((self.dome_motion_s \
           == 'TRACKING') or (self.dome_motion_s == 'ANTICIPATE') or \
           (self.dome_motion_s == 'GO_TO'))): 
               reply = self.open_dome(wait=1)
               if re.search(self.fail, reply):
                    self.log_write('science_obs: %s: %s' % (self.fail,reply))
                    reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail
          
          try:
       	       # Prepare to move telescope: First get desired coordinates
               point = astrocoords(ra*15.0,dec,equinox=equinox)
	       point.shift([raoff,decoff],arcsec=1,arcmin=0)
	       pointeph = ephem.readdb("%s,f|S,%s,%s,15.0,%.2f" % ('None', \
	        point.sxg()[0], point.sxg()[1], equinox))
	       self.palomar.date = ephem.now()
	       pointeph.compute(self.palomar)
	       airmass = self.tnow.airmass(self.palomar,pointeph)

               # Now Telescope Coordinates
	       telescope = astrocoords(self.telescope_ra_s*15.0,\
	        self.telescope_dec_s)
               teleph = ephem.readdb("%s,f|S,%s,%s,15.0,%.2f" % ('None', \
                telescope.sxg()[0], telescope.sxg()[1], 2000.0))
               teleph.compute(self.palomar)

	       # calculate offset
               offsetrad = ephem.separation(pointeph, teleph)
               offseta = float(offsetrad) * 180.0 / math.pi * 3600.0

          except (TypeError,ValueError,LookupError,AttributeError):
               print 'Error with Coordinates: %f %f' % (ra, dec) 
               self.log_write('science_obs: %s: Error with coordinates: ' \
                '%f %f' % (self.fail, ra, dec))
               return self.fail

          # Move the Telescope 
	  if (not name == self.object_s) or (abs(offseta) > \
           self.telescope_error): 
	       reply = self.move_tel(name, point.radeg()/15.0,\
	                             point.dcdeg(),equinox=equinox,wait=wait)
               if re.search(self.paramrange, reply):
                    print 'Target inaccessible'
                    return self.paramrange
               elif not re.search(self.success, reply):
                    self.log_write('science_obs: %s: %s' % (self.fail,reply))
                    reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Update Object Information
	  self.object_s = name
	  self.object_ra_s = ra
	  self.object_dec_s = dec

          # Move Filter
          if not (self.filterlist.has_key(filter)):
               print 'No Filter %s in Wheel' % filter
               self.log_write('science_obs: %s: %s: %s' % (self.fail, \
                self.nofilt, filter))
               return self.nofilt
          elif not ((self.filter_s == self.filterlist[filter]) and \
           (self.filter_motion_s == 'STATIONARY')):
               reply = self.move_filter(filter, wait=wait)
               if not re.search(self.success, reply):
                    self.log_write('science_obs: ' + reply)
                    reply2 = self.handle_error(level=1,errstr='FLT: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Move Focus
          focus = self.get_focus(filter, airmass=airmass)
          if (focus == -1.0):
               self.log_write('science_obs: %s: Error Determining Best Focus'\
                % self.fail)
               return self.fail
          if (self.defocus_s == 1):
               focus += self.defocus
          if (self.focus_position_s == -1.0) or (abs(self.focus_position_s\
           - focus) > self.focus_error) or (not self.focus_motion_s == \
           'STATIONARY'): 
               reply = self.move_focus(focus, wait=wait)
               if not re.search(self.success, reply):
                    self.log_write('science_obs: %s: %s' % (self.fail,reply))
                    reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail
          
          # Wait for Telescope Motion to Stop
	  while (self.telescope_motion_s == 'SLEWING') or \
	   (self.telescope_motion_s == 'MOVING'):
               reply = self.async_recv(source='TCS',command='GOPOS')

               if (re.search(self.fail, reply)):
                    print 'No Valid Reply from Move'
                    self.log_write('science_obs: %s: No Valid Return from \
                     GOPOS' % self.fail)
                    reply2 = self.handle_error(level=1, errstr='TCS: %s' \
                     % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail
               elif not self.check_coords():
	            print 'Telescope Move Unsuccessful'
		    self.log_write('science_obs: %s: Telescope Move\
		     Unsuccessful' % self.fail)
		    return self.fail

          # Wait for Filter motion to stop
          while (not self.filter_s == self.filterlist[filter]) or \
           (not self.filter_motion_s == 'STATIONARY'):
               reply = self.async_recv(source='FLT', command='MOVE')
               if re.search(self.fail, reply):
                    print 'No Valid Reply from filter wheel'
                    self.log_write('science_obs: %s: No Valid reply from \
                     filter wheel' % self.fail)
                    reply2 = self.handle_error(level=1,errstr='FLT: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Wait for Focus Motion to Stop
          while (not self.focus_motion_s == 'STATIONARY') or \
           (self.focus_position_s == -1.0) or (abs(self.focus_position_s - \
           focus) > self.focus_error):
               reply = self.async_recv(source='TCS',command='GOFOCUS')
               if re.search(self.fail, reply):
                    print 'No Valid Reply from Focus'
                    self.log_write('science_obs: %s: No Valid reply from \
                     focus' % self.fail)
                    reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # If CCD status unknown, try to deal intelligently with this
          if (self.ccd_s == 'UNKNOWN'):
               print 'CCD Controller in Unknown State, Trying to Reinit'
               reply = self.reset_controller()
               if re.search(self.fail, reply):
                    self.log_write('science_obs: %s: Could not reinit \
                     controller' % self.fail)
                    reply2 = self.handle_error(level=1,errstr='DHE:%s' % reply) 
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Wait for CCD to finish reading out
          while not (self.ccd_s == 'IDLE'):
               reply = self.async_recv(source='DHE',command='EXPOSE')
               if re.search(self.fail,reply):
                    print 'No Valid Reply from CCD'
                    self.log_write('science_obs: %s: No Vaild reply from CCD'\
                     % self.fail)
                    reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Grab any old images
          self.ftp_files()

          # Set Binning and ROI
	  if not ([xbin, ybin] == self.ccd_bin_s):
               reply = self.set_binning(xbin, ybin)
               if (re.search(self.fail, reply)):
                    self.log_write('science_obs: %s: %s' % (self.fail,reply))
                    reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

	  if not ([x1, y1, xl, yl] == self.ccd_roi_s):
               reply = self.set_roi(x1, y1, xl, yl)
               if (re.search(self.fail, reply)):
                    self.log_write('science_obs: %s: %s' % (self.fail,reply))
                    reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Take exposure
          reply = self.take_images(1, exptime, name, wait=wait)
          if not (re.search(self.success,reply)):
               self.log_write('science_obs: %s: %s' % (self.fail,reply))
               return reply 
          else:
               print 'Done with Science Observation'
               self.log_write('science_obs: %s' % self.success)
               return self.success

     #########################################################################
     # afternoon(self, biaslist, flatlist)
     #########################################################################

     def afternoon(self, biaslist, flatlist):

          print 'Executing Afternoon Calibrations'

          self.get_full_status()
          self.write_full_status()

	  # Make sure dome closed
	  if not (self.dome_shutters_s == 'CLOSED'):
	       reply = self.close_dome(wait=1)
               if (re.search(self.fail, reply)):
	            self.log_write('afternoon: %s: %s' % (self.fail, reply)) 
		    reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
		    if re.search(self.fail,reply2):
		         print 'Could not recover from Error'
			 return self.goodnight
	            else:
	                 print 'Recovered from Error'
			 return self.fail

	  # Take Bias Frames
          for entry in biaslist:
               reply = self.take_bias(number=entry[7], xbin=entry[0], \
                ybin=entry[1], x1=entry[2], y1=entry[3], xl=entry[4], \
                yl=entry[5], binroi=entry[6])
               if not (re.search(self.success, reply)):
                    self.log_write('afternoon: %s: %s' % (self.fail, reply))
                    return reply

          # Reset binning/roi settings to default
          reply = self.set_binning()
          if (re.search(self.fail, reply)):
               self.log_write('afternoon: %s: %s' % (self.fail, reply))
               reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
               if re.search(self.fail, reply2):
                    print 'Could not Recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error' 
                    return self.fail

          reply = self.set_roi()
          if (re.search(self.fail, reply)):
               self.log_write('afternoon: %s: %s' % (self.fail, reply))
               reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
               if re.search(self.fail, reply2):
                    print 'Could not Recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail

          self.binroi_s = 'A'

          for entry in flatlist:

               # Check to see if filter in wheel 
               if not (self.filterlist.has_key(entry[0])):
                    print 'No Filter %s in Wheel' % entry[0] 
                    self.log_write('afternoon: %s: %s: %s' % (self.fail, \
                     self.nofilt, entry[0]))
                    pass
               else:
                    # Execute Dome Flat
                    reply = self.dome_flat(entry[0],entry[1],entry[2])
                    if not (re.search(self.success, reply)):
                         self.log_write('afternoon: %s: %s' % (self.fail, \
                          reply))
			 if not (self.lamp_s == 'OFF'):
			      reply = self.lamp_off()
                         return reply
            
          # Turn Lamp_off
          if (self.lamp_s == 'ON'):
               reply = self.lamp_off()

          # Begin Processing
          if (self.PROC == 1) and (self.p60redux_s == 0):
               time.sleep(300)
               self.spawn_redux()

          print 'Done with Afternoon Calibrations'
          return self.success

     #########################################################################
     # morning(self, biaslist, flatlist)
     #########################################################################

     def morning(self, biaslist, flatlist):

          print 'Executing Morning Calibrations'

          self.homedir = 'morning'
          self.get_full_status()
          self.write_full_status()

          # Remove files from previous night
          os.system('rm -f %s%s/raw/*' % (self.rootdir, self.homedir))
          os.system('rm -f %s%s/misc/*' % (self.rootdir, self.homedir))
          os.system('rm -f %s%s/proc/20*r.fits' % (self.rootdir, \
           self.homedir))
          os.system('rm -f %s%s/proc/Flat*' % (self.rootdir, self.homedir)) 
          os.system('ssh p60@%s \'rm -f %s%s/*\'' % (self.aview_ip, \
           self.baugiroot, self.homedir))

	  # Make Sure dome is closed
          if not (self.dome_shutters_s == 'CLOSED'):
	       reply = self.close_dome(wait=1)
	       if (re.search(self.fail, reply)):
	            self.log_write('morning: %s: %s' % (self.fail, reply))
	            reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
	            if re.search(self.fail,reply2):
	                 print 'Could not recover from Error'
	                 return self.goodnight
	            else:
	                 print 'Recovered from Error'
	                 return self.fail

          # Take all the images
          for entry in biaslist:
               reply = self.take_bias(number=entry[7], xbin=entry[0], \
                ybin=entry[1], x1=entry[2], y1=entry[3], xl=entry[4], \
                yl=entry[5], binroi=entry[6])
               if not (re.search(self.success, reply)):
                    self.log_write('morning: %s: %s' % (self.fail, reply))
                    return reply

          # Reset binning/roi settings to default
          reply = self.set_binning()
          if (re.search(self.fail, reply)):
               self.log_write('morning: %s: %s' % (self.fail, reply))
               reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
               if re.search(self.fail, reply2):
                    print 'Could not Recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail

          reply = self.set_roi()
          if (re.search(self.fail, reply)):
               self.log_write('morning: %s: %s' % (self.fail, reply))
               reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
               if re.search(self.fail, reply2):
                    print 'Could not Recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail

          self.binroi_s = 'A'

          # Loop over entries in flatlist
          for entry in flatlist:

               # Make Sure Filter is in Wheel
               if not (self.filterlist.has_key(entry[0])):
                    print 'No Filter %s in Wheel' % entry[0]
                    self.log_write('morning: %s: %s: %s' % (self.fail, \
                     self.nofilt, entry[0]))
                    pass
               else:
                    reply = self.dome_flat(entry[0],entry[1],entry[2])
                    if not (re.search(self.success, reply)):
			 if not (self.lamp_s == 'OFF'):
			      reply = self.lamp_off()
                         self.log_write('morning: %s: %s' % (self.fail, reply))
                         pass

          # Turn Lamp_off
          if (self.lamp_s == 'ON'):
               reply = self.lamp_off()

          # Stow Telescope in Safety Position
          if not (self.check_stow(dome=self.stow_dome_safe)):
               reply = self.stow_tel(hourang=self.stow_ha, dec=self.stow_dec, \
                dome=self.stow_dome_safe, wait=1)

          # Look for bad flats
          #print 'Searching for Bad Flats ... ',
          #os.system('mv %s%s/raw/Flat*.fits %s%s/flats' % (self.rootdir, \
           #self.homedir, self.rootdir, self.homedir))
          #os.chdir('%s%s/flats' % (self.rootdir, self.homedir))
          #dir = os.listdir('.')
          #for flat in dir:
               #result = badflat(flat)
               #if not (result == [0]):
                    #os.remove(flat)
                    #root, null = flat.split('.')
                    #os.system('rm %s%s/raw/%s*' % (self.rootdir, self.homedir,\
                     #root))
          #os.system('mv %s%s/flats/*.fits %s%s/raw' % (self.rootdir, \
           #self.homedir, self.rootdir, self.homedir))
          #print 'Done'

          # Apply Keywords
          print 'Applying Keywords ... ',
          os.chdir('%s%s/raw' % (self.rootdir, self.homedir))
          os.system(self.p60keys1_x)
          print 'Done'

          # Process Files
          print 'Processing Files'
          os.chdir('%s%s/proc' % (self.rootdir, self.homedir))
          os.system(self.p60redux1_x)

          os.chdir(self.rootdir)
          print 'Done with Morning Calibrations'
          return self.success

     #########################################################################
     # standard_obs(self, object, obslist, ra, dec, equinox=2000.0, 
     #  wait=0)
     #########################################################################

     def standard_obs(self, name, obslist, ra, dec, equinox=2000.0, \
      wait=0):

          print 'Performing Photometric Standard Observations'

          # Update Status and ftp files
          self.get_full_status()
          self.write_full_status()
          self.stdloop += 1
          self.imgtype_s = 'PHOTOSTD'
          self.p60prnm_s = 'CALIBRATION'
          self.p60prtm_s = 0
          self.p60prpi_s = 'CALIBRATION'
          self.p60prid_s = 'CALIBRATION'
          self.mscsize_s = '1,1'
          self.mscposn_s = '0,0'
          self.objtype_s = 'standard'
          self.binroi_s = 'A'

          if (not self.remote_close_s == 'OKAY') or (not self.weather_s == \
           'OKAY') or (not self.ready_s == 'READY'): 
               print 'Telescope Not Ready to Observe'
               self.log_write('standard_obs: %s: Telescope Not Ready to \
                Observe' % self.fail)
               return self.fail

          # Spawn Reduction if necessary
          if (self.PROC == 1) and (self.p60redux_s == 0):
               self.spawn_redux()

          # Open Dome
          if not ( (self.dome_shutters_s == 'OPEN') and ( (self.dome_motion_s \
           == 'TRACKING') or (self.dome_motion_s == 'ANTICIPATE')  \
           or (self.dome_motion_s == 'GO_TO') ) ):
               reply = self.open_dome(wait=1)
               if (re.search(self.fail, reply)):
                    self.log_write('standard_obs: %s: %s' % (self.fail, \
                     reply))
                    reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
                    if re.search(self.fail,reply2):
                         print 'Could not recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          try:
               # Prepare to move telescope: First get desired coordinates
               point = astrocoords(ra*15.0,dec,equinox=equinox)
               pointeph = ephem.readdb("%s,f|S,%s,%s,15.0,%.2f" % ('None', \
                point.sxg()[0], point.sxg()[1], equinox))
               self.palomar.date = ephem.now()
               pointeph.compute(self.palomar)
               airmass = self.tnow.airmass(self.palomar,pointeph)
	  
               # Now Telescope Coordinates
               telescope = astrocoords(self.telescope_ra_s*15.0,\
                self.telescope_dec_s)
               # calculate offset
               offset = point.diff(telescope,arcsec=1,arcmin=0,degree=0)

          except (TypeError,ValueError,LookupError,AttributeError):
               print 'Error with Coordinates: %f %f ' % (ra, dec) 
               self.log_write('set_coords: %s: Error with coordinates ' \
                '%f %f' % (self.fail, ra, dec))
               return self.fail

          # Move the Telescope 
          if (not name == self.object_s) or (abs(offset[0]) \
            > self.telescope_error) or (abs(offset[1]) > self.telescope_error):
               reply = self.move_tel(name, point.radeg()/15.0,\
                point.dcdeg(),equinox=equinox,wait=wait)
               if re.search(self.paramrange, reply):
                    print 'Target inaccessible'
                    return self.paramrange
               elif not re.search(self.success, reply):
                    self.log_write('standard_obs: %s: %s' % (self.fail,reply))
                    reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Update Object Information
          self.object_s = name
          self.object_ra_s = ra
          self.object_dec_s = dec

          # Loop over each entry in the list
          for entry in obslist:

               filter = entry[0]
               exptime = entry[1]

               # Move the filter
               if not (self.filterlist.has_key(filter)):
                    print 'No Filter %s in Wheel' % filter
                    self.log_write('standard_obs: %s: %s: %s' % (self.fail, \
                     self.nofilt, filter))
                    pass
               elif not ( (self.filter_s == self.filterlist[filter]) and \
                (self.filter_motion_s == 'STATIONARY') ):
                    reply = self.move_filter(filter, wait=wait)
                    if (re.search(self.fail, reply)):
                         self.log_write('standard_obs: %s: %s' % (self.fail,\
                          reply))
                         reply2 = self.handle_error(level=1,errstr='FLT: %s' \
                          % reply)
                         if re.search(self.fail,reply2):
                              print 'Could not recover from Error'
                              return self.goodnight
                         else:
                              print 'Recovered from Error'
                              return self.fail

               # Move the focus
               focus = self.get_focus(filter, airmass=airmass)
               if (focus == -1.0):
                    self.log_write('standard_obs: %s: Bad Focus Value' % \
                     self.fail)
                    return self.fail
               elif (self.focus_position_s == -1.0) or (abs(focus - \
                self.focus_position_s) > self.focus_error) or (not \
                self.focus_motion_s == 'STATIONARY'):
                    reply = self.move_focus(focus, wait=wait)
                    if (re.search(self.fail, reply)):
                         self.log_write('standard_obs: %s: %s' % (self.fail,\
                          reply))
                         reply2 = self.handle_error(level=1,errstr='TCS: %s' \
                          % reply)
                         if re.search(self.fail,reply2):
                              print 'Could not recover from Error'
                              return self.goodnight
                         else:
                              print 'Recovered from Error'
                              return self.fail

               # Wait for Telescope Motion to Stop
               while (self.telescope_motion_s == 'SLEWING') or \
                (self.telescope_motion_s == 'MOVING'):
                    reply = self.async_recv(source='TCS',command='GOPOS')
                    if (re.search(self.fail, reply)):
                         print 'No Valid Reply from Move'
                         self.log_write('science_obs: %s: No Valid Return from \
                          GOPOS' % self.fail)
                         reply2 = self.handle_error(level=1, errstr='TCS: %s' \
                          % reply)
                         if re.search(self.fail, reply2):
                              print 'Could not Recover from Error'
                              return self.goodnight
                         else:
                              print 'Recovered from Error'
                              return self.fail
                    elif not self.check_coords():
                         print 'Telescope Move Unsuccessful'
                         self.log_write('science_obs: %s: Telescope Move\
                          Unsuccessful' % self.fail)
                         return self.fail

               # Wait for Focus to Stop
               while (self.focus_position_s == -1.0) or (abs(focus - \
                self.focus_position_s) > self.focus_error) or (not \
                self.focus_motion_s == 'STATIONARY'):
                    reply = self.async_recv(source='TCS',command='GOFOCUS')
                    if re.search(self.fail,reply):
                         print 'No Valid reply from Focus'
                         self.log_write('standard_obs: %s: No Vaild reply \
                          from Focus' % self.fail)
                         reply2 = self.handle_error(level=1,errstr='TCS: %s' \
                          % reply)
                         if re.search(self.fail,reply2):
                              print 'Could not recover from Error'
                              return self.goodnight
                         else:
                              print 'Recovered from Error'
                              return self.fail

               # Wait for filter to stop
               while (not self.filter_s == self.filterlist[filter]) or \
                (not self.filter_motion_s == 'STATIONARY'): 
                    reply = self.async_recv(source='FLT',command='MOVE')
                    if re.search(self.fail, reply):
                         print 'No Valid reply from Filter'
                         self.log_write('standard_obs: %s: No valid reply \
                          from filter wheel' % self.fail)
                         reply2 = self.handle_error(level=1,errstr='FLT: %s' \
                          % reply)
                         if re.search(self.fail,reply2):
                              print 'Could not recover from Error'
                              return self.goodnight
                         else:
                              print 'Recovered from Error'
                              return self.fail

               # Wait for CCD to finish exposing
               while not (self.ccd_s == 'IDLE'):
                    reply = self.async_recv(source='DHE',command='EXPOSE')
                    if re.search(self.fail, reply):
                         print 'No Valid reply from CCD'
                         self.log_write('standard_obs: %s: No valid reply \
                          from CCD' % self.fail)
                         reply2 = self.handle_error(level=1,errstr='DHE: %s' \
                          % reply)
                         if re.search(self.fail,reply2):
                              print 'Could not recover from Error'
                              return self.goodnight
                         else:
                              print 'Recovered from Error'
                              return self.fail

               # Grab any old images
               self.ftp_files()

               # Set Binning and ROI
	       if not (self.ccd_bin_s == [1,1]):
                    reply = self.set_binning()
                    if (re.search(self.fail, reply)):
                         self.log_write('standard_obs: %s: %s' % \
			  (self.fail,reply))
                         reply2 = self.handle_error(level=1,\
			  errstr='DHE: %s' % reply)
                         if re.search(self.fail,reply2):
                              print 'Could not recover from Error'
                              return self.goodnight
                         else:
                              print 'Recovered from Error'
                              return self.fail

	       if not (self.ccd_roi_s == [1, 1, 2048, 2048]):
                    reply = self.set_roi()
                    if (re.search(self.fail, reply)):
                         self.log_write('standard_obs: %s: %s' % \
			  (self.fail,reply))
                         reply2 = self.handle_error(level=1,\
			  errstr='DHE: %s' % reply)
                         if re.search(self.fail,reply2):
                              print 'Could not recover from Error'
                              return self.goodnight
                         else:
                              print 'Recovered from Error'
                              return self.fail

               # Take images
               reply = self.take_images(1, exptime, name, wait=wait) 
               if not (re.search(self.success, reply)):
                    self.log_write('standard_obs: %s: %s' % (self.fail, \
                     reply))
                    return reply

          print 'Done with Photometric Standards'
          self.log_write('standard_obs: %s' % self.success)
          return self.success

     ########################################################################
     # do_fast_focus(self, name, filter, median, number, step, exptime, ra,
     #  dec, equinox=2000.0)
     ########################################################################

     def do_fast_focus(self, name, filter, median, number, step, exptime, ra, \
      dec, equinox=2000.0):

          print 'Performing Fast Focus Loop'

          # Refresh status info
          self.get_full_status()
          self.write_full_status()

          # Check ready
          if (not self.remote_close_s == 'OKAY') or (not self.weather_s == \
           'OKAY') or (not self.ready_s == 'READY'):
               print 'Telescope Not Ready to Observe'
               self.log_write('do_fast_focus: %s: Telescope Not Ready to \
                Observe' % self.fail)
               return self.fail

          # Spawn reduction if necessary
          if (self.p60fastfoc_s == 0):
               self.spawn_fastfoc()

          # Check Dome
          if not ( (self.dome_shutters_s == 'OPEN') and ((self.dome_motion_s \
           == 'TRACKING') or (self.dome_motion_s == 'ANTICIPATE') or \
           (self.dome_motion_s == 'GO_TO') ) ):
               reply = self.open_dome(wait=1)
               if (re.search(self.fail, reply)):
                    self.log_write('do_fast_focus: %s: %s' % (self.fail,reply))
                    reply2 = self.handle_error(level=1,errstr='TCS: %s' %reply)
                    if re.search(self.fail, reply2):
                         print 'Could not Recover from error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Find Nearest Hipparcos Star
          reply = self.find_hip(ra=ra, dec=dec)
          if (reply == -1):
               self.log_write('do_fast_focus: %s: Could not find nearby HIP \
                star' % self.fail)
               return self.fail
          else:
               hipra = reply[0] 
               hipdec = reply[1] 

          # Go to nearest Hipparcos Star (including offset)
	  hipcoo = astrocoords(hipra*15.0,hipdec)
	  hipcoo.shift([90.0, -90.0],arcsec=1,arcmin=0)
          hipeph = ephem.readdb("%s,f|S,%s,%s,15.0,%.2f" % ('None', \
           hipcoo.sxg()[0], hipcoo.sxg()[1], equinox))
          self.palomar.date = ephem.now()
          hipeph.compute(self.palomar)
          airmass = self.tnow.airmass(self.palomar,hipeph)

          reply = self.move_tel('Hipparcos', hipcoo.radeg()/15.0, \
	   hipcoo.dcdeg(), wait=0)
          if (re.search(self.fail, reply)):
               self.log_write('do_fast_focus: %s: %s' % (self.fail,reply))
               reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
               if re.search(self.fail, reply2):
                    print 'Could not Recover from error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail
          else:
	       self.object_s = 'Hipparcos'
	       self.object_ra_s = hipra
	       self.object_dec_s = hipdec

          # Check Filter
          if not (self.filterlist.has_key(filter)):
               print 'No Filter %s in Wheel' % filter
               self.log_write('do_fast_focus: %s: %s: %s' % (self.fail, \
                self.nofilt, filter))
               return self.nofilt

          # Move Filter
          if not ( (self.filter_s == self.filterlist[filter]) and \
           (self.filter_motion_s == 'STATIONARY') ):
               reply = self.move_filter(filter, wait=0)
               if (re.search(self.fail, reply)):
                    self.log_write('do_fast_focus: %s: %s' % (self.fail, \
                     reply))
                    reply2 = self.handle_error(level=1,errstr='FLT: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Move Focus
          if (median == -1.0):
               median = self.get_focus(filter, airmass=airmass)
               if (median == -1.0):
                    print 'Error Finding Best Focus Value'
                    self.log_write('do_fast_focus: %s: Error Finding Focus \
                     Value' % self.fail)
                    return self.fail
          lowfocus = median - ((number - 1) * step / 2.0) 
          #reply = self.move_focus(lowfocus, wait=0)
          reply = self.move_focus(lowfocus - step, wait=0)
          if (re.search(self.fail, reply)):
               self.log_write('do_fast_focus: %s: %s' % (self.fail,reply))
               reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
               if re.search(self.fail, reply2):
                    print 'Could not recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail

          # Wait for Telescope Motion to Stop
          reply = self.async_recv(source='TCS', command='GOPOS')
          if (re.search(self.fail, reply)):
               print 'No Valid Reply from GOPOS'
               self.log_write('do_fast_focus: %s: No Valid Return from \
                GOPOS' % self.fail)
               reply2 = self.handle_error(level=1, errstr='TCS: %s' \
                % reply)
               if re.search(self.fail, reply2):
                    print 'Could not Recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail
          elif not self.check_coords():
	       print 'Error Moving Telescope'
	       self.log_write('do_fast_focus: %s: Error Moving Telescope' \
	        % self.fail)
	       return self.fail

          # Wait for filter to settle
          while not ( (self.filter_s == self.filterlist[filter]) and \
           (self.filter_motion_s == 'STATIONARY') ):
               reply = self.async_recv(source='FLT', command='MOVE')
               if (re.search(self.fail,reply)):
                    print 'No Valid Reply from Filter Wheel'
                    self.log_write('do_fast_focus: %s: No reply from Filter' \
                     % self.fail)
                    reply2 = self.handle_error(level=1,errstr='FLT: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Wait for focus to settle
          while (not self.focus_motion_s == 'STATIONARY') or \
           (self.focus_position_s == -1.0) or (abs((lowfocus - step) - \
           self.focus_position_s) > self.focus_error):
               reply = self.async_recv(source='TCS',command='GOFOCUS')
               if (re.search(self.fail,reply)):
                    print 'No Valid Reply from Focus'
                    self.log_write('do_fast_focus: %s: No reply from Focus' \
                     % self.fail)
                    reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Wait for ccd to return
          while not (self.ccd_s == 'IDLE'):
               reply = self.async_recv(source='DHE',command='EXPOSE')
               if (re.search(self.fail,reply)):
                    print 'No Valid Reply from CCD'
                    self.log_write('do_fast_focus: %s: No reply from CCD' \
                     % self.fail)
                    reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Grab any old files
          self.ftp_files()

          # Set binning / roi
	  if not (self.ccd_bin_s == [1, 1]):
               reply = self.set_binning()
               if (re.search(self.fail, reply)):
                    self.log_write('do_fast_focus: %s: %s' % (self.fail, reply))
                    reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

	  if not (self.ccd_roi_s == [1, 1, 2048, 2048]):
               reply = self.set_roi()
               if (re.search(self.fail, reply)):
                    self.log_write('do_fast_focus: %s: %s' % (self.fail, reply))
                    reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Turn off idle mode and autoclear mode
          reply = self.set_idle('OFF')
          if (re.search(self.fail, reply)):
               self.log_write('do_fast_focus: %s: %s' % (self.fail, reply))
               reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
               if re.search(self.fail, reply2):
                    print 'Could not recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail

          reply = self.set_autoclear('OFF')
          if (re.search(self.fail, reply)):
               self.log_write('do_fast_focus: %s: %s' % (self.fail, reply))
               reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
               if re.search(self.fail, reply2):
                    print 'Could not recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail

          # New Change for bad first focus step 8/17/2007
          # Make initial focus step
          reply = self.adj_focus(step, wait=1)
          if (re.search(self.fail, reply)):
               print 'Failure Adjusting Focus'
               self.log_write('do_fast_focus: %s: %s' % (self.fail, reply))
               reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
               if re.search(self.fail, reply2):
                    print 'Could not recover from Error'
                    return self.goodnight
               else:
                    print 'Recovered from Error'
                    return self.fail

          # Loop over focus values
          for i in range(1, number):

               # Open Shutter
               reply = self.set_shutter('open')
               if (re.search(self.fail, reply)):
                    print 'Failure Opening Shutter'
                    self.log_write('do_fast_focus: %s: %s' % (self.fail, \
                     reply))
                    reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

               time.sleep(exptime)

               # Close Shutter
               reply = self.set_shutter('close')
               if (re.search(self.fail, reply)):
                    print 'Failure Closing Shutter'
                    self.log_write('do_fast_focus: %s: %s' % (self.fail, \
                     reply))
                    reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

               # Parallel Shift
               reply = self.parallel_shift()
               if (re.search(self.fail, reply)):
                    print 'Failure with Parallel Shift'
                    self.log_write('do_fast_focus: %s: %s' % (self.fail, \
                     reply))
                    reply2 = self.handle_error(level=1,errstr='DHE: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

               # Adjust focus
               reply = self.adj_focus(step, wait=1)
               if (re.search(self.fail, reply)):
                    print 'Failure Adjusting Focus'
                    self.log_write('do_fast_focus: %s: %s' % (self.fail, \
                     reply))
                    reply2 = self.handle_error(level=1,errstr='TCS: %s' % reply)
                    if re.search(self.fail, reply2):
                         print 'Could not recover from Error'
                         return self.goodnight
                    else:
                         print 'Recovered from Error'
                         return self.fail

          # Prepare for final image
          self.imgtype_s = 'SAOFOCUS'
          self.p60prnm_s = 'CALIBRATION'
          self.p60prtm_s = 0
          self.p60prpi_s = 'CALIBRATION'
          self.p60prid_s = 'CALIBRATION'
          self.mscsize_s = '1,1'
          self.mscposn_s = '0,0'
          self.objtype_s = 'sao'
          self.lowfocus = lowfocus
          self.focusstep = step
          self.focusnum = number
          self.object_s = 'Hipfocus'
	  self.binroi_s = 'A'
          self.loopnum += 1

          # Take Image
          reply = self.take_images(1, exptime, self.object_s, foc=1, wait=0)
          if not (re.search(self.success, reply)):
               print 'Failure Taking Image'
               self.log_write('do_fast_focus: %s: %s' % (self.fail, reply))
               return self.fail 

          print 'Done Executing Fast Focus Loop'
          return self.success 

     #########################################################################
     # end_night(self, datatran=1, number=1, wait=1)
     #########################################################################

     def end_night(self, datatran=1, number=1, wait=1):

          print 'Executing End of Night Tasks'

          # Reset variables
          self.get_full_status()
          self.write_full_status()

          # Wait for ccd status to become idle
          while not (self.ccd_s == 'IDLE'):
               reply = self.async_recv()
               if re.search(self.fail,reply):
                    self.log_write('wrap_up_night: %s: No Valid Reply from \
                     CCD' % self.fail) 
                    break

          # Ftp any remaining files
          self.ftp_files()

          # Close Dome
          if not (self.dome_shutters_s == 'CLOSED'):
               reply = self.close_dome(wait=1)
               self.log_write('close_dome: %s' % reply)

          # Stow Telescope
          if not (self.check_stow(dome=self.stow_dome_safe)):
               reply = self.stow_tel(dome=self.stow_dome_safe, wait=1)
               self.log_write('stow_tel: %s' % reply)

          # Spawn datatran
          if (datatran == 1):
               self.datatran = os.popen('%s %s%s %d %d >> /dev/null' % \
                (self.datatran_x, self.rootdir, self.homedir, number, wait)) 
               self.datatran_s = 1
          
          return self.success

     ########################################################################
     # clean_up(self, dir, number=1)
     ########################################################################

     def clean_up(self, dir, number=1, wait=0):

          print 'Transferring Data'

          # Hold off for processing to finish?
	  if wait==1:
               time.sleep(3600.0)

          # Extract homedir and rootdir
          try:
               temp = dir.split('/')
               homedir = temp[len(temp)-1]
               rootdir = '/'
               for i in range(1,len(temp)-1):
                    rootdir += '%s/' % temp[i]
          except:
               print 'Improperly Specified Directory %s' % dir
               return 0
 
          # Move Calibration Files around
          for filter in self.filterlist:
               if os.path.exists('%s/proc/Flat-%s.fits' % (dir, filter)):
                    shutil.copyfile('%s/proc/Flat-%s.fits' % (dir, filter), \
                     '%scalib/Flat-%s.fits' % (rootdir, filter))

          # Move Status_Science.txt
          if os.path.exists('%s/proc/Status_Science.txt' % dir):
               shutil.copyfile('%s/proc/Status_Science.txt' % dir, \
                '%ssci_logs/%s-sci.log' % (rootdir, homedir))

          # Clean Up Proc Directory
          os.chdir('%s/proc' % dir)
               
          # Create stop file for PTF data transfer
          if not os.path.exists("PTF%s.stop" % homedir):
               os.system("touch PTF%s.stop" % homedir)
               if (wait==1):     
                    time.sleep(600.0)
          
          dirlist = os.listdir('.')

          for file in dirlist:
               if re.search('b20\d*r.fits',file):
                    print 'Pipeline Processing Error: %s' % file
                    self.log_write('clean_up: Pipeline Processing Error')
                    return 0
               elif re.search('20\d*r.fits', file) or re.search('.cat', file)\
                or re.search('.reg', file) or re.search('.lis', file):
                    os.remove(file)
               elif re.search('mask', file) or re.search('stars', file):
                    os.rename(file, '../misc/%s' % file)
               elif re.search('BPM', file) or re.search('Bias', file) or \
                re.search('Flat', file):
                    os.rename(file, '../calib/%s' % file)
               elif re.search(self.image_status, file):
                    os.rename(file, '../%s' % file)
               elif re.match('^PTF2[0-9]{7}.', file):
                    os.rename(file, '../misc/%s' % file)
               elif not re.search('20\d*p.fits', file):
                    os.remove(file) 

          # Look for Keyword errors
          os.chdir('%s/proc' % dir)
          contents = os.listdir('.')
          for file in contents:
               test = check_proc_keywords(file)
               if (test == 0):
                    print 'Keyword Error: %s' % file
                    return 0
          os.chdir('%s/raw' % dir)
          contents = os.listdir('.')
          for file in contents:
               test = check_raw_keywords(file)
               if (test == 0):
                    print 'Keyword Error: %s' % file 
                    return 0

          # Make sure all files processed
          for file in contents:
               [imgtype] = get_head(file, ['IMGTYPE'])
               if (imgtype == 'SCIENCE'):
                    procfile = file[:14] + 'p.fits'
                    if not os.path.exists('../proc/%s' % procfile):
                         print 'Some Files Not Processed: %s' % file
                         return 0
          os.chdir(dir)

          # Generate webpage with Observation logs
          os.system('/home/p60/bin/summary.csh %s' % homedir)
	  #os.system('/home/60inch/bin/summary.csh %s' % homedir)
          os.chdir(rootdir)
	  os.system('/home/p60/bin/autolog2.pl %s' % homedir)
          #os.system('/home/60inch/bin/autolog2.pl %s' % homedir)
          os.chdir(dir)

          # Create md5 checksum
          md5cmd = 'md5sum calib/* proc/* raw/* > %s.md5' % homedir 
          os.system(md5cmd)

          # create tarball
          os.chdir('../')
          tarcmd =  'tar cMfffffffffffffffffff ' + homedir + '.tar.1 ' + homedir
          tarcmd += '.tar.2 ' + homedir + '.tar.3 ' + homedir + '.tar.4 '
          tarcmd += homedir + '.tar.5 ' + homedir + '.tar.6 ' + homedir + \
             '.tar.7 ' + homedir + '.tar.8 ' + homedir + '.tar.9 ' + homedir \
             + '.tar.10 ' + homedir + '.tar.11 ' + homedir + '.tar.12 ' \
             + homedir + '.tar.13 ' + homedir + '.tar.14 ' \
             + homedir + '.tar.15 ' + homedir + '.tar.16 ' \
             + homedir + '.tar.17 ' \
             + homedir + '.tar.18 ' + homedir + '.tar.19 -L 1000000 ' \
             + homedir + '/calib ' \
             + homedir + '/raw ' + homedir + '/proc ' + homedir + '/' \
             + homedir + '.log ' + homedir + '/' + homedir + '.md5'
          os.system(tarcmd)

          # Rename tar files to signify how many there are
          temp = os.popen('ls ' + homedir + '.tar.* | wc', 'r', -1)
          temp2 = temp.readlines()
          temp3 = temp2[0].split()
          numfiles = int(temp3[0])

          for i in range(numfiles):
               mvcmd = 'mv ' + homedir + '.tar.' + str(i + 1) + ' ' \
                  + homedir + '.' + 'R0' + str(number) + '.tar.' + str(i+1) + \
                  'of' + str(numfiles)
               os.system(mvcmd)

          # gzip tar files
          zipcmd = 'gzip ' + homedir + '*.tar*'
          os.system(zipcmd)

          # Ftp files to IPAC
          for i in range(numfiles):
               ftpcmd = 'ncftpput ftp.astro.caltech.edu /incoming/grb/60inch '\
                 + homedir + '.' + 'R0' + str(number) + '.tar.' + str(i+1) + \
                 'of' + str(numfiles) + '.gz'
               os.system(ftpcmd)

          # Return to /dst directory
          os.chdir(self.rootdir)

          return 0

     ########################################################################
     # parse_afternoon(self, data):
     ########################################################################

     def parse_afternoon(self, data):

          biaslist = []
          flatlist = []

          try:
               lines = data.splitlines()
               for i in range(1, len(lines), 2):
                    if re.search('BIAS', lines[i]):
                         [null, xbin, ybin, x1, y1, xl, yl, number] = \
                          lines[i].split()
                         binroi = lines[i+1].split()[2]
                         [xbin, ybin, x1, y1, xl, yl, number] = self.set_type(\
                          [xbin, ybin, x1, y1, xl, yl, number], IntType) 
                         [binroi] = self.set_type([binroi], StringType)
                         if (xbin == -1) or (ybin == -1) or (x1 == -1) or \
                          (y1 == -1) or (xl == -1) or (yl == -1) or \
                          (number == -1) or (binroi == 'UNKNOWN'):
                              self.log_write('parse_afternoon: Error from \
                               data\n%s' % data)
                              return [ [], [] ]
                         else:
                              biaslist.append((xbin, ybin, x1, y1, xl, yl, \
                               binroi, number))
                    else:
                         [filter, number, exptime] = lines[i].split()
                         [filter] = self.set_type([filter], StringType)
                         [number] = self.set_type([number], IntType)
                         [exptime] = self.set_type([exptime], FloatType)
                         if (filter == 'UNKNOWN') or (number == -1) or \
                          (exptime == -1.0):
                              self.log_write('parse_afternoon: Error from \
                               data\n%s' % data)
                              return [ [], [] ]
                         else:
                              flatlist.append((filter, number, exptime))
               
          except (LookupError,ValueError,TypeError,AttributeError):
               self.log_write('parse_afternoon: Error from data\n%s' % data)
               return [ [], [] ]

          return [biaslist, flatlist]

     ########################################################################
     # parse_morning(self, data):
     ########################################################################

     def parse_morning(self, data):
     
          reply = self.parse_afternoon(data)
          return reply

     ########################################################################
     # parse_standard_obs(self, data):
     ########################################################################

     def parse_standard_obs(self, data):

          obslist = []

          try:
               lines = data.splitlines()
               name = lines[1]
               [name] = self.set_type([name], StringType)
               [ra, dec, equinox, null, null] = lines[2].split()
               [ra, dec, equinox] = self.set_type([ra, dec, equinox], \
                FloatType)
               if (name == 'UNKNOWN') or (ra == -1.0) or (dec == -1.0) or \
                (equinox == -1.0):
                    self.log_write('parse_standard: Error from data\n%s' \
                     % data)
                    return ['', [], -1, -1, -1]

               for i in range(3, len(lines)-1):
                    [filt, exptime] = lines[i].split()
                    [filt] = self.set_type([filt], StringType)
                    [exptime] = self.set_type([exptime], FloatType)
                    if (filt == 'UNKNOWN') or (exptime == -1.0):
                         self.log_write('parse_standard: Error from data\n%s'\
                          % data)
                         return ['', [], -1, -1, -1]
                    else:
                         obslist.append((filt, exptime))
                    
               return [name, obslist, ra, dec, equinox]

          except (LookupError,ValueError,TypeError,AttributeError):
               self.log_write('parse_standard: Error from data\n%s' % data)
               return ['', [], -1, -1, -1]

     ########################################################################
     # parse_do_sao_focus(self, data):
     #########################################################################

     def parse_do_sao_focus(self, data):

          try:
               lines = data.splitlines()
               name = lines[1]
               [ra, dec, equinox, null, null] = lines[2].split()
               [filt, median, number, step, exptime] = lines[3].split()
               [name, filt] = self.set_type([name, filt], StringType)
               [ra, dec, equinox, step, exptime] = self.set_type(\
                [ra, dec, equinox, step, exptime], FloatType)
               [number] = self.set_type([number], IntType)
               if (name == 'UNKNOWN') or (filt == 'UNKNOWN') or (ra == -1.0)\
                or (dec == -1.0) or (equinox == -1.0) or (step == 1.0) or \
                (exptime == -1.0) or (number == -1):
                    self.log_write('parse_do_sao_focus: Error from data\n%s' \
                     % data)
                    return ['', '', -1, -1, -1, -1, -1, -1, -1]

               if (median == 'CURRENT'):
                    median = -1
               else:
                    [median] = self.set_type([median], FloatType)
                    if (median == -1.0):
                         self.log_write('parse_do_sao_focus: Error from data\
                          \n%s' % data)
                         return ['', '', -1, -1, -1, -1, -1, -1, -1]

               return [name, filt, median, number, step, exptime, ra, dec,\
                equinox]

          except (IndexError,ValueError):
               self.log_write('parse_do_sao_focus: Error from data\n%s' % data)
               return ['', '', -1, -1, -1, -1, -1, -1, -1]

     #######################################################################
     # parse_science_obs(self, data):
     #######################################################################

     def parse_science_obs(self, data):

          try:
               lines = data.splitlines()
               name = lines[1]
               [ra, dec, equinox, raoff, decoff] = lines[2].split()
               filt = lines[3]
               exptime = lines[4]
               [xbin, ybin] = lines[5].split()
               [x1, y1, xl, yl] = lines[6].split()
               self.set_keywords(lines[7:])
            
               [name, filt] = self.set_type([name, filt], StringType)
               [ra, dec, equinox, raoff, decoff, exptime] = self.set_type(\
                [ra, dec, equinox, raoff, decoff, exptime], FloatType)
               [xbin, ybin, x1, y1, xl, yl] = self.set_type(\
                [xbin, ybin, x1, y1, xl, yl], IntType)
               
               if (name == 'UNKNOWN') or (filt == 'UNKNOWN') or (ra == -1.0)\
                or (dec == -1.0) or (equinox == -1.0) or (exptime == -1.0) \
                or (xbin == -1) or \
                (ybin == -1) or (x1 == -1) or (y1 == -1) or (xl == -1) or \
                (yl == -1):
                    self.log_write('parse_science_obs: Error from data\n%s' \
                     % data)
                    return ['', '', -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, \
                     -1, -1]
               else:
                    return [name, filt, exptime, ra, dec, equinox, raoff, \
                     decoff, xbin, ybin, x1, y1, xl, yl]

          except(IndexError,ValueError):
               self.log_write('parse_science_obs: Error from data\n%s' % data)
               return ['', '', -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]

     ########################################################################
     # set_keywords(self, lines)
     #######################################################################

     def set_keywords(self, lines):

          self.defocus_s = 0

          try:
               for line in lines:

                    rexp = re.search('(KEYWORD)\s*(\w*)\s*(\".*\")', line)
                    rexp2 = re.search('(KEYWORD)\s*(\w*)\s*(\'.*\')', line)
                    if rexp:
                         [null, keyword, value] = \
                          [rexp.group(1), rexp.group(2), rexp.group(3)]
                    elif rexp2:
                         [null, keyword, value] = \
                          [rexp2.group(1), rexp2.group(2), rexp2.group(3)]
                    else:
                         [null, keyword, value] = line.split()     

                    if re.search('P60PRID',keyword):
                         self.p60prid_s = value
                    elif re.search('P60PRTM',keyword):
                         self.p60prtm_s = int(value)
                    elif re.search('P60PRPI',keyword):
                         self.p60prpi_s = value
                    elif re.search('P60PRNM',keyword):
                         self.p60prnm_s = value
                    elif re.search('MSCSIZE',keyword):
                         self.mscsize_s = value
                    elif re.search('MSCPOSN',keyword):
                         self.mscposn_s = value
                    elif re.search('OBJTYPE',keyword):
                         self.objtype_s = value
                    elif re.search('BINROI',keyword):
                         self.binroi_s = value
                    elif re.search('DEFOCUS',keyword):
                         self.defocus_s = value
                    else:
                         print 'Unknown Keyword: %s' % keyword
                         self.generic_keys_s = '%s     %s\n' % (keyword, value)
               return

          except (LookupError, TypeError):
               print 'Bad Keyword Entry'
               self.p60prid_s = 'UNKNOWN'
               self.p60prtm_s = 0
               self.p60prpi_s = 'UNKNOWN'
               self.p60prnm_s = 'UNKNOWN'
               self.mscsize_s = '1,1'
               self.mscposn_s = '0,0'
               self.binroi_s = 'A'

     #########################################################################
     # def wrap_up_night(self)
     #########################################################################

     def wrap_up_night(self):

          # Wait for ccd status to become idle
          while not (self.ccd_s == 'IDLE'):
               reply = self.async_recv()
               if re.search(self.fail,reply):
                    self.log_write('wrap_up_night: %s: No Valid Reply from \
                     CCD' % self.fail) 
                    break

          # Ftp any remaining files
          self.ftp_files()

          # Close OSS
          self.oss_close()

          # Close Logfile
          self.close_logfile()

          # Close reduction pipeline
          self.close_redux()

          # Move Facility Status
          self.get_homedir()
          os.system('mv %s %s%s/' % (self.facility_status, \
           self.rootdir, self.homedir))

          # Get Rid of raw files in proc directory
          os.system('rm %s%s/proc/20*r*' % (self.rootdir, self.homedir))
          
          # Get Rid of hipfocus files
          os.system('rm %s%s/hipfocus/*' % (self.rootdir, self.homedir))
          #os.system('mv %s%s/hipfocus/* %spfocus/' % 
           #(self.rootdir, self.homedir, self.rootdir))

          return self.success
     
     ########################################################################
     # handle_error(level=level, errstr=errstr):
     ########################################################################

     def handle_error(self, level=1, errstr=null):

          print 'Handling Previous Error'
          print 'Level = %i, Error String = %s' % (level, errstr)
          
          # Make sure we are not still exposing
          if (self.ccd_s == 'EXPOSE') or (self.ccd_s == 'READOUT'):
               self.async_recv(source='DHE', timeout=240.0)

          # TCS errors
          if (level == 1) and re.search('TCS', errstr):
          
               if re.search('BUSY', errstr):
                    if re.search('TELINIT', errstr):
                         reply = self.direct_cmd('TCS STOPINIT', short=1)
                    elif re.search('OPEN',errstr) or re.search('CLOSE',errstr):
                         reply = self.direct_cmd('TCS STOPSHTRS', short=1)
                    elif re.search('GOPOS',errstr) or re.search('STOW',errstr)\
                     or re.search('GOREF',errstr):
                         reply = self.direct_cmd('TCS STOPTEL', short=1)
                    elif re.search('MOVE', errstr):
                         reply = self.direct_cmd('TCS STOPMOVE', short=1)
                    elif re.search('GOFOCUS', errstr) or \
                     re.search('ADJFOCUS', errstr):
                         reply = self.direct_cmd('TCS STOPFOCUS', short=1)
                    else:
                         reply = self.direct_cmd('TCS STOPINIT', short=1)
                         reply = self.direct_cmd('TCS STOPSHTRS', short=1)
                         reply = self.direct_cmd('TCS STOPSHTRS', short=1)
                         reply = self.direct_cmd('TCS STOPMOVE', short=1)
                         reply = self.direct_cmd('TCS STOPFOCUS', short=1)
          
                    if re.search(self.done, reply):
                         return self.fixed
                    else:
                         reply = self.handle_error(level=2)
                         return reply
               elif re.search('NOCONTROL', errstr):
                    reply = self.take_control()
                    if re.search(self.success, reply):
                         return self.fixed
                    else:
                         reply = self.handle_error(level=2)
                         return reply
               elif re.search('NOTREADY', errstr):
                    reply = self.check_ready()
                    if re.search(self.success, reply):
                         return self.fixed
                    elif (self.weather_s == 'NOT_OKAY') or \
                     (self.remote_close_s == 'NOT_OKAY'):
                         return self.ok
                    else:
                         reply = self.handle_error(level=2)
                         return reply
               elif re.search('PARAMRANGE',errstr) and \
                re.search('GOPOS',errstr):
                    return self.ok
               # Unexpected TCS Error
               else:
                    reply = self.handle_error(level=2)
                    return reply 
                    
          # Filter Wheel Errors
          elif (level==1) and (re.search('FLT', errstr)): 
               self.recycle_power(dev='flt')
               reply = self.flt_init(wait=1)
               if not re.search(self.fail, reply):
                    return self.fixed
               else:
                    reply = self.handle_error(level=2)
                    return reply

          # DHE Errors
          elif (level == 1) and re.search('DHE', errstr):
               self.recycle_power(dev='dhe')
               reply = self.reset_controller()
               if not re.search(self.fail, reply):
                    return self.fixed
               else:
                    reply = self.handle_error(level=2)
                    return reply 

          # Level 2 (or unknown Level 1)
          elif (level == 1) or (level == 2):
               time.sleep(10.0)
               reply = self.direct_cmd('TCS STOPFOCUS', short=1)
               reply = self.comreply_send('TCS ADJFOCUS 0.05',short=0,wait=1,\
                source='TCS', asynccommand='ADJFOCUS')
               if not re.search(self.fail, reply):
                    return self.ok
               else:
                    self.arcview_disconnect()
                    reply = self.arcview_connect()
                    if re.search(self.success, reply):
                         return self.fixed
                    else:
                         reply = self.handle_error(level=3)
                         return reply
               #reply = self.handle_error(level=3)
               return reply

          # Level 3
          elif (level == 3):
               self.email_notify('Warning: Restarting ArcVIEW!')
               self.arcview_close()
               reply = self.arcview_spawn()
               if re.search(self.success, reply):
                    reply2 = self.initialize()
                    if re.search(self.success, reply2):
                         return self.fixed
                    else:
                         return self.fail
               else:
                    return self.fail
                   
##############################################################################
 
