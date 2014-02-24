#! /usr/bin/env python

from ocs import *
from oscar_session import *

###########################################################################

def main():

     # Full shutdown if system experience exception
     shutdown = 'y'

     # Parse Command line
     try:
          opts, args = getopt.getopt( sys.argv[1:], 
                       "haiovqdbtnwrf:s:", 
                       ["help", "auto", "interactive", "observe", "verbose", \
                        "quiet", "day", "blocking", "test", "noproc", \
                        'wait', "rawscores", "fwait", "swait"])
     except getopt.GetoptError:
          usage()
          sys.exit(2)

     ######################
     
     # Defaults
     AUTO = 1
     SKIPSTART = 1
     VERBOSE = 0
     DAY = 0
     ASYNC = 1
     PROC = 1
     TEST = 0
     WAIT = 0
     RAW = 0
     FWAIT = 4.0
     SWAIT = 1000.0

     # Options parsing
     for opt, val in opts:
          # Usage info only
          if opt in ('-h', '--help'):
               usage()
               sys.exit(1)
          # Run in automated mode
          elif opt in ('-a', '--auto'):
               AUTO = 1
          # Run in interactive mode
          elif opt in ('-i', '--interactive'):
               AUTO = 0
          # Skip early night standard field
          elif opt in ('-o', '--observe'):
               SKIPSTART = 1
          # Verbose operation
          elif opt in ('-v', '--verbose'):
               VERBOSE = 1
          # Quiet Operation
          elif opt in ('-q', '--quiet'):
               VERBOSE = 0
          # Day Mode
          elif opt in ('-d', '--day'):
               DAY = 1
          # Blocking mode
          elif opt in ('-b', '--blocking'):
               ASYNC = 0
          # Test mode
          elif opt in ('-n', '--noproc'):
               PROC = 0
          # Wait mode
          elif opt in ('-w', '--wait'):
               WAIT = 1
          # Use raw scores from scheduler
          elif opt in ('-r', '--rawscores'):
               RAW = 1
          # Non-default focus wait time
          elif opt in ('-f', '--fwait'):
               FWAIT = float(val)
          # Non-default standard wait time
          elif opt in ('-s', '--swait'):
               SWAIT = float(val)
          else:
               sys.stderr.write("Unmatched option %s\n" % opt)

     ########################

     # Define the ocs object
     ocs = ocs_ooriented(TEST=TEST, PROC=PROC, ASYNC=ASYNC, DAY=DAY, \
      VERBOSE=VERBOSE, SKIPSTART=SKIPSTART, AUTO=1)
     oss = oscar_session()
  
     try:

          try: 

               # Now time for the big daily loop
               while (1):
  
                    # Wait until midnight (UT) to begin (if necessary)
                    while ((time.gmtime()[3] > 15) or (time.gmtime()[3] < 1)):
                         print 'Waiting until afternoon'
                         time.sleep(600)

                    # Take a break for the electronics (if possible)
                    #if (oss.sunalt() > 5):
                    #     ocs.take_break() 

                    # Initialize the system
                    reply = ocs.initialize(aclose=0)
                    if re.search(ocs.fail, reply):
                         print 'Error Initializing System'
                         sys.exit(3)

                    # NEW: Start Scheduler After 5 PM Local Time
                    while (time.localtime()[3] > 17) and (time.localtime()[3] \
 	             < 18):
                         print 'Waiting to read target list'
                         time.sleep(300)

                    # Spawn Scheduler
                    if (ocs.DAY == 0) and (ocs.TEST == 0):
                         reply = ocs.oss_start(raw=RAW,fwait=FWAIT,swait=SWAIT)
                         if re.search(ocs.fail, reply):
                              print 'Error Spawning Scheduler'
                              sys.exit(4)

                    # Get First Request from OSS
                    data = ocs.osssock_recv()
 
                    # Main Nightly Loop
                    while (1):
 
                         # Record Request from oss
                         print '%s\n%s' % (data, time.asctime())
                         ocs.log_write(data)
 
                         # Afternoon
                         if ( data[:9] == 'AFTERNOON' ):
                              [biaslist, flatlist] = ocs.parse_afternoon(data)
                              reply = ocs.afternoon(biaslist, flatlist)

                         # Morning
                         elif ( data[:7] == 'MORNING' ):
                              [biaslist, flatlist] = ocs.parse_morning(data)
                              reply = ocs.morning(biaslist, flatlist)

                         # Photo Standard
                         elif ( data[:8] == 'STANDARD' ):
                              [name, obslist, ra, dec, equinox] = \
                               ocs.parse_standard_obs(data)
                              reply = ocs.standard_obs(name, obslist, ra, dec,\
                               equinox=equinox, wait=WAIT)

                         # Sao Focus Loop
                         elif ( data[:10] == 'DOSAOFOCUS' ):
                              [name, filter, median, number, step, exptime, \
                               ra,dec,equinox] = ocs.parse_do_sao_focus(data)
                              reply = ocs.do_fast_focus(name, filter, median, \
                               number, step, exptime, ra, dec, equinox=equinox)

                         # Science Observations
                         elif ( data[:7] == 'SCIENCE' ):
                              [name, filter, exptime, ra, dec, equinox, raoff,\
                               decoff, xbin, ybin, x1, y1, xl, yl] = \
                               ocs.parse_science_obs(data)
                              reply = ocs.science_obs(name, filter, exptime, \
                               ra, dec, equinox=equinox, raoff=raoff, \
                               decoff=decoff, xbin=xbin, ybin=ybin, x1=x1, \
                               y1=y1, xl=xl, yl=yl, wait=WAIT) 

                         # End Night
                         elif ( data[:8] == 'ENDNIGHT' ):
                              reply = ocs.end_night(datatran=1,number=1,wait=1)
  
                         # Open Dome
                         elif ( data[:8] == 'OPENDOME' ):
                              reply = ocs.open_dome(wait=1)

                         # Check ready
                         elif ( data[:10] == 'CHECKREADY' ):
                              reply = ocs.check_ready()

                         # Disconnect OSS
                         elif ( data[:4] == 'DONE' ):
                              reply = ocs.wrap_up_night()

                         # Otherwise Unrecognized
                         else:
                              print 'Unknown Routine: %s' % data
                              reply = ocs.fail

                         # If ocs tells us to shutdown
                         if re.search(ocs.goodnight, reply):
                              print 'Shutting Down System'
                              sys.exit(4)

                         # If we are done for the night break
                         elif ( data[:4] == 'DONE' ):
                              print 'Goodnight'
                              break
                         # Otherwise send reply and grab a new command
                         else:
                              ocs.osssock_send(reply)
			      # Check to see if we have a manual command
			      if os.path.exists(p60manual):
			           print 'Accessing Manual Command'
				   infile = open(p60manual)
				   lines = infile.readlines()
				   data = ''
				   for line in lines:
				        data += line
			           os.remove(p60manual)
			      # If not, get next command from scheduler
			      else:
                                   print 'Waiting for data from OSS'
                                   data = ocs.osssock_recv()

                    # Going to sleep
                    print 'Going to Sleep'
                    time.sleep(14400)

          except (KeyboardInterrupt), e:
               print 'Keyboard Interrupt Received'
               shutdown = raw_input('Full Hardware Shutdown (y/n)?:')
               while (shutdown != 'y') and (shutdown != 'n'):
                    shutdown = raw_input('Full Hardware Shutdown (y/n)?:')

          except Exception, e:
               print 'Unknown Exception: %s' % e
               print 'OCS Closing Unexpectedly'

     finally:

          if not (ocs.osssock_s == 0):
               ocs.osssock_send('GOODNIGHT')

          ocs.wrap_up_night()

          if (shutdown == 'n'):
               sys.exit(5)
          else:
               ocs.email_notify('P60 Closing Unexpectedly')
               ocs.end_night(datatran=0)

#############################################################################

def usage():

     (xdir, xname) = os.path.split(sys.argv[0])
     print "Usage: " + xname +  " [-h] [--auto] [--interactive] [--observe]" \
      + " [--verbose] [--quiet] [--day] [--blocking] [--noproc] [--test]" + \
      " [--wait]"

#############################################################################

if __name__=="__main__":

     main()


               


          

            
     
