#! /usr/local/bin/python

import re

##########################################################################
# Global Variables for Use by the P60 OCS
##########################################################################

# Rootdir (needed for later variables)
rootdir = '/home/p60/'		     # Root directory for ocs operations

# OSS Variables
oss_ip = '198.202.125.198' 	     # IP address of OSS host machine
oss_port = 1326                      # Port for OSS Connection
oss_x = 'oscar_auto.py'              # OSS Executable filename

# Arcview Variables
#aview_ip = '198.202.125.195'	     # IP address of Arcview host machine (baugi)
aview_ip = '198.202.125.202'	     # IP address of Arcview host machine (baugi2)
comreply_port = 1320		     # Port for Command/Reply Channel
async_port = 1325		     # Port for Asynchronous Channel
comreply_short_tout = 15	     # Timeout for comreply short commands
comreply_long_tout = 300	     # timeout for comreply long commands
async_tout = 300		     # Default async timeout
longreply = 'DONE'	             # Expected initial reply from long command

# Filter List    		     # Real Filter name - Arcview name mapping
filterlist = {		   	     # *ONLY FILTERS CURRENTLY IN WHEEL*
	'H6602' : 'H6602',
	'B' : 'B',
	'V' : 'V',
	'R' : 'R',
	'I' : 'I',
	'g' : 'g',
	'H6614' : 'H6614',
	'ip' : 'ip',
	'zp' : 'zp',
	'Ha20' : 'Ha20',
	'Ho20' : 'Ho20',
	'r' : 'Cr',
        'upr' : 'upr',
        'gpr' : 'gpr',
        'rpr' : 'rpr',
        'ipr' : 'ipr',
        'zpr' : 'zpr'}
rfilterlist = {
	'H6602' : 'H6602',
	'B' : 'B',
	'V' : 'V',
	'R' : 'R',
	'I' : 'I',
	'g' : 'g',
	'H6614' : 'H6614',
	'ip' : 'ip',
	'zp' : 'zp',
	'Ha20' : 'Ha20',
        'Ho20' : 'Ho20',
        'Cr' : 'r',
        'upr' : 'upr',
        'gpr' : 'gpr',
        'rpr' : 'rpr',
        'ipr' : 'ipr',
        'zpr' : 'zpr'}

home_filter = 'B'

# String Variables
done = 'DONE'
success = 'SUCCESS'
ok = 'OK'
fail = 'FAILURE'
error = 'ERROR'
NOT = 'NOT'
null = 'NULL'
nofilt = 'NOFILT'
redo = 'REDO'
wetness = 'WETNESS'
paramrange = 'PARAMRANGE'
illegal = 'ILLEGAL'
fixed = 'FIXED'
goodnight = 'GOODNIGHT'
rawtxt = re.compile(r'.fits')

# Email Notification Variables
email_subject = 'P60 Shutting Down!'
email_from = 'p60@mimir.palomar.caltech.edu'
email_to = '5105088220@txt.att.net'

# Standard Dome Stow Position
stow_ha = (11.0 / 3.0)
stow_dec =  50.0
stow_dome = 40.0
stow_dome_safe = 220.0
stow_error = 10.0                        # Like wise for stow pointing (arcmin)
dome_error = 1.0

# Executables
p60keys_x = 'p60keywords.py >> /dev/null'	# Keywords 
p60keys1_x = 'p60keywords_0.py >> /dev/null'  # Keywords (1x thru)
p60redux_x = 'iqp60 -w 15:00 \"20*r.fits\" >> /dev/null' # Reduction
p60redux1_x = 'iqp60 -o \"20*r.fits\" >> /dev/null'     # Redux (1x thru)
p60saofoc_x = 'saofocus -w 14:00 \"20*r.fits\" >> /dev/null' # Saofocus
p60fastfoc_x = 'fastfocus -w 14:00 \"20*r.fits\" >> /dev/null' # Fast Focus
p60copyfoc_x = 'p60copyfoc >> /dev/null'	# Copy focus files
datatran_x = 'clean_up.py'			# Data transfer
p60transient_x = 'ptf2transient.py'

# Misc Global Variables
p60manual = rootdir + 'p60manual.txt'  	  # Manual Commands file
foctable = rootdir + 'foctable/foctable.txt'	  # Current Focus Table
baugiroot = '/disk/images/'			# Rootname for files on baugi
targetlist = rootdir + 'targets/p60-targets.lis'	# Target list
standardlist = rootdir + 'targets/p60-standard.lis'	# Standards list
telescope_error = 1.0			# Allowable offset from nominal
                                        # telescope position (arcsec)
focus_error = 0.02			# Allow focus error (mm)
offset_error = 1.0			# Allow Offset Error (arcsec)
sao_bright = 8.5			# Brightest SAO star we'll obseve
sao_dim = 10.0				# Dimmest SAO star we'll observe

image_status = 'Status_Science.txt'
facility_status = rootdir + 'Status_Facility.txt'
defocus = 0.1				# Offset for defocussed images
Qsize = 7				# Size of Pointing Update queue
min_point = 10.0			# Min Pointing Correction to apply (")
max_point = 180.0			# Max pointing correction to apply (")

