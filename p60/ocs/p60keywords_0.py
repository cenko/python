#! /usr/local/bin/python

import sys, math
import time, mx.DateTime
import os
import re
import string
import shutil
from iqutils import *

# Global Variables
pixscale = 0.3787
pa = -0.74 * math.pi / 180.0 

###############################################################################

def get_obsdate(file):

     if check_head(file, 'UTSHUT'):
          utshut = get_head(file, 'UTSHUT')
     else:
          print 'No UTSHUT Keyword in %s!' % file
          return -1

     try:
          temp = utshut.split('-')
          obsdate = temp[0] + temp[1]
          temp2 = temp[2].split('T')
          obsdate += temp2[0]
          temp3 = temp2[1].split(':')
          obsdate += temp3[0] + temp3[1]
          if (temp3[2][1] == '.'):
               obsdate += '0' + temp3[2][0]
          else:
               obsdate += temp3[2][:2]
          return obsdate
     except LookupError:
          print 'Error Parsing UTSHUT from %s' % file
          return -1

##############################################################################

def write_new_keys(oldroot):

     global cd1_1, cd1_2, cd2_1, cd2_2, pixscale

     # Open keysfile
     keysfile = oldroot + '.xkeys'
     keys = open(keysfile, 'a', -1)

     # Translate UTSHUT to mx.DateTime object
     if check_head(oldroot + '.fits', 'UTSHUT'):
          utshut = get_head(oldroot + '.fits', 'UTSHUT')
     else:
          print 'No UTSHUT Keyword in %s!' % oldroot + '.fits' 
          return -1
     year = int(utshut[:4])
     month = int(utshut[5:7])
     day = int(utshut[8:10])
     hour = int(utshut[11:13])
     minute = int(utshut[14:16])
     second = float(utshut[17:])
     uttime = mx.DateTime.DateTime(year, month, day, hour, minute, second)
    
     # Write UT Time keywords
     keys.write('OBSDATE     ' + str(uttime.date) + '\n')
     keys.write('OBSTIME     ')
     if (uttime.hour < 10.0):
          keys.write('0' + str(uttime.hour) + ':')
     else:
          keys.write(str(uttime.hour) + ':')
     if (uttime.minute < 10.0):
          keys.write('0' + str(uttime.minute) + ':')
     else:
          keys.write(str(uttime.minute) + ':')
     if (uttime.second < 10.0):
          keys.write('0%.3f' % uttime.second + '\n')
     else:
          keys.write('%.3f' % uttime.second + '\n')
     keys.write('OBSJD       ' + str(uttime.jdn) + '\n')
     keys.write('OBSMJD      ' + str(uttime.mjd) + '\n')
 
     # Translate to Local Time
     ltime = uttime.localtime()
     keys.write('OBSCVDT     ' + str(ltime.date) + '\n')
     keys.write('OBSCVTM     ')
     if (ltime.hour < 10.0):
          keys.write('0' + str(ltime.hour) + ':')
     else:
          keys.write(str(ltime.hour) + ':')
     if (ltime.minute < 10.0):
          keys.write('0' + str(ltime.minute) + ':')
     else:
          keys.write(str(ltime.minute) + ':')
     if (ltime.second < 10.0):
          keys.write('0%.3f' % ltime.second + '\n')
     else:
          keys.write('%.3f' % ltime.second + '\n')

     # Write original name keyword
     obsdate = get_obsdate(oldroot + '.fits')
     keys.write('ORIGNAME    ' + obsdate + 'r.fits' + '\n')

     # Get relevant keywords from file
     (telra, teldec, ccdsum) = get_head(oldroot \
        + '.fits',  ['TELRA','TELDEC','CCDSUM'])
     detsec = get_head(oldroot + '.fits', 'DETSEC', extn=1)
     if (telra == '') or (teldec == '') or (detsec == ''):
          return -1          

     # Simple ones first
     keys.write('CTYPE1     RA---TAN\n')
     keys.write('CTYPE2     DEC--TAN\n')

     # Get binning info
     if (ccdsum == ''):
          ccdsum = '1 1'
          keys.write('CCDSUM      1 1\n')
     binning = ccdsum.split()[0]

     # Pixel Scale Info
     pixs = pixscale * int(binning)
     keys.write('PIXSCALE     ' + str(pixs) + '\n')
     keys.write('PIXSCAL1     ' + str(pixs) + '\n') 
     keys.write('PIXSCAL2     ' + str(pixs) + '\n')

     # CD Matrix (from add_wcs)
     cd1_1 = -pixs * math.sin(pa) / 3600
     cd2_2 = pixs * math.sin(pa) / 3600
     cd1_2 = pixs * math.cos(pa) / 3600
     cd2_1 = pixs * math.cos(pa) / 3600
     keys.write('CD1_1     ' + str(cd1_1) + '\n')
     keys.write('CD1_2     ' + str(cd1_2) + '\n')
     keys.write('CD2_1     ' + str(cd2_1) + '\n')
     keys.write('CD2_2     ' + str(cd2_2) + '\n')

     # Now to get ROI info
     xrange, yrange = detsec.split(',')
     x1, x2  = xrange[1:len(xrange)].split(':')
     y1, y2  = yrange[:len(yrange)-1].split(':')

     # Write IAXIS keywords
     iaxis1 = int(x2) - int(x1) + 1
     iaxis2 = (int(y2) - int(y1) + 1) * 2
     update_head(oldroot + '.fits', 'IAXIS', '2', extn=0)
     update_head(oldroot + '.fits', 'IAXIS1', str(iaxis1), extn=0)
     update_head(oldroot + '.fits', 'IAXIS2', str(iaxis2), extn=0)

     # CRPIX Keywords
     crpix1 = int(int(x1) + (iaxis1 / 2))
     keys.write('CRPIX1     ' + str(crpix1) + '\n')
     update_head(oldroot + '.fits', 'CRPIX2', '1024', extn=0)
     update_head(oldroot + '.fits', 'CRPIX2', '1024', extn=1)
     update_head(oldroot + '.fits', 'CRPIX2', '0', extn=2)

     # RA and DEC Keywords
     target = astrocoords(telra, teldec)
     deltax = crpix1 - 1024
     deltara = deltax * cd1_1
     deltadec = deltax * cd2_1
     target.shift([deltara, deltadec])
     targetra, targetdec = target.sxg()
     targetra_deg, targetdec_deg = target.deg()
     keys.write('RA     ' + str(targetra) + '\n')
     keys.write('DEC     ' + str(targetdec) + '\n')
     keys.write('CRVAL1     ' + str(targetra_deg) + '\n')
     keys.write('CRVAL2     ' + str(targetdec_deg) + '\n')

     return 1

############################################################################
     
filedict = {}
rawtxt = re.compile(r'.fits')
#intyear = time.gmtime()[0]
year = re.compile(r'20\d{12}r')

# Look for new fits files and process them
raw = os.listdir('.')
for i in range(len(raw)):

     try:

               print raw[i] 
               if (rawtxt.search(raw[i])):

                    # Skip Bias & Flats
                    if year.match(raw[i]):

                         filedict[raw[i]] = 1
 
                    elif not (filedict.has_key(raw[i])):

                         # Wait to make sure full image has been transferred
                         time.sleep(10)

                         (oldroot, ext) = raw[i].split('.')

                         # Write Additional keywords to xkeys file
                         reply = write_new_keys(oldroot)
                         if (reply == -1):
                              filedict[raw[i]] = 1
                              pass

                         # Add additiona keywords to .fits file 
                         addkeycmd = 'pyaddkeys -f ' + oldroot + '.xkeys ' \
                                     + oldroot + '.fits'
                         os.system(addkeycmd)
                         time.sleep(1)

                         # Gets obsdate
                         obsdate = get_obsdate(raw[i])

                         # For now, copy fits file (later rename)
                         shutil.copy(oldroot + '.fits', obsdate + 'r.fits')
                         shutil.copy(obsdate + 'r.fits', '../proc/' + \
                                     obsdate + 'r.fits')

                         # For now, move old files to ../misc 
                         tempcmd = 'mv ' + oldroot + '* ../misc'
                         os.system(tempcmd)
                         time.sleep(1)
                         
                         # Add file to dictionary
		         filedict[raw[i]] = 1
                         newfile = obsdate + 'r.fits'
                         filedict[newfile] = 1

     except Exception, e:

          print >>sys.stderr, e




