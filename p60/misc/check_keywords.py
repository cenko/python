#! /usr/local/bin/python

from iqutils import *

def check_proc_keywords(file):

     keywords = ['CTYPE1', 'CTYPE2', 'CRPIX1', 'CRPIX2', \
               'CRVAL1', 'CRVAL2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', \
               'EQUINOX', 'P60PRID', 'P60PRNM', \
               'P60PRPI', 'P60PRTM', 'OBJECT', 'OBJTYPE', 'OBSDATE', \
               'OBSTIME', 'OBSMJD', 'OBSLST', 'OBSCVDT', 'OBSCVTM', 'DARK', \
               'IMGTYPE', 'FILTER', 'EXPTIME', 'ALTITUDE', 'AZIMUTH', \
               'AIRMASS', 'HOURANG', 'MOONDEG', 'AUTOTRIG', \
               'PROCESSD', 'PROCVER', 'SEEING', 'SKYBKG', 'NSTARS', \
               'EXTINCT', 'PROCPROB']
     results = get_head(file, keywords)

     fail = 0
     for i in range(len(results)):
          if (results[i] == ''):
               print 'Keyword Error: No Keyword ' + keywords[i]
               fail = 1
          else:
               pass

     if (fail == 0):
          return 1
     else:
          return 0

def check_raw_keywords(file):

     keywords = ['CTYPE1', 'CTYPE2', 'CRPIX1', 'CRPIX2', \
               'CRVAL1', 'CRVAL2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', \
               'EQUINOX', 'P60PRID', 'P60PRNM', \
               'P60PRPI', 'P60PRTM', 'OBJECT', 'OBJTYPE', 'OBSDATE', \
               'OBSTIME', 'OBSMJD', 'OBSLST', 'OBSCVDT', 'OBSCVTM', 'DARK', \
               'IMGTYPE', 'FILTER', 'EXPTIME', 'ALTITUDE', 'AZIMUTH', \
               'AIRMASS', 'HOURANG', 'MOONDEG', 'AUTOTRIG']
     results = get_head(file, keywords)

     fail = 0
     for i in range(len(results)):
          if (results[i] == ''):
               print 'Keyword Error: No Keyword ' + keywords[i]
               fail = 1
          else:
               pass

     if (fail == 0):
          return 1
     else:
          return 0


