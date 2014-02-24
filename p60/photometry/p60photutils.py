from iqutils import *
from pyraf import iraf
from iqpkg import *
from mx.DateTime import *

iraf.digiphot()
iraf.apphot()
phot = iraf.digiphot.apphot.phot

def p60zeropt(ifile, refstars, filter):

	infile = open(ifile)
	line = infile.readline().rstrip()
	refstars.set_mag(filter)

	while not (line == ''):
     		[naxis1, naxis2] = get_head(line, ['NAXIS1','NAXIS2'])
		refstars.wcs2pix(line)
     		image = line.split('.')
		root = ''
                for i in range(0,len(image)-1):
                    root+='%s.' % image[i]
     		xyfile = open(root + 'xy', 'w')
     		for star in refstars:
			if ( (star.xval > 1) and (star.xval < naxis1) and \
			 (star.yval > 1) and (star.yval < naxis2) ):
               		     xyfile.write('%10.3f%10.3f\n' % (star.xval, star.yval))
     		xyfile.close()
     		phot.coords = root + 'xy'
     		phot.output = root + 'mag'
		phot.verbose = no
                phot.calgorithm = "centroid"
     		phot(line,interactive=no)
     		stars = Starlist(root + 'mag')
     		zp, zpu = stars.zeropt(refstars,method="mean",rejout=no,
                                       useflags=no)
#    		print '%25s%10.3f%10.3f' % (line, zp, zpu)
     		line = infile.readline().rstrip()
     
def p60agmag(ifile, ot):

	infile = open(ifile)
	line = infile.readline().rstrip()
	while not (line == ''):
     		ot.wcs2pix(line)
     		image = line.split('.')
		root = ''
                for i in range(0,len(image)-1):
                    root+='%s.' % image[i]
      		coofile = open(root + 'coo', 'w')
      		for star in ot:
           		coofile.write('%10.3f%10.3f\n' % (star.xval, star.yval))
      		coofile.close()
     		phot.coords = root + 'coo'
     		phot.output = root + 'grb'
#		phot.verbose = yes
                #phot.calgorithm = "none"
     		phot(line, interactive=no)
     		line = infile.readline().rstrip()
 
def p60doall(ifile, ofile, refstars, filter, ot, t0):

     # Do the Zeropoint first
     p60zeropt(ifile, refstars, filter)

     # Then get the ot mags
     p60agmag(ifile, ot)

     # Now we need to format everything nicely
     inlist = iraffiles('@'+ifile)
     outfile = open(ofile, 'w')

     for image in inlist:
     
          im=image.split('.')
	  root = ''
          for i in range(0,len(im)-1):
              root+='%s.' % im[i]
          
          # Get the time
	  [utshut,exptime] = get_head(image, ['UTSHUT','EXPTIME'])
	  (utdate,uttime) = utshut.split('T')
	  year = int(utdate.split('-')[0])
	  month = int(utdate.split('-')[1])
	  day = int(utdate.split('-')[2])
	  hour = int(uttime.split(':')[0])
	  minute = int(uttime.split(':')[1])
	  second = float(uttime.split(':')[2])
	  t1 = DateTime(year, month, day, hour, minute,second)
	  dt = (t1 - t0).seconds + (exptime / 2.0)

	  # Get the Magnitude and Statistical error
          stars = Starlist(root + 'grb')
	  smag = stars[0].mags[filter]
	  smagu = stars[0].magus[filter + 'U']

	  # Get the Photometric Zeropoint error
	  refstars.wcs2pix(image)
	  refstars.set_mag(filter)
	  stars = Starlist(root + 'mag')
	  zpmag, zpmagu = stars.zeropt(refstars,method="mean",rejout=yes,
                                       useflags=no)

	  # Total Magnitude & Error
	  magu = sqrt( pow(smagu,2) + pow(zpmagu,2) )
	  mag = smag + zpmag

	  # Write everything out
	  outfile.write('%s\t%.2f\t%.2f\t%.3f\t%.3f\t%.3f\t%.3f\n' % \
	   (utshut, dt, exptime, mag, magu, smagu, zpmagu))

     outfile.close()

##############################################################################

def p60sexdoall(ilis, ofile, refstars, filter, ot, t0):

     inlist = iraffiles(ilis)
     refstars.set_mag(filter)
     outfile = open(ofile, 'w')

     for image in inlist:

          # Get zeropoint conversion and magnitude of OT
          root = image.split('.')[0]
	  iqobjs(image, 3.0, 50000.0, skyval="!SKYBKG", wtimage='', fwhm=1.5, pix=0.378, gain=2.3)
	  stars = Starlist('%s.reg' % root)
	  refstars.wcs2pix(image)
	  ot.wcs2pix(image)
          zp, zpu = stars.zeropt(refstars, method='mean',rejout=no,
                                 useflags=no)
          a, b = stars.match(ot)
	  otzmag = a.mags()[0]
          otmagu = a.magus()[0]
	  otmag = float(otzmag) + zp

	  # Now we need the date / time elapsed
	  [utshut,exptime] = get_head(image, ['UTSHUT','EXPTIME'])
	  (utdate,uttime) = utshut.split('T')
	  year = int(utdate.split('-')[0])
	  month = int(utdate.split('-')[1])
	  day = int(utdate.split('-')[2])
	  hour = int(uttime.split(':')[0])
	  minute = int(uttime.split(':')[1])
	  second = float(uttime.split(':')[2])
	  t1 = DateTime(year, month, day, hour, minute,second)
	  dt = (t1 - t0).seconds + (exptime / 2.0)

          # Write everything out
	  outfile.write('%s\t%10.2f%10.2f%10.3f%10.3f%10.3f\n' % \
	   (utshut, dt, exptime, otmag, otmagu, zpu))

     outfile.close()

#############################################################################
# get time midpoint for coadded exposure
#############################################################################

def gettmid(filelist):
   
   inlis = iraffiles('@'+filelist)
   im1 = inlis[0]
   ut1 = get_head(im1,'UTSHUT')
   t1 = DateTime(int(ut1[:4]), int(ut1[5:7]), int(ut1[8:10]), int(ut1[11:13]), int(ut1[14:16]), float(ut1[17:]))
   expt1 = get_head(im1,'EXPTIME')
   t1toadd = float(expt1) / 2.0
   t1 += RelativeDateTime(seconds=t1toadd)
   counter = 0
   for i in range(1,len(inlis)):
      imi = inlis[i] 
      uti = get_head(imi, 'UTSHUT')
      ti = DateTime(int(uti[:4]), int(uti[5:7]), int(uti[8:10]), int(uti[11:13]), int(uti[14:16]), float(uti[17:]))
      expti = get_head(imi, 'EXPTIME')
      titoadd = float(expti) / 2.0
      ti += RelativeDateTime(seconds=titoadd)
      counter += (ti-t1).days
   midt = t1 + (counter / len(inlis))
   return midt

##############################################################################
# xrtfl2fnu: Converts XRT light curves from Swift XRT Light Curve
# Repository to units of flux density
##############################################################################

def xrtfl2fnu(infile, outfile, ulfile, beta, nu0=2, nu1=0.3, nu2=10.0):

   keV2Hz = 2.24e17
   inf = open(infile)
   outf = open(outfile, 'a+')
   ulf = open(ulfile, 'a+')
   lines = inf.readlines()
   for line in lines:
      temp = line.rstrip().split()
      t = float(temp[0])
      dt = float(temp[1]) - float(temp[2])
      fl = float(temp[3])
      fnu = fl * 1e29 * (1-beta) * pow((keV2Hz * nu0), -beta) / ( pow((keV2Hz * nu2), 1-beta) - pow((keV2Hz * nu1), 1-beta) )
      dfl = float(temp[4])
      if (dfl == 0):
         dfnu = fnu
         ulf.write('%20g%20g%20g%20g\n' % (t, dt, fnu, dfnu))
      else:
         dfnu = dfl * 1e29 * (1-beta) * pow((keV2Hz * nu0), -beta) / ( pow((keV2Hz * nu2), 1-beta) - pow((keV2Hz * nu1), 1-beta) )
         outf.write('%20g%20g%20g%20g\n' % (t, dt, fnu, dfnu))
   outf.close()
   inf.close()

def dat2fl(infile, outfile, zp, A):

   inf = open(infile)
   outf = open(outfile, 'a+')
   lines = inf.readlines()
   for line in lines:
      temp = line.rstrip().split()
      td = float(temp[1])
      dt = float(temp[2])
      fl = zp * pow(10, -(float(temp[3]) - A) / 2.5)
      dfl = zp * pow(10, -(float(temp[3]) - float(temp[4]) - A) / 2.5) - fl
      #dfl2 = zp * pow(10, -(float(temp[3]) - float(temp[5]) - A) / 2.5) - fl
      outf.write('%15g%15g%15g%15g\n' % (td, dt, fl, dfl))
   outf.close()
   inf.close()

def ul2fl(infile, outfile, zp, A):

   inf = open(infile)
   outf = open(outfile, 'a+')
   lines = inf.readlines()
   for line in lines:
      temp = line.rstrip().split()
      td = float(temp[0])
      dt = float(temp[1])
      fl = zp * pow(10, -(float(temp[2]) - A) / 2.5)
      outf.write('%15g%15g%15g\n' % (td, dt, fl))
   outf.close()
   inf.close()

def dat2textab(infile, filter, mon):
   
   inf = open(infile)
   lines = inf.readlines()
   for line in lines:
      temp = line.rstrip().split()
      year = int(temp[0][:4])
      month = int(temp[0][5:7])
      day = int(temp[0][8:10])
      hour = int(temp[0][11:13])
      minute = int(temp[0][14:16])
      second = float(temp[0][17:])
      dt = float(temp[1])
      exptime = float(temp[2])
      mag = float(temp[3])
      dmag = float(temp[4])
      t0 = DateTime(year, month, day, hour, minute, second)
      td = t0.day + (((t0.second / 60.0) + t0.minute) / 60.0 + t0.hour) / 24.0
      if (dt < 1e4):
         print '$\cdots$ & %i %s %.4f & %.1f & %s & %.1f & $%.2f \pm %.2f$ \\\\' % (t0.year, mon, td, dt, filter, exptime, mag, dmag)
      elif (dt < 1e5):
         print '$\cdots$ & %i %s %.4f & $%.3f \\times 10^{4}$ & %s & %.1f & $%.2f \pm %.2f$ \\\\' % (t0.year, mon, td, dt / 1e4, filter, exptime, mag, dmag)
      elif (dt < 1e6):
         print '$\cdots$ & %i %s %.4f & $%.3f \\times 10^{5}$ & %s & %.1f & $%.2f \pm %.2f$ \\\\' % (t0.year, mon, td, dt / 1e5, filter, exptime, mag, dmag)
      else:
         print '$\cdots$ & %i %s %.4f & $%.3f \\times 10^{6}$ & %s & %.1f & $%.2f \pm %.2f$ \\\\' % (t0.year, mon, td, dt / 1e6, filter, exptime, mag, dmag)


def dat2fl2(inf, outf, zp, A, nu):
   infile=open(inf)
   outfile=open(outf,'w')
   lines=infile.readlines()
   for line in lines:
      temp=line.split()
      t=float(temp[0])
      dt=float(temp[2])
      fl=zp*pow(10, -(float(temp[3])-A)/2.5)
      dfl=zp*pow(10, -(float(temp[3])-A-float(temp[4]))/2.5)-fl
      outfile.write('%15.1f%15.3f%10.3g%15.2f%15.2f\n' % (t, dt, nu, fl, dfl))
   infile.close()
   outfile.close()
   return

###########################################################################

def natxrtbin(infile, outfile, t1, t2):

   inf=open(infile)
   outf=open(outfile, 'a+')
   lines=inf.readlines()
   sum_crate=0; sum_weights=0;

   for line in lines:
      if not re.search('#', line):
         temp=line.rstrip().split()
         if (float(temp[0])*1000.0 > t1) and (float(temp[1])*1000.0 < t2):
            crate=float(temp[5])/float(temp[7])/float(temp[2])/1000.0
            crateerr=float(temp[6])/float(temp[7])/float(temp[2])/1000.0 
	    sum_crate+=crate/pow(crateerr,2)
            sum_weights+=1/pow(crateerr,2)

   outf.write('%12.3f%12.3f%15.8f%15.8f\n' % ((t1+t2)/2.0, (t2-t1), sum_crate/sum_weights, sqrt(1/sum_weights)))
   outf.close()
   inf.close()

############################################################################

def p60_anet(imfile):

   inf = open(imfile)
   lines = inf.readlines()
   for line in lines:
      image = line.rstrip()
      #root='ccd.'+image.split('.')[1]+'.0'
      root = image.split('.')[0]
      anetcmd = "anet_p60.py %s > %s.wcs" % (image, root)
      print anetcmd
      success = 0
      while not success:
         os.system(anetcmd)
         outf = open("%s.wcs" % root)
         outlines = outf.readlines()
         outf.close()
         if re.search("Solved", outlines[len(outlines)-2]):
            update_head(image, 'IQWCS', 1)
            success=1
         elif re.search("Failed", outlines[len(outlines)-2]) or \
              re.search("Running", outlines[len(outlines)-2]):
            update_head(image, 'IQWCS', 0)
            success=1
         else:
            success=0
