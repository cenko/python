
"""OSCAR Defaults, Constants, and Classes"""

######################################################################
# 
# To-Do and Thoughts
#
#  1 - MySQL interaction via methods of the target list...
#
# Simple stuff:
# 
#  1 - Delete hourang weighting (and all traces thereof)
#  2 - Where is the "keyword binroi" coming from?
#      * Added by oscar_session, I think
#  
# $Id: oscar_defs.py,v 1.10 2007/06/20 17:45:45 derekfox Exp derekfox $
# 
######################################################################

import copy, os, shutil, glob, sys, string, re, math, time, random

import numpy 
import ephem
import socket
import signal
import getopt

from types import *
from math import *
from time_now import *

######################
## Global Variables ##
######################

# The only global variables we allow are "macros" / constants that
# NEVER CHANGE. 

# Convert degrees to radians
DTOR=0.0174532925199433
# Seconds in a day
SECSPERDAY=86400.0
# Self-explanatory
TWOPI=6.283185307179586
# Radians per hour
RADSPERHOUR=TWOPI/24.0
# Length of sidereal day in solar days
SIDEREALDAY=1.00273790935

# MJD for January 1, 1970 (the Unix "epoch")
MJDUNIXEP=40587.0

# MJD for noon on Dec 31, 1899 (the PyEphem "epoch")
MJDEPHEP=15019.5

# Instrument field of view in arcseconds
FOV=774.0

# Approximate detector readout time in seconds
TREADOUT=30.0    # new CCD
# TREADOUT=160.0 # old CCD

# Name of science status file
SCISTATFILE="Status_Science.txt"

# Name of new target file
NEWTGTFILE="oscar_new_target.lis"
TMPTGTFILE="_tmp_target.lis"

# Name of new priority file
NEWPRIFILE="oscar_new_priority.lis"
TMPPRIFILE="_tmp_priority.lis"

# Coordinates of Palomar Observatory
palomar=ephem.Observer()

# Give:  Longitude (deg E), Latitude (deg N), Elevation (m)
palomar.long, palomar.lat, palomar.elev = \
  '-116.8630',  '33.3560',      1706.0
# could also define:  mean temperature (C), barometric pressure
# (could set to zero since I think atmospheric refraction is accounted
#  for by the TCS)

# Operational restrictions for the P60
#LIM_HOURANG_EAST=-6.4
LIM_HOURANG_EAST=-6.2
#LIM_HOURANG_WEST=+6.4
LIM_HOURANG_WEST=+6.2
LIM_DEC_SOUTH=-41.8
LIM_DEC_NORTH=90.0
LIM_ELEVATION=10.0

# Just a useful default
MAX_AIRMASS=20.0

# Important astronomical objects..
palomar.date=ephem.now()
sun=ephem.Sun()
sun.compute(palomar)
moon=ephem.Moon()
moon.compute(palomar)
razero=ephem.readdb("Meridian,f|S,0.0,33.3560,15.0,2000.00")
razero.compute(palomar)

# Weight functions
weight_keys=['seeing','airmass','extinction','sky','moondeg']
weight_lims={'seeing':[0.5,5.0],
             'airmass':[1.0,5.75],
             'extinction':[0.0,10.0],
             'sky':[1.0,20.0],
             'moondeg':[0.0,180.0],
             'night':[0.0,24.0]}

# Status variables
oscar_values={'seeing':2.0,
              'extinction':1.0,
              'sky':1.0}

# Dome lamp exposures
dome_exposures={'U':300.0,
                'B':40.0,
                'V':15.0,
                'R':4.0,
                'I':2.0,
                'g':30.0,
                'r':4.0,
                'ip':4.0,
                'zp':4.0,
                'zs':5.0,
                'zl':15.0,
                'Ha20':120.0,
                'Ho20':120.0,
                'Ha100':40.0,
                'H6571':120.0,
                'H6593':120.0,
                'H6602':120.0,
                'H6614':120.0,
                'H6625':120.0,
                'H6640':120.0,
                'upr':300.0,
                'gpr':50.0,
                'rpr':10.0,
                'ipr':4.0,
                'zpr':4.0}

# Exposures for standard-field calibrations
standard_exposures={'B':30.0,
                    'V':30.0,
                    'R':30.0,
                    'I':30.0,
                    'g':30.0,
                    'r':30.0,
                    'ip':30.0,
                    'zp':30.0,
                    'zs':30.0,
                    'zl':30.0,
                    'gpr':30.0,
                    'rpr':30.0,
                    'ipr':30.0,
                    'zpr':30.0}

# Standard settings for the BINROI keyword
# This is every setting that is NOT full-chip, single-pixel
binroi_defaults=[[[1,1],[1025,1,1024,2048]],    # half-chip = B
                 [[1,1],[1025,513,1024,1024]]]  #  1/4-chip = C

# Standard defocus setting
defocus_default=0.10 # mm defocus

######################################################################

def binroi2str(binroi):

    [bin,roi]=binroi
    out="%d,%d,%d,%d,%d,%d" % (bin[0],bin[1],roi[0],roi[1],
                               roi[2],roi[3])
    return out

####################

def getlines(filename):

    infile=open(filename)
    lines=infile.readlines()
    infile.close()
    lines2=[]
    for i in range(len(lines)):
        if (re.search("\S+",lines[i])):
            lines2.append(lines[i].rstrip())
    return lines2

####################

def backupfile(filename,n=1):

    backfile=filename+(".%d" % n)
    while os.path.exists(backfile):
        n+=1
        backfile=filename+(".%d" % n)
    
    return backfile

####################

def spiral_step(jmax,kmax,jzero,kzero,nstep):

    """ determine next step in a CCW spiral mosaic 
        total grid is jmax x kmax in size;
        starting point (0-indexed) is (jzero,kzero);
        and we want the step #nstep.                  """

    # Size of the array in FOVs
    nr=(jmax-1)/2
    nd=(kmax-1)/2

    # There are two basic spirals:  One starts out with one step to
    # the right; the other starts out with one step up.  
    dmos=(jzero>nr)

    # Error checking
    if nstep>=(nr+1)*(nd+1):
        sys.stderr.write("Asked for more steps than in spiral")
        return [0,0]

    # Execute the spiral
    istep=0; ii=0
    [jmos,kmos]=[jzero,kzero]
    [jtgt,ktgt]=[jzero,kzero]

    while istep<nstep:

        # Directions are 0=right, 1=up, 2=left, 3=down

        # Target square for next move
        if dmos==0:
            jtgt+=2
        elif dmos==1:
            ktgt+=2
        elif dmos==2:
            jtgt-=2
        elif dmos==3:
            ktgt-=2
        else:
            # Error
            pass

        # Update settings for legal step
        if jtgt>=0 and jtgt<jmax and ktgt>=0 and ktgt<kmax:
            istep+=1
            [jmos,kmos]=[jtgt,ktgt]

        # Check if we need to change direction...
        ii+=1
        isqrt=int(math.sqrt(ii))

        # This happens when we complete a square or 
        #    one side of a new square.   
        if ii==isqrt*isqrt or ii==isqrt*(isqrt+1):
            dmos = (dmos+1) % 4
            dstep=0

    # The answer
    return [jmos,kmos]

####################

def spiral_mosaic(jmax,kmax,nstep):

    # Total number of steps in the mosaic
    ntot=jmax*kmax

    # Size of the array in FOVs
    nr=(jmax-1)/2
    nd=(kmax-1)/2

    # Determine oddness/evenness
    jodd = (nr % 2)
    kodd = (nd % 2)

    # Number of steps in the four interleaved spirals
    jmax1=nr;    kmax1=nd
    jmax2=nr+1;  kmax2=nd+1
    jmax3=nr+1;  kmax3=nd
    jmax4=nr;    kmax4=nd+1

    # Total number of steps for each
    ntot1=jmax1*kmax1
    ntot2=jmax2*kmax2
    ntot3=jmax3*kmax3
    ntot4=jmax4*kmax4

    # Starting points for sub-mosaics
    if jodd and kodd:
        jzero1=nr;   kzero1=nd
        jzero2=nr+1; kzero2=nd-1
        jzero3=nr+1; kzero3=nd
        jzero4=nr;   kzero4=nd-1
    elif jodd and not kodd:
        jzero1=nr;   kzero1=nd-1
        jzero2=nr+1; kzero2=nd
        jzero3=nr+1; kzero3=nd-1
        jzero4=nr;   kzero4=nd
    elif kodd and not jodd:
        jzero1=nr+1; kzero1=nd
        jzero2=nr;   kzero2=nd-1
        jzero3=nr;   kzero3=nd
        jzero4=nr+1; kzero4=nd-1
    elif not jodd and not kodd:     
        jzero1=nr+1; kzero1=nd-1
        jzero2=nr;   kzero2=nd
        jzero3=nr;   kzero3=nd-1
        jzero4=nr+1; kzero4=nd
    else:
        # problem?
        pass

    if nstep<ntot1:
        # First spiral
        nstep1=nstep
        [jmos,kmos]=spiral_step(jmax,kmax,jzero1,kzero1,nstep1)
    elif nstep<ntot1+ntot2:
        # Second spiral
        nstep2=nstep-ntot1
        [jmos,kmos]=spiral_step(jmax,kmax,jzero2,kzero2,nstep2)
    elif nstep<ntot1+ntot2+ntot3:
        # Third spiral
        nstep3=nstep-ntot1-ntot2
        [jmos,kmos]=spiral_step(jmax,kmax,jzero3,kzero3,nstep3)
    else:
        # Fourth spiral
        nstep4=nstep-ntot1-ntot2-ntot3
        [jmos,kmos]=spiral_step(jmax,kmax,jzero4,kzero4,nstep4)

    return [jmos,kmos]

##############################

def sphdist(ra1,dc1,ra2,dc2,arcsec=1,arcmin=0,degree=0):

    """ Use Haversine's formula to calculate distance on a sphere """

    ddtor=0.0174532925199433

    ra1r=ra1*ddtor
    dc1r=dc1*ddtor
    ra2r=ra2*ddtor
    dc2r=dc2*ddtor
    
    dra=ra2r-ra1r
    ddc=dc2r-dc1r
    a=sin(ddc/2)**2 + cos(dc1r)*cos(dc2r)*sin(dra/2)**2
    cdeg=2*asin(min(1,sqrt(a)))/ddtor

    if arcsec and not arcmin and not degree:
        cans=3600*cdeg
    elif arcmin:
        cans=60*cdeg
    else:
        cans=cdeg

    return cans

##############################

def sphdist2(ra1,dc1,ra2,dc2,arcsec=1,arcmin=0,degree=0):

    """ Use Law of Cosines to calculate distance on a sphere """

    ddtor=0.0174532925199433
    
    dra=(ra2-ra1)*ddtor

    cd = sin(dc1*ddtor)*sin(dc2*ddtor) + \
         cos(dc1*ddtor)*cos(dc2*ddtor)*cos(dra)
    c=acos(cd)/ddtor

    if arcsec:
        cans=3600*c
    elif arcmin:
        cans=60*c
    else:
        cans=c

    return cans

##############################

def radecdiff(ra1,dc1,ra2,dc2,arcsec=1,arcmin=0,degree=0):

    # RA Distance
    radiff=sphdist(ra1,dc1,ra2,dc1,
                   arcsec=arcsec,arcmin=arcmin,degree=degree)
    # RA Sign
    rasign=+1
    if ra1<ra2:
        rasign=-1
    # This will be wrong if the two RA's straddle 24h...
    if abs(ra1-ra2)>180:
	# ...so correct it
	rasign=-rasign
    # RA Answer
    radiff *= rasign

    # Dec Distance
    dcdiff=sphdist(ra2,dc1,ra2,dc2,
                   arcsec=arcsec,arcmin=arcmin,degree=degree)
    # Dec Sign
    dcsign=+1
    if dc1<dc2:
        dcsign=-1
    # Dec Answer
    dcdiff *= dcsign

    # The End
    return [radiff,dcdiff]

######################################################################

class astrocoords:

    def __init__(self,ra,dec,equinox=2000.0,system="fk5"):

        self.set(ra,dec,equinox=equinox,system=system)

    def set(self,ra,dec,equinox=2000.0,system="fk5"):
            
        # Parse RA
        if type(ra)==StringType:
            ra_reg1=re.search("(\d+)[:\s](\d+)[:\s]([\d\.]+)",ra)
            ra_reg2=re.search("([\d\.]+)",ra)

            if ra_reg1:
                radeg = 15*(float(ra_reg1.group(1)) + \
                            (float(ra_reg1.group(2)) + \
                             float(ra_reg1.group(3))/60.0)/60.0)
            elif ra_reg2:
                radeg = float(ra_reg2.group(1))
            else:
                raise ValueError, "Failed to parse input RA-string \"%s\" as RA"
        elif type(ra)==FloatType:
            radeg = ra
        else:
            try:
                radeg=float(ra)
            except:
                raise ValueError, "Failed to parse input RA as float"

        # Check range limits on RA
        if radeg<0 or radeg>360:
            raise ValueError, "Parsed RA as illegal value %.5f" % radeg

        # Parse Dec
        if type(dec)==StringType:
            dc_reg1=re.search("([\+\-]?\d+)[:\s](\d+)[:\s]([\d\.]+)",dec)
            dc_reg2=re.search("([\+\-]?[\d\.]+)",dec)

            if dc_reg1:
                sgndec=+1
                if dec[0]=='-':
                    sgndec=-1
                absdcd=abs(float(dc_reg1.group(1)))
                dcdeg = sgndec*(absdcd + \
                                (float(dc_reg1.group(2)) + \
                                 float(dc_reg1.group(3))/60.0)/60.0)
            elif dc_reg2:
                dcdeg = float(dc_reg2.group(1))
            else:
                raise ValueError, "Failed to parse input Dec-string \"%s\" as Dec"
        elif type(dec)==FloatType:
            dcdeg = dec
        else:
            try:
                dcdeg=float(dec)
            except:
                raise ValueError, "Failed to parse input Dec as float"

        # Check range limits on Dec
        if abs(dcdeg)>90:
            raise ValueError, "Parsed Dec as illegal value %.5f" % dcdeg

        # Grab the equinox
        xeqx=float(equinox)

        # Check range limits on Equinox
        if xeqx<0 or xeqx>3000:
            raise ValueError, "Illegal Equinox value %.5f" % xeqx

        # Save everything
        self.__radeg=radeg
        self.__dcdeg=dcdeg

        # Direct properties
        self.equinox=xeqx
        self.system=system

    def sxg(self):
        radeg=self.__radeg
        absdec=abs(self.__dcdeg)
        sgndec=cmp(self.__dcdeg,0)
        
        rah=int(radeg/15.0)
        ram=int(60*(radeg/15.0-rah))
        ras=60*(60*(radeg/15.0-rah)-ram)

        dcsgn="+"
        if sgndec<0:
            dcsgn="-"

        dcd=int(absdec)
        dcm=int(60*(absdec-abs(dcd)))
        dcs=60*(60*(absdec-abs(dcd))-dcm)

        return ["%02i:%02i:%07.4f" % (rah,ram,ras),
                dcsgn+"%i:%02i:%06.3f" % (dcd,dcm,dcs)]

    def deg(self):
        return [self.__radeg,self.__dcdeg]

    def rasxg(self):
        rasxg,dcsxg=self.sxg()
        return rasxg

    def dcsxg(self):
        rasxg,dcsxg=self.sxg()
        return self.dcsxg

    def radeg(self):
        return self.__radeg

    def dcdeg(self):
        return self.__dcdeg

    def set_radeg(self,radeg):
        try:
            xradeg=float(radeg)
        except:
            print "radeg must be floating-point RA in degrees"
            return
        if xradeg<0 or xradge>=360:
            print "radeg must be floating-point RA in degrees"
            return
        self.__radeg=xradeg

    def set_dcdeg(self,dcdeg):
        try:
            xdcdeg=float(dcdeg)
        except:
            print "dcdeg must be floating-point Declination in degrees"
            return
        if xdcdeg<-90 or xdcdge>90:
            print "dcdeg must be floating-point Declination in degrees"
            return
        self.__dcdeg=xdcdeg

    def set_deg(self,radeg,dcdeg):
        self.set_radeg(radeg)
        self.set_dcdeg(dcdeg)

    def rad(self):
        dtor=180.0/math.pi
        rarad=self.__radeg/dtor
        dcrad=self.__dcdeg/dtor
        return [rarad,dcrad]

    def shift(self,shift,arcsec=0,arcmin=0):
        dtor=180.0/math.pi

        if arcsec:
            # Shift is quoted in arcsec
            ashift=[shift[0]/3600.0,shift[1]/3600.0]
        elif arcmin:
            # Shift is quoted in arcmin
            ashift=[shift[0]/60.0,shift[1]/60.0]
        else:
            # Shift is quoted in degrees (default)
            ashift=shift

        # New RA
        if abs(self.__dcdeg)<90:
            radeg=self.__radeg+ashift[0]/math.cos(self.__dcdeg/dtor)
        else:
            # Degenerate case
            radeg=self.__radeg

        # Correct for RA negative or >360
        while radeg<0:
            radeg += 360
        while radeg>=360:
            radeg -= 360

        # New Dec
        dcdeg=self.__dcdeg+ashift[1]
        if abs(dcdeg)>90:
            # Oops, over the pole
            sgndec=cmp(dcdeg,0)
            dcdeg=sgndec*(180-abs(dcdeg))

        # Save 'em
        self.__radeg=radeg
        self.__dcdeg=dcdeg

    def diff(self,other,arcsec=0,arcmin=0,degree=1):
        [ra1,dc1]=self.deg()
        [ra2,dc2]=other.deg()
        coodiff=radecdiff(ra1,dc1,ra2,dc2,arcsec=arcsec,
                          arcmin=arcmin,degree=degree)
        return coodiff

######################################################################

class moving_ephem:

    """ Class for specifying moving objects without current orbital
        approximations """

    def __init__(self,moving,name="",observer=None,
                 tnow=None,fake=0,time=0.0,equinox=2000.0):

        # Target properties
        self.name=name
        self.equinox=equinox
        self.coords=astrocoords(0.0,0.0)

        #####################
        # 
        # Format of the "moving" list:
        # 
        # Each element is a list giving:
        # 
        #   [MJD, RA deg, Dec deg]  --OR--
        #   [MJD, RA deg, Dec deg, RA arcsec/hour, Dec arcsec/hour]
        #
        # If MJD<MJD[0]  then we use first position (+ optional motion)
        # If MJD>MJD[-1] then we use last position (+ optional motion)
        # 
        #####################

        # Be careful in parsing the "moving" list
        moveout=[]
        for i in xrange(len(moving)):
            movel=moving[i]
            coo=astrocoords(movel[1],movel[2])
            if len(movel)>=5:
                try:
                    moveout.append([float(movel[0]),
                                    coo.radeg(),
                                    coo.dcdeg(),
                                    float(movel[3]),
                                    float(movel[4])])
                except:
                    print "Failed to parse a moving element"
            elif len(movel)>=3:
                try:
                    moveout.append([float(movel[0]),
                                    coo.radeg(),
                                    coo.dcdeg()])
                except:
                    print "Failed to parse a moving element"
            else:
                print "Moving element had insufficient elements"
            
        # This is what we keep
        self.moving=moveout
        self.sort_moving()

        # Establish our tnow object
        if type(tnow)==NoneType:
            self.tnow=time_now(fake=fake,time=time)
        else:
            self.tnow=tnow

        # Establish our coords & ephemeris objects
        self.set_eph()

        # Compute ourselves
        if type(observer)!=NoneType:
            self.compute(observer)

    def sort_moving(self):

        """ Sort moving ephemeris entries into chronological order """

        mjds1=[]
        for move in self.moving:
            mjds1.append(move[0])
        mjds2=copy.deepcopy(mjds1)
        mjds2.sort()

        newmove=[]
        for mjd in mjds2:
            ix=mjds1.index(mjd)
            if not (self.moving[ix] in newmove):
                newmove.append(self.moving[ix])

        self.moving=newmove

    def set_coords(self):

        # Number of steps in ephemeris
        nsteps=len(self.moving)

        # Determine which entries to interpolate between
        i=0
        mjdnow=self.tnow.mjd()
        if mjdnow>self.moving[0][0]:
            for i in xrange(nsteps-1):
                if mjdnow>=self.moving[i][0] and \
                   mjdnow<self.moving[i+1][0]:
                    break

        # Use the given ephemeris to get our position
        dthours=24.0*(mjdnow-self.moving[i][0])
        if len(self.moving[i])==5:
            # Parse the moviing list
            rabase,ramove=self.moving[i][1],self.moving[i][3]
            dcbase,dcmove=self.moving[i][2],self.moving[i][4]
            # Convert RA angular motion to coordinate motion
            racmove=ramove/math.cos(dcbase*DTOR)
            # Calculate moving position
            radeg=(rabase + dthours*racmove/3600.0) % 360.0
            dcdeg=dcbase + dthours*dcmove/3600.0
            # "over the pole"
            if abs(dcdeg)>90:
                if dcdeg>0:
                    dcdeg = 180.0 - dcdeg
                else:
                    dcdeg = -180.0 - dcdeg
        else:
            radeg,dcdeg=self.moving[i][1],self.moving[i][2]

        # Set our coords object
        self.coords.set(radeg,dcdeg,equinox=self.equinox)

        # Set additional properties
        self.radeg=radeg
        self.dcdeg=dcdeg

    def set_eph(self):

        self.set_coords()
        [rasxg,dcsxg]=self.coords.sxg()
        self.eph=ephem.readdb("%s,f|S,%s,%s,15.0,%.2f" % \
                              (self.name,rasxg,dcsxg,self.equinox))

    def compute(self,observer):

        self.set_eph()
        self.eph.compute(observer)

        # Set ephem-type properties
        self.ra=self.eph.a_ra
        self.a_ra=self.eph.a_ra 
        self.g_ra=self.eph.g_ra
        self.dec=self.eph.a_dec
        self.a_dec=self.eph.a_dec
        self.g_dec=self.eph.g_dec
        self.alt=self.eph.alt
        self.az=self.eph.az
        self.radius=self.eph.radius
        
        # These will not be very correct but oh well...
        #self.eph.compute(observer)
        #self.rise_time=observer.next_rising(self.eph)
        #self.rise_az=self.eph.az
        #self.eph.compute(observer)
        #self.transit_time=observer.next_transit(self.eph)
        #self.transit_alt=self.eph.alt
        #self.set_time=self.eph.set_time
        #self.set_az=self.eph.set_az
        #self.eph.compute(observer)
       

######################################################################

class oscar_exposure:

    def __init__(self,filter,exptime,imult,jmax=1,kmax=1,
                 obsq=[],edither=[],obsi=[],target="Parent",
                 treadout=TREADOUT):

        # Check for reasonable filter
        if type(filter)!=StringType:
            raise TypeError, "Exposure filter must be passed as a string"

        # Check for reasonable exposure time
        xexptime=float(exptime)
        if xexptime<0 or xexptime>10000:
            raise ValueError, "Exposure time out of range: %f" % xexptime

        # Check for reasonable exposure multiple
        iimult=int(imult)
        if iimult<0:
            raise ValueError, "Exposure multiple out of range: %i" % iimult

        # Check mosaic setting
        if jmax<1 or kmax<1:
            sys.stderr.write("Illegal mosaic settings in "+ \
                             "oscar_exposure definition\n")
            jmax=1
            kmax=1

        # Check for reasonable exposure dithers
        if type(edither)!=ListType:
            raise TypeError, "Exposure dither settings must "+ \
                             "be passed as a string list"

        # Define the exposure
        self.filter=filter
        self.exptime=xexptime
        self.treadout=treadout
        self.imult=iimult
        self.jmax=jmax
        self.kmax=kmax
        self.edither=copy.deepcopy(edither)
        self.target=target
        self.done=0
        self.skipped=0

        # Partial mosaic grid
        blank=numpy.array(0,'UInt8')
        self.obsq=numpy.resize(blank,(self.jmax,self.kmax))
        if len(obsq)>=1:
            ii=0
            for j in xrange(self.jmax):
                for k in xrange(self.kmax):
                    self.obsq[j,k]=int(obsq[ii])>0
                    ii += 1

        # Exposure completion tracking
        if len(obsi)<1:
            self.obsi=numpy.resize(blank,(self.jmax,self.kmax))
        else:
            self.obsi=obsi

        # Check for doneness
        self.check_done()

    ####################
        
    def square(self):

        mindone=self.mindone()
        if mindone==self.imult:
            return [-1,-1]
        
        for nstep in xrange(0,self.jmax*self.kmax):
            [j,k]=[nstep % self.jmax, floor(nstep/self.jmax)]
            if self.obsq[j,k]:
                continue
            if self.obsi[j,k]==mindone:
                break

        return [j,k]

    def spiral(self):

        mindone=self.mindone()
        if mindone==self.imult:
            return [-1,-1]
        
        for nstep in xrange(0,self.jmax*self.kmax):
            [j,k]=spiral_mosaic(self.jmax,self.kmax,nstep)
            if self.obsq[j,k]:
                continue
            if self.obsi[j,k]==mindone:
                break

        return [j,k]

    def update(self,j=0,k=0):

        self.obsi[j,k] += 1
        self.check_done()

    def mindone(self):

        # This puts "imult" in the obsi matrix everywhere we're not
        # supposed to observe (that is, where obsq[ii]==1)
        iactive=numpy.choose(self.obsq,(self.obsi,self.imult))
        return iactive.min()

    def check_done(self):

        self.done=(self.mindone()>=self.imult)

    def duration(self,treadout=-1):

        """ Total exposure + readout time, in seconds """

        if treadout<0:
            treadout=self.treadout
        return (self.exptime+treadout) * self.imult * \
               (self.jmax*self.kmax - self.obsq.sum())

    def time_left(self,treadout=-1):

        """ Exposure + readout time remaining, in seconds """
        
        if treadout<0:
            treadout=self.treadout

        tleft = (self.exptime+treadout) * \
                (self.imult*self.jmax*self.kmax
                 - self.obsi.sum()
                 - self.imult*self.obsq.sum())

        if tleft<0:
            tleft=0.0

        return tleft

    def reset(self,zero=0):

        self.obsi[:]=zero
        self.skipped=0
        self.check_done()

    def skip(self):

        # Legal positions get imult, illegal positions left at zero
        fullobsi=numpy.choose(self.obsq,(self.imult,0))

        self.obsi[:]=fullobsi[:]
        self.skipped=1
        self.done=1

    def __str__(self):

        strout="%dx %s for %.2f s" % \
                (self.imult,self.filter,self.exptime)

        if type(self.target)!=StringType:
            strout += " (Target %s)" % self.target.name
        return strout
    
######################################################################

class oscar_target:

    def __init__(self,name,ra,dec,exposures,priority=1,
                 offset=[],nr=0,nd=0,mscq=[],dither=[],weight={},
                 binning=[],roi=[],timing={},mjdstart=0.0,
                 repeat=1,obsn=0,obst=0.0,tnow=None,fake=0,time=0.0,
                 observer=palomar,keyword={},property={},equinox=2000.0,
                 system="fk5",pmra=0.0,pmdec=0.0,sunlim=-12.0,
                 twilight=None,subtargets=[],badfilt=[],
                 nonsidereal="",planet="",moving=[]):

        # Defaults
        name_maxlen=19
        def_offset=[0.0,0.0]
        def_dither=[60,60]
        def_binning=[1,1]
        def_roi=[1,1,2048,2048]
        def_weight={'airmass':[2.5,0.25],
                    'moondeg':[150.0,0.1],
                    'night':[24.0,0.5]}
        def_timing={'log':['1.0']}
        def_keyword={'P60PRID':'TEST',
                     'P60PRNM':'"SCIENCE TEST"',
                     'P60PRPI':'FOX_CENKO',
                     'P60PRTM':'24'}

        # Establish tnow object
        if type(tnow)==NoneType:
            self.tnow=time_now(fake=fake,time=time)
        else:
            self.tnow=tnow

        # Establish observer object (default is J2000)
        self.observer=observer

        # Object name
        if len(name)>name_maxlen:
            print "Target name '%s' is more than %d characters, please shorten." % (name,name_maxlen)
            sys.exit(3)
        self.name=name.replace('-','_') # no hyphens in target names

        # Ephemeris object and coordinate-type settings
        self.equinox=equinox
        self.system=system
        self.nonsidereal=nonsidereal
        self.planet=planet
        self.moving=moving
        if len(nonsidereal)>0 or len(planet)>0 or len(moving)>0:

            # Nonsidereal conic-section orbits
            if len(nonsidereal)>0:
                try:
                    self.eph=ephem.readdb("%s,%s" % (self.name,nonsidereal))
                except:
                    print "Error interpreting nonsidereal ephemeris "+\
                          ("'%s'" % nonsidereal)
                    sys.exit(3)

            # Common planets
            elif len(planet)>0:
                if re.search('mercury',planet,re.I):
                    self.eph=ephem.Mercury()
                elif re.search('venus',planet,re.I):
                    self.eph=ephem.Venus()
                elif re.search('moon',planet,re.I):
                    self.eph=ephem.Moon()
                elif re.search('mars',planet,re.I):
                    self.eph=ephem.Mars()
                elif re.search('jupiter',planet,re.I):
                    self.eph=ephem.Jupiter()
                elif re.search('saturn',planet,re.I):
                    self.eph=ephem.Saturn()
                elif re.search('uranus',planet,re.I):
                    self.eph=ephem.Uranus()
                elif re.search('neptune',planet,re.I):
                    self.eph=ephem.Neptune()
                elif re.search('pluto',planet,re.I):
                    self.eph=ephem.Pluto()
                else:
                    print "Unrecognized planet '%s'" % planet
                    sys.exit(3)

            # Moving targets
            elif len(moving)>0:
                # Created our own superclass
                self.eph=moving_ephem(moving)

            # Create ephemeris & coordinates objects
            self.eph.compute(self.observer)
            self.coords=astrocoords(str(self.eph.ra),str(self.eph.dec),
                                    equinox,system)

        else:
            # Regular Sidereal targets with RA/Dec positions and proper motions
            self.coords=astrocoords(ra,dec,equinox,system)
            self.pmra=pmra
            self.pmdec=pmdec
            [rasxg,dcsxg]=self.coords.sxg()
            self.eph=ephem.readdb("%s,f|S,%s,%s,15.0,%.2f" % \
                                  (self.name,rasxg,dcsxg,equinox))
            self.eph.compute(self.observer)

        # Basic properties
        self.priority=float(priority)
        self.badfilt=badfilt
        self.sunlim=sunlim
        self.property=property
        self.lastscore=0.0
        self.done=0

        # Record of observations
        self.repeat=repeat
        self.obsn=obsn
        self.obst=obst

        # Shape of mosaic grid
        self.nr=nr
        self.nd=nd
        # Partial mosaic specification, if given
        self.mscq=mscq

        # MJD Start (for entry into queue)
        mjdtime=self.mjdnow()
        if mjdstart>50000.0:
            self.mjdstart=mjdstart
        else:
            self.mjdstart=mjdtime

        # Offset
        if len(offset)>1:
            self.offset=offset[0:2]
        else:
            self.offset=def_offset

	# Dither
	if len(dither)>0:
	    self.dither=dither
	else:
	    self.dither=def_dither

        # Binning
        if len(binning)==1:
            if type(binning)==ListType:
                self.binning=[binning[0],binning[0]]
            else:
                self.binning=[binning,binning]
        elif len(binning)==2:
            if binning[0]!=binning[1]:
                print "Anisotropic binning not allowed"
            self.binning=[binning[0],binning[0]]
        else:
            self.binning=def_binning

        # ROI
        if len(roi)==4:
            self.roi=roi
        else:
            self.roi=def_roi

        # Additional header keywords
        if len(keyword)>0:
            self.keyword=keyword
        else:
            self.keyword=def_keyword

        # Weight functions
        self.weight=def_weight
        for key in weight.keys():
            self.weight[key]=weight[key]

        # Timing weight function
        if len(timing)>0:
            self.timing=timing
        else:
            self.timing=def_timing

        # Characteristics of requested exposures
        self.exposures=exposures
        self.nexp=len(exposures)
        self.subtargets=subtargets

        # Special source properties
        [radeg,decdeg]=self.coodeg()
        minairmass=1/math.cos(decdeg*DTOR-self.observer.lat)
        self.minairmass=minairmass
        absmaxz=weight_lims['airmass'][1]
        maxtransitd=180-decdeg-self.observer.lat/DTOR
        if math.acos(1/absmaxz)/DTOR>maxtransitd:
            maxairmass=1/math.cos(maxtransitd*DTOR)
        else:
            maxairmass=absmaxz
        self.maxairmass=maxairmass

        # Targets that we are happy to observe "partially"
        self.partialok=self.property.has_key('PARTIALOK') and \
                            self.property['PARTIALOK']

        # Set non-visible targets to "done"
        if minairmass>self.weight['airmass'][0]:
            self.done=1
            sys.stderr.write("Target "+name+" will never " +
                             "satisfy airmass constraint\n")

        # Check against Palomar limits
        if self.dcdeg()<LIM_DEC_SOUTH:
            self.done=1
            sys.stderr.write("Target "+name+" violates the southern " +
                             "pointing restriction\n")

        # Score tracking
        sckeys=self.weight.keys()
        for tm_key in self.timing.keys():
            sckeys.append('time_'+tm_key)
        self.scores={}
        for sckey in sckeys:
            self.scores[sckey]=[0.0,1.0]

        # Set twilight times
        if type(twilight)==NoneType:
            self.set_twilight()
        else:
            self.twilight=twilight

        # Calculate target duration property
        self.set_duration()

        # Determine alive times
        [up1,up2]=self.set_alive_times()

        # Make sure target is fully observable
        if type(up1)==NoneType:
            print "Target '%s' will not satisfy its constraints tonight" % \
                  self.name
        else:
            uphours=24*(up2-up1)
            if uphours<=0:
                print ("Target '%s' will not be up long enough "+ \
                       "to be observed") % self.name
                uphours=0.1 # nonsense value
            if self.weight.has_key('night'):
                # Hours of observable time for target
                self.weight['night'][0]=uphours

        # Check if we're done
        self.check_done()

    ##############################

    def hourang_from_airmass(self,airmass):

        # Observer latitude in radians
        latrad=float(self.observer.lat)
        # Relevant coordinate pole
        if latrad>=0:
            prad=0.5*math.pi
        else:
            prad=-0.5*math.pi
        # Source elevation for given airmass, in radians
        altrad=math.asin(1.0/airmass)
        # Source RA & declination, in radians
        [rarad,dcrad]=self.coorad()
        # Law of cosines gives hour angle in radians
        acosarg=(math.cos(0.5*math.pi-altrad) - \
                 math.cos(prad-dcrad)*math.cos(prad-latrad)) / \
                 (math.sin(prad-dcrad)*math.sin(prad-latrad))
        # cos(a) = cos(b)*cos(c) + cos(A)*sin(b)*sin(c)
        #       -->
        # cos(A) = (cos(a)-cos(b)*cos(c))/(sin(b)*sin(c))
        if abs(acosarg)<1:
            hourrad=math.acos(acosarg)
        else:
            print "Airmass out of range in hourang_from_airmass"
            hourrad=12.0*15.0*DTOR

        # Convert to hours
        hourang=hourrad/DTOR/15.0

        return hourang

    def airmass_from_hourang(self,hourang):

        # Hour angle in radians
        hourrad=hourang*15*DTOR
        # Observer latitude in radians
        latrad=float(self.observer.lat)
        # Relevant coordinate pole
        if latrad>=0:
            prad=0.5*math.pi
        else:
            prad=-0.5*math.pi
        # Source RA & declination, in radians
        [rarad,dcrad]=self.coorad()
        # Law of cosines gives altitude in radians
        altrad=0.5*math.pi - \
               math.acos(math.cos(prad-dcrad)*math.cos(prad-latrad) + \
                         math.cos(hourrad)*math.sin(prad-dcrad)*\
                                           math.sin(prad-latrad))
        # cos(a) = cos(b)*cos(c) + cos(A)*sin(b)*sin(c)
        if altrad<LIM_ELEVATION*DTOR:
            # Impose a maximum limit
            airmass=MAX_AIRMASS
        else:
            # Airmass from altitude
            airmass=1.0/math.sin(altrad)

        return airmass

    ##############################

    def now(self):
        return self.tnow.now()

    def date(self):
        return self.tnow.date()

    def mjdnow(self):
        return self.tnow.mjd()

    def lstnow(self):
        return self.tnow.lst(self.observer)

    def lstrange(self):
        return self.tnow.lstrange(self.observer)

    def lstmid(self):
        return self.tnow.lstmid(self.observer)

    def coosxg(self):
        return self.tnow.objsxg(self.observer,self.eph)

    def coodeg(self):
        return self.tnow.objdeg(self.observer,self.eph)

    def coorad(self):
        return self.tnow.objrad(self.observer,self.eph)

    def radeg(self):
        return self.tnow.objradeg(self.observer,self.eph)

    def dcdeg(self):
        return self.tnow.objdecdeg(self.observer,self.eph)

    def altitude(self):
        return self.tnow.altitude(self.observer,self.eph)

    def azimuth(self):
        return self.tnow.azimuth(self.observer,self.eph)

    def airmass(self):
        return self.tnow.airmass(self.observer,self.eph)

    def moondeg(self):
        return self.tnow.moondeg(self.observer,self.eph)

    def hourang(self):
        return self.tnow.hourang(self.observer,self.eph)

    def transit_time(self):
        return self.tnow.transit_time(self.observer,self.eph)

    def in_sky_times(self):
        return self.tnow.objinsky(self.observer,self.eph)

    def alt_times(self,alt):
        return self.tnow.alt_times(self.observer,self.eph,alt)

    def up_times(self,alt,sunalt=-12.0):
        return self.tnow.objup(self.observer,self.eph,alt=alt,sunalt=sunalt)

    def up_hours(self,alt,sunalt=-12.0):
        [u1,u2]=self.up_times(alt,sunalt=sunalt)
        return 24*(u2-u1)

    ##############################

    def time_left(self):

        """Duration of all exposures left to be done for this
           target, including readout times, in hours"""

        # "PARTIALOK" targets are always "almost done"
        if self.partialok:
            return 0.01

        # Ordinary targets have their time left calculated
        tleft=0.0
        for xpr in self.exposures:
            # Remaining time for exposure, in seconds
            tleft += xpr.time_left()

        return tleft/3600.0

    ##############################

    def set_duration(self):

        """Total duration of all exposures in the target, including
           readout times, in hours"""

        ttot=0.0
        
        for xpr in self.exposures:
            # Total time for exposure
            ttot += xpr.duration()

        self.duration=ttot/3600.0
        return self.duration

    def set_twilight(self):
        self.twilight=self.tnow.twilight(self.observer,sunalt=self.sunlim)
        return self.twilight

    def set_alive_times(self):

        # We consider only the hour angle and airmass cutoffs,
        # plus morning/evening twilight

        # Hour-angle cutoffs
        [hacut1,hacut2]=[-12,12]
        if self.weight.has_key('hourang'):
            hacutoff=self.weight['hourang'][0]
            hacut1=-hacutoff
            hacut2=+hacutoff

        hacut1=max([hacut1,LIM_HOURANG_EAST])
        hacut2=min([hacut2,LIM_HOURANG_WEST])

        # Airmass cutoff corresponding to the hour angle
        [amcut1,amcut2]=[self.airmass_from_hourang(hacut1),
                         self.airmass_from_hourang(hacut2)]
        
        # Times at airmass cutoff
        if self.weight.has_key('airmass'):
            # Target's airmass cutoff
            amcutoff=self.weight['airmass'][0]
        else:
            # Maximum airmass for the facility
            amcutoff=1/math.sin(LIM_ELEVATION*DTOR)

        amcut1=min([amcut1,amcutoff])
        amcut2=min([amcut2,amcutoff])

        # Convert to altitude
        [altcut1,altcut2]=[math.asin(1.0/amcut1)/DTOR,
                           math.asin(1.0/amcut2)/DTOR]

        # Calculate the uptime for that altitude
        [up1,up2]=self.up_times(altcut1)
        if amcut1!=amcut2:
            [junk,up2]=self.up_times(altcut2)

        # Set property
        self.alive_times=[up1,up2]
        # And return
        return [up1,up2]

    ##############################

    def night_value(self):

        """ How many hours until the next source-set (accounting for
            duration of the target)?
            
            - Negative if we think source-set is closer than time_left()
            - A large number less than 24.0 if source has set for tonight
            - Not computable or error:  Returns 25.0 """

        # Uptimes (in days, via ephem)
        [up1,up2]=self.alive_times
        if type(up1)==NoneType:
            return 25.0

        # Adjust set-time for target duration
        # up2 -= self.time_left()/24.0

        # Night hours range
        uphours=24*(up2-up1)

        # Check for possibility
        if up2<up1:
            return 25.0

        # Compare current time to night hours range
        now=self.now()
        if now>=up1 and now<=up2:
            # Source is up
            ntval=24*(up2-now)-self.time_left()
        elif now<up1:
            # Source is not up yet
            ntval=24*(up2-now)-self.time_left()
        else:
            # Source has already set
            ntval=24*(1-(now-up2))

        return ntval

    ##############################

    def alive_now(self):

        """ Are we currently observable? """

        # Uptimes (in days, via ephem)
        [up1,up2]=self.alive_times
        if type(up1)==NoneType:
            return 0

        # Compare current time to uptimes
        now=self.now()
        if now>=up1 and now<=up2:
            alive_now=1
        else:
            alive_now=0

        return alive_now

    ##############################

    def weight_function(self,weight_key,key_val,lims):

        # Impose absolute limits
        value=key_val
        if value<lims[0]:
            value=lims[0]
        elif value>lims[1]:
            value=lims[1]

        # Interpret the weighting function
        [wt_zero,wt_slope]=self.weight[weight_key]

        if value>wt_zero:
            return 0

        vlim=wt_zero
        if vlim>lims[1]:
            vlim=lims[1]

        vslope=wt_slope
        if abs(vslope)>1:
            vslope=cmp(vslope,0)*1.0

        # Apply the weighting function
        if value>vlim:
            score=0
        else:
            wt_mean=(lims[1]-lims[0])/(vlim-lims[0])
            maxslope=2*wt_mean/(vlim-lims[0])
            xslope=vslope*maxslope
            vmean=0.5*(vlim-lims[0])
            score=wt_mean - xslope*(value-vmean)
            # This last check should not be necessary
            if score<=0: score=0

        # Done
        return score

    ##############################

    def timing_periodic(self,mjdnow,meanep,coeffs,slop):

        mjdlast=self.obst
        xslop=float(slop)

        tdays=mjdnow-float(meanep)
        tlast=mjdlast-float(meanep)
        phase=float(coeffs[0])
        phlast=float(coeffs[0])
        for i in range(0,len(coeffs)):
            phase += float(coeffs[i])*tdays**i
            phlast += float(coeffs[i])*tlast**i

        if round(phase)==round(phlast):
            # In the same cycle as our last set of observations!
            tscore=0
        else:
            # In a new cycle -- try to hit phase==0
            phrem=(phase % 1.0)
            if phrem>0.5:
                phrem -= 1
            aphase=abs(phrem)
            if aphase>xslop:
                tscore=0
            else:
                tscore=(xslop-aphase)/(xslop*xslop)

        return tscore

    ##############################
    
    def timing_function(self,time_key,mjdnow):

        score=1
        mjdstart=self.mjdstart
        mjdlast=self.obst

        if time_key=='none':
            # No time dependence: No effect
            score=1
	elif time_key=='window':
            # Window:  Even weighting throughout a specified window
            mjdwdw=self.timing[time_key][0:2]
            [mjd1,mjd2]=[float(mjdwdw[0]),float(mjdwdw[1])]
            if mjdnow>mjd1 and mjdnow<mjd2:
                score=1
            else:
                score=0
        elif time_key=='monitor':
            # Monitoring: Want to get (one iteration of) all exposures
            # in at approximately same time
            speed=float(self.timing[time_key][0])
            if mjdlast==0:
                score=1+math.log(1+speed*(mjdnow-mjdstart))
            elif mjdnow<mjdlast+0.5/speed:
                # Not too close together...
                score=0
            elif mjdnow>=mjdlast+0.5/speed:
                score=1+math.log(speed*(mjdnow-mjdlast))
            else:
                score=1
        elif time_key=='log':
            # Log-timing:  Relative to MJDstart time
            speed=float(self.timing[time_key][0])
            if mjdnow>mjdstart:
                score=1+math.log(1+speed*(mjdnow-mjdstart))
            else:
                score=1
        elif time_key=='phase':
            # Phase-timing:  Try to get close to zero-phase
            [meanep,coeffs,slop]=self.timing[time_key][0:3]
            coeffx=coeffs.split(',')
            # Do the calculation (see routine for details)
            score=self.timing_periodic(mjdnow,meanep,coeffx,slop)
        elif time_key=='repeat':
            # Repeat:  Take observations at a fixed cadence
            # (really just a special case of periodic)
            [meanep,wait,slop]=self.timing[time_key][0:3]
            coeffx=[0.0,1.0/float(wait)]
            if mjdlast==0:
                # Never observed before
                score=1+math.log(1+(mjdnow-mjdstart))
            else:
                # Just a periodic function, now
                score=self.timing_periodic(mjdnow,meanep,coeffx,slop)
        else:
            # Default has no effect
            score=1

        return score

    ##############################
    
    def score(self,badfilt=[],boostactive=2.0):

        # Update bad-filter list, if necessary
        for ff in badfilt:
            if ff not in self.badfilt:
                self.badfilt.append(ff)

        # Grab the values & limits dictionaries
        values=copy.deepcopy(oscar_values)
        limits=copy.deepcopy(weight_lims)

        # Add source-specific values
        values['moondeg']=self.moondeg()
        values['airmass']=self.airmass()
        values['night']=self.night_value()

        # Airmass range
        if self.minairmass<limits['airmass'][1]:
            limits['airmass'][0]=self.minairmass
        if self.maxairmass<limits['airmass'][1]:
            limits['airmass'][1]=self.maxairmass

        # Nighttime range
        [twi1,twi2]=self.twilight
        limits['night'][1]=24*(twi2-twi1)

        # Start with the raw priority
        score=self.priority

        # Check for violation of hard limits
        hourang=self.hourang()
        if hourang<LIM_HOURANG_EAST or hourang>LIM_HOURANG_WEST:
            score=0

        # Check if we are done
        if self.done:
            score=0

        # Check if we are illegal
        if len(self.property)>0:
            if self.property.has_key('ILLEGAL') and \
                   self.property['ILLEGAL']:
                score=0

        # Check if we have bad-filter exposures
        if self.has_badfilt():
            score=0

        # Partially complete targets have "frozen scores"...
        if self.check_partdone() and score>0:
            if self.lastscore>0:
                score = boostactive*self.lastscore

            # Check only for violation of limits
            for weight_key in self.weight.keys():

                # If we don't have a value and limits, we can't score
                if not values.has_key(weight_key):
                    continue
                if not limits.has_key(weight_key):
                    continue
                scmult = self.weight_function(weight_key,
                                              values[weight_key],
                                              limits[weight_key])
                # Check if multiplier is zero
                if scmult>0:
                    pass
                else:
                    score=0

                # Maintain score record
                self.scores[weight_key]=[values[weight_key],scmult]

        # Full calculation for the rest
        else:

            # Adjust for the various weights
            for weight_key in self.weight.keys():

                # If we don't have a value and limits, we can't score
                if not values.has_key(weight_key):
                    continue
                if not limits.has_key(weight_key):
                    continue
                scmult = self.weight_function(weight_key,
                                              values[weight_key],
                                              limits[weight_key])
                # Apply the multiplier
                if score>0:
                    score *= scmult

                # Maintain full score record
                self.scores[weight_key]=[values[weight_key],scmult]

            # Adjust for timing
            mjdtime=self.mjdnow()
            for time_key in self.timing.keys():
                if time_key=='none':
                    scmult=1
                else:
                    scmult = self.timing_function(time_key,mjdtime)
                if score>0:
                    score *= scmult
                self.scores['time_'+time_key]=[mjdtime,scmult]

            # Save this
            self.lastscore=score

        # Report the score
        return score

    ##############################

    def choose_exposure(self, badfilt=[]):

        # Add specified bad filters to our list, if necessary
        for bf in badfilt:
            if bf not in self.badfilt:
                self.badfilt.append(bf)

        # Default return value -> error condition
        [doi,doj,dok]=[-1,-1,-1]
        
        # Refuse to choose if we're done
        if self.done:
            return [doi,doj,dok]

        # Find the next exposure
        exps=self.exposures
        for i in xrange(len(exps)):
            # Pass over if exposure is done
            if exps[i].done:
                continue
            # Skip if exposure filter is a bad filter
            if exps[i].filter in self.badfilt:
                exps[i].skip()
                continue
            # Skip if exposure subtarget is not visible
            if type(exps[i].target)!=StringType:
                # Subtarget exposure: check for visibility
                if not exps[i].target.alive_now():
                    exps[i].skip()
                    continue
            # Exposure is acceptable - do it!
            doi=i
            break

        if doi<0:
            # No exposure found - probably we are done
            self.check_done()
            return [doi,doj,dok]

        if doi>=0:
            # Get the [j,k] from this exposure
            [doj,dok]=exps[doi].spiral()
        
        if doi<0 or doj<0 or dok<0:
            # Something wrong, couldn't find the right exposure/visit
            sys.stderr.write("Had trouble choosing an exposure\n")
            for xpr in exps:
                sys.stderr.write("   "+xpr.__str__()+"\n")

        return [doi,doj,dok]

    ##############################

    def update(self,i,j,k):

        self.exposures[i].update(j,k)
        self.check_done()

    ##############################

    def has_badfilt(self,badfilt=[]):

        # Update bad-filter list, if necessary
        for bf in badfilt:
            if bf not in self.badfilt:
                self.badfilt.append(bf)

        # Do we have any exposures in any bad filter?
        havebad=0

        if len(self.badfilt)>0:
            for xpr in self.exposures:
                if xpr.filter in self.badfilt:
                    havebad=1
                    break

        return havebad

    ##############################

    def check_partdone(self):

        # Are exposures complete?
        somedone=0
        someundone=0
        for xpr in self.exposures:
            if xpr.done:
                somedone=1
            else:
                someundone=1

        # Return answer
        partdone=(somedone and someundone)
        return partdone

    ##############################

    def donefrac(self):

        # Fraction of our total time done
        fracdone=1.0 - self.time_left()/self.duration
        return fracdone

    ##############################

    def check_done(self):

        # Is obsn>=repeat?
        done=1
        if self.obsn>=self.repeat:
            self.done=done
            return done

        # Are all exposures complete?
        for xpr in self.exposures:
            if not xpr.done:
                done=0
                break

        # If so, we are "done" with one iteration
        if done:
            # Reset individual exposures
            for xpr in self.exposures:
                xpr.reset()
            # Increment our iteration count
            self.obsn += 1
            # Update the last iteration time
            self.obst = self.mjdnow()
            # Check if we are completely done
            done=(self.obsn>=self.repeat)

        # Set property and return answer
        self.done=done
        return done

    ##############################

    def __str__(self):

        strout  = "oscar_target\n"
        strout += "   name:   %s\n" % self.name
        strout += "   priority:  %.2f\n" % self.priority
        strout += "   repeat:    %d\n" % self.repeat
        sxg=self.coosxg()
        strout += "   coords: %s, %s\n" % (sxg[0],sxg[1])
        strout += "   offset: %.1f, %.1f\n" % (self.offset[0],
                                               self.offset[1])
        strout += "   binning:  %d x %d\n" % (self.binning[0],
                                              self.binning[1])
        strout += "   ROI:   %d,%d  %d,%d\n" % (self.roi[0],self.roi[1],
                                                self.roi[2],self.roi[3])
        strout += "   exposures:\n"
        for exp in self.exposures:
            strout += "     %dx %s for %.2fs" % \
                      (exp.imult,exp.filter,exp.exptime)
            if type(exp.target)!=StringType:
                strout += " (Target %s)" % exp.target.name
            strout += "\n"
        strout += "   weights:\n"
        for wt_key in self.weight.keys():
            [cutoff,slope]=self.weight[wt_key]
            strout += "     %s:  %.2f, %.2f\n" % (wt_key,cutoff,slope)
        strout += "   time weights:\n"
        for tm_key in self.timing.keys():
            strout += "     %s:  " % tm_key
            strout += str(self.timing[tm_key])+"\n"
        strout += "   scores:\n"
        strout += "     total:  %.2f\n" % self.score()
        for sc_key in self.scores.keys():
            [scval,scmult]=self.scores[sc_key]
            strout += "     %s:  %.2f -> %.2f\n" % \
                      (sc_key,scval,scmult)

        return strout

######################################################################

class oscar_target_list:

    def __init__(self,catalog="",tnow=None,fake=0,time=0.0,
                 badfilt=[],sunlim=-12.0,rawranks=0,verbose=0,
                 observer=None):

        # Properties
        self.filters=[]
        self.badfilt=[]
        self.sunlim=sunlim
        self.rawranks=rawranks
        self.verbose=verbose

        # Establish tnow object
        if type(tnow)==NoneType:
            self.tnow=time_now(fake=fake,time=time)
        else:
            self.tnow=tnow

        # Oscar target definition defaults
        self.def_priority=1.0
        self.def_repeat=1
        self.def_mjdstart=0.0
        self.def_offset=[]
        self.def_dither=[]
        self.def_binning=[]
        self.def_roi=[]
        self.def_exps=[oscar_exposure('R',60,1)]
        self.def_weight={}
        self.def_timing={}
        self.def_mscnr=0
        self.def_mscnd=0
        self.def_keyword={}
        self.def_property={}

        # Establish observer object
        if type(observer)==NoneType:
            self.observer=palomar
        else:
            self.observer=observer

        # Define the twilight range for tonight
        self.set_twilight()

        # Read in the target catalog, if specified
        self.targets=[]
        self.ranking=[]
        if len(catalog)>0:
            self.read_targets(catalog)

    ##############################

    def parse_target(self,lines,def_priority=None,def_repeat=None,
            def_mjdstart=None,def_mscnr=None,def_mscnd=None,
            def_exps=None,def_offset=None,def_dither=None,
            def_binning=None,def_roi=None,def_weight=None,
            def_timing=None,def_keyword=None,def_property=None):

        """Parse the list of text lines that define a single target;
           return the corresponding target"""

        # Default settings
        is_default=0
        [ra,dec,eqnx,pmra,pmdec,nonsidereal,planet,moving] = \
              [0.0,0.0,2000.0,0.0,0.0,"","",[]]
        owntime=0

        # Arguments
        if type(def_priority)==NoneType:
            def_priority=self.def_priority
        if type(def_repeat)==NoneType:
            def_repeat=self.def_repeat
        if type(def_mjdstart)==NoneType:
            def_mjdstart=self.def_mjdstart
        if type(def_mscnr)==NoneType:
            def_mscnr=self.def_mscnr
        if type(def_mscnd)==NoneType:
            def_mscnd=self.def_mscnd
        if type(def_exps)==NoneType:
            def_exps=self.def_exps
        if type(def_offset)==NoneType:
            def_offset=self.def_offset
        if type(def_dither)==NoneType:
            def_dither=self.def_dither
        if type(def_binning)==NoneType:
            def_binning=self.def_binning
        if type(def_roi)==NoneType:
            def_roi=self.def_roi
        if type(def_weight)==NoneType:
            def_weight=self.def_weight
        if type(def_timing)==NoneType:
            def_timing=self.def_timing
        if type(def_keyword)==NoneType:
            def_keyword=self.def_keyword
        if type(def_property)==NoneType:
            def_property=self.def_property

        # Initialize target parameters
        priority=def_priority
        mjdstart=def_mjdstart
        repeat=def_repeat
        mscnr=def_mscnr
        mscnd=def_mscnd
        mscq=[]

        # Arrays / dictionaries require a deepcopy or new version
        exps=[]
        offset=copy.deepcopy(def_offset)
        dither=copy.deepcopy(def_dither)
        binning=copy.deepcopy(def_binning)
        roi=copy.deepcopy(def_roi)
        weight=copy.deepcopy(def_weight)
        timing=copy.deepcopy(def_timing)
        keyword=copy.deepcopy(def_keyword)
        property=copy.deepcopy(def_property)
        subtargets=[]
        obsn=0
        obst=0.0

        ####################

        # Target name/coords
        line=lines.pop(0)
        els=line.split()

        # Default specification
        if els[0].upper()=='DEFAULT':
            name='DEFAULT'
            is_default=1

        # Nonsidereal target specification
        elif els[0].upper()=='NONSIDEREAL':
            name=els[1]
            rens=re.search('^\S+\s+\S+\s+((\"[^\"]+\")|(\'[^\']+\'))',
                           line,re.I)
            nonsidereal=rens.group(1)[1:-1] # exclude the quotes

        # Planetary target specification
        elif els[0].upper()=='PLANET':
            name=els[1]
            planet=name

        # Moving target specification
        elif els[0].upper()=='MOVING':
            name=els[1]
            if len(els)>4:
                [ra,dec,eqnxs]=els[2:5]
                try:
                    eqnx=float(eqnxs)
                except:
                    sys.stderr.write("Failed to parse Equinox '%s' as float, using default %.2f" % (eqnxs,eqnx))

        # Standard Fixed Target -- Coordinates + Equinox
        else:
            if len(els)>=4:
                [name,ra,dec,eqnxs]=els[0:4]
                try:
                    eqnx=float(eqnxs)
                except:
                    sys.stderr.write("Failed to parse Equinox '%s' as float, using default %.2f" % (eqnxs,eqnx))

            # Proper motions in arcsec / year
            if len(els)>=6:
                [pmra,pmdec]=[float(els[4]),float(els[5])]
            else:
                [pmra,pmdec]=[0.0,0.0]

        # Target parameters
        while lines and re.search('^\s+',lines[0]):

            # Grab a line
            line=lines.pop(0)
            # List of nonwhitespace elements
            els=line.split()
            # Beginning whitespace
            resp=re.search('^(\s+)\S',line)
            spaces=resp.group(1)
            
            if re.search('^Priority',els[0],re.I):
                # Target priority
                priority=float(els[1])

            elif re.search('^Ephem',els[0],re.I):
                # Ephemeris entry for moving target
                moving.append(els[1:])

            elif re.search('^Equi',els[0],re.I):
                # Direct setting of equinox for coordinates
                try:
                    eqnx=float(els[1])
                except:
                    sys.stderr.write("Failed to parse Equinox '%s' as float, using default %.2f" % (els[1],eqnx))
                
            elif re.search('^Offset',els[0],re.I):
                # Offset from center of full CCD
                if len(els)>2:
                    offset=[float(els[1]),float(els[2])]
                else:
                    sys.stderr.write("Offset requested but not specified for %s\n" % name)

            elif re.search('^Dither',els[0],re.I):
                # Size of dither box
                if len(els)>2:
                    dither=[float(els[1]),float(els[2])]
                else:
                    sys.stderr.write("Dither box requested but not specified for %s\n" % name)

            elif re.search('^Bin',els[0],re.I):
                # On-chip binning
                if len(els)>=2:
                    binning=[int(els[1]),int(els[1])]
                else:
                    sys.stderr.write("Binning requested but not legally specified for %s\n" % name)

            elif re.search('^ROI',els[0],re.I):
                # Region of interest to read out
                if len(els)>4:
                    roi=[int(els[1]),int(els[2]),
                         int(els[3]),int(els[4])]
                else:
                    sys.stderr.write("ROI requested but not legally specified for %s\n" % name)

            elif re.search('^Exp',els[0],re.I):
                # Exposure specification
                filter=els[1]
                exptime=float(els[2])
                imult=1
                edither=[]
                # Multiplier
                if len(els)>3:
                    imult=int(els[3])
                # User-defined dither (edither)
                if len(els)>4:
                    for elx in els[4:]:
                        edither.append(elx)
                # Accumulate list of exposures
                exps.append(oscar_exposure(filter,exptime,imult,
                                  jmax=2*mscnr+1,kmax=2*mscnd+1,
                                  obsq=mscq,edither=edither))

            elif re.search('^Target',els[0],re.I):
                # Identify subtarget definition lines
                rest=re.search('^\s+Target\s+(\S.+)$' ,line,re.I)
                tlines=[rest.group(1)]
                while lines and re.search('^%s\s+' % spaces,lines[0]):
                    tlines.append(lines.pop(0))
                # Define the new target object
                subtarget=self.parse_target(tlines,def_priority=priority,
                               def_offset=offset,def_dither=dither,
                               def_binning=binning,def_roi=roi,
                               def_keyword=keyword,def_property=property)
                subtargets.append(subtarget)
                # Convert its exposures to our target's
                for subexp in subtarget.exposures:
                    sub2=copy.deepcopy(subexp)
                    sub2.target=subtarget
                    exps.append(sub2)

            elif re.search('^Weight',els[0],re.I):
                # Parametric weighting specification
                wtkey=els[1].lower()
                if len(els)>3:
                    [wtzero,wtslope]=els[2:4]
                    weight[wtkey]=[float(wtzero),float(wtslope)]
                else:
                    sys.stderr.write("Weight type %s requested but not specified for %s\n" % (wtkey,name))

            elif re.search('^Tim',els[0],re.I):
                # Timing specification
                tmwt=els[1]
                tmwt=tmwt.lower()
                tmvals=[]
                for elx in els[2:]:
                    tmvals.append(elx)
                if not owntime:
                    # The first timing weight-setting replaces the
                    # default time weightings
                    timing={}
                    owntime=1
                timing[tmwt]=tmvals

            elif re.search('^Mosaic',els[0],re.I):
                # Mosaic specification
                if len(els)>2:
                    mscnr=int(els[1])
                    mscnd=int(els[2])
                    if len(els)>3:
                        mscq=els[3].split(',')
                else:
                    sys.stderr.write("Mosaic requested but not specified for %s\n" % name)

            elif re.search('^Key',els[0],re.I):
                # Header keyword fields
                if re.search('^[\"\']',els[2],re.I):
                    # Multi-word keyword value
                    rek=re.search('Key\S+\s+\S+\s+((\"[^\"]+\")|(\'[^\']+\'))',line,re.I)
                    if rek:
                        keyword[els[1]]=rek.group(1)
                    else:
                        print "Bad multiword keyword for %s" % els[1]
                else:
                    # Remaining properties are single words
                    keyword[els[1]]=els[2]

            elif re.search('^Prop',els[0],re.I):
                # Arbitrary integer properties
                if len(els)>2:
                    property[els[1].upper()]=int(els[2])
                else:
                    sys.stderr.write("Property %s designated but not set for %s\n" % (els[1],name))

            elif re.search('^MJDSta',els[0],re.I):
                # MJD of target entry into the queue
                mjdstart=float(els[1])

            elif re.search('^Repeat',els[0],re.I):
                # Requested number of iterations for the target
                repeat=int(els[1])

            elif re.search('^Obsn',els[0],re.I):
                # Number of iterations done so far
                obsn=int(els[1])

            elif re.search('^Obst',els[0],re.I):
                # MJD of the last iteration
                obst=float(els[1])

            else:
                sys.stderr.write("Skipped line for target %s:  '%s'\n" %
                                 (name,line))
                pass

        ####################

        # No null exposures (breaks stuff)
        if len(exps)==0:
            exps=copy.deepcopy(def_exps)

        # Reset DEFAULT parameters for default targets
        if is_default:
            self.def_priority=priority
            self.def_mjdstart=mjdstart
            self.def_repeat=repeat
            self.def_mscnr=mscnr
            self.def_mscnd=mscnd

            # Arrays / dictionaries require a deepcopy
            self.def_offset=copy.deepcopy(offset)
            self.def_dither=copy.deepcopy(dither)
            self.def_binning=copy.deepcopy(binning)
            self.def_roi=copy.deepcopy(roi)
            self.def_exps=copy.deepcopy(exps)
            self.def_weight=copy.deepcopy(weight)
            self.def_timing=copy.deepcopy(timing)
            self.def_keyword=copy.deepcopy(keyword)
            self.def_property=copy.deepcopy(property)

            target="Default"

        else:

            target=oscar_target(name,ra,dec,exps,
                       subtargets=subtargets,equinox=eqnx,
                       pmra=pmra,pmdec=pmdec,priority=priority,
                       offset=offset,dither=dither,binning=binning,
                       roi=roi,nr=mscnr,nd=mscnd,mscq=mscq,weight=weight,
                       timing=timing,keyword=keyword,property=property,
                       mjdstart=mjdstart,repeat=repeat,obsn=obsn,obst=obst,
                       nonsidereal=nonsidereal,planet=planet,moving=moving,
                       tnow=self.tnow,sunlim=self.sunlim,
                       twilight=self.twilight,badfilt=self.badfilt,
                       observer=self.observer)

        return target

    ##############################

    def read_targets(self,catalog,reload=0,debug=0):

        if not os.path.exists(catalog):
            print "Failed to find catalog file %s for reading." % catalog
            return

        if reload:
            self.targets=[]
        
        lines=getlines(catalog)
        while lines:
            line=lines.pop(0)
            # Comment line
            if re.search('^#',line):
                continue
            # Blank line
            if not re.search('\S+',line):
                continue
            # Beginning a target sequence
            tlines=[line]
            while lines and re.search('^\s+',lines[0]):
                tlines.append(lines.pop(0))
            # Define the new target object
            if debug:
                #print tlines[0]
                print tlines
            target=self.parse_target(tlines)
            # Skip default targets, but keep the rest
            if type(target)!=StringType:
                self.targets.append(target)

        # Sort in order of decreasing priority
        self.sort_targets()

        # Fill the target rankings
        self.rank_targets()

        # Update filter list
        self.update_filters()

    ##############################

    def express_target(self,target,indent="",parent=None):

        """Convert target to a target specification"""

        # Useful shorthands
        t=target
        p=parent
        isparent=(type(p)==NoneType)

        # Accumulate the specification as a list of strings
        lines=[]

        # Name + Nonsidereal specification
        if not isparent:
            namepfx=indent+"Target  "
        else:
            namepfx=indent
            
        if len(t.nonsidereal)>0:
            # Nonsidereal target
            lines.append(namepfx+'Nonsidereal %s  "%s"\n' % \
                         (t.name,t.nonsidereal))
        elif len(t.planet)>0:
            # Planetary target
            lines.append(namepfx+"Planet %s\n" % t.name)
        elif len(t.moving)>0:
            # Moving target
            [tradeg,tdcdeg]=t.coodeg()
            teqnx=t.equinox
            lines.append(namepfx+"Moving %s  %.6f %.6f  %.2f\n" % \
                         (t.name,tradeg,tdcdeg,teqnx))
        else:
            # Name RA Dec Equinox (PMRA PMDEC)
            [tradeg,tdcdeg]=t.coodeg()
            teqnx=t.equinox
            line=namepfx + "%s  %.6f %.6f  %.2f" % \
                          (t.name,tradeg,tdcdeg,teqnx)
            if abs(t.pmra)>0 or abs(t.pmdec)>0:
                line += "%.5f %.5f" % (t.pmra,t.pmdec)
            lines.append(line+"\n")

        # Indent string
        indent += "    "

        # Ephemeris info for moving targets
        if len(t.moving)>0:
            for move in t.moving:
                line=indent+"Ephemeris %s %s %s" % \
                             (move[0],move[1],move[2])
                if len(move)>=5:
                    line += (" %s %s" % (move[3],move[4]))
                lines.append(line+"\n")

        # Target basics
        if isparent:

            # Priority & Repeat count
            lines.append(indent+"Priority %.2f\n" % t.priority)
            lines.append(indent+"Repeat %d\n" % t.repeat)

            # Observation tracking
            lines.append(indent+"Obsnum %d\n" % t.obsn)
            lines.append(indent+"Obstime %.6f\n" % t.obst)

            # MJD Entry into queue
            lines.append(indent+"MJDStart %.6f\n" % t.mjdstart)

            # Weighting criteria
            for wt in t.weight.keys():
                lines.append(indent+"Weight %s %.2f %.2f\n" %
                             (wt,t.weight[wt][0],t.weight[wt][1]))

            # Time-weighting
            for wt in t.timing.keys():
                line=indent+"Timing %s" % wt
                for x in t.timing[wt]:
                    line += " %s" % x
                lines.append(line+"\n")

        # Mosaicking (not inherited)
        if t.nr>0 or t.nd>0:
            line="Mosaic %i %i" % (t.nr,t.nd)
            if len(t.mscq)>0:
                line += "  "+','.join(t.mscq)
            lines.append(indent+line+"\n")

        # Inherited properties
        if isparent or t.offset!=p.offset:
            lines.append(indent+"Offset %.1f %.1f\n" %
                         (t.offset[0],t.offset[1]))
        if isparent or t.dither!=p.dither:
            lines.append(indent+"Dither %.1f %.1f\n" %
                         (t.dither[0],t.dither[1]))
        if isparent or t.binning!=p.binning:
            lines.append(indent+"Binning %d %d\n" %
                         (t.binning[0],t.binning[1]))
        if isparent or t.roi!=p.roi:
            lines.append(indent+"ROI %d %d %d %d\n" %
                         (t.roi[0],t.roi[1],t.roi[2],t.roi[3]))

        # Header keywords
        if len(t.keyword.keys())>0:
            for k in t.keyword.keys():
                if isparent or k not in p.keyword.keys() or \
                   t.keyword[k]!=p.keyword[k]:
                    lines.append(indent+"Keyword %s %s\n" % (k,t.keyword[k]))

        # Additional properties
        if len(t.property.keys())>0:
            for k in t.property.keys():
                if isparent or k not in p.property.keys() \
                   or t.property[k]!=p.property[k]:
                    lines.append(indent+"Property %s %d\n" % (k,t.property[k]))

        # Exposures
        newt=None
        for exp in t.exposures:
            if type(exp.target)!=StringType:
                # This exposure belongs to a subtarget
                if type(newt)==NoneType or exp.target!=newt:
                    # We've got a new subtarget
                    newt=exp.target
                    newlines=self.express_target(newt,indent=indent,
                                                 parent=t)
                    for newline in newlines:
                        lines.append(newline)
                else:
                    # This exposure is associated with a subtarget
                    # that has already been processed, so skip it  
                    continue
            else:
                # This exposure belongs to the parent target
                line=indent+"Exposure %s %.2f %d" % \
                          (exp.filter,exp.exptime,exp.imult)
                if len(exp.edither)>0:
                    for edith in exp.edither:
                        line += " %s" % edith
                lines.append(line+"\n")
                # Clear the new-target setting
                newt=None

        return lines

    ##############################

    def write_targets(self,catalog):

        # Sort in order of decreasing priority
        self.sort_targets()

        if os.path.exists(catalog):
            if self.verbose:
                print "Target catalog %s will be overwritten" % catalog
            os.remove(catalog)

        catout=open(catalog,"w")

        for t in self.targets:
            tlines=self.express_target(t)
            for line in tlines:
                catout.write(line)

        catout.close()

    ##############################

    def sort_targets(self):

        dctgt={}
        for t in self.targets:
            prio=t.priority
            if prio in dctgt.keys():
                dctgt[prio].append(t)
            else:
                dctgt[prio]=[t]

        dkeys=dctgt.keys()
        dkeys.sort()
        dkeys.reverse()

        newt=[]
        for prio in dkeys:
            for t in dctgt[prio]:
                newt.append(t)

        self.targets=newt

    ##############################

    def rank_targets(self,boostactive=2.0):

        # Fix boostactive=1.0 for "rawranks"
        if self.rawranks:
            boostactive=1.0

        sctgt=[]
        for i in xrange(len(self.targets)):
            t=self.targets[i]
            tscore=t.score(boostactive=boostactive)
            # Here is where rawranks is applied
            if self.rawranks and tscore>0:
                tscore=t.priority
            sctgt.append([tscore,i])

        # Out of targets!
        if len(sctgt)<1:
            return -1

        # Sorted list
        sctgt.sort()
        sctgt.reverse()

        # Save the ranked list
        self.ranking=sctgt

        # Uh-oh, we didn't find any nonzero scores
        if sctgt[0][0]==0:
            return -1

        # Return index of the top-ranked target
        topix=sctgt[0][1]
        return topix

    ##############################

    def set_twilight(self):
        self.twilight=self.tnow.twilight(self.observer,
                                         sunalt=self.sunlim)

    def update_filters(self):

        # Gather list of filters from target list
        filts={}
        for t in self.targets:
            if not t.done:
                for exp in t.exposures:
                    filts[exp.filter]=1

        # Condense into a list
        self.filters=filts.keys()

    def list_of_filters(self):
        return self.filters

    def update(self):
        self.sort_targets()
        self.update_filters()

    ##############################

    # Give the object list-like behavior

    def index(self,value):
        return self.targets.index(value)

    def count(self,value):
        return self.targets.count(value)

    def append(self,value):
        self.targets.append(value)
        self.update()

    def __len__(self):
        return len(self.targets)

    def __getitem__(self,ix):
        return self.targets[ix]

    def __setitem__(self,ix,value):
        self.targets[ix]=value
        self.update()

    def __delitem__(self,ix):
        del self.targets[ix]
        self.update()

    def __getslice__(self,ix,jx):
        return self.targets[ix:jx]

    def __setslice__(self,ix,jx,values):
        self.targets[ix:jx]=values
        self.update()

    def __delslice__(self,ix,jx):
        self.targets[ix:jx]=[]

    def __concat__(self,s2):
        self.targets=self.targets.append(s2)
        self.update()

    def __contains__(self,element):
        return (element in self.targets)

