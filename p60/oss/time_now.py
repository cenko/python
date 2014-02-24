
import math
import ephem
from types import *

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

# Just a useful default
MAX_AIRMASS=20.0

######################################################################

class time_now:

    def __init__(self,fake=0,time=0.0):

        self.mjdepoch=MJDEPHEP
        self.unixepoch=MJDUNIXEP
        self.twopi=TWOPI
        self.dtor=DTOR
        self.radsperhour=RADSPERHOUR
        self.secsperday=SECSPERDAY
        self.siderealday=SIDEREALDAY
        self.fake=(fake!=0)
        self.fakenow=0.0
        self.sun=ephem.Sun()
        self.moon=ephem.Moon()
        self.razero=ephem.readdb("Meridian,f|S,0.0,33.3560,15.0,2000.00")
        
        if self.fake:
            if time==0.0:
                self.fakenow=ephem.now()
            else:
                self.fakenow=time

    ##############################

    def now(self):

        if self.fake:
            now=self.fakenow
        else:
            now=ephem.now()

        return now

    def makefake(self,now=None):

        if type(now)==NoneType:
            now=self.now()

        self.fakenow=now
        self.fake=1

    def makereal(self):

        self.fake=0

    ##############################

    def ephem(self):

        return self.now()

    def date(self):

        return ephem.date(self.now())

    def mjd(self,now=None):

        if type(now)==NoneType:
            now=self.now()

        return self.mjdepoch+self.now()

    def unix(self,now=None):

        if type(now)==NoneType:
            now=self.now()

        return self.secsperday * \
               (self.mjdepoch+self.now()-self.unixepoch)

    ##############################

    def lst(self,observer,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        lst=float(observer.sidereal_time())/self.radsperhour
        return lst

    def sunpos(self,observer,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        self.sun.compute(observer)
        return [self.sun.az/self.dtor,self.sun.alt/self.dtor]

    def sunaz(self,observer,now=None):
        [sunaz,sunalt]=self.sunpos(observer,now=now)
        return sunaz

    def sunalt(self,observer,now=None):
        [sunaz,sunalt]=self.sunpos(observer,now=now)
        return sunalt

    def sunradec(self,observer,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        self.sun.compute(observer)
        return [str(self.sun.a_ra),str(self.sun.a_dec)]

    ##############################

    def tonight(self,observer,now=None,verbose=1):

        if type(now)==NoneType:
            now=self.now()

        try:
            # Sunset first -- If sun is down, want previous setting
            if self.sunalt(observer,now=now)<0:
                sunset=observer.previous_setting(self.sun, start=now)
            # Otherwise want next setting
            else:
                sunset=observer.next_setting(self.sun, start=now)

            # Sunrise calc -- Always want next rising
            sunrise=observer.next_rising(self.sun, start=now)

            # Done
            return [sunset,sunrise]
        
        except ephem.CircumpolarError:
            if verbose:
                print "No sunset or sunrise at this time of year"
            return [None, None]

    ##############################

    def twilight(self,observer,sunalt=-12.0,now=None,verbose=1):

        if type(now)==NoneType:
            now=self.now()

        # We want the times between evening and morning twilight
        # corresponding to the specified sun-altitude

        try:

            # Modify the horizon setting
            observer.horizon=sunalt*self.dtor

            # If the sun is down, we want previous setting time
            if self.sunalt(observer,now=now)<sunalt:
                twilight1=observer.previous_setting(self.sun, start=now)
            # Otherwise want next setting
            else:
                twilight1=observer.next_setting(self.sun, start=now)

            # End of twilight calc -- Always want next rising
            twilight2=observer.next_rising(self.sun, start=now)

        except ephem.CircumpolarError:
            if verbose:
                print "No sunset or sunrise at this time of year"
            [twilight1,twilight2] = [None, None]
        
        # Reset the observer horizon
        observer.horizon=0

        # Done
        return [twilight1,twilight2]

    ##############################

    def moonpos(self,observer,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        self.moon.compute(observer)
        return [self.moon.az/self.dtor,self.moon.alt/self.dtor]

    def moonradec(self,observer,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        self.moon.compute(observer)
        return [str(self.moon.a_ra),str(self.moon.a_dec)]

    def moonphase(self,observer,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        self.moon.compute(observer)
        return self.moon.phase

    def moonisup(self,observer,now=None):

        [az,alt]=self.moonpos(observer,now=now)
        return alt>0

    def moonup(self,observer,now=None):

        if type(now)==NoneType:
            now=self.now()

        # We're interested in the time range when the moon is up
        # during "tonight"
        [sunset,sunrise]=self.tonight(observer,now=now)
        if sunset==None or sunrise==None:
            print "Error: Cannot calculate sun rise or set"
            return [None, None]
 
        try:
            # If moon is up now, want previous rise time
            if self.moonisup(observer,now=now):
                moonstart=observer.previous_rising(self.moon,start=now)
            # Otherwise want next rise time
            else:
                moonstart=observer.next_rising(self.moon,start=now)

            # Always want next setting
            moonstop=observer.next_setting(self.moon,start=now)

            # Adjust based on sun
            if moonstop<sunset or moonstart>sunrise:
                print "Moon rises only during sunlight hours"
                return [None, None]
            if moonstart<sunset:
                moonrise=sunset
            else:
                moonrise=moonstart
            if moonstop>sunrise:
                moonset=sunrise
            else:
                moonset=moonstop

            return [moonrise, moonset]

        except ephem.CircumpolarError:
            print "Error: Moon never rises or sets"
            return [None, None]

    ##############################

    def lstrange(self,observer,sunalt=0,now=None):

        """What is the LST at sunset/twilight and sunrise/twilight?"""

        if sunalt==0:
            [sunset,sunrise]=self.tonight(observer,now=now)
        else:
            [sunset,sunrise]=self.twilight(observer,sunalt=sunalt,now=now)
            
        lstset=self.lst(observer,now=sunset)
        lstrise=self.lst(observer,now=sunrise)

        return [lstset,lstrise]

    def lstmid(self,observer,sunalt=0,now=None):

        """What is the LST at the middle of the night?"""

        [lst1,lst2]=self.lstrange(observer,sunalt=sunalt,now=now)
        lstmid=ephem.hours(0.5*(lst1+lst2)*RADSPERHOUR)

        return lstmid

    ##############################

    def alt_times(self,observer,object,alt,
                  now=None,antipodes=0,verbose=0):

        """When does object pass through a critical altitude?"""

        if type(now)==NoneType:
            now=self.now()

        # Calculate the two times bracketing the current transit when
        # object will pass through the altitude alt (in degrees)
        try:

            # Set horizon
            observer.horizon=alt*self.dtor

            # If object is up, use previous rising
            if self.objisup(observer,object,alt=alt,now=now):
                t1=observer.previous_rising(object,start=now)
            # Otherwise use next rising
            else:
                t1=observer.next_rising(object,start=now)

            # Always want next setting
            t2=observer.next_setting(object,start=now)

        except ephem.CircumpolarError:

            print "Object will not reach that altitude"
            [t1,t2]=[None,None]

        # Fix horizon and return
        observer.horizon=0
        return [t1,t2]

    ####################################

    def objisup(self,observer,object,alt=0.0,now=None):

        """Is object above a certain altitude?"""

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        object.compute(observer)
        
        return object.alt/self.dtor>alt

    def objinsky(self,observer,object,now=None,verbose=0):

        """When is object above the horizon?"""

        if type(now)==NoneType:
            now=self.now()

        [objrise,objset]=self.alt_times(observer,object,0,now=None)

        # Done
        return [objrise,objset]

    def objup(self,observer,object,alt=0,sunalt=-12.0,now=None,verbose=0):

        """When is object up with a night sky?"""

        if type(now)==NoneType:
            now=self.now()

        # Determine dark time
        [twi1,twi2]=self.twilight(observer,sunalt=sunalt,now=now)

        # Error check
        if type(twi1)==NoneType:
            if verbose:
                print "Sun doesn't reach requested altitude, sorry"
            return [None,None]

        # Determine when object is in sky
        [pass1,pass2]=self.alt_times(observer,object,alt,now=twi1)

        # Error check
        if type(pass1)==NoneType:
            if verbose:
                print "Object never reaches critical altitude, sorry"
            return [None,None]

        # If there's still no overlap, source is not visible
        if pass2<twi1 or pass1>twi2:
            if verbose:
                print "Condition not satisfied during night hours"
            return [None,None]

        # Get the overlap times
        if pass1>twi1:
            objstart=pass1
        else:
            objstart=twi1

        if pass2<twi2:
            objstop=pass2
        else:
            objstop=twi2

        # Done
        return [objstart,objstop]

    ####################

    def objsxg(self,observer,object,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        object.compute(observer)

        return [str(object.a_ra),str(object.a_dec)]

    def objdeg(self,observer,object,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        object.compute(observer)

        return [object.a_ra/self.dtor,
                object.a_dec/self.dtor]

    def objrad(self,observer,object,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        object.compute(observer)

        return [object.a_ra*1.0,object.a_dec*1.0]

    def objradeg(self,observer,object,now=None):
        
        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        object.compute(observer)

        return object.a_ra/self.dtor

    def objdecdeg(self,observer,object,now=None):
        
        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        object.compute(observer)

        return object.a_dec/self.dtor

    def altitude(self,observer,object,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        object.compute(observer)

        return object.alt/self.dtor

    def azimuth(self,observer,object,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        object.compute(observer)

        return object.az/self.dtor

    def airmass(self,observer,object,now=None,max=MAX_AIRMASS):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        object.compute(observer)

        if object.alt<0:
            airmass=max
        elif object.alt<10.0*self.dtor:
            airmass=6.0
        else:
            airmass=1/math.sin(object.alt)

        return airmass

    def moondeg(self,observer,object,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        self.moon.compute(observer)
        object.compute(observer)
        moondeg=180.0-ephem.separation([self.moon.a_ra,self.moon.a_dec],
                                       [object.a_ra,object.a_dec])/self.dtor

        return moondeg

    def hourang(self,observer,object,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        object.compute(observer)

        hourang=24*(observer.date - 
                    observer.next_transit(object,start=now))/self.siderealday
        if abs(hourang)>12:
            hourang -= 24*cmp(hourang,0)

        return hourang

    def transit_time(self,observer,object,now=None):

        if type(now)==NoneType:
            now=self.now()

        observer.date=now
        object.compute(observer)
        transit_time=observer.next_transit(object,start=now)

        return transit_time

    ##############################

    def __str__(self):

        s = "%.5f" % (self.mjdepoch+self.now())
        return s

    def __repr__(self):

        if self.fake:
            s="Fake"
        else:
            s="Real"

        s += " time_now object at " + str(ephem.date(self.now())) + " UT"

        return s

    def __add__(self,other):

        if self.fake:
            self.fakenow += other

        return self

    def __sub__(self,other):

        if self.fake:
            self.fakenow -= other

        return self

    def __iadd__(self,other):

        if self.fake:
            self.fakenow += other

        return self

    def __isub__(self,other):

        if self.fake:
            self.fakenow -= other

        return self

    def __float__(self):

        return self.now()

    def __int__(self):

        return int(self.now())

    def __long__(self):

        return long(self.now())

######################################################################

