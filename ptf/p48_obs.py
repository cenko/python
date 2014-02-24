#! /usr/bin/env python

import ephem, sys
from mx.DateTime import *
from iqutils import *

def main(ra, dec):

	targ = astrocoords(float(ra), float(dec))

	palomar=ephem.Observer()
	palomar.long, palomar.lat, palomar.elev = '-116.8630', '33.3560', 1706.0
	palomar.date = "2013/12/1"

	try:
	 	sun=ephem.Sun()
		sun.compute(palomar)
		palomar.horizon = "-12:00"
		tm = palomar.next_rising(sun)
		te = palomar.next_setting(sun)

		palomar.horizon = "30:00"
		etarg = ephem.readdb("Null,f|S,%s,%s,15.0,2000.0" % (targ.sxg()[0], targ.sxg()[1]))
		etarg.compute(palomar)
		tr = palomar.next_rising(etarg)
		ts = palomar.next_setting(etarg)

		if te < tr < ts < tm:
			dt1 = ts - tr
		elif (te < ts < tm < tr) or (tr < te < ts < tm):
			dt1 = ts - te
		elif (te < tr < tm < ts) or (ts < te < tr < tm):
			dt1 = tm - tr
		elif (tr < te < tm < ts) or (te < tm < ts < tr) or (ts < tr < te < tm):
			dt1 = tm - te
		elif (te < ts < tr < tm):
			dt1 = (ts - te) + (tm - tr)
		else:
			dt1 = 0

		palomar.date = "2014/1/1"
		etarg.compute(palomar)
		tr = palomar.next_rising(etarg)
		ts = palomar.next_setting(etarg)
		sun.compute(palomar)
		palomar.horizon = "-12:00"
		tm = palomar.next_rising(sun)
		te = palomar.next_setting(sun)

		if te < tr < ts < tm:
			dt2 = ts - tr
		elif (te < ts < tm < tr) or (tr < te < ts < tm):
			dt2 = ts - te
		elif (te < tr < tm < ts) or (ts < te < tr < tm):
			dt2 = tm - tr
		elif (tr < te < tm < ts) or (te < tm < ts < tr) or (ts < tr < te < tm):
			dt2 = tm - te
		elif (te < ts < tr < tm):
			dt2 = (ts - te) + (tm - tr)
		else:
			dt2 = 0
			
	except (ephem.NeverUpError,ephem.AlwaysUpError) as e:
		return [0.0, 0.0]
		
	#print "%10.4f%10.4f%10.2f%10.2f" % (float(ra), float(dec), dt1*24.0,dt2*24.0)
	return [dt1, dt2]
	
if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])
