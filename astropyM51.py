from astropy import units as u
from astropy.coordinates import SkyCoord
"""
http://astropy.readthedocs.org/en/latest/coordinates/index.html

estimate the linear distance from the center of M51 
to the center of region in Corder et al 2008
"""


gc = SkyCoord('13h29m56.2s', '+47d13m50s', frame='icrs')
r1 = SkyCoord('13h29m50.381s','+47d12m02.28s', frame='icrs')

distance = gc.separation(r1).arcsec 

print "separation in arcsec = ", distance
print "separation in kpc = ", distance*0.046