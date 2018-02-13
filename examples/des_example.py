#!/usr/bin/env python
'''
Example program for use with value-added DES astrometry solutions
'''
from __future__ import print_function
import os
import numpy as np
from pixmappy import DESMaps, Gnomonic

# Show that we are using CAL_PATH to find files
try:
    print('Reading files on path', os.environ['CAL_PATH'])
except KeyError:
    print('Set CAL_PATH environment variable to include directories with astrometric solutions')

# Read solution set from default files.
# Optional argument inhibits use/creation of python pickle files.
maps = DESMaps(use_pkl=False)

# Acquire the WCS for desired exposure, CCD combination
# ValueError will be raised if there is no astrometric solution for this combination.
expnum = 514157
detpos = 'N8'
ccdnum = 39
wcs = maps.getDESWCS(expnum, detpos)

# Map a single pair of pixel coordinates to RA, Dec, specifying object color
x = 1088.4
y = 3115.2
c = 1.4
ra,dec = wcs.toSky(x,y,c)

print('x/y',x,y,'map to RA/Dec:',ra,dec)

# We can also request a WCS by the CCDNUM of the device
wcs2 = maps.getDESWCS(expnum, ccdnum)
print('With CCDNUM:')
ra,dec = wcs2.toSky(x,y,c)
print('x/y',x,y,'map to RA/Dec:',ra,dec)

# Map these back
print('Inversion yields', wcs.toPix(ra,dec,c))

# Map a set of coordinates with varying colors.
# c could be a scalar and would be broadcast
x = np.arange(900.,1500.,100.)
y = x + 1000.
c = np.arange(0.9,1.5,0.1)
ra,dec = wcs.toSky(x,y,c)
print('Array outputs:')
print('RA:',ra)
print('Dec:',dec)

# Now we're going to ask the WCS to reproject coordinates
# into a gnomonic projection about a new point, in
# units of degrees

ra0 = 37.
dec0 = -25.
wcs.reprojectTo( Gnomonic(ra0,dec0) )

# And we just use the wcs as a function object to map pixel
# coordinates into projection-plane coordinates:
x = 555.
y = 4321.
c = 0.5
print('xi,eta are', wcs(x,y,c))

# We can get the local Jacobian of the map from pixel
# to project coordinates
print('Jacobian:',wcs.jacobian(x,y,c))

# Now get the estimated covariance matrix of astrometric
# errors for this exposure.
cov = maps.getCovariance(expnum)
print('Covariance matrix:',cov)
# Is there a warning about quality of this covariance?
print('Warning for covariance?',maps.covarianceWarning(expnum))

