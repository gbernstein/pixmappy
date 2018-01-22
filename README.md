# pixmappy
Python interface to gbdes pixel map (astrometry) solutions

`PixelMapCollection` class can read the YAML format astrometry solutions produced by Gary's `WcsFit` program (in gbdes repo).  This class can then issue a `PixelMap` instance, which is a map from one 2d coordinate system ("pixel") to another ("world") 2d system.  A `PixelMap` instance can be used as a function mapping one (or many) coordinate pairs.  An `inverse` method does reverse mapping, and the local `jacobian` of the map is available also.

A `WCS` is a `PixelMap` that additionally specifies a (de)projection from the world 2d system onto the celestial sphere.  It uses `astropy.coordinates.SkyCoord` objects to represent celestial positions, and has `toSky` and `toPix` methods to map between the pixel system and `SkyCoord`s.  

The YAML deserialization is the slow part, using `PyYAML`.  `PixelMapCollection` will first try to use the C-based YAML loader, which requires libyaml to be available.  If this fails it falls back to the 20x slower all-Python loader.

The type of mapping that can be expressed is very flexible, and `PixelMaps` can be compounded into chains of tranformations.  See the `gbdes` source files for more documentation on the types of transformations available and their YAML encodings.  One type of map is a `TemplateMap`, which uses a lookup table stored in some other YAML file.  Two such template files used for DECam solutions are part of this repository and installed by setup.py.  The setup area will be searched automatically for requested template files.  Additional paths to search can be specified by the `CAL_PATH` environment variable, using the usual colon-separated list format.

## DES Astrometric Solutions

The `DESMaps` class derives from `PixelMapCollection` and is specialized to read astrometric solutions derived for all of the useful exposures in the Y4A1 internal data release.  Upon creation of an instance of this class, some YAML and FITS files containing WCS parameters for all these exposures are read.  The user can then request production of a `WCS` appropriate to any combination of exposure number and focal-plane detector.  Quick instructions for doing so are as follows:

* Acquire this repository and run `python setup.py install`
* Acquire and unpack the `y4a1_astrometry.tar.gz` file containing the solution information.
* Make sure that the environment variable `CAT_PATH` contains the directory into which these data were placed.
* Run your python code!
Here is an example program:
```python
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
wcs = maps.getDESWCS(expnum, detpos)

# Map a single pair of pixel coordinates to RA, Dec, specifying object color
x = 1088.4
y = 3115.2
c = 1.4
ra,dec = wcs.toSky(x,y,c)

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
```
    
