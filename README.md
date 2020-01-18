# pixmappy
Python interface to gbdes pixel map (astrometry) solutions

`PixelMapCollection` class can read the YAML format astrometry solutions produced by Gary's `WcsFit` program (in gbdes repo).  This class can then issue a `PixelMap` instance, which is a map from one 2d coordinate system ("pixel") to another ("world") 2d system.  A `PixelMap` instance can be used as a function mapping one (or many) coordinate pairs.  An `inverse` method does reverse mapping, and the local `jacobian` of the map is available also.

A `WCS` is a `PixelMap` that additionally specifies a (de)projection from the world 2d system onto the celestial sphere.  It uses `astropy.coordinates.SkyCoord` objects to represent celestial positions, and has `toSky` and `toPix` methods to map between the pixel system and `SkyCoord`s.  

The YAML deserialization is the slow part, using `PyYAML`.  `PixelMapCollection` will first try to use the C-based YAML loader, which requires libyaml to be available.  If this fails it falls back to the 20x slower all-Python loader.

The type of mapping that can be expressed is very flexible, and `PixelMaps` can be compounded into chains of tranformations.  See the `gbdes` source files for more documentation on the types of transformations available and their YAML encodings.  One type of map is a `TemplateMap`, which uses a lookup table stored in some other YAML file.  Two such template files used for DECam solutions are part of this repository and installed by setup.py.  The setup area will be searched automatically for requested template files.  Additional paths to search can be specified by the `CAL_PATH` environment variable, using the usual colon-separated list format.

## DES Astrometric Solutions

The `DESMaps` class derives from `PixelMapCollection` and is specialized to read astrometric solutions derived for all of the useful exposures in the Y6A1 internal data release.  Upon creation of an instance of this class, some YAML and FITS files containing WCS parameters for all these exposures are read.  The user can then request production of a `WCS` appropriate to any combination of exposure number and focal-plane detector.  Quick instructions for doing so are as follows:

* Acquire this repository and run `python setup.py install`
* The Y6A1 astrometric solutions are included in the data directory of this repo
  and will be accessed by default by `DESMaps` (If you are a DES member you can 
  also find these files at https://cdcvs.fnal.gov/redmine/projects/des-y6/wiki/Y6A1_Astrometric_Solutions).
* If you are using some other set of solutions, 
  make sure that the environment variable `CAL_PATH` contains the
  directory into which these data were placed (_e.g._ /xxx/ALTERNATE_ASTROMETRY).
* Run your python code!
* Note that the color argument `c` is assumed to be _g-i_ for the DES
  data.  A value is required, since the solutions include differential
  chromatic refraction in the atmosphere and (for _gr_ bands) lateral
  color in the corrector.  Use a value of `c=0.61` if you don't know
  your true color and want something that is not crazy.
* If you request a solution for an exposure/CCD pair that is not in
  the solution set, a `ValueError` exception will be raised.

### Use in GalSim

The `GalSimWCS` class allows `pixmappy` astrometric solutions to be used 
within `GalSim`. `GalSimWCS` will read a solution any one of three ways:
* Giving a `PixelMapCollection` and the name of a WCS within that collection.
* Giving a `yaml_file` which encodes a `PixelMapCollection`, plus the name of a WCS
  within that collection.
* Specifying `use_DESMaps=True`, which triggers use of the `DESMaps` class 
  and file-access methods described above, plus the exposure number and ccd number of
  the desired astrometric map.  Choosing this option will by default
  access the Y6A1_ASTROMETRY solutions included in the repo, and is the simplest
  way to use the class.  

### Astrometric error estimation

The astrometric solutions above reduce uncertainties in the
instrumental map to <3 mas RMS.  The dominant errors that remain are
(1) errors in the object centroid due to image noise, and (2)
stochastic astrometric distortions due to atmospheric turbulence.
The measurement errors (1) are typically taken from the
`ERRAWIN_IMAGE` column of `SExtractor` catalogs, since all of the
astrometric solutions are referenced to `[XY]WIN_IMAGE` centroids.

The atmospheric errors (2) are anisotropic and differ from exposure to
exposure as weather conditions vary.  We have estimated the covariance
matrix of these efforts from the residuals to the fits used to
establish the solutions.
_[In detail: these are estimated from the 2-pt correlation function of the astrometric residuals, which eliminates terms such as measurement noise which should not correlate between distinct stars.]_
A call to `DESMaps.getCovariance` will return the 2x2 covariance
matrix estimated for a given exposure.  You should add to this the
(usually diagonal) covariance matrix for measurement noise to obtain
an estimate of the total error on a given source's position.  Some of
these covariance estimates are questionable (usually due to too few
stars used to estimate them).  These can be identified by calling
`DESMaps.covarianceWarning().`

A caveat is that these covariance matrices are estimates, and noisy
ones at that.  Some may be defective but not flagged by the warning,
so be cautious in your reliance on their accuracy.

### Example program

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

```
    
