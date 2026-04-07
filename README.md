# pixmappy
Python interface to gbdes pixel map (astrometry) solutions

`PixelMapCollection` class can read the YAML format astrometry solutions produced by Gary's `WcsFit` program (in gbdes repo).  This class can then issue a `PixelMap` instance, which is a map from one 2d coordinate system ("pixel") to another ("world") 2d system.  A `PixelMap` instance can be used as a function mapping one (or many) coordinate pairs.  An `inverse` method does reverse mapping, and the local `jacobian` of the map is available also.

A `WCS` is a `PixelMap` that additionally specifies a (de)projection from the world 2d system onto the celestial sphere.  It uses `astropy.coordinates.SkyCoord` objects to represent celestial positions, and has `toSky` and `toPix` methods to map between the pixel system and `SkyCoord`s.  

The YAML deserialization is the slow part, using `PyYAML`.  `PixelMapCollection` will first try to use the C-based YAML loader, which requires libyaml to be available.  If this fails it falls back to the 20x slower all-Python loader.

The type of mapping that can be expressed is very flexible, and `PixelMaps` can be compounded into chains of tranformations.  See the `gbdes` source files for more documentation on the types of transformations available and their YAML encodings.  One type of map is a `TemplateMap`, which uses a lookup table stored in some other YAML file.  Various template files required for DECam solutions are part of this repository.  The setup area will be searched automatically for requested template files.  Additional paths to search can be specified by the `CAL_PATH` environment variable, using the usual colon-separated list format.

## DELVE Astrometric solutions
As of v1.2.0, there is a `DelveMap` class (derived from `PixelMapCollection` that will return a `WCS` or `PixelMap` instance for any combination of DECam exposure number and CCD identifier. The repository includes the relevant calibration files, which will automatically be found by the code.  _v1.2.0 does not yet include some final details of the DECam solutions nor cover all the exposures - these will be added in later versions._ 

The basic steps to using these solutions are:
* Acquire this repository and run `pip install .` from the top level directory of the repo, with your target Python environment activated.
* In your code, have the following steps:
  ```python
  import pixmappy as pm
  dmc = pm.DelveMaps()
  mywcs = dmc.getDelveWCS(442144, 23)  #Your expnum, ccdnum in the arguments
  ra, dec = mywcs.toSky( x, y, color)  # Transform pixel coord to sky coords

  # Or you can transform to the local gnomonic projection
  mymap = dmc.getDelveMap(442144, 23)
  u,v = mymap(x,y,c)
  ```
* Note that the color argument `c` is assumed to be _g-i_ for the DECam
  data.  A value is required, since the solutions include differential
  chromatic refraction in the atmosphere and lateral
  color in the corrector.  Use a value of `c=0.61` if you don't know
  your true color and want something that is not crazy.
* If you request a solution for an exposure/CCD pair that is not in
  the solution set, a `ValueError` exception will be raised.

### Additional exposure information
The file `delveExposures.hdf5` included the repo's `data/` directory includes additional information about each exposure that could be of use.  This table can be read as an `astropy.table.Table` and includes the following columns:

* __expnum__: DECam exposure number
* __band__: Filter used, one of _griz_.
* __mjdmid__: MJD of the temporal midpoint of the exposure
* __exptime__: exposure time, in seconds
* __t_eff__: "effective exposure time," giving fraction of the actual exposure time that would have yielded the same point-source depth under some nominal, nearly-ideal conditions of seeing, background, and clouds.
* __nite__: Date on the start of the night of the exposure, in format like 20121102.
* __pole__: The ICRS (RA, Dec) (in degrees) of the approximate pointing axis of the exposure.  This is also the pole of the gnomonic projection used for this exposure.
* __hpix__: The healpixel (NSIDE=32, RING ordering) containing the pole of the exposure.
* __airmass, ha__: The airmass and hour angle (in degrees) of the exposure.
* __parallactic__: Parallactic angle, in _radians_ from North through East.
* __obsicrs__: The Cartesian coordinates of the observatory, in barycentric ICRS coordinates, at the midpoint of the exposure (in AU).
* __nCCD30__: Number of CCDs that contain >=30 matches to Gaia stars.  Numbers below 60 indicate problems (e.g. bright stars) with some CCDs.
* __5p2uv__: This is a 2x5 matrix that maps a Gaia 5-parameter solution (u0,v0,pm_ra,pm_dec,parallax) into the two positions (u,v) that the star should have during this exposure.  The u's and v's are the RA and DEC projected around the __pole__ of the exposure, in degrees.
  


## DES Astrometric Solutions

The current versions retain the `DESMaps` class (derived from `PixelMapCollection`) that instantiates the astrometric solutions derived for the DES Y6A1 internal data release.  It works the same way as `DelveMaps` but its solutions should be superceded by `DelveMaps`.

The two capabilities listed below are implemented for `DESMaps` but not yet implemented or tested for `DelveMaps`:

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
    
