# pixmappy
Python interface to gbdes pixel map (astrometry) solutions

`PixelMapCollection` class can read the YAML format astrometry solutions produced by Gary's `WcsFit` program (in gbdes repo).  This class can then issue a `PixelMap` instance, which is a map from one 2d coordinate system ("pixel") to another ("world") 2d system.  A `PixelMap` instance can be used as a function mapping one (or many) coordinate pairs.  An `inverse` method does reverse mapping, and the local `jacobian` of the map is available also.

A `WCS` is a `PixelMap` that additionally specifies a (de)projection from the world 2d system onto the celestial sphere.  It uses `astropy.coordinates.SkyCoord` objects to represent celestial positions, and has `toSky` and `toPix` methods to map between the pixel system and `SkyCoord`s.  

The YAML deserialization is the slow part, using `PyYAML`.  `PixelMapCollection` will first try to use the C-based YAML loader, which requires libyaml to be available.  If this fails it falls back to the 20x slower all-Python loader.

The type of mapping that can be expressed is very flexible, and `PixelMaps` can be compounded into chains of tranformations.  See the `gbdes` source files for more documentation on the types of transformations available and their YAML encodings.  One type of map is a `TemplateMap`, which uses a lookup table stored in some other YAML file.  Two such template files used for DECam solutions are part of this repository and installed by setup.py.  The setup area will be searched automatically for requested template files.  Additional paths to search can be specified by the `CAL_PATH` environment variable, using the usual colon-separated list format.
