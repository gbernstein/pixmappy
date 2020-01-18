
try:
    import galsim
except ImportError:
    class GalSimWCS(object):
        def __init__(self, *args, **kwargs):
            raise NotImplementedError(
                'Unable to import galsim. The GalSimWCS interface is not available.')
else:

    from .decaminfo import DECamInfo
    from .PixelMapCollection import PixelMapCollection
    from .DESMaps import DESMaps
    import os
    import astropy.coordinates
    import numpy as np
    from . import files

    class GalSimWCS(galsim.wcs.CelestialWCS):
        """A wrapper of the `PixelMapCollection` class that can be used as a galsim `WCS` type.

        This can be constructed using either:
        * an existing `PixelMapCollection` object (supplied as `pmc`), or
        * a generic serialized PixelMapCollection in a YAML file, 
          (supplied by `yaml_file`), or
        * The `DESMaps` class, which accesses files in the format of
          the DES Y6A1_ASTROMETRY release, which one selects via
          `useDESMaps=True`.  In this case one can supply names for
          the `exposure_file,guts_file,resids_file,affine_file` that
          the `DESMaps` needs, or any of these will use Y6A1 default
          names if they are supplied as `None` (This is default
          behavior).

        Exactly one of these must be true: `pmc is not None`;
        `yaml_file is not None`; `useDESMaps`.

        The WCS within the collection is selected by name.  One can
        either specify the name directly with `wcs_name`; or a name
        can be constructed using DES conventions from an exposure
        number `exp` and `ccdnum`.

        :param pmc:         An existing pixmappy.PixelMapCollection instance [default: None]
        :param yaml_file: The yaml file with the PixelMapCollection
                            description. [default: None]
        :param useDESMaps:  If `True`, use `DESMaps` to construct WCS. [default: False]

        :param dir:         Optional directory to prepend to all filename 
                            arguments. [default: None]

        :param wcs_name: The name of the WCS within the
                            PixelMapCollection to use.  [default:
                            None; either `wcs_name` or (`exp` and
                            `ccdnum`) is required.  DESMaps require
                            the exp+ccdnum.]
        :param exp:         The exposure number of the desired WCS. [default: None]
        :param ccdnum:      The CCD number of the desired WCS. [default: None]

        :param exposure_file: FITS file holding binary table of DES
                            per-exposure info for DESMaps.  [default:
                            None; if `useDESMaps` then the file in the
                            Y6A1_ASTRONOMY release will be used in
                            this default case.]  
        :param guts_file:   YAML file holding static DECam distortions
                            for DESMaps.  [default: None; (same
                            behavior as above).  
        :param resids_file: FITS file holding 2d residual adjustment
                            maps for DESMaps [default: None; (same
                            behavior as above).  
        :param affine_file: FITS file holding time-dependent DECam CCD
                            affine tweaks for DESMaps [default: None;
                            (same behavior as above).

        :param origin:      Optional origin position for the image coordinate system.
                            If provided, it should be a PositionD or
                            PositionI.  [default: None]
        :param cache:       Cache this file's PixelMapCollection in the GalSimWCS.cache dict?
                            [default: True] 
        :param default_color: The default color to use if this WCS involves color terms and
                            `wcs.toWorld` or similar methods do not pass in a color value.
                            [default: None, which means an exception will be raised if no
                            color term is in the map and no color value is provided]
        """
        _opt_params = { "origin" : galsim.PositionD, "ccdnum": int,
                        "dir" : str, "guts_file" : str,
                        "exposure_file" : str, "resids_file" : str, "affine_file" : str,
                        "default_color" : float}
        _single_params = [ { "wcs_name" : str, "exp" : int },
                               {"yaml_file" : str, "use_DESMaps":bool}]
        _takes_rng = False
        
        info = DECamInfo()
        cache = dict()

        def __init__(self, pmc=None, yaml_file=None, use_DESMaps=False, dir=None,
                     wcs_name=None, exp=None, ccdnum=None,
                     exposure_file=None, guts_file=None, resids_file=None, affine_file=None,
                     origin=None, cache=True, default_color=None):
            self._color = default_color

            # Make sure only one method is in use:
            count = int(pmc is not None) + int(yaml_file is not None) + int(use_DESMaps)
            if count!=1:
                raise TypeError("Must provide exactly one of yaml_file, pmc, or use_DESMaps")
            
            if pmc is not None:
                self._pmc = pmc
                self._tag = 'pmc=%r'%pmc  # Used in __repr__

            if yaml_file is not None:
                if dir is not None:
                    yaml_file = os.path.join(dir,yaml_file)
                if yaml_file in self.cache:
                    pmc = self.cache[yaml_file]
                else:
                    pmc = PixelMapCollection(yaml_file)
                    if cache:
                        self.cache[yaml_name] = pmc
                self._tag = 'yaml_file=%r'%yaml_file
                self._pmc = pmc
                
            if use_DESMaps:
                if exp is None or ccdnum is None:
                    raise TypeError("exp and ccdnum must be provided when using DESMaps")

                self._tag = 'use_DESMaps=True'
                
                if exposure_file is None:
                    exposure_file = files.default_exposure_file
                else:
                    self._tag = self._tag + ', exposure_file=%s'%exposure_file

                if guts_file is None:
                    guts_file = files.default_guts_file
                else:
                    self._tag = self._tag + ', guts_file=%s'%guts_file

                if resids_file is None:
                    resids_file = files.default_resids_file
                else:
                    self._tag = self._tag + ', resids_file=%s'%resids_file

                if affine_file is None:
                    affine_file = files.default_affine_file
                else:
                    self._tag = self._tag + ', affine_file=%s'%affine_file
                    
                if dir is not None:
                    exposure_file = os.path.join(dir,exposure_file)
                    guts_file = os.path.join(dir,guts_file)
                    resids_file = os.path.join(dir,resids_file)
                    affine_file = os.path.join(dir,affine_file)
                    self._tag = self._tag + ', dir=%s'%dir

                # We'll cache the DESMaps object by the exposure_file name
                if exposure_file in self.cache:
                    pmc = self.cache[exposure_file]
                else:
                    pmc = DESMaps(guts_file=guts_file,
                                  exposure_file=exposure_file,
                                  resids_file=resids_file,
                                  affine_file=affine_file)
                    if cache:
                        self.cache[exposure_file] = pmc
                self._pmc = pmc

            # Now extract the desired WCS from our PixelMapCollection or DESMaps
            if use_DESMaps:
                if exp is None or ccdnum is None:
                    raise TypeError("DESMaps require exp,ccdnum")
                ccdname = self.info.ccddict[ccdnum]
                self._wcs_name = 'D%s/%s'%(exp, ccdname)  #Used by __eq__
                self._wcs = pmc.getDESWCS(exp, ccdname)
                self._tag = self._tag + ', exp=%r, ccdnum=%r'%(exp,ccdnum)
                
            else:
                if wcs_name is None:
                    if exp is None or ccdnum is None:
                        raise TypeError("Must provide either wcs_name or (exp,ccdnum)")
                    ccdname = self.info.ccddict[ccdnum]
                    self._wcs_name = 'D%s/%s'%(exp, ccdname)
                elif exp is not None or ccdnum is not None:
                    raise TypeError("Cannot provide both wcs_name and (exp,ccdnum)")
                else:
                    self._wcs_name = wcs_name

                self._wcs = pmc.getWCS(self._wcs_name)
                self._tag = self._tag + ', wcs_name=%r'%self._wcs_name

            # Set origin, if any
            if origin is None:
                self._origin = galsim.PositionD(0,0)
            else:
                if isinstance(origin, galsim.PositionI):
                    origin = galsim.PositionD(origin.x, origin.y)
                elif not isinstance(origin, galsim.PositionD):
                    raise TypeError("origin must be a PositionD or PositionI argument")
                self._origin = origin


        @property
        def pmc(self): return self._pmc

        @property
        def wcs_name(self):
            # Note: older versions of pixmappy didn't set _wcs_name if it was None.
            # This is fixed above, but to accommodate reading in older serialized pixmappy
            # objects, we check for the attribute existing here.
            return self._wcs_name if hasattr(self,'_wcs_name') else None
                
        @property
        def origin(self): return self._origin

        @classmethod
        def clear_cache(cls):
            """Clear the cache of PixelMapCollections.
            
            The PixelMapCollection objects that are read in from a file are often rather large,
            and a typical use case involves getting many wcs objects from the same file.
            So the GalSimWCS class caches them avoid needing to read in the file repeatedly for
            each WCS you want to extract from it.

            However, if you are done with a particular input file, you might not want to keep
            it around anymore.  So ``pixmappy.GalSimWCS.clear_cache()`` will release the memory
            currently being used by the cache.

            You can also modify the cache yourself if you want (say to remove a particular element
            rather than all objects).  It is a dict indexed by the the yaml filename or the exposure
            filename for DESMaps.
            """
            cls.cache.clear()

        def _radec(self, x, y, c=None):
            ra, dec = self._wcs.toSky(x, y, c=c )
            ra *= galsim.degrees / galsim.radians
            dec *= galsim.degrees / galsim.radians
            return ra, dec

        def _xy(self, ra, dec, c=None):
            ra *= galsim.radians / galsim.degrees
            dec *= galsim.radians / galsim.degrees
            x, y = self._wcs.toPix(ra, dec, c=c)
            return x, y

        def _newOrigin(self, origin):
            ret = self.copy()
            ret._origin = origin
            return ret

        def _writeHeader(self, header, bounds):
            raise NotImplementedError("Cannot write PixelMap to a fits header")

        @staticmethod
        def _readHeader(header):
            raise NotImplementedError("Cannot read PixelMap from a fits header")

        def copy(self):
            # The copy module version of copying the dict works fine here.
            import copy
            return copy.copy(self)

        def __eq__(self, other):
            return (isinstance(other, GalSimWCS) and
                    self._tag == other._tag and
                    self.wcs_name == other.wcs_name and
                    self.origin == other.origin ) 

        def __repr__(self):
            # Should eval back into itself
            s = "pixmappy.GalSimWCS(%s, origin=%r"%(self._tag,  self.origin)
            if self._color is not None:
                s += ', default_color=%r'%self._color
            s += ')'
            return s

        def __hash__(self): return hash(repr(self))

        def __getstate__(self):
            # The naive pickling works, but it includes _pmc, which is huge, and not actually
            # necessary for the functioning of the object.  (It's just for information purposes
            # really.)  So remove it from the dict to be pickled.
            d = self.__dict__.copy()
            d['_pmc'] = None
            return d
