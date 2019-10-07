
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
    import os
    import astropy.coordinates
    import numpy as np

    class GalSimWCS(galsim.wcs.CelestialWCS):
        """A wrapper of the PixelMapCollection class that can be used as a galsim WCS type.

        This can be constructed using either a `file_name` or an existing PixelMapCollection
        object (as `pmc`).

        :param file_name:   The yaml file with the PixelMapCollection description.
        :param dir:         Optional directory to prepend to `file_name`. [default: None]
        :param pmc:         An existing pixmappy.PixelMapCollection instance [default: None]
        :param wcs_name:    The name of the WCS within the PixelMapCollection to use.
                            [default: None; either wcs_name or (exp and ccdnum) is required.]
        :param exposure_file:  FITS file holding binary table of DES per-exposure info.
                            [default: None; if present, this triggers using DESMaps rather than
                            a regular PixelMapCollection.]
        :param resids_file: FITS file holding 2d residual adjustment maps for DECam devices
                            [default: None]
        :param affine_file: FITS file holding time-dependent DECam CCD affine tweaks.
                            [default: None
        :param exp:         The exposure number of the desired WCS. [default: None]
        :param ccdnum:      The CCD number of the desired WCS. [default: None]
        :param origin:      Optional origin position for the image coordinate system.
                            If provided, it should be a PositionD or PositionI.
                            [default: None]
        :param cache:       Cache this file's PixelMapCollection in the GalSimWCS.cache dict?
                            [default: True]
        :param default_color:   The default color to use if this WCS involves color terms and
                            `wcs.toWorld` or similar methods to not pass in a color term.
                            [default: None, which means an exception will be raised if no color
                            term is provided]
        """
        _req_params = { "file_name" : str }
        _opt_params = { "origin" : galsim.PositionD, "ccdnum": int,
                        "exposure_file" : str, "resids_file" : str, "affine_file" : str }
        _single_params = [ { "wcs_name" : str, "exp" : int } ]
        _takes_rng = False
        
        info = DECamInfo()
        cache = dict()

        def __init__(self, file_name=None, dir=None, pmc=None, wcs_name=None,
                     exposure_file=None, resids_file=None, affine_file=None,
                     exp=None, ccdnum=None, origin=None, cache=True, default_color=None):
            self._color = default_color
            if file_name is not None:
                if dir is not None:
                    file_name = os.path.join(dir,file_name)
                    exposure_file = os.path.join(dir,exposure_file) if exposure_file else None
                    resids_file = os.path.join(dir,resids_file) if resids_file else None
                    affine_file = os.path.join(dir,affine_file) if affine_file else None
                if pmc is not None:
                    raise TypeError("Cannot provide both file_name and pmc")
                if file_name in self.cache:
                    pmc = self.cache[file_name]
                else:
                    if exposure_file is None:
                        pmc = PixelMapCollection(file_name)
                    else:
                        pmc = DESMaps(guts_file=file_name,
                                      exposure_file=exposure_file,
                                      resids_file=resids_file,
                                      affine_file=affine_file)
                    if cache:
                        self.cache[file_name] = pmc
                self._tag = 'file_name=%r'%file_name
            else:
                self._tag = 'pmc=%r'%pmc
            
            if pmc is None:
                raise TypeError("Must provide either file_name or pmc")
            self._pmc = pmc
            if wcs_name is not None:
                if exp is not None or ccdnum is not None:
                    raise TypeError("Cannot provide both wcs_name and (exp,ccdnum)")
                self._wcs_name = wcs_name
            else:
                if exp is None or ccdnum is None:
                    raise TypeError("Must provide either wcs_name or (exp,ccdnum)")
                self.exp = exp
                self.ccdnum = ccdnum
                self.ccdname = self.info.ccddict[ccdnum]
                self._wcs_name = 'D%s/%s'%(self.exp, self.ccdname)
            self._wcs = pmc.getWCS(self._wcs_name)

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
        def wcs_name(self): return self._wcs_name
                
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
            rather than all objects).  It is a dict indexed by the the file_name.
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
                    self._wcs_name == other._wcs_name and
                    self.origin == other.origin )

        def __repr__(self):
            s = "pixmappy.GalSimWCS(%s, wcs_name=%r, origin=%r"%(
                    self._tag, self._wcs_name, self.origin)
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

