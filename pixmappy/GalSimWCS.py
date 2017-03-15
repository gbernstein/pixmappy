
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
        :param exp:         The exposure number of the desired WCS. [default: None]
        :param ccdnum:      The CCD number of the desired WCS. [default: None]
        :param origin:      Optional origin position for the image coordinate system.
                            If provided, it should be a PositionD or PositionI.
                            [default: None]
        :param cache:       Cache this file's PixelMapCollection in the GalSimWCS.cache dict?
                            [default: True]
        """
        _req_params = { "file_name" : str }
        _opt_params = { "origin" : galsim.PositionD, "ccdnum": int }
        _single_params = [ { "wcs_name" : str, "exp" : int } ]
        _takes_rng = False
        
        info = DECamInfo()
        cache = dict()

        def __init__(self, file_name=None, dir=None, pmc=None, wcs_name=None,
                     exp=None, ccdnum=None, origin=None, cache=True):
            self._color = None
            if file_name is not None:
                if dir is not None:
                    file_name = os.path.join(dir,file_name)
                if pmc is not None:
                    raise TypeError("Cannot provide both file_name and pmc")
                if file_name in self.cache:
                    pmc = self.cache[file_name]
                else:
                    pmc = PixelMapCollection(file_name)
                    if cache:
                        self.cache[file_name] = pmc
                self._tag = 'file_name = ' + file_name
            else:
                self._tag = 'pmc = '+repr(pmc)
            
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
            coord = self._wcs.toSky(x, y, c=c )
            ra = coord.icrs.ra.rad
            dec = coord.icrs.dec.rad
            return ra, dec

        def _xy(self, ra, dec, c=None):
            sky_coord = astropy.coordinates.SkyCoord(ra,dec,unit='rad')
            x, y = self._wcs.toPix(sky_coord, c=c)
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
            return "galsim.GalSimWCS(%s, wcs_name=%s, origin=%r)"%(
                    self._tag, self._wcs_name, self.origin)

        def __hash__(self): return hash(repr(self))


