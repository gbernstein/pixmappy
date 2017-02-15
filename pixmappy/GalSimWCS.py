
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
        """
        _req_params = { "file_name" : str }
        _opt_params = { "origin" : galsim.PositionD, "ccdnum": int }
        _single_params = [ { "wcs_name" : str, "exp" : int } ]
        _takes_rng = False
        
        info = DECamInfo()
        cache = dict()

        def __init__(self, file_name=None, dir=None, pmc=None, wcs_name=None,
                     exp=None, ccdnum=None, origin=None):
            if file_name is not None:
                if dir is not None:
                    file_name = os.path.join(dir,file_name)
                if pmc is not None:
                    raise TypeError("Cannot provide both file_name and pmc")
                if file_name in self.cache:
                    pmc = self.cache[file_name]
                else:
                    pmc = PixelMapCollection(file_name)
                    self.cache[file_name] = pmc
                self._tag = 'file_name = ' + file_name
            else:
                self._tag = 'pmc = '+repr(pmc)
            
            if pmc is None:
                raise TypeError("Must provide either file_name or pmc")
            self._pmc = pmc
            self._origin = origin
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
            print('radec for x,y = ',x,',',y)
            coord = self._wcs.toSky( (x,y), c=c )
            print('coord = ',coord)
            return coord.ra, coord.dec

        def _xy(self, ra, dec, c=None):
            print('xy for ra,dec = ',ra,',',dec)
            xy = self._wcs.toPix( astropy.coordinates.SkyCoord(ra,dec), c=c )
            print('xy = ',xy)
            return xy[0], xy[1]

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


