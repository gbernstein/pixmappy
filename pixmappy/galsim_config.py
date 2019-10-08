
# This class makes the pixmappy.GalSimWCS usable as a GalSim config wcs type.
# To use it, just add pixmappy to the modules section of a config file.
# Then in the config file, the wcs type is called 'Pixmappy'.

try:
    import galsim
    from .GalSimWCS import GalSimWCS

    class PixmappyBuilder(galsim.config.WCSBuilder):

        def buildWCS(self, config, base, logger):

            req = { "file_name" : str,
                    "exp" : str,
                    "ccdnum" : int
                  }
            opt = { "dir" : str,
                    "exposure_file" : str, "resids_file" : str, "affine_file" : str
                  }
            kwargs, safe = galsim.config.GetAllParams(config, base, req=req, opt=opt)

            logger.info('Loading WCS for %s ccd %s',kwargs['exp'],kwargs['ccdnum'])
            wcs = GalSimWCS(**kwargs)
            logger.info('Done loading pixmappy WCS')

            wcs._color = 0  # For now.  Maybe make this settable somehow.
            return wcs

    galsim.config.RegisterWCSType('Pixmappy', PixmappyBuilder())

except:
    pass  # Don't fail if galsim isn't available.  Only fail it they try to use
          # the PixmappyBuilder and galsim isn't available.
