# Paths to standard files for DECam maps

from importlib.resources import files as irfiles
data_dir = irfiles('pixmappy').joinpath('data')

# Default names for the DES pixel maps solutions
default_guts_file='y6a1.guts.astro'
default_exposure_file='y6a1.exposureinfo.fits'
default_resids_file='y6a1.astroresids.fits'
default_affine_file='y6a1.affine.fits'

default_delve_affine = 'epochShifts.fits'
default_delve_exposures = 'delveExposures.hdf5'
default_delve_tweaks = 'delveTweaks_f4.hdf5'
default_delve_guts = 'delve.guts.astro'

# Default path will be data directory given in pyproject
default_cal_path = str(data_dir)
