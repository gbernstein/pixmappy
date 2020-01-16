import os
root_dir = os.path.dirname(__file__)
data_dir = os.path.join(root_dir, 'data')

# Default names for the DES pixel maps solutions
default_guts_file='y6a1.guts.astro'
default_exposure_file='y6a1.exposureinfo.fits'
default_resids_file='y6a1.astroresids.fits'
default_affine_file='y6a1.affine.fits'

# Default path will be this data directory
default_cal_path = data_dir
