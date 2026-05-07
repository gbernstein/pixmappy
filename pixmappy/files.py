# Paths to standard files for DECam maps

import os
from importlib.resources import files as irfiles

data_dir = irfiles('pixmappy').joinpath('data')

# Default names for the DES pixel maps solutions
default_guts_file='y6a1.guts.astro'
default_exposure_file='y6a1.exposureinfo.fits'
default_resids_file='y6a1.astroresids.fits'
default_affine_file='y6a1.affine.fits'

default_delve_affine = 'epochAffine5D.hdf5'
default_delve_exposures = 'delveExposures.hdf5'
default_delve_tweaks = 'delveTweaks5D.hdf5'
default_delve_guts = 'delve.guts.astro'
default_delve_traps = 'trapTable.hdf5'

# Default path will be data directory given in pyproject
default_cal_path = str(data_dir)

def findOnPath(filename, envPathName='CAL_PATH'):
    '''Look for existing file with the name <filename> using the paths
    in the colon-separated list stored in environment variable
    with name <envPathName>.  Searches current directory first.
    If filename is an absolute path, just tries that.

    :param filename: full absolute path to file, or a relative path
    :param envPathName: environment variable which optionally stores
        a list of paths to search, in order of decreasing priority, 
        for a relative pathname. [default=`CAL_PATH`].  Current directory
        is always appended to the path.
    :returns: full path to file
    :raises: IOError if file is not found on any path.
    '''
    if os.path.isabs(filename):
        if os.path.isfile(filename):
            return filename
        else:
            raise IOError('Absolute path <' + filename + '> is non-existent file')
    else:
        paths = []
        if envPathName in os.environ:
            paths += os.environ[envPathName].split(':')
            pathFound = True
        else:
            # Use the default path if none is in the environment
            paths += default_cal_path.split(':')
            pathFound = False
        # And the current directory is always searched last
        paths.append('')
        
        for p in paths:
            path = os.path.join(p,filename)
            if os.path.isfile(path):
                return path

        # If we get here we have failed
        if pathFound:
            raise IOError('Cannot find file ' + filename + ' in path ' + envPathName)
        else:
            raise IOError('Cannot find file ' + filename + ' in default path ' + default_cal_path)
