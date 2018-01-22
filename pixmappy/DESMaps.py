import os
import astropy.io.fits as pf
import numpy as np

from . import PixelMapCollection, Identity, Constant, ColorTerm, Polynomial, Composite, WCS

def findOnPath(filename, envPathName='CAL_PATH'):
    '''Look for existing file with the name <filename> using the paths
    in the colon-separated list stored in environment variable
    with name <envPathName>.  Searches current directory first.
    If filename is an absolute path, just tries that.
    '''
    if os.path.isabs(filename):
        if os.path.isfile(filename):
            return filename
        else:
            raise IOError('Absolute path <' + filename + ' is non-existent file')
    else:
        paths = []
        if envPathName in os.environ:
            paths += os.environ[envPathName].split(':')
        # And last resort will be current directory
            paths.append('')
        for p in paths:
            path = os.path.join(p,filename)
            if os.path.isfile(path):
                return path
        raise IOError('Can not find file ' + filename + ' in path ' + envPathName)

class DESMaps(PixelMapCollection):
    '''DESMaps is an extension of PixelMapCollection that allows the
    user to build WCS/PixelMaps for DES survey exposures by extracting
    exposure-specific information from custom FITS tables or the DESDM
    database.  The user must have a local copy of the YAML file
    specifying the PixelMaps for the "guts" of astrometric solution -
    the time-independent specifications of camera distortions, and the
    small tweaks for different observing epochs, as well as local
    copies of the templates for the tree-ring and edge distortions.
    The CAL_PATH will be used to search for these, as for
    PixelMapCollection.  

    The user must also have access either to the
    DESDM database (by providing an easyaccess connection) or to a
    FITS file containing the tabulated astrotrmetric info for each
    exposure/ccd.
    '''
    exposureName = 'D{:06d}Z{:03d}'  # String to format to get exposure name
    wcsName = 'D{:06d}/{:s}'         # String to format to get WCS name for expo/detpos pair
    basemapName = 'D{:06d}/{:s}/base' # String to format for PixelMap name
    
    def __init__(self, conn=None,
                     guts_file='y4a1.guts.astro',
                     exposure_file='y4a1.expastro.fits', ccd_file='y4a1.ccdastro.fits',
                     **kwargs):
        '''Create PixelMapCollection that can create new entries for specified DES
        exposure number / CCD combinations using stored astrometric solutions.  These
        will be created using information in DESDM database if a connection is given.
        Otherwise these will be sought in two FITS binary tables. 

        guts_file: locally available YAML file with time-invariant portions of solution.
        conn:  easyaccess connection to dessci database.  
        exposure_file:  FITS file holding binary table of DES per-exposure info
        ccd_file: FITS file holding binary table of DES per-ccd info
        Other kwargs passed to PixelMapCollection
        '''

        # Find the guts_file and initialize with it
        path = findOnPath(guts_file)
        super(DESMaps, self).__init__(filename=path, **kwargs)

        self.conn = conn
        if self.conn is None:
            # Read in the tabular information from FITS files
            path = findOnPath(ccd_file)
            self.ccdtab = pf.getdata(path,1)

            path = findOnPath(exposure_file)
            self.exptab = pf.getdata(path,1)

        return

    def getDESMap(self, expnum, detpos):
        '''Acquire PixelMap for specified exposure number / CCD combination
        '''
        name = self.basemapName.format(expnum,detpos)
        if not self.hasMap(name):
            self._acquireWCS(expnum,detpos)
        return self.getMap(name)

    def getDESWCS(self, expnum, detpos):
        '''Acquire WCS for specified exposure number / CCD combination
        '''
        name = self.wcsName.format(expnum,detpos)
        if not self.hasWCS(name):
            self._acquireWCS(expnum,detpos)
        return self.getWCS(name)

    def _acquireWCS(self, expnum, detpos):
        '''Acquire info on exposure/detpos combo from database/files and 
        add it to the PixelMapCollection.
        '''

        if self.conn is not None:
            raise NotImplementedError('Database access to astrometry solutions not ready yet')

        # Find row of the ccdtable containing this expnum/detpos
        i0 = np.searchsorted(self.ccdtab['expnum'],expnum)
        # Desired detpos might be anywhere within next 62 rows of first expnum
        tmp = np.where( np.logical_and( self.ccdtab['expnum'][i0:i0+62]==expnum,
                                        self.ccdtab['detpos'][i0:i0+62]==detpos))[0]
        if len(tmp)<1:
            raise ValueError('No CCD solution found for ' + self.wcsName.format(expnum,detpos))
        elif len(tmp)>1:
            raise ValueError('Error, multiple CCD solutions for ' + self.wcsName.format(expnum,detpos))

        ccd_row = i0 + tmp[0]
        zone = self.ccdtab['zone'][ccd_row]

        # Find the row of exposure table corresponding to this expnum / zone
        exp_row = np.searchsorted(self.exptab['expnum'],expnum)
        # Search for match
        while True:
            if exp_row > len(self.exptab) or self.exptab['expnum'][exp_row]!=expnum:
                raise ValueError('No  solution found for expnum/zone {:06d}/{:03d}'.format(expnum,zone))
            if self.exptab['zone'][exp_row]==zone:
                break
            exp_row = exp_row + 1
            
        # Make a dictionary that we'll add to the PixelMapCollection
        pixmaps = {}
        # Make WCS dictionary entry
        basemap = self.basemapName.format(expnum,detpos)
        wcs = {'Type':'WCS',
               'MapName':basemap,
               'Projection':{'Type':'Gnomonic',
                             'Xi':0.,
                             'Eta':0.,
                             'Orientation':{'RA':self.exptab['ra'][exp_row],
                                            'Dec':self.exptab['dec'][exp_row],
                                            'PA':0.}},
               'Scale':0.0174532925199433}
        # Add this WCS spec to the dictionary
        pixmaps['WCS'] = {self.wcsName.format(expnum,detpos):wcs}

        # Build the PixelMap elements of this map:
        # Start with the instrumental solution, already in PixelMapCollection:
        elements = ['{:s}{:s}/{:s}'.format(self.exptab['filter'][exp_row],
                                           self.exptab['epoch'][exp_row],
                                           detpos)]

        # Then DCR map, if this exposure needs one
        dcr_map = 'D{:06d}/dcr'.format(expnum)
        if np.any(self.exptab['dcr'][exp_row]):
           elements.append(dcr_map)

           # Make the DCR map if we don't have it
           if not self.hasMap(dcr_map):
               d = self.exptab['dcr'][exp_row]
               dcr = {'Type':'Color',
                      'Reference': d[2],
                      'Function':{'Type':'Constant',
                                  'Parameters':d[:2].tolist()}}
               pixmaps[dcr_map] = dcr

        # Then the exposure solution
        expo_map = self.exposureName.format(expnum,zone)
        elements.append(expo_map)

        # Add the composite to the new pixmaps
        pixmaps[basemap] = {'Type':'Composite',
                            'Elements':elements}
        
        # Now create the polynomial exposure solution if we don't have it already
        if not self.hasMap(expo_map):
            poly = {'Type':'Poly',
                'XMin': -1,
                'XMax': 1,
                'YMin': -1,
                'YMax': 1,
                'Tolerance': 2.778e-07,
                'XPoly':{'SumOrder':True,
                        'OrderX': 3,
                        'Coefficients': self.exptab['xpoly'][exp_row].tolist()},
                'YPoly':{'SumOrder':True,
                        'OrderX': 3,
                        'Coefficients': self.exptab['ypoly'][exp_row].tolist()}}
            pixmaps[expo_map] = poly


        # Add new pixmaps to the PixelMapCollection
        self.update(pixmaps)


        

