import os
import astropy.io.fits as pf
import numpy as np
from scipy.interpolate import RectBivariateSpline

from . import PixelMap,PixelMapCollection, Identity, Constant, ColorTerm, Polynomial, Composite, WCS
from . import files

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
            paths += files.default_cal_path.split(':')
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
            raise IOError('Cannot find file ' + filename + ' in default path ' + files.default_cal_path)

def arg2detpos(arg_in):
    '''Return a DECam DETPOS specification given either a DETPOS string
    or CCDNUM integer.
    :raises: `ValueError` if input does not correspond to a DECam CCD'''
    if type(arg_in)==str:
        return arg_in
    elif type(arg_in)==int:
        ccdnum2detpos = {1:'S29',  2:'S30',  3:'S31',  4:'S25',  5:'S26',  6:'S27',
                 7:'S28',  8:'S20',  9:'S21', 10:'S22', 11:'S23', 12:'S24',
                13:'S14', 14:'S15', 15:'S16', 16:'S17', 17:'S18', 18:'S19',
                19:'S8',  20:'S9',  21:'S10', 22:'S11', 23:'S12', 24:'S13',
                25:'S1',  26:'S2',  27:'S3',  28:'S4',  29:'S5',  30:'S6',
                31:'S7',  32:'N1',  33:'N2',  34:'N3',  35:'N4',  36:'N5',
                37:'N6',  38:'N7',  39:'N8',  40:'N9', 41:'N10', 42:'N11',
                43:'N12', 44:'N13', 45:'N14', 46:'N15', 47:'N16', 48:'N17',
                49:'N18', 50:'N19', 51:'N20', 52:'N21', 53:'N22', 54:'N23',
                55:'N24', 56:'N25', 57:'N26', 58:'N27', 59:'N28', 60:'N29',
                61:'N30', 62:'N31'}
        if arg_in not in ccdnum2detpos:
            raise ValueError('Invalid DECam CCD number {:d}'.format(arg_in))
        return ccdnum2detpos[arg_in]
    else:
        raise ValueError('DECam CCD number must be str or int')

class DECamTweak():
    '''DECamTweak applies a 2d lookup table of astrometric corrections to
    measured pixel positions, based on gridded mean astrometric residuals
    for each CCD, and time-specific affine transforms per CCD.'''
    def __init__(self, resids_file=files.default_resids_file,
                     affine_file=files.default_affine_file):
        '''
        DECam tweak specifications will be read from two files, searching
        through the CAL_PATH specified by environment (or in pixmappy/data
        directory if no CAL_PATH is given), followed by current dir.

        :param resids_file: FITS file holding gridded astrometric residuals
          for each CCD.
        :param affine_file: FITS file giving affine shifts per CCD and epoch
        '''

        # Open the file of tweaks and create spline lookup tables for each
        # device
        if resids_file is None:
            self.tweaks = None
        else:
            ff = pf.open(findOnPath(resids_file))
            self.tweaks = {}
            for hdu in ff[1:]:
                detpos = hdu.header['EXTNAME']
                binpix = hdu.header['BINPIX']
                nx = hdu.data.shape[2]
                ny = hdu.data.shape[1]
                # Locate grid points, in 1-indexed pixel system
                xvals = binpix * np.arange(nx) + 0.5*binpix + 1
                yvals = binpix * np.arange(ny) + 0.5*binpix + 1
                bbox = [1, nx*binpix+1, 1, ny*binpix+1]
                # Create linear spline for x and y components
                # Note that data array comes in with (y,x) indexing
                self.tweaks[detpos] = (RectBivariateSpline(xvals, yvals, hdu.data[0].transpose(),
                                                           bbox=bbox, kx=1, ky=1),
                                   RectBivariateSpline(xvals, yvals, hdu.data[1].transpose(),
                                                           bbox=bbox, kx=1, ky=1))
            ff.close()

        if affine_file is None:
            self.affine = None
        else:
            bigtab = pf.getdata(findOnPath(affine_file),1)
            # Split the table up into little tables for each detpos
            dps = np.unique(bigtab['detpos'])
            self.affine = {}
            for dp in dps:
                self.affine[dp] = bigtab[bigtab['detpos']==dp]
        return

    def getDataFor(self, detpos, mjd):
        """Get just the pieces of the big tables that are relevant for a single ccd and time
        :returns: spline interpolator, affine table for this detector at this time.
        """
        dp = arg2detpos(detpos)
        if self.tweaks is not None:
            if dp not in self.tweaks:
                raise IndexError('No 2d tweaks available for detpos',dp)
            spline = self.tweaks[dp]
        else:
            spline = None

        if self.affine is not None:
            if dp not in self.affine:
                raise IndexError('No affine tweaks available for detpos',dp)
            # Find the row for this MJD
            iRow = np.searchsorted(self.affine[dp]['mjd'], mjd, side='right')-1
            rr = self.affine[dp][iRow]
        else:
            rr = None

        return (spline, rr)

    @staticmethod
    def tweakFromData(data, xpos, ypos):
        """Tweak using the data for a particular detpos and mjd.
        :param data: the spline,tweak tuple returned by `getDataFor()`
        :param xpos,ypos: input pixel position arrays 
        :returns: x,y output (tweaked) pixel positions.
        """
        spline, rr = data

        if spline is not None:
            xpos, ypos = xpos-spline[0](xpos,ypos,grid=False), \
                         ypos-spline[1](xpos,ypos,grid=False)

        if rr is not None:
            xx = np.copy(xpos) # temporary copy
            xpos = xpos - rr['x0'] + (rr['mag']+rr['e1'])*(xpos-1024.5) \
                          + (rr['e2']+rr['rot'])*(ypos-2048.5)
            ypos = ypos - rr['y0'] + (rr['e2']-rr['rot'])*(xx-1024.5) \
                          + (rr['mag']-rr['e1'])*(ypos-2048.5)

        return xpos,ypos

    def tweak(self, detpos, mjd, xpos, ypos):
        '''Apply tweak to data
        :param detpos: which CCD the positions are from
        :param mjd: MJD of exposure 
        :param xpos,ypos: input pixel positions
        :returns: x,y output (tweaked) pixel positions
        '''
        return self.tweakFromData(self.getDataFor(detpos, mjd), xpos, ypos)


    def tweakTable(self, tab, detpos, mjd, xkey='xpix', ykey='ypix'):
        ''' Tweak the contents of the two columns of the table 
        giving pixel positions of objects.
        :param tab: the table holding pixel coordinates
        :param detpos: CCD of the data
        :param mjd: MJD of the exposure
        :param xkey, ykey: column names for the pixel positions
             [default:'xpix','ypix']
        :returns: nothing.  Pixel positions are tweaked in-place in
             the table.
        '''
        xx, yy = self.tweak(detpos, mjd, tab[xkey],tab[ykey])
        tab[xkey] = xx
        tab[ykey] = yy
        return

class Tweak(PixelMap):
    '''PixelMap that implements the small time- and detpos-dependent
    adjustments to pixel positions described by the DECamTweak class above.
    '''
    @staticmethod
    def type():
        return 'Tweak'

    # A class variable contains the DECamTweak instance and the files
    # from which it came.  Mixing tweaks in the same run will be an error.

    tweaker = None
    residsFile = None
    affineFile = None
    
    def __init__(self, name, **kwargs):
        '''PixelMap that implements DECam astrometric tweaks and affine shifts.
        :param name: name to be given to this `PixelMap`
        :param Detpos: CCD that positions are on
        :param MJD: time of exposure that positions are from
        :param ResidsFile: name of FITS file holding DECam tweaks, if any
        :param AffineFile: name of FITS file holding DECam CCD shifts, if any
        '''
        super(Tweak,self).__init__(name)
        
        if self.tweaker is None:
            # Need to read in a tweak file.
            # Note: Setting class variables requires using class name, not self.
            if 'ResidsFile' in kwargs:
                Tweak.residsFile = kwargs['ResidsFile']
            if 'AffineFile' in kwargs:
                Tweak.affineFile = kwargs['AffineFile']
            Tweak.tweaker = DECamTweak(resids_file = self.residsFile, affine_file = self.affineFile)
        else:
            # Check that any requested files agree with the one we have loaded already
            if ('residsFile' in kwargs and residsFile != kwargs['residsFile']) or \
               ('affineFile' in kwargs and affineFile != kwargs['affineFile']):
                raise ValueError('Tweak maps use inconsistent lookup table files')
               
        # Now save away the mjd and detpos for this instance.
        if 'Detpos' not in kwargs or 'MJD' not in kwargs:
            raise ValueError('Missing Detpos or MJD in Tweak PixelMap')
        self.dp = kwargs['Detpos']
        self.mjd = kwargs['MJD']
        self.tweak_data = self.tweaker.getDataFor(self.dp, self.mjd)

    def __call__(self, x, y, c=None):
        '''Apply tweaks to DECam pixel positions
        :param x,y: pixel coordinate arrays or scalars
        :param c:   color of source(s).  Not used.
        :returns: x, y tweaked pixel positions.
        '''
        return DECamTweak.tweakFromData(self.tweak_data, x, y)


class DESMaps(PixelMapCollection):
    '''DESMaps is an extension of PixelMapCollection that allows the
    user to build WCS/PixelMaps for DES survey exposures by extracting
    exposure-specific information from custom FITS tables.
    The user must also have a local copy of the YAML file
    specifying the PixelMaps for the "guts" of astrometric solution -
    the time-independent specifications of camera distortions, and the
    small tweaks for different observing epochs, as well as local
    copies of the templates for the tree-ring and edge distortions.
    Environment variable CAL_PATH gives the path to search for these
    files; this module's data directory will be searched if no CAL_PATH
    is in environment.  Current directory is searched last.
    '''
    exposureName = 'D{:06d}'  # String to format to get exposure name
    wcsName = 'D{:06d}/{:s}'         # String to format to get WCS name for expo/detpos pair
    basemapName = 'D{:06d}/{:s}/base' # String to format for PixelMap name
    tweakName = 'D{:06d}/{:s}/twk'  # String to create Tweak map name
    
    def __init__(self,
                 guts_file=files.default_guts_file,
                 exposure_file=files.default_exposure_file,
                 resids_file=files.default_resids_file,
                 affine_file=files.default_affine_file,
                 **kwargs):
        '''Create PixelMapCollection that can create new entries for specified DES
        exposure number / CCD combinations using stored astrometric solutions.  These
        will be sought in local files.  Defaults for these file names are those from
        Y6A1 astrometry release.  An argument of `None` indicates that the file is absent.

        :param guts_file:   locally available YAML file with time-invariant portions of solution.
        :param exposure_file: FITS file holding binary table of DES per-exposure info
        :param resids_file: FITS file holding 2d residual adjustment maps for DECam devices (None to skip)
        :param affine_file: FITS file holding time-dependent DECam CCD affine tweaks (None to skip)

        Other kwargs are passed to PixelMapCollection
        '''

        # Add the tweaker to PixelMapCollection atoms
        PixelMapCollection.addAtom(Tweak)
        
        # Find the guts_file and initialize with it
        path = findOnPath(guts_file)
        super(DESMaps, self).__init__(filename=path, **kwargs)

        # Read in the tabular information from FITS files
        path = findOnPath(exposure_file)
        self.exptab = pf.getdata(path,1)
        self.residsFile = resids_file
        self.affineFile = affine_file
        return

    def getDESMap(self, expnum, detpos):
        '''Acquire PixelMap for specified exposure number / CCD combination.

        :param expnum:  exposure number for the desired `PixelMap`
        :param detpos:  CCD number or detpos string for desired `PixelMap`
        :returns: A valid `PixelMap` for this exposure/CCD
        '''
        detpos = arg2detpos(detpos)
        name = self.basemapName.format(expnum,detpos)
        if not self.hasMap(name):
            self._acquireWCS(expnum,detpos)
        return self.getMap(name)

    def getDESWCS(self, expnum, detpos):
        '''Acquire WCS for specified exposure number / CCD combination

        :param expnum:  exposure number for the desired `WCS`
        :param detpos:  CCD number or detpos string for desired `WCS`
        :returns: A valid `WCS` for this exposure/CCD
        '''

        detpos = arg2detpos(detpos)
        name = self.wcsName.format(expnum,detpos)
        if not self.hasWCS(name):
            self._acquireWCS(expnum,detpos)
        return self.getWCS(name)

    def getCovariance(self, expnum, defaultError=10.):
        '''Return the estimated covariance matrix for atmospheric
        astrometric errors in the selected exposure.  Units are
        in mas^2.  [0] axis points east, [1] axis north. Call
        covarianceWarning()  to check for potentially invalid matrix.
        A circular error of defaultError radius is returned if
        there is no valid matrix for this expnum.

        :param expnum:  exposure number for desired atmospheric turbulence
        :param defaultError: turbulence error (in mas) to assign if a value
               is not available in the file. [default: 10 mas]
        :returns: A 2x2 covariance matrix for astrometric turbulence in
               this exposure. (units of mas^2)
        '''
        # Find the row of exposure table corresponding to this expnum 
        #(by default, searchsorted returns matching row if one is equal)
        exp_row = np.searchsorted(self.exptab['expnum'],expnum)
        cov = self.exptab['cov'][exp_row]
        out = np.zeros( (2,2), dtype=float)
        if cov[0]<=0.:
            out[0,0] = defaultError*defaultError
            out[1,1] = out[0,0]
        else:
            out[0,0] = cov[0]
            out[1,1] = cov[1]
            out[0,1] = cov[2]
            out[1,0] = cov[2]
        return out

    def covarianceWarning(self,expnum):
        '''Returns True if the estimated covariance matrix for expnum
        is suspicious because of negative or too-small eigenvalues.

        :param expnum:  exposure number for desired atmospheric turbulence
        '''
        # Find the row of exposure table corresponding to this expnum 
        exp_row = np.searchsorted(self.exptab['expnum'],expnum)
        return self.exptab['cov'][exp_row][0] <= 0

    def _acquireWCS(self, expnum, detpos):
        '''Acquire info on exposure/detpos combo from files and 
        add it to the PixelMapCollection.
        '''

        # Find the row of exposure table corresponding to this expnum 
        exp_row = np.searchsorted(self.exptab['expnum'],expnum)
        if exp_row > len(self.exptab) or self.exptab['expnum'][exp_row]!=expnum:
            raise ValueError('No solution found for expnum {:06d}'.format(expnum))
            
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
        elements = []
        
        # Start with DECam tweaks, if they are in use:
        if self.residsFile is not None or self.affineFile is not None:
            twk = self.tweakName.format(expnum,detpos)
            elements.append(twk)
            if not self.hasMap(twk):
                # Need to create the Tweak map
                pixmaps[twk] = {'Type':'Tweak',
                     'ResidsFile':self.residsFile,
                     'AffineFile':self.affineFile,
                     'Detpos':detpos,
                      'MJD':self.exptab['mjd'][exp_row]}
                
        # Next the instrumental solution, already in PixelMapCollection:
        elements.append('{:s}{:s}/{:s}'.format(self.exptab['band'][exp_row],
                                           self.exptab['epoch'][exp_row],
                                           detpos))

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
        expo_map = self.exposureName.format(expnum)
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
