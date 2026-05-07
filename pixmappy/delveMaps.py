#!/usr/bin/env python
# New pixmaps defined using full DES+DELVE

from astropy.table import Table,vstack
import numpy as np
from astropy.time import Time
from . import PixelMap,PixelMapCollection
from .decaminfo import detpos2ccdnum,ccdnum2detpos, arg2detpos, nite2day, day2nite, DECAM_MJD0,REF_COLOR
from . import files
from scipy.interpolate import interp1d


class ShiftFinder:
    def __init__(self, affine_file=files.default_delve_affine):
        '''Class that will generate Linear PixelMap appropriate to
        a given DECam CCD on a given night of observing. Input
        path to the shifts file created from DES+DELVE data on
        construction.  Then call this object with the (MJD or NITE)
        of the observation and the CCDNUM.

        The returned Linear instance will have a name like
        `20121101/S14` with the indicated nite being
        the night of observations.
        
        Also returns the epoch value and the polynomial set to use.'''
        
        self.tab = Table.read(files.findOnPath(affine_file))
        self.starts = self.tab['startDay']
    def __call__(self, mjd, ccdnum):
        if mjd > 20000000:
            # This is a NITE, not an MJD
            nite = int(mjd)
            # Convert NITE notation to the DECam day count
            day =  nite2day(nite)
        elif mjd<10000:
            # This is a day count
            day = mjd
            nite =  day2nte(day)
        else:
            day = np.floor(mjd - DECAM_MJD0 - 0.7)  # Obs of a given NITE are ~0.9-1.5 days past 00:00 UT of the NITE
            nite =  day2nite(day)
        index = np.searchsorted(self.starts, day, side='right') - 1
        if index<0:
            raise(ValueError('Requested MJD '+str(mjd) + ' is before epochs begin'))
        startDay = self.starts[index]
        epoch = self.tab['epoch'][index]
        detpos = ccdnum2detpos[ccdnum]
        name = '{:08d}/{:s}'.format(nite,detpos)

        dt = day - self.tab['d0'][index]
        affine = self.tab['affine'][index,ccdnum,:,:] @ np.array([1,dt])
            
        # PixMappy class wants flattened 2x3
        mm = {'Type':'Linear', 'name':name, 'Coefficients':affine.T.flatten()}
        poly = 'Y1'
        if nite >=  20180619:  # start of Y6
            poly = 'Y6'
        
        return mm, epoch, poly
    

# New lateral color map
class LC(PixelMap):

    @staticmethod
    def type():
        return 'LC'

    def __init__(self, name, **kwargs):

        '''PixelMap that makes the recomputed lateral color correction for
        g or r band (25 Mar 2026)'''
        # These are tabulations of the corrections.
        # Locations of reference points for new shifts
        super(LC,self).__init__(name)

        self.old_poly = np.array([-9.16033066e-06,  1.92890317e-05, -2.19422989e-05])

        MAS = 1 / 3600. / 1000.   # in degrees
        band = kwargs['band']
        if band=='g':
            self.coeffs = np.array( [[-11.2*MAS, +29.5*MAS + self.old_poly[0]],
                                     [-22.4*MAS, -65.7*MAS + self.old_poly[1]],
                                     [ +6.4*MAS, +60.8*MAS + self.old_poly[2]]])
        elif band=='r':
            self.coeffs = np.array( [[ -1.6, +1.3],
                                     [ 10.6, -0.8],
                                     [-10.4, -5.5]]) * MAS
        elif band=='i':
            self.coeffs = np.array( [[ -6.0, +1.7],
                                     [ 15.5, -2.0],
                                     [ -6.8, -2.0]]) * MAS
        elif band=='z':
            self.coeffs = np.array( [[ -4.6, +0.7],
                                     [ 11.8, +0.7],
                                     [ -6.8, -1.7]]) * MAS
        return
    
    def _dr(self, rsq, c):
        cc = c - REF_COLOR
        return (self.coeffs[0,0] + cc*self.coeffs[0,1]) \
             + (self.coeffs[1,0] + cc*self.coeffs[1,1]) * rsq \
             + (self.coeffs[2,0] + cc*self.coeffs[2,1]) * rsq *rsq
             
    def __call__(self, u, v, c):
        '''Apply tweaks to DECam pixel positions
        :param x,y: uv plane coordinate arrays about optic axis, degrees
        :param c:   g-i color of source(s).
        :returns: u, v tweaked positions.
        '''
        if np.array(c).ndim>0:
            uv = np.stack([u,v],axis=-1)
            rsq = np.sum(uv*uv, axis=-1)
            dr = self._dr(rsq,c)
            uv = uv* (1 + dr)[:,np.newaxis]
            return uv[:,0], uv[:,1]
        else:
            # work properly for single scalar input
            rsq = u*u+v*v
            dr = self._dr(rsq,c)

            return u*(1 + dr),v*(1+dr)

# New DCR map
class DCR(PixelMap):
    @staticmethod
    def type():
        return 'DCR'

    def __init__(self, name, **kwargs):
        '''Tranformation of uv coordinates with new DCRs.
        Arguments:
        `band`: 'g', 'r', 'i', or 'z'
        `airmass`:   sec(z)
        `parallactic`:  parallactic angle, radians'''

        super(DCR,self).__init__(name)

        if 'band' not in kwargs or 'airmass' not in kwargs or 'parallactic' not in kwargs:
            raise ValueError('Missing arguments for DCR PixelMap')
        
        airmass = kwargs['airmass']
        parallactic = kwargs['parallactic']
        band = kwargs['band']
        
        # These are tabulations of the corrections.
        # Locations of reference points for new shifts
        self.uv0 = np.sqrt(airmass*airmass-1) * np.array([np.sin(parallactic), np.cos(parallactic)])
        
        self.includeOld = True  # False would omit the older corrections
        
        # DCR corrections derived in notebook
        # Now color factors
        cmid = np.arange(0.1,4,0.2)
        dcr_c = {'g': np.array([-0.00165, -0.00669, -0.00142, +0.00818, +0.01672, 
                                +0.02357, +0.02939, +0.03532, +0.04135, +0.04725, 
                                +0.05266, +0.05775, +0.06256, +0.06744, +0.07255, 
                                +0.07770, +0.08277, +0.08762, +0.09054, +0.09125]),
             'r': np.array([-0.00395, -0.00230, -0.00079, +0.00065, +0.00226, 
                            +0.00387, +0.00529, +0.00677, +0.00843, +0.01007, 
                            +0.01202, +0.01454, +0.01730, +0.02028, +0.02312, 
                            +0.02561, +0.02805, +0.03056, +0.03292, +0.03521]),
             'i': np.array([-0.00118, -0.00073, -0.00022, +0.00018, +0.00067, 
                            +0.00124, +0.00175, +0.00240, +0.00313, +0.00398, 
                            +0.00503, +0.00638, +0.00768, +0.00891, +0.00994, 
                            +0.01125, +0.01267, +0.01392, +0.01510, +0.01619]),
             'z': np.array([-0.00056, -0.00031, -0.00014, +0.00011, +0.00034, 
                            +0.00053, +0.00072, +0.00089, +0.00102, +0.00118, 
                            +0.00139, +0.00164, +0.00190, +0.00219, +0.00247, 
                            +0.00287, +0.00337, +0.00385, +0.00433, +0.00481]) }

        if band in 'griz':
            dcr = dcr_c[band] / 3600.  # Convert from arcsec to degrees
        else:
            # Y band has no corrections
            dcr = np.zeros_like(dcr_c['z'])
        if band in 'Y':
            # Add back in previous slopes (which were in mas/mag):
            dcrConstant = {'g':45.0, 'r':8.4, 'i':3.2, 'z':1.4, 'Y':1.1}  
            # Convert mas/mag to degree/mag
            oldFactor = dcrConstant[band] / 3600. / 1000.
            dcr += (cmid-REF_COLOR) * oldFactor

        # Shift cfunc to be zero at g-i=REF_COLOR
        tmp = interp1d(cmid,dcr,kind='linear')
        self.cfunc = interp1d(cmid, dcr - tmp(REF_COLOR),
                              kind='linear', bounds_error=False, fill_value='extrapolate')

    def __call__(self, u, v, c):
        '''Apply DCR correction to (u,v) coordinates (in degrees)'''
        
        if np.array(c).ndim>0:
            duv = self.cfunc(c)[:,np.newaxis] * self.uv0
            return u+duv[:,0], v+duv[:,1]
        else:
            # Work properly for scalar inputs
            duv = self.cfunc(c) * self.uv0
            return u+duv[0], v+duv[1]
 
class Trap(PixelMap):
    @staticmethod
    def type():
        return 'Trap'

    # Static table with trap information
    tab = None
    trapfile = files.default_delve_traps
    
    def __init__(self, name, **kwargs):
        '''This PixelMap applies nominal serial trap correction
        to the x pixel coordinate based on lookup tables.  It does
        *not* apply a correction to the y coordinate of objects
        behind the parallel traps on CCD #57.  

        The `mas` method will return an evaluation of the nominal
        size of the x trap, in milliarcsec.  A value of -100 will
        be returned for objects subject to the CCD57 parallel trap,
        indicating that the y (=u) coordinate will be unreliable.

        Arguments:
        `ccdnum`: 
        `band`: Must be in `griz`
        `nite`: Evening of observation, e.g. 20180422
        `trapfile`: Any desired alternative trap table to load.
           (Only works for first instance created, since table is static)
        '''

        super(Trap,self).__init__(name)
        
        if Trap.tab is None:
            # Need to load the table
            if 'trapfile' in kwargs:
                # Use given file instead of default
                Trap.trapfile = kwargs['trapfile']
            Trap.tab = Table.read(files.findOnPath(Trap.trapfile))
        elif 'trapfile' in kwargs:
            if kwargs['trapfile']!=Trap.trapfile:
                raise ValueError('Attempt to change trap table after opening')
        
        self.ccdnum = kwargs['ccdnum']
        use = Trap.tab['ccdnum']==self.ccdnum
        if not np.any(use):
            # No traps for this CCD
            self.dx = None
            return
        self.xStart = Trap.tab['xStart'][use]
        self.xStop = Trap.tab['xStop'][use]
        dv = []
        # Convert NITE notation to the DECam day count from DECAM_MJD0
        nite = kwargs['nite']
        day =  np.floor(Time('{:04d}-{:02d}-{:02d}'.format(nite//10000, (nite//100)%100, nite%100)).mjd - 56232)
        for kv in Trap.tab['knots_'+kwargs['band']][use]:
            # Interpolate amplitude of shift to given MJD
            nk = np.count_nonzero(kv[0])  # How many knots
            if day<=kv[0,0]:
                # Const before first knot
                dv.append(kv[1,0])
            elif day>=kv[0,nk-1]:
                #...or after last knot
                dv.append(kv[1,nk-1])
            else:
                for j in range(0,nk-1):
                    f = (day-kv[0,j])/(kv[0,j+1]-kv[0,j])
                    if f>0:
                        dv.append((1-f)*kv[1,j] + f*kv[1,j+1])
                        break
        self.dx = -np.array(dv) / 264.   # Change from mas of v to pixels of x
        return
    def __get_dx(self,x):
        '''Evaluate the shift for an array of x values'''
        dx = np.zeros_like(x,dtype=float)
        if self.dx is None:
            # No traps
            pass
        else:
            for i in range(len(self.dx)):
                dx += self.dx[i] * np.logical_and(x>self.xStart[i], x<self.xStop[i])
        return dx
    
    def __call__(self, x, y, c=None):
        '''Apply shifts to x, if any'''
        xx = np.array(x)
        if self.dx is not None:
            xx += self.__get_dx(xx)
        return xx, np.array(y)

    def mas(self,x,y):
        '''Return the nominal trap shift in mas for sources at (x,y).
        Postiive values are size of the serial (x, v direction) trap.
        Negative value indicates inside the parallel trap zone on CCD57,
        size not estimated.'''
        out = np.abs(self.__get_dx(x)) * 264  # Convert to mas
        if self.ccdnum==57:
            # Mark the zone of CCD 57 having bad parallel trap
            bad = np.all(np.stack([x>=193,
                                   x<641,
                                   y>=3073],axis=0),axis=0)
            out[bad] = -100  # Mark as -100 mas trap.
        return out
       
class DelvePoly(PixelMap):
    @staticmethod
    def type():
        return 'DelvePoly'

    # Static tables for the class that hold the calibration information
    y1Table = None
    y6Table = None

    def __init__(self, name, **kwargs):
        '''Tranformation of xy to uv coordinates, including delve 2d tweaks.
        Arguments:
        `band`: 'g', 'r', 'i', or 'z'
        `ccdnum`:   
        `poly`: 'Y1' or 'Y6'
        '''

        super(DelvePoly,self).__init__(name)

        if 'path' in kwargs:
            self.dtpath = kwargs['path']
        else:
            self.dtpath = files.default_delve_tweaks
            
        self.dtpath = files.findOnPath(self.dtpath)

        if 'band' not in kwargs or 'ccdnum' not in kwargs or 'poly' not in kwargs:
            raise ValueError('Missing arguments for DelvePoly PixelMap')


        if kwargs['poly']=='Y1':
            # Load the tweak table if needed
            if DelvePoly.y1Table is None:
                DelvePoly.y1Table = Table.read(self.dtpath,path='Y1')
            tab = DelvePoly.y1Table
        elif kwargs['poly']=='Y6':
            # Load the tweak table if needed
            if DelvePoly.y6Table is None:
                DelvePoly.y6Table = Table.read(self.dtpath,path='Y6')
            tab = DelvePoly.y6Table
        else:
            raise ValueError('DelvePoly given invalid poly epoch ' + kwargs['poly'])

        b = kwargs['band']
        ccdnum = kwargs['ccdnum']
        j = np.where(tab['ccdnum']==ccdnum)[0][0]
        self.ucoeff = tab['coeff_u'+b][j]
        self.vcoeff = tab['coeff_v'+b][j]
        self.u2d = tab['map_u'][j]
        self.v2d = tab['map_v'][j]

        # z band - rescale the map
        zFactor = 0.8
        if b=='z':
            self.u2d = self.u2d * zFactor
            self.v2d = self.v2d * zFactor

    def _poly(self,x,y,cc):
        # Hard-coded quartic polynomial
        xx = x/1024.-1  # Scale between -1 and 1 for 0 -2048 values
        yy = y/2048.-1
        out = cc[0] + xx*(cc[1] + xx*(cc[3] + xx*(cc[6] + xx*cc[10])))
        out += yy*( cc[2] + xx*(cc[4] + xx*(cc[7] + xx*cc[11])))
        y2 = yy*yy  #y^2 now
        out += y2*(cc[5] + xx*(cc[8] + xx*cc[12]))
        y2 *= yy  #y^3 now
        out += y2*(cc[9] + xx*cc[13])
        out += cc[14]*y2*yy  # y^4

        return out

    def __call__(self, x, y, c=None):
        '''Map (x,y) pixel coords (1-indexed) into (u,v)'''

        # Apply polynomial
        u = self._poly(x,y,self.ucoeff)
        v = self._poly(x,y,self.vcoeff)

        ### apply 2d map adjustment, no interpolation
        xbin = np.floor(x/16).astype(np.int16)
        xbin = np.clip(xbin, 0, 127)
        ybin = np.floor(y/16).astype(np.int16)
        ybin = np.clip(ybin, 0, 255)

        u += self.u2d[ybin,xbin]
        v += self.v2d[ybin,xbin]

        return u, v
        
class DelveMaps(PixelMapCollection):
    '''DelveMaps is an extension of PixelMapCollection that allows the
    user to build WCS/PixelMaps for all DECam exposures in DELVE DR3 by 
    extracting exposure-specific information from custom tables.
    The user must also have a local copy of the YAML file
    specifying the PixelMaps for the "guts" of astrometric solution -
    the time-independent specifications of certain camera distortions.
    Environment variable CAL_PATH gives the path to search for these
    files; this module's data directory will be searched if no CAL_PATH
    is in environment.  Current directory is searched last.
    '''
    wcsName = 'D{:08d}/{:s}'         # String to format to get WCS name for expo/detpos pair
    basemapName = 'D{:08d}/{:s}/base' # String to format for PixelMap name

    
    def __init__(self, exposures_file=files.default_delve_exposures,
                 **kwargs):
        '''Create PixelMapCollection that can create new entries for specified DES
        exposure number / CCD combinations using stored astrometric solutions.  These
        will be sought in local files.  
        Arguments:
        `exposure_file`: File holding table of WCS information for all exposures. Use
        default if absent.

        Other kwargs are passed to PixelMapCollection
        '''

        # Make these file names fixed...
        guts_file=files.default_delve_guts
        affine_file=files.default_delve_affine

        # Add the tweaker to PixelMapCollection atoms
        PixelMapCollection.addAtom(DelvePoly)
        PixelMapCollection.addAtom(DCR)
        PixelMapCollection.addAtom(LC)
        PixelMapCollection.addAtom(Trap)
        
        # Find the guts_file and initialize with it
        path = files.findOnPath(guts_file)
        super(DelveMaps, self).__init__(filename=path, **kwargs)

        # Read in the tabular exposure information
        path = files.findOnPath(exposures_file)
        self.exptab = Table.read(path)
        self.exptab.sort('expnum')   ### Could be pre-sorted ????

        # Get affine shifts
        self.sf = ShiftFinder(affine_file=affine_file)

        return


    def trapMapFor(self,expnum, detpos):
        '''Returns the name of the `Trap` map for this exposure/detector pair.
        Gives the PixelMapCollection its specifications if it does not exist yet.
        Arguments:
        `expnum`:
        `detpos`: detpos or ccdnum.
        '''
        detpos = arg2detpos(detpos)
        exp_row = np.searchsorted(self.exptab['expnum'],expnum)
        if exp_row > len(self.exptab) or self.exptab['expnum'][exp_row]!=expnum:
            raise ValueError('No solution found for expnum {:08d}'.format(expnum))
        nite = self.exptab['nite'][exp_row]
        band = self.exptab['band'][exp_row]
        name = '{:s}/{:s}/{:08d}/trap'.format(band,detpos,nite)

        # Register the map if it doesn't exist
        if not self.hasMap(name):
            self.update( {name:{'Type':'Trap',
                                'band':band,
                                'ccdnum':detpos2ccdnum[detpos],
                                'nite':nite}})
        return name         

    def getDelveMap(self, expnum, detpos):
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

    def getDelveWCS(self, expnum, detpos):
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

    ### Not implemented yet def getCovariance(self, expnum, defaultError=10.):

    def _acquireWCS(self, expnum, detpos):
        '''Acquire info on exposure/detpos combo from files and 
        add it to the PixelMapCollection.
        '''

        # Find the row of exposure table corresponding to this expnum 
        exp_row = np.searchsorted(self.exptab['expnum'],expnum)
        if exp_row > len(self.exptab) or self.exptab['expnum'][exp_row]!=expnum:
            raise ValueError('No solution found for expnum {:08d}'.format(expnum))
        row = self.exptab[exp_row]
        # Make a dictionary that we'll add to the PixelMapCollection
        pixmaps = {}
        # Make WCS dictionary entry
        ccdnum = detpos2ccdnum[detpos]
        basemap = self.basemapName.format(expnum,detpos)
        wcs = {'Type':'WCS',
               'MapName':basemap,
               'Projection':{'Type':'Gnomonic',
                             'Xi':0.,
                             'Eta':0.,
                             'Orientation':{'RA':row['pole'][0],
                                            'Dec':row['pole'][1],
                                            'PA':0.}},
               'Scale':0.0174532925199433}
        # Add this WCS spec to the dictionary
        pixmaps['WCS'] = {self.wcsName.format(expnum,detpos):wcs}

        # Build the PixelMap elements of this map:
        elements = []
        
        if self.hasMap(basemap):
            # Return existing map
            return self.getMap(basemap)            

        # Define the dcr PixelMap for this exposure
        dcrname = 'D{:07d}/dcr'.format(expnum)
        band = row['band']
        if not self.hasMap(dcrname):
            dcrmap = {'Type':'DCR','band':band,
                                   'airmass':row['airmass'],
                                   'parallactic':row['parallactic']}
            pixmaps[dcrname] = dcrmap

        # Aquire ccdshift and build exposure solution,
        # including new shift, DCR, and traps.
        shift, epoch, poly = self.sf(row['mjdmid'], ccdnum)
        polyname = '{:s}/{:s}/{:s}/dpoly'.format(band,poly,detpos),
        lcname = '{:s}/color2'.format(band),
        expo_map = 'D{:08d}/expo'.format(expnum)
        trapname = self.trapMapFor(expnum,detpos)

        elements = [trapname,
                    '{:s}/{:s}/rings'.format(band,detpos),
                    '{:s}/{:s}/lowedge'.format(band,detpos),
                    '{:s}/{:s}/highedge'.format(band,detpos),
                    polyname,
                    shift['name'],
                    lcname,
                    dcrname,
                    expo_map]
            
        if not self.hasMap(lcname):
            pixmaps[lcname] = {'Type':'LC','band':band}
        if not self.hasMap(polyname):
            pixmaps[polyname] = {'Type':'DelvePoly',
                                    'band':band,'poly':poly,
                                    'ccdnum':ccdnum}
        if not self.hasMap(shift['name']):
            # Add the shift map to the collection if it's new
            nn = shift.pop('name')
            pixmaps[nn] = shift

        # Now the polynomial exposure solution if we don't have it already
        if not self.hasMap(expo_map):
            coeffs = row['coeffs'].copy()
            # Coefficients stored are for the delta; add in the original value too.
            coeffs[1,0] += 1.
            coeffs[2,1] += 1.
            poly = {'Type':'Poly',
                'XMin': -1,
                'XMax': 1,
                'YMin': -1,
                'YMax': 1,
                'Tolerance': 2.778e-07,
                'XPoly':{'SumOrder':True,
                        'OrderX': 3,
                        'Coefficients': coeffs[:,0].tolist()},
                'YPoly':{'SumOrder':True,
                        'OrderX': 3,
                        'Coefficients': coeffs[:,1].tolist()}}
            pixmaps[expo_map] = poly

            
        # Add the composite to the new pixmaps
        pixmaps[basemap] = {'Type':'Composite',
                            'Elements':elements}

        # Add new pixmaps to the PixelMapCollection
        self.update(pixmaps)
