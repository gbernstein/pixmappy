#!/usr/bin/env python
# New pixmaps defined using full DES+DELVE

from astropy.table import Table,vstack
import numpy as np
from astropy.time import Time
from . import PixelMap,PixelMapCollection
from .decaminfo import detpos2ccdnum,ccdnum2detpos
from .DESMaps import findOnPath, arg2detpos
from . import files

from scipy.interpolate import interp1d

class ColorConverter:
    def __init__(self):
        # Build the LUTs

        self.lut = {0: lambda c:c}  # Identify for c=g-i
        # c=1 for r-z
        cgi = np.array([ [+0.213, +0.62974],[+0.281, +0.76145],[+0.350, +0.89118],[+0.419, +1.03354],[+0.488, +1.16993],
                         [+0.556, +1.31010],[+0.625, +1.43978],[+0.694, +1.55955],[+0.762, +1.69363],[+0.831, +1.78694],
                         [+0.900, +1.87093],[+0.968, +1.95833],[+1.037, +2.04732],[+1.106, +2.12380],[+1.175, +2.19398],
                         [+1.243, +2.26417],[+1.312, +2.33342],[+1.381, +2.39121],[+1.449, +2.44899],[+1.518, +2.50711],
                         [+1.587, +2.56154],[+1.655, +2.61501],[+1.724, +2.66849],[+1.793, +2.72167],[+1.862, +2.77487],
                         [+1.930, +2.82807],[+1.999, +2.88029],[+2.068, +2.92684],[+2.136, +2.97458],[+2.205, +3.02232],
                         [+2.274, +3.03318],[+2.342, +3.03152],[+2.411, +3.02985],[+2.480, +2.97022],[+2.549, +2.81817],
                         [+2.617, +2.67112],[+2.686, +2.52818]])
        self.lut[1] = interp1d(*cgi.T, kind='linear', bounds_error=False, fill_value='extrapolate')

        # c=2 for r-i
        cgi = np.array([ [+0.168, +0.69850],[+0.215, +0.84963],[+0.263, +1.00121],[+0.311, +1.14671],[+0.358, +1.29264],
                         [+0.406, +1.42588],[+0.453, +1.55002],[+0.501, +1.69622],[+0.548, +1.77039],[+0.596, +1.84978],
                         [+0.643, +1.93751],[+0.691, +2.02799],[+0.738, +2.11881],[+0.786, +2.18230],[+0.833, +2.24579],
                         [+0.881, +2.31035],[+0.929, +2.37114],[+0.976, +2.42488],[+1.024, +2.47862],[+1.071, +2.53351],
                         [+1.119, +2.58704],[+1.166, +2.64057],[+1.214, +2.69410],[+1.261, +2.74842],[+1.309, +2.80261],
                         [+1.356, +2.85680],[+1.404, +2.90949],[+1.451, +2.95961],[+1.499, +3.00974],[+1.547, +3.05986],
                         [+1.594, +3.06845],[+1.642, +3.07045],[+1.689, +3.07246],[+1.737, +3.01998],[+1.784, +2.89600],
                         [+1.832, +2.77881],[+1.879, +2.66295]])
        self.lut[2] = interp1d(*cgi.T, kind='linear', bounds_error=False, fill_value='extrapolate')

        # c=3 for g-r
        cgi = np.array([ [+0.473, +0.61569],[+0.514, +0.67445],[+0.554, +0.73321],[+0.594, +0.79132],[+0.635, +0.84942],
                         [+0.675, +0.90939],[+0.716, +0.96742],[+0.756, +1.02515],[+0.797, +1.08288],[+0.837, +1.14404],
                         [+0.878, +1.20176],[+0.918, +1.25948],[+0.959, +1.31670],[+0.999, +1.38729],[+1.039, +1.45931],
                         [+1.080, +1.52914],[+1.120, +1.59506],[+1.161, +1.65923],[+1.201, +1.50902],[+1.242, +1.70546],
                         [+1.282, +1.84609],[+1.323, +1.95275],[+1.363, +2.08349],[+1.403, +2.23949],[+1.444, +2.38948],
                         [+1.484, +2.52209],[+1.525, +2.62990],[+1.565, +2.71299],[+1.606, +2.77229],[+1.646, +2.79889],
                         [+1.687, +2.79145]])
        self.lut[3] = interp1d(*cgi.T, kind='linear', bounds_error=False, fill_value='extrapolate')

        #c=4 for Gaia bp-rp
        cgi = np.array([ [+0.829, +0.63470],[+0.904, +0.71517],[+0.979, +0.80412],[+1.053, +0.89107],[+1.128, +0.99223],
                         [+1.203, +1.08915],[+1.278, +1.18611],[+1.353, +1.28663],[+1.428, +1.38774],[+1.503, +1.48995],
                         [+1.578, +1.59112],[+1.652, +1.69052],[+1.727, +1.78321],[+1.802, +1.87394],[+1.877, +1.95449],
                         [+1.952, +2.03678],[+2.027, +2.11540],[+2.102, +2.18226],[+2.176, +2.24625],[+2.251, +2.30908],
                         [+2.326, +2.37816],[+2.401, +2.44290],[+2.476, +2.50697],[+2.551, +2.57278],[+2.626, +2.61841],
                         [+2.701, +2.66405],[+2.775, +2.70088],[+2.850, +2.71635],[+2.925, +2.73183],[+3.000, +2.73308],
                         [+3.075, +2.72145],[+3.150, +2.70983],[+3.225, +2.71390],[+3.299, +2.72286],[+3.374, +2.73182]])
        self.lut[4] = interp1d(*cgi.T, kind='linear', bounds_error=False, fill_value='extrapolate')

    def __call__(self, vals, c):
        '''Convert colors other than g-i to the g-i that gives equivalent
        mean DCR in g band.  If input color is >10, output will be 99.,
        as a no-data flag value.  

        First argument is an array of color values.
        The last argument is the input color system:
        0:  g-i (identity transformation in this case)
        1:  r-z
        2:  r-i
        3:  g-r (avoid this)
        4:  Gaia bp-rp
        '''
        return np.where(vals < 10, self.lut[c](vals), 99.)

class ShiftFinder:
    def __init__(self, path='epochShifts.fits'):
        '''Class that will generate Linear PixelMap appropriate to
        a given DECam CCD on a given night of observing. Input
        path to the shifts file created from DES+DELVE data on
        construction.  Then call this object with the (MJD or NITE)
        of the observation and the CCDNUM.

        The returned Linear instance will have a name like
        `20121101/S14` with the indicated nite being
        the night of observations.
        
        Also returns the epoch value and the polynomial set to use.'''
        
        self.tab = Table.read(path)
        self.DECAM_MJD0 = self.tab.meta['MJD0']
        self.starts = self.tab['startDay']
    def __call__(self, mjd, ccdnum):
        if mjd > 20000000:
            # This is a NITE, not an MJD
            nite = int(mjd)
            # Convert NITE notation to the DECam day count
            day =  np.floor(Time('{:04d}-{:02d}-{:02d}'.format(nite//10000, (nite//100)%100, nite%100)).mjd - self.DECAM_MJD0)
        elif mjd<10000:
            # This is a day count
            day = mjd
            ymd = Time(day+self.DECAM_MJD0, format='mjd').ymdhms
            nite =  int(ymd[0]*10000 + ymd[1]*100 + ymd[2])            
        else:
            day = np.floor(mjd - self.DECAM_MJD0 - 0.7)  # Obs of a given NITE are ~0.9-1.5 days past 00:00 UT of the NITE
            ymd = Time(day+self.DECAM_MJD0, format='mjd').ymdhms
            nite =  int(ymd[0]*10000 + ymd[1]*100 + ymd[2])
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

        '''self.includeOld = True   # Combine new tweaks with old polynomial
        rmid = np.array([0.14265489, 0.24708551, 0.31898602, 0.37742935, 0.42796466,
           0.47313273, 0.51434951, 0.5525    , 0.58818117, 0.62181823,
           0.65372682, 0.6841488 , 0.71327443, 0.74125653, 0.76822007,
           0.79426879, 0.81948993, 0.84395769, 0.8677358 , 0.88515029,
           0.89657206, 0.90785014, 0.91898982, 0.92999608, 0.9408736 ,
           0.95162679, 0.96225982, 0.97277663, 0.98318096, 0.99347632,
           1.00366609, 1.01375344, 1.0237414 , 1.03363285, 1.04343054,
           1.05313708, 1.06275497, 1.07228659, 1.08173424, 1.09110008,
           1.1003862 ])
        # Values of new radial shifts (converted to degrees)
        dr_new = np.array([0.00029848, 0.00045575, 0.00099359, 0.00064308, 0.0008012 ,
           0.00122591, 0.00104458, 0.00095653, 0.00086855, 0.00112881,
           0.00113355, 0.00101216, 0.0009162 , 0.0010048 , 0.00117232,
           0.00154478, 0.00178469, 0.00159779, 0.00146685, 0.00177162,
           0.00207413, 0.00256249, 0.00298204, 0.00349814, 0.00386177,
           0.00444939, 0.00509988, 0.00570929, 0.00680259, 0.00760426,
           0.00807254, 0.00848433, 0.00938779, 0.01006415, 0.01018218,
           0.01041629, 0.01055002, 0.00978171, 0.00851019, 0.00698859,
           0.00297378]) / 3600.
        # Coefficients of old correction polynomial, (r, r^3, r^5)
        # (with r and resid both in degrees)
        # Build LUT that uses r^2 as x and dr/r as y
        x = np.concatenate([[0,], rmid*rmid])
        y = np.concatenate([[0,], dr_new / rmid])
        self.rfunc = interp1d(x,y, kind='linear', bounds_error=False,
                            fill_value=(y[0], y[-1]))

        # Now color factors
        cmid = np.arange(0.1,4,0.2)
        y = np.array([-1.78058871e-01, -7.88957863e-02, -2.67864600e-02, -6.70515556e-03,
       -8.94105499e-04,  1.22985517e-04, -7.61278786e-04,  1.03024975e-03,
        3.97318282e-03,  1.11527501e-02,  2.63817935e-02,  4.91845151e-02,
        1.06382047e-01,  2.08982573e-01,  3.65628275e-01,  5.62582647e-01,
        7.90394250e-01,  1.01409756e+00,  1.20539754e+00,  1.33343773e+00])

        self.cfunc = interp1d(cmid, y, kind='linear', bounds_error=False,
                            fill_value=(y[0], y[-1]))
        '''
        self.old_poly = np.array([-9.16033066e-06,  1.92890317e-05, -2.19422989e-05])
        self.oldRef = 0.61  # Reference color

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
        cc = c - self.oldRef
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
            ### dr = self.cfunc(c) * self.rfunc(rsq)
            if False: ##self.includeOld:
                dr += (c-self.oldRef)*(self.old_poly[0] + self.old_poly[1]*rsq + self.old_poly[2]*rsq*rsq)
            else:
                dr = self._dr(rsq,c)
            uv = uv* (1 + dr)[:,np.newaxis]
            return uv[:,0], uv[:,1]
        else:
            # work properly for single scalar input
            rsq = u*u+v*v
            ### dr = self.cfunc(c) * self.rfunc(rsq)
            if False: ##self.includeOld:
                dr += (c-self.oldRef)*(self.old_poly[0] + self.old_poly[1]*rsq + self.old_poly[2]*rsq*rsq)
            else:
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

        oldRef = 0.61
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
            dcr += (cmid-oldRef) * oldFactor

        # Shift cfunc to be zero at g-i=0.61 (oldRef)
        tmp = interp1d(cmid,dcr,kind='linear')
        self.cfunc = interp1d(cmid, dcr - tmp(0.61),
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
 
class DelvePoly(PixelMap):
    @staticmethod
    def type():
        return 'DelvePoly'

    # Static tables for the class that hold the calibration information
    y1Table = None
    y6Table = None

    # And the chip-specific jump locations.  These are in units of 16 pix first
    xJump = np.zeros((63,2),dtype=float)
    # CCDs with a midline jump:
    for detpos in ('S30','S27','S20','S24','S15','S19','S5','N3','N14','N16','N17'):
        xJump[detpos2ccdnum[detpos],1] = 64
    i=detpos2ccdnum['S25']
    xJump[i] = np.array([64,76])
    i=detpos2ccdnum['S29']
    xJump[i] = np.array([64,90])
    i=detpos2ccdnum['S6']
    xJump[i] = np.array([64,98])
    i=detpos2ccdnum['N1']
    xJump[i] = np.array([28,64])
    i=detpos2ccdnum['N6']
    xJump[i] = np.array([25,64])
    i=detpos2ccdnum['N7']
    xJump[i] = np.array([12,64])
    i=detpos2ccdnum['N22']
    xJump[i] = np.array([49,64])  # Just the worse of 2 adjoining regions.


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
            
        self.dtpath = findOnPath(self.dtpath)

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
        self.jumpBegin = DelvePoly.xJump[ccdnum,0]
        self.jumpEnd = DelvePoly.xJump[ccdnum,1]

        # z band - rescale the map
        zFactor = 0.8
        if b=='z':
            self.u2d = self.u2d * zFactor
            self.v2d = self.v2d * zFactor

    def _poly(self,x,y,cc):
        # Hard-coded quartic polynomial and "jump"
        xx = x/1024.-1  # Scale between -1 and 1 for 0 -2048 values
        yy = y/2048.-1
        out = cc[0] + xx*(cc[1] + xx*(cc[3] + xx*(cc[6] + xx*cc[10])))
        out += yy*( cc[2] + xx*(cc[4] + xx*(cc[7] + xx*cc[11])))
        y2 = yy*yy  #y^2 now
        out += y2*(cc[5] + xx*(cc[8] + xx*cc[12]))
        y2 *= yy  #y^3 now
        out += y2*(cc[9] + xx*cc[13])
        out += cc[14]*y2*yy  # y^4

        # Now the "jump", if any
        if cc[-1]!=0.:
            # The jump limits are in binned 16x16 1-indexed pixels
            xx = x/16
            out += cc[-1]*np.logical_and(xx>=self.jumpBegin, xx<self.jumpEnd)
        return out

    def __call__(self, x, y, c=None):
        '''Map (x,y) pixel coords (1-indexed) into (u,v)'''

        # Apply polynomial and any v jump
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

    
    def __init__(self,
                 **kwargs):
        '''Create PixelMapCollection that can create new entries for specified DES
        exposure number / CCD combinations using stored astrometric solutions.  These
        will be sought in local files.  Defaults for these file names are those from
        Y6A1 astrometry release.  An argument of `None` indicates that the file is absent.


        Other kwargs are passed to PixelMapCollection
        '''

        # Make these file names fixed...
        guts_file=files.default_delve_guts
        exposures_file=files.default_delve_exposures
        affine_file=files.default_delve_affine

        # Add the tweaker to PixelMapCollection atoms
        PixelMapCollection.addAtom(DelvePoly)
        PixelMapCollection.addAtom(DCR)
        PixelMapCollection.addAtom(LC)
        
        # Find the guts_file and initialize with it
        path = findOnPath(guts_file)
        super(DelveMaps, self).__init__(filename=path, **kwargs)

        # Read in the tabular exposure information
        path = findOnPath(exposures_file)
        self.exptab = Table.read(path)

        # Get affine shifts
        path = findOnPath(affine_file)
        self.sf = ShiftFinder(path=path)

        return

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
        # including new shift and DCR.
        shift, epoch, poly = self.sf(row['mjdmid'], ccdnum)
        polyname = '{:s}/{:s}/{:s}/dpoly'.format(band,poly,detpos),
        lcname = '{:s}/color2'.format(band),
        expo_map = 'D{:08d}/expo'.format(expnum)

        elements = ['{:s}/{:s}/rings'.format(band,detpos),
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
