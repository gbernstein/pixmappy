'''
Implementation of gbtools/astrometry PixelMap classes in Python.
This is meant as a read-only version of the astrometry library.
Capabilities include:
* deserialize the YAML files specifying PixelMaps and WCS's
* Execute forward (pixel to world) transformations with any map/wcs
* Calculate local derivatives of maps
* Execute inverse transformations using generic solver method

All transformations can be done on arrays of coordinates as well as scalars.
That is, the input (x,y) or (ra,dec) arguments may be either float values or
numpy arrays.  In the latter case, the return values will also be numpy arrays.
Sky positions, RA, Dec, are in degrees.

Template PixelMaps access their templates from YAML-formatted files of their own.
The CAL_PATH environment variable can be used to give a list of paths to search
for these template files.
'''

# ??? How does invert tolerance work, which dimension??
# ?? Could speed up interpolation of tables if desired as they're equal-spaced.

import numpy as np
import astropy.coordinates as co
import future
import yaml
from scipy.optimize import root
import os
import sys
try:
    import cPickle as pickle
except ImportError:
    import pickle

from . import files

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from .files import data_dir

class PixelMap(object):
    ''' Base class for transformations from one 2d system ("pixel") to another ("world").
    Each derived class must implement __call__ to execute this transform on array of
    shape (2) or (N,2).  Base class implements some routines such as taking local derivatives
    and solving for map inverse.
    c is always a color (or array of them)
    '''
    def __init__(self, name):
        self.name = name

    def jacobian(self, x, y, c=None, step=1.):
        '''Use finite differences with indicated step size to get Jacobian derivative
        matrix at each point x, y.  If x, y are arrays then the returned array has
        shape (2,2,N) with [i,j,n] giving dx_i/dx_j at nth point.
        '''
        dx = np.array(self(x+step, y, c))
        dx -= self(x-step, y, c)
        dy = np.array(self(x, y+step, c))
        dy -= self(x, y-step, c)
        dx /= 2.*step
        dy /= 2.*step
        out = np.array([dx, dy]).T
        return out

    def inverse(self, xw, yw, xp, yp, c=None, tol=0.0001):
        '''
        Fill the arrays xp, yp with the solutions to PixelMap(xp, yp)=xw, yw.  The input values
        of xp, yp are the initial guess.  Tolerance can be altered from default, which
        is given in the pixel space ???
        '''
        # Need to call the solver row by row, will be slow.
        class resid(object):
            ''' Callable giving deviation of output from target
            '''
            def __init__(self, pmap, targetx, targety, c=None):
                self.target = np.array([targetx, targety])
                self.pmap = pmap
                self.c = c

            def __call__(self, xyp):
                xp, yp = xyp
                return np.array(self.pmap(xp, yp, self.c)) - self.target

        try:
            npts = len(xw)
        except TypeError:
            xyp = np.array([xp, yp])
            xp, yp = root(resid(self.__call__, xw, yw, c), xyp, tol=tol).x
            return xp, yp
        else:
            if len(yp) != npts or len(yp) != npts or len(yw) != npts:
                raise ValueError('Mismatch of xw, yw, xp, yp point counts in PixelMap.inverse')

            for i in range(npts):
                xp[i], yp[i] = self.inverse(xw[i], yw[i], xp[i], yp[i], c, tol)
            return xp, yp
    
'''
Following are various derived classes of PixelMaps.  They each have a class attribute
field "type" which corresponds to the "Type" field in the YAML node.
'''
    
class Identity(PixelMap):
    '''Identity map
    '''
    @staticmethod
    def type():
        return 'Identity'

    def __init__(self, name, **kwargs):
        super(Identity,self).__init__(name)

    def __call__(self, x, y, c=None):
        return x, y

class Constant(PixelMap):
    '''Constant shift pixel map
    '''
    @staticmethod
    def type():
        return 'Constant'

    def __init__(self, name, **kwargs):
        super(Constant,self).__init__(name)
        if 'Parameters' not in kwargs:
            raise TypeError('Missing Parameters in Constant PixelMap spec')
        if len(kwargs['Parameters']) !=2:
            raise TypeError('Wrong # of parameters in Constant PixelMap spec')
        self.shift = np.array(kwargs['Parameters'],dtype=float)

    def __call__(self, x, y, c=None):
        x += self.shift[0]
        y += self.shift[1]
        return x, y

class Linear(PixelMap):
    '''Affine transformation pixel map
    '''
    @staticmethod
    def type():
        return 'Linear'

    def __init__(self, name, **kwargs):
        super(Linear,self).__init__(name)
        if 'Coefficients' not in kwargs:
            raise TypeError('Missing Coefficients in Linear PixelMap spec')
        if len(kwargs['Coefficients']) !=6:
            raise TypeError('Wrong # of coefficients in Linear PixelMap spec')
        p = np.array(kwargs['Coefficients'],dtype=float).reshape((2,3))
        self.shift = p[:,0]
        self.m = p[:,1:]

    def __call__(self, x, y, c=None):
        # Note: with only a 2x2 matrix, it is faster to do by hand rather than using np.dot
        xw = self.m[0,0] * x
        yw = self.m[1,1] * y
        #xw += self.m[0,1] * y
        y *= self.m[0,1]
        xw += y
        #yw += self.m[1,0] * x
        x *= self.m[1,0]
        yw += x
        xw += self.shift[0]
        yw += self.shift[1]
        return xw, yw

class Polynomial(PixelMap):
    '''2d polynomial pixel map
    '''
    @staticmethod
    def type():
        return 'Poly'

    def __init__(self, name, **kwargs):
        super(Polynomial,self).__init__(name)
        # Read the scaling bounds that will map interval into [-1,+1]
        xmin = kwargs['XMin']
        xmax = kwargs['XMax']        
        ymin = kwargs['YMin']
        ymax = kwargs['YMax']        
        self.shift = np.array( [0.5*(xmin+xmax), 0.5*(ymin+ymax)])
        self.scale = np.array( [2./(xmax-xmin), 2./(ymax-ymin)])
        # Read coefficients, x then y
        self.coeffs = []
        for d in (kwargs['XPoly'], kwargs['YPoly']):
            sumOrder = bool(d['SumOrder'])
            xOrder = int(d['OrderX'])
            if sumOrder:
                yOrder = xOrder
            else:
                yOrder = int(d['OrderY'])
            c = np.zeros( (xOrder+1, yOrder+1), dtype=float)
            v = np.array(d['Coefficients'])
            # Here is where we are very specific to the ordering of coeffs
            # in the gbtools/utilities/Poly2d.cpp code:
            if sumOrder:
                # Ordering is 1, x, y, x^2, xy, y^2, x^3, xy^2
                ix = 0
                iy = 0
                N = 0
                for cc in v:
                    c[ix,iy] = cc
                    if ix>0:
                        iy = iy+1
                        ix = ix-1
                    else:
                        N = N+1
                        ix = N
                        iy = 0
            else:
                # ordering is 1., y, y^2, y^3 ... , x, xy, xy^2,
                c = v.reshape(xOrder, yOrder)
            self.coeffs.append(np.array(c,copy=False))

    def __call__(self, x, y, c=None):
        x -= self.shift[0]
        y -= self.shift[1]
        x *= self.scale[0]
        y *= self.scale[1]
        try:
            from galsim.utilities import horner2d
            xw = horner2d(x, y, self.coeffs[0])
            yw = horner2d(x, y, self.coeffs[1])
        except ImportError:
            xw = np.polynomial.polynomial.polyval2d(x, y, self.coeffs[0])
            yw = np.polynomial.polynomial.polyval2d(x, y, self.coeffs[1])
        return xw, yw

class Template(PixelMap):
    '''Mapping defined by scaling of a 1-dimensional lookup table
    '''
    @staticmethod
    def type():
        return 'Template'

    # A class variable contains all the templates read so far, since they
    # are very slow to get from YAML.
    libraries = {}
    
    def __init__(self, name, **kwargs):
        super(Template,self).__init__(name)
        if kwargs['HasSplit']:
            raise ValueError('Template pixel map not coded for split templates')
        fname = kwargs['Filename']
        if fname not in self.libraries:
            # Read in a new library file
            path = None
            if os.path.isabs(fname):
                # Just attempt read from an absolute path or if there's no path
                path = fname
            elif os.path.isfile(os.path.join(data_dir, fname)):
                # If the file is in our data directory, use that.
                path = os.path.join(data_dir, fname)
            else:
                # Search along path for the file
                path1 = [files.data_dir]
                if 'CAL_PATH' in os.environ:
                    path1 += os.environ['CAL_PATH'].split(':')
                # And last resort will be current directory
                path1.append('')
                for p in path1:
                    if os.path.isfile(os.path.join(p,fname)):
                        path = os.path.join(p,fname)
                        break
            if path is None:
                raise IOError('Can not find template library file ' + fname)
            with open(path) as f:
                self.libraries[fname] = yaml.load(f,Loader=Loader)

        # Now find the desired template
        if kwargs['LowTable'] not in self.libraries[fname]:
            raise RuntimeError('Did not find map ' + kwargs['LowTable'] + ' in file ' + fname)
        self.scale = float(kwargs['Parameter'])
        tab = self.libraries[fname][kwargs['LowTable']]
        
        self.axis = tab['Axis']  # 'X', 'Y', 'R"
        if self.axis=='R':
            self.center = np.array([tab['XCenter'],tab['YCenter']], dtype=float)
        self.vals = np.array(tab['Values'])
        self.args = float(tab['ArgStart']) + float(tab['ArgStep']) * \
          np.arange(self.vals.shape[0])

    def __call__(self, x, y, c=None):
        # Note that np.interp does not exploit the equal spacing of arguments???
        # C++ code Lookup1d.cpp uses endpoints when beyond table bounds, as does
        # np.interp() by default.
        if self.axis=='X':
            x1 = np.interp(x, self.args, self.vals)
            x1 *= self.scale
            x += x1
        elif self.axis=='Y':
            y1 = np.interp(y, self.args, self.vals)
            y1 *= self.scale
            y += y1
        elif self.axis=='R':
            xc = x - self.center[0]
            yc = y - self.center[1]
            rad = np.hypot(xc,yc)
            dr = np.interp(rad, self.args, self.vals)
            dr *= self.scale
            dr /= rad
            xc *= dr
            yc *= dr
            x += xc
            y += yc
        else:
            raise ValueError('Unknown Template axis type ' + self.axis)
        return x, y

class Piecewise(PixelMap):
    '''Mapping defined by piecewise-linear function
    '''
    @staticmethod
    def type():
        return 'Piecewise'

    def __init__(self, name, **kwargs):
        super(Piecewise,self).__init__(name)
        self.axis = kwargs['Axis']  # 'X', 'Y', 'R"
        if self.axis=='R':
            self.center = np.array([kwargs['XCenter'],kwargs['YCenter']], dtype=float)
        self.vals = np.array(kwargs['Parameters'])
        # First and last elements are always set to zero in C++:
        self.vals[0] = 0.
        self.vals[-1] = 0.
        self.args = float(kwargs['ArgStart']) + float(kwargs['ArgStep']) * \
          np.arange(self.vals.shape[0])

    def __call__(self, x, y, c=None):
        # Note that np.interp does not exploit the equal spacing of arguments???
        if self.axis=='X':
            x += np.interp(x, self.args, self.vals)
        elif self.axis=='Y':
            y += np.interp(y, self.args, self.vals)
        elif self.axis=='R':
            xc = x - self.center[0]
            yc = y - self.center[1]
            rad = np.hypot(xc,yc)
            dr = np.interp(rad, self.args, self.vals)
            dr /= rad
            xc *= dr
            yc *= dr
            x += xc
            y += yc
        else:
            raise ValueError('Unknown Piecewise axis type ' + self.axis)
        return x, y

class ColorTerm(PixelMap):
    '''Pixel map that is color times deviation of some other PixelMap
    '''
    @staticmethod
    def type():
        return 'Color'

    def __init__(self, name, pmap=None, **kwargs):
        '''Need to provide the modified atomic PixelMap as map argument.
        '''
        super(ColorTerm,self).__init__(name)
        if pmap is None:
            raise ValueError('No modified map specified for ColorTerm')
        self.pmap = pmap
        self.reference = float(kwargs['Reference'])

    def __call__(self, x, y, c=None):
        if c is None:
            raise ValueError('ColorTerm requires non-null c argument')
        x1 = np.array(x)
        y1 = np.array(y)
        xw, yw = self.pmap(x1, y1, c)
        xw -= x
        yw -= y
        xw *= (c-self.reference)
        yw *= (c-self.reference)
        xw += x
        yw += y
        return xw, yw

class Composite(PixelMap):
    '''Pixel map defined as composition of other PixelMaps
    '''
    @staticmethod
    def type():
        return 'Composite'

    def __init__(self,name,elements):
        '''Create a compound PixelMap from a list of PixelMap elements
        '''
        super(Composite,self).__init__(name)
        self.elements = elements

    def __call__(self, x, y, c=None):
        # call all elements in succession
        x,y = self.elements[0](x,y,c)
        for el in self.elements[1:]:
            x,y = el(x,y,c)
        return x,y

class WCS(PixelMap):
    '''Pixel map that is augmented by a projection from the world coordinates
    to definitive sky position, and toSky() will return RA, Dec arrays
    corresponding to the inputs.
    Also can alternatively supply a second projection such that __call__ method 
    returns the coordinates in the second projection's system.
    A projection is anything that has toSky() and toPix() methods going to/from
    RA, Dec.
    Scale factor is multiplied into world coordinates to turn them into degrees.
    '''
    def __init__(self, name, pmap, projection, scale=1.):
        super(WCS,self).__init__(name)
        self.pmap = pmap
        self.projection = projection
        self.scale = scale
        self.dest = None

    def reprojectTo(self, projection):
        self.dest = projection

    def toSky(self, x, y, c=None):
        ''' Return sky coordinates corresponding to array
        '''
        # This is the front end user interface.  Let's not modify the user's input x,y values.
        # Everything else is allowed to modify the inputs to make the outputs more efficiently.
        x1 = np.array(x, dtype=float)
        y1 = np.array(y, dtype=float)
        xw, yw = self.pmap(x1,y1,c)
        xw *= self.scale
        yw *= self.scale
        return self.projection.toSky(xw, yw)

    def toPix(self, ra, dec, c=None, guess=np.array([0.,0.])):
        ''' Return pixel coordinates corresponding to input RA, Dec.
        guess is a starting guess for solver, which can be either a single
        coordinate pair or an array matching coords length.
        '''
        xw, yw = self.projection.toXY(ra, dec)
        xw /= self.scale
        yw /= self.scale
        xp = np.empty_like(xw)
        yp = np.empty_like(yw)
        xp.fill(guess[0])
        yp.fill(guess[1])
        xp, yp = self.pmap.inverse(xw, yw, xp, yp, c)  # ??? tolerance???
        return xp, yp

    def __call__(self, x, y, c=None):
        '''Map the input coordinates to the coordinates either in the original
        projection, or in the reprojected system if one has been given.
        '''
        # This is the front end user interface.  Let's not modify the user's input x,y values.
        # Everything else is allowed to modify the inputs to make the outputs more efficiently.
        x1 = np.array(x, dtype=float)
        y1 = np.array(y, dtype=float)
        xw,yw = self.pmap(x1,y1,c)
        if self.dest is not None:
            # Execute reprojection
            ra, dec = self.projection.toSky(xw*self.scale, yw*self.scale)
            xw, yw = self.dest.toXY(ra, dec)
            xw /= self.scale
            yw /= self.scale
        return xw, yw

class PixelMapCollection(object):
    '''Class that holds a library of PixelMap/WCS specifications deserialized from
    a YAML file.  One can then request any PixelMap or WCS by name and be given a
    functional realization of map with that name.  Realizations are cached so that they
    are not remade every time.
    '''
    def __init__(self, filename=None, use_pkl=True):
        '''Create PixelMapCollection from the named YAML file

        If use_pkl=False, this will always read in the given `filename` YAML file.

        If use_pkl=True (the default), look for a pickle file named `filename + '.pkl'`,
        to read in instead, which tends to be much faster.  If it is not there, it
        will read in the YAML file and then write out the pkl file as a pickled version
        of the PixelMapCollection to significantly speed up subsequent I/O operations
        on this file.
        '''

        if filename is None:
            # Start with empty dictionary
            self.root = {}
        else:
            pkl_filename = filename + '.pkl'
            if use_pkl and  os.path.isfile(pkl_filename):
                with open(pkl_filename) as f:
                    self.root = pickle.load(f)
            else:
                with open(filename) as f:
                    self.root = yaml.load(f,Loader=Loader)
                if use_pkl:
                    with open(pkl_filename, 'wb') as f:
                        pickle.dump(self.root, f)

        # Extract the WCS specifications into their own dict
        if 'WCS' in self.root:
            self.wcs = self.root.pop('WCS')
        else:
            self.wcs = {}  # Empty dictionary
        self.realizedMap = {}
        self.realizedWCS = {}

    # Build a static dictionary of atomic PixelMap types
    atoms = {t.type():t for t in (Identity, Constant, Linear,
                                  Polynomial, Template, Piecewise)}

    def update(self, d):
        '''Read new map/wcs specs from a supplied dictionary.  Replaces
        any duplicate names. Exception is generated if a new map or wcs
        replaces one that has already been realized.  Dictionary will be altered
        '''
        # First check for overwrite of realized PixelMap:
        intersect = filter(self.realizedMap.has_key, d.keys())
        if intersect:
            raise ValueError('attempt to update already-realized PixelMaps: '+str(intersect))
        # And WCS:
        if 'WCS' in d:
            intersect = filter(self.realizedWCS.has_key, d['WCS'].keys())
            if intersect:
                raise ValueError('attempt to update already-realized WCS: '+str(intersect))
        # Proceed with updates
        if 'WCS' in d:
            self.wcs.update(d.pop('WCS'))
        self.root.update(d)
        
    def hasMap(self, name):
        return name in self.root

    def hasWCS(self, name):
        return name in self.wcs

    def parseAtom(self, name, **kwargs):
        '''Build a PixelMap realization for one of the atomic types
        from the specs in the kwargs dictionary
        '''
        if kwargs['Type'] not in self.atoms:
            raise ValueError('Unknown PixelMap atomic type ' + kwargs['Type'])
        return self.atoms[kwargs['Type']](name,**kwargs)
    
    def getMap(self, name):
        ''' Return realization of PixelMap with given name
        '''
        if name in self.realizedMap:
            return self.realizedMap[name]
        if name not in self.root:
            raise ValueError('No PixelMap with name ' + name)

        specs = self.root[name]  # The dict with specifications of this map
        if specs['Type']==ColorTerm.type():
            # Special procedure for color terms since we'll parse its captive
            pmap = ColorTerm(name, pmap=self.parseAtom('captive',**specs['Function']),
                                **specs)
        elif specs['Type']==Composite.type():
            # Another special procedure for compound maps, read all its members
            elements = [self.getMap(el) for el in specs['Elements']]
            pmap = Composite(name, elements)
        else:
            # Map is atomic
            pmap = self.parseAtom(name, **specs)
        # Add new map to those already realized
        self.realizedMap[name] = pmap
        return pmap            

    def getWCS(self, name):
        ''' Return realization of WCS with the given name
        '''
        if name in self.realizedWCS:
            return self.realizedWCS[name]
        if name not in self.wcs:
            raise ValueError('No WCS with name ' + name)

        specs = self.wcs[name]  # The dict with specifications of this map
        # First get the PixelMap that it modifies
        pmap = self.getMap(specs['MapName'])
        # Scale to turn world coords into degrees.  Note that C++ uses radians
        # for celestial coordinates and here we are in degrees by default, so
        # adjust the scale for this
        scale = specs['Scale'] * 180. / np.pi
        # Read the projection.  We recognize only 2 kinds right now:
        pspecs = specs['Projection']
        if pspecs['Type']==ICRS.type():
            projection = ICRS()
        elif pspecs['Type']==Gnomonic.type():
            ra = float(pspecs['Orientation']['RA'])
            dec = float(pspecs['Orientation']['Dec'])
            pa = float(pspecs['Orientation']['PA'])
            projection = Gnomonic(ra,dec,rotation=pa)
        else:
            raise ValueError('Do not know how to parse projection of type ',pspecs['Type'])
        # Make the WCS 
        w = WCS(name, pmap=pmap, projection=projection, scale=scale)
        # Cache it
        self.realizedWCS[name] = w
        return w
    
########################################################
### Projections: we have just two.
########################################################

class ICRS(object):
    ''' Class giving the (trivial) projection from ICRS to ICRS coordinates, i.e.
    "pixel" coordinates are just the ICRS RA and Dec in degrees.
    '''
    @staticmethod
    def type():
        return 'ICRS'

    def __init__(self):
        pass

    def toSky(self, x, y):
        return x, y

    def toXY(self, ra, dec):
        return ra, dec
    
class Gnomonic(object):
    ''' Class representing a gnomonic projection about some point on
    the sky.  Can be used to go between xi,eta coordinates and ra,dec.
    All xy units are assumed to be in degrees as are the ra, dec, and PA of
    the projection pole.
    '''
    @staticmethod
    def type():
        return 'Gnomonic'

    def __init__(self, ra, dec, rotation=0.):
        '''
        Create a Gnomonic transformation by specifying the position of the
        pole (in ICRS degrees) and rotation angle of the axes relative
        to ICRS north.
        '''
        self.pole_ra = ra
        self.pole_dec = dec
        self.rotation = rotation
        self.frame = None

    def _set_frame(self):
        pole = co.SkyCoord(self.pole_ra, self.pole_dec, unit='deg',frame='icrs')
        self.frame = pole.skyoffset_frame(rotation=co.Angle(self.rotation,unit='deg'))

    def toSky(self, x, y):
        '''
        Convert xy coordinates in the gnomonic project (in degrees) into ra, dec.
        '''
        try:
            import galsim
            pole = galsim.CelestialCoord(self.pole_ra * galsim.degrees,
                                         self.pole_dec * galsim.degrees)
            x *= -3600.  # GalSim wants these in arcsec, not degrees
            y *= 3600.   # Also, a - sign for x, since astropy uses +ra as +x direction.
            # apply rotation
            if self.rotation != 0.:
                # TODO: I'm not sure if I have the sense of the rotation correct here.
                #       The "complex wcs" test has PA = 0, so I wasn't able to test it.
                #       There may be a sign error on the s terms.
                s, c = (self.rotation * galsim.degrees).sincos()
                x, y = x*c - y*s, x*s + y*c
            # apply projection
            ra, dec = pole.deproject_rad(x, y, projection='gnomonic')
            ra *= galsim.radians / galsim.degrees
            dec *= galsim.radians / galsim.degrees
            return ra, dec

        except ImportError:
            if self.frame is None: self._set_frame()

            # Get the y and z components of unit-sphere coords, x on pole axis
            y, z = x, y
            y *= np.pi / 180.
            z *= np.pi / 180.
            temp = np.sqrt(1 + y*y + z*z)
            y /= temp
            z /= temp
            dec = np.arcsin(z)
            ra = np.arcsin(y / np.cos(dec))
            coord = co.SkyCoord(ra, dec, unit='rad', frame=self.frame)
            return coord.icrs.ra.deg, coord.icrs.dec.deg

    def toXY(self, ra, dec):
        '''
        Convert RA, Dec into xy values in the gnomonic projection, in degrees
        '''
        try:
            import galsim
            pole = galsim.CelestialCoord(self.pole_ra * galsim.degrees,
                                         self.pole_dec * galsim.degrees)
            ra *= galsim.degrees / galsim.radians
            dec *= galsim.degrees / galsim.radians
            # apply projection
            x, y = pole.project_rad(ra, dec, projection='gnomonic')
            x /= -3600.
            y /= 3600.
            # apply rotation
            if self.rotation != 0.:
                s, c = (self.rotation * galsim.degrees).sincos()
                x, y = x*c + y*s, -x*s + y*c
            return x, y

        except ImportError:
            if self.frame is None: self._set_frame()

            coord = co.SkyCoord(ra, dec, unit='deg')
            s = coord.transform_to(self.frame)
            # Get 3 components on unit sphere
            x = np.cos(s.lat.radian)*np.cos(s.lon.radian)
            y = np.cos(s.lat.radian)*np.sin(s.lon.radian)
            z = np.sin(s.lat.radian)
            out_x = y/x * (180. / np.pi)
            out_y = z/x * (180. / np.pi)
            return out_x, out_y
