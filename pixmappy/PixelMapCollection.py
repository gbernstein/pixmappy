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

class PixelMap(object):
    ''' Base class for transformations from one 2d system ("pixel") to another ("world").
    Each derived class must implement __call__ to execute this transform on array of
    shape (2) or (N,2).  This base class implements some routines such as taking local derivatives
    and solving for map inverse.
    c is always a color (or array (N) of them)
    '''
    def __init__(self, name):
        self.name = name

    def jacobian(self, x, y, c=None, step=1.):
        '''Use finite differences with indicated step size to get Jacobian derivative
        matrix at each point x, y.  If x, y are arrays then the returned array has
        shape (2,2,N) with [i,j,n] giving dx_i/dx_j at nth point.

        :param x,y: pixel coordinate arrays or scalars
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], any color terms in the map will raise an exception.
        :param step: Step size (in pixel units) for finite difference calculation.
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
        ''' Fill the arrays xp, yp with the solutions to PixelMap(xp, yp)=xw, yw.  

        :param xw,yw: input world coordinate arrays or scalars
        :param xp,yp: output pixel coordinates from inverse maps.  The input values are taken
                    as initial guesses for the root-finding routine.  
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], any color terms in the map will raise an exception.
        :param tol: tolerance for termination of solution (in pixel space??).  [default: 1e-4]
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
                raise ValueError('Mismatch of xw, yw, xp, yp, c point counts in PixelMap.inverse')

            # Allow for color to either be a scalar or a vector:
            if  c is None or type(c)==float or np.ndim(c)==0:
                for i in range(npts):
                    # c is None or a scalar:
                    xp[i], yp[i] = self.inverse(xw[i], yw[i], xp[i], yp[i], c, tol)
            else:
                # c should be an array matching xy's:
                if len(c) != npts:
                    raise ValueError('Mismatch of c array size to xy in PixelMap.inverse')
                for i in range(npts):
                    xp[i], yp[i] = self.inverse(xw[i], yw[i], xp[i], yp[i], c[i], tol)
            return xp, yp
    
'''
Following are various derived classes of PixelMaps.  They each have a class attribute
field "type" which corresponds to the "Type" field in the YAML node.
'''
    
class Identity(PixelMap):
    @staticmethod
    def type():
        return 'Identity'

    def __init__(self, name, **kwargs):
        '''Identity map - no free parameters
        :param name:  map name
        '''
        super(Identity,self).__init__(name)

    def __call__(self, x, y, c=None):
        '''Apply identity transformation to pixel coordinates.

        :param x,y: pixel coordinate arrays or scalars
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], any color terms in the map will raise an exception.
        :returns: xw, yw world coordinate arrays or scalars
        '''
        return x,y

class Constant(PixelMap):
    @staticmethod
    def type():
        return 'Constant'

    def __init__(self, name, **kwargs):
        '''
        Constant-shift pixel map.
        :param Parameters: 2-element array of constant offsets
        :param name:  map name
        '''
        super(Constant,self).__init__(name)
        if 'Parameters' not in kwargs:
            raise TypeError('Missing Parameters in Constant PixelMap spec')
        if len(kwargs['Parameters']) !=2:
            raise TypeError('Wrong # of parameters in Constant PixelMap spec')
        self.shift = np.array(kwargs['Parameters'],dtype=float)
        
    def __call__(self, x, y, c=None):
        '''Apply shift transformation to pixel coordinates.

        :param x,y: pixel coordinate arrays or scalars
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], any color terms in the map will raise an exception.
        :returns: xw, yw world coordinate arrays or scalars
        '''
        x += self.shift[0]
        y += self.shift[1]
        return x, y

class Linear(PixelMap):
    @staticmethod
    def type():
        return 'Linear'

    def __init__(self, name, **kwargs):
        '''Affine transformation pixel map
        xw = a_xx * xp + a_xy * yp + b_x,
        yw = a_yx * xp + a_yy * yp + b_y

        :param Coefficients: 6-element array [b_x, a_xx, a_xy, b_y, a_yx, a_yy]
        :param name:  map name
        '''
        super(Linear,self).__init__(name)
        if 'Coefficients' not in kwargs:
            raise TypeError('Missing Coefficients in Linear PixelMap spec')
        if len(kwargs['Coefficients']) !=6:
            raise TypeError('Wrong # of coefficients in Linear PixelMap spec')
        p = np.array(kwargs['Coefficients'],dtype=float).reshape((2,3))
        self.shift = p[:,0]
        self.m = p[:,1:]

    def __call__(self, x, y, c=None):
        '''Apply linear transformation to pixel coordinates.

        :param x,y: pixel coordinate arrays or scalars
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], any color terms in the map will raise an exception.
        :returns: xw, yw world coordinate arrays or scalars
        '''
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
        '''Polynomial transformation pixel map
        xp, yp are first linearly transformed so that
        (XMin,XMax) -> (-1,+1)
        (YMin,YMax) -> (-1,+1)
        Then xw and yw are produced as polynomial functions of the rescaled x,y.

        :param name:  map name
        :param XMin,XMax: x values that are rescaled to -1, +1
        :param YMin,YMax: y values that are rescaled to -1, +1
        :param XPoly: dictionary describing xw polynomial
        :param YPoly: dictionary describing yw polynomial

        The polynomial dictionaries have the following keys:
        OrderX:   Order of the polynomial in the rescaled x parameter
        OrderY:   Order of the polynomial in the rescaled y parameter.  Defaults to OrderX.
        SumOrder: boolean which, if true, limits polynomial terms to those
                  with sum of powers in x and y to be <=OrderX.
        Coefficients: array of polynomial coefficients, which are the free parameters
                  of the map.  Ordering is 1, x, y, x^2, xy, y^2, x^3, xy^2, ...,
                  i.e. first sort key is total order, second sort key is x order.
        '''
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
        '''Polynomial transformation of pixel coordinates.

        :param x,y: pixel coordinate arrays or scalars
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], any color terms in the map will raise an exception.
        :returns: xw, yw world coordinate arrays or scalars
        '''
        x =  x-self.shift[0]
        y =  y-self.shift[1]
        x =  x*self.scale[0]
        y =  y*self.scale[1]
        try:
            from galsim.utilities import horner2d
            xw = horner2d(x, y, self.coeffs[0])
            yw = horner2d(x, y, self.coeffs[1])
            # horner2d returns 0-d array if x,y are scalars.  Get them back
            if np.ndim(xw)==0:
                xw = xw.item()
                yw = yw.item()
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
        '''Transform pixel coordinates by adding multiple of a piecewise-linear
        function of either xp, yp, or radius from some center.  The displacement
        is also applied in either the x, y, or radial direction according to
        whether the argument of the LUT is x, y, or radius.

        Arguments outside the range of the LUT are assigned the values of the first
        or last node of the LUT.

        The map has one free parameter, which is a multiplicative factor to be
        applied to this lookup-table function.

        :param name:  map name
        :param Parameter: the multiplicative scaling to apply (free parameter)
        :param Filename:  the YAML file that contains library template specification
        :param LowTable:  the name (key) of the desired template in the YAML library

        The template library YAML file should be structured as a dictionary, with
        each entry also being a dictionary.  The dict for a single template should
        have these members:
        `Axis`:  one of ['X','Y','R'] specifying whether the argument of the lookup
                 function should be xp, yp, or `hypot(xp-x0,yp-y0)`.
        `XCenter`,`YCenter`: the values of x0,y0 used if `Axis='R'`
        `ArgStart`: the argument of the first node in the lookup table.
        `ArgStep`:  the spacing between nodes of the lookup table
        `Values`:   array of the values of the function at each node.
        '''

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
            elif os.path.isfile(os.path.join(files.data_dir, fname)):
                # If the file is in our data directory, use that.
                path = os.path.join(files.data_dir, fname)
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
        '''Add multiple of template to pixel coordinates.

        :param x,y: pixel coordinate arrays or scalars
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], any color terms in the map will raise an exception.
        :returns: xw, yw world coordinate arrays or scalars
        '''
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
        '''Transformation that adds a piecewise-linear
        function of either xp, yp, or radius from some center.  The displacement
        is also applied in either the x, y, or radial direction according to
        whether the argument of the LUT is x, y, or radius.

        The values at the first and last node of the LUT are fixed at zero, and
        arguments outside the range of the LUT are assigned zero displacement as
        well.  
        The free parameters of this map are the values at each node of the LUT
        (except for the first and last ones).

        :param name:  map name
        :param Parameter: the multiplicative scaling to apply (free parameter)
        :param Axis:  one of ['X','Y','R'] specifying whether the argument of the lookup
                 function should be xp, yp, or `hypot(xp-x0,yp-y0)`.
        `XCenter`,`YCenter`: the values of x0,y0 used if `Axis='R'`
        `ArgStart`: the argument of the first node in the lookup table.
        `ArgStep`:  the spacing between nodes of the lookup table
        `Parameters`:   array of the values of the function at each node.  The first
                    and last members will be set to zero.  The remaining members
                    will be free parameters of this map.
        '''
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
        '''Add piecewise-linear function to pixel coordinates.

        :param x,y: pixel coordinate arrays or scalars
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], any color terms in the map will raise an exception.
        :returns: xw, yw world coordinate arrays or scalars
        '''
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
        '''A PixelMap that is defined as the color of the source times
        the deviation of some other PixelMap,
        xw(x,y,c) = x + (c-Reference)*(pmap(x,y,c).x - x)
        and similarly for y.

        :param name: name for this transformation
        :param Reference: the reference color at which the shift is zero
        :param pmap: the `PixelMap` to scale by color
        '''

        super(ColorTerm,self).__init__(name)
        if pmap is None:
            raise ValueError('No modified map specified for ColorTerm')
        self.pmap = pmap
        self.reference = float(kwargs['Reference'])

    def __call__(self, x, y, c=None):
        '''Apply color term to pixel coordinates.

        :param x,y: pixel coordinate arrays or scalars
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], an ValueError will be raised.
        :returns: xw, yw world coordinate arrays or scalars
        '''
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
    @staticmethod
    def type():
        return 'Composite'

    def __init__(self,name,elements):
        '''PixelMap defined as composition of other PixelMaps

        :param name: the name of the compounded map
        :param elements: a list of `PixelMap` objects in the order that they will
        be applied to the pixel coordinates.
        '''
        super(Composite,self).__init__(name)
        self.elements = elements

    def __call__(self, x, y, c=None):
        '''Apply compound transformation to pixel coordinates.

        :param x,y: pixel coordinate arrays or scalars
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], any color terms in the map will raise an exception.
        :returns: xw, yw world coordinate arrays or scalars
        '''
        # call all elements in succession
        x,y = self.elements[0](x,y,c)
        for el in self.elements[1:]:
            x,y = el(x,y,c)
        return x,y

class WCS(PixelMap):
    '''A WCS is a PixelMap that is augmented by a projection from the world coordinates
    to definitive sky position (ra/dec). The method toSky() will return RA, Dec arrays
    corresponding to the input pixel coordinates.
    Also can alternatively supply a second projection such that __call__ method 
    returns the coordinates in the second projection's system.
    A projection is anything that has toSky() and toPix() methods going to/from
    RA, Dec.  The projections `ICRS` and `Gnomonic` are defined in this file.
    The `scale` factor is multiplied into world coordinates to turn them into degrees.
    '''
    def __init__(self, name, pmap, projection, scale=1.):
        '''Create a WCS with
        :param name: the WCS name.  The `WCS` name can duplicate a `PixelMap` name since
        these are maintained in separate namespaces.
        :param pmap: a `PixelMap` giving the transformation from pixel to world coordinates
        :param projection: the projection from world coordinates of the `PixelMap` into RA/Dec
        :param scale: factor to apply to world coordinates to turn them into degrees, applied
           before executing projection.
        '''
        super(WCS,self).__init__(name)
        self.pmap = pmap
        self.projection = projection
        self.scale = scale
        self.dest = None

    def reprojectTo(self, projection):
        '''(Re)define the WCS to yield world coordinates in the provided projection.
        Typically this will be a Gnomonic object with a new RA,Dec as projection origin.
        :param projection: The Projection used to map RA,Dec back into a world coordinate system.
        '''
        self.dest = projection

    def toSky(self, x, y, c=None):
        ''' Return ICRS sky coordinates corresponding to x,y pixel coordinates (arrays ok)
        :param x,y: pixel coordinate arrays or scalars
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], any color terms in the map will raise an exception.
        :returns ra,dec: sky coordinates (in degrees)
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
        ''' Return pixel coordinates x,y corresponding to input RA, Dec.x
        guess is a starting guess for solver, which can be either a single
        coordinate pair or an array matching coords length.
        :param ra,dec: sky coordinate arrays or scalars (assumed in degrees)
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], any color terms in the map will raise an exception.
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
        '''Map the input pixel coordinates to the world coordinates of the original
        `PixelMap`, or to the world coordinates of the reprojected system if one has been given.
        :param x,y: pixel coordinate arrays or scalars
        :param c:   color of source(s).  A scalar will be broadcast if x,y are arrays.
                    If `c=None` [default], any color terms in the map will raise an exception.
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
    are not remade every time.  The cache can be cleared to keep it from getting large.
    '''
    def __init__(self, filename=None, use_pkl=False):
        '''Create PixelMapCollection from the named YAML file

        If use_pkl=False, this will always read in the given `filename` YAML file.

        If use_pkl=True (the default), look for a pickle file named `filename + '.pkl'`
        to read instead, which tends to be much faster.  If `use_pkl=True` and
        no pkl file exists, it will read in the YAML file and then write out the pkl
        file as a pickled version of the PixelMapCollection to significantly speed up 
        subsequent I/O operations on this file.

        Specifications for new `PixelMap` or `WCS`'s can be added to this collection
        by passing a dictionary to the `update()` method.

        :param filename: path to YAML file specifying all `PixelMap`s and `WCS`s.  If `None`, begin
        with an empty collection (which is the default).
        :param use_pkl: `True` to read/write from pickled version of `filename`. [default=`False`]
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

    # Give other code the ability to add new atomic types:
    @classmethod
    def addAtom(cls, newAtom):
        cls.atoms[newAtom.type()] = newAtom
        
    def update(self, d):
        '''Read new map/wcs specs from a supplied dictionary.  Replaces
        any duplicate names. Exception is generated if a new map or wcs
        replaces one that has already been realized.
        :param d: the dictionary of new map/wcs specs, which will be altered.
        '''
        # First check for overwrite of realized PixelMap:
        intersect = set(self.realizedMap.keys()).intersection(set(d.keys()))
        if intersect:
            raise ValueError('attempt to update already-realized PixelMaps: '+str(intersect))
        # And WCS:
        if 'WCS' in d:
            intersect = set(self.realizedMap.keys()).intersection(set(d['WCS'].keys()))
            if intersect:
                raise ValueError('attempt to update already-realized WCS: '+str(intersect))
        # Proceed with updates
        if 'WCS' in d:
            self.wcs.update(d.pop('WCS'))
        self.root.update(d)
        
    def hasMap(self, name):
        '''Return `True` iff a `PixelMap` of the given name is in this collection
        '''
        return name in self.root

    def hasWCS(self, name):
        '''Return `True` iff a `WCS` of the given name is in this collection
        '''
        return name in self.wcs

    def parseAtom(self, name, **kwargs):
        '''Build a PixelMap realization with the given name for one of the atomic types
        from the specs in the kwargs dictionary.  Returns the `PixelMap` instance created.
        '''
        if kwargs['Type'] not in self.atoms:
            raise ValueError('Unknown PixelMap atomic type ' + kwargs['Type'])
        return self.atoms[kwargs['Type']](name,**kwargs)
    
    def getMap(self, name):
        ''' Return a realization of PixelMap with given name.  Build a realization
        and add it to the cache if it was not already there.
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
        ''' Return realization of WCS with the given name.  Build a realization
        and add it to the cache if it was not already there.
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

    def clearCache(self):
        '''Clear `PixelMap` and `WCS` caches.  Python garbage collection will get rid
        of the ones that are not currently in use.
        '''
        self.realizedWCS = {}
        self.realizedMap = {}
        return
    
########################################################
### Projections: we have just two.
########################################################

class ICRS(object):
    ''' Class giving the (trivial) projection from ICRS to ICRS coordinates, i.e.
    "world" coordinates are just the ICRS RA and Dec in degrees.
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
    the projection pole.  Uses astropy.coordinates.
    '''
    @staticmethod
    def type():
        return 'Gnomonic'

    def __init__(self, ra, dec, rotation=0.):
        '''
        Create a Gnomonic transformation by specifying the position of the
        pole (in ICRS degrees) and rotation angle of the axes relative
        to ICRS north.

        :param ra,dec: ICRS RA and Declination of the pole of the projection.
        :param rotation: position angle (in degrees) of the projection axes.
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
            import coord
            pole = coord.CelestialCoord(self.pole_ra * coord.degrees,
                                        self.pole_dec * coord.degrees)
            deg_per_radian = coord.radians / coord.degrees
            # Coord wants these in radians, not degrees
            # Also, a - sign for x, since astropy uses +ra as +x direction.
            x /= -deg_per_radian
            y /= deg_per_radian
            # apply rotation
            if self.rotation != 0.:
                # TODO: I'm not sure if I have the sense of the rotation correct here.
                #       The "complex wcs" test has PA = 0, so I wasn't able to test it.
                #       There may be a sign error on the s terms.
                s, c = (self.rotation * coord.degrees).sincos()
                x, y = x*c - y*s, x*s + y*c
            # apply projection
            ra, dec = pole.deproject_rad(x, y, projection='gnomonic')
            ra *= deg_per_radian
            dec *= deg_per_radian
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
            import coord
            pole = coord.CelestialCoord(self.pole_ra * coord.degrees,
                                        self.pole_dec * coord.degrees)
            deg_per_radian = coord.radians / coord.degrees
            ra /= deg_per_radian
            dec /= deg_per_radian
            # apply projection
            x, y = pole.project_rad(ra, dec, projection='gnomonic')
            x *= -deg_per_radian
            y *= deg_per_radian
            # apply rotation
            if self.rotation != 0.:
                s, c = (self.rotation * coord.degrees).sincos()
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
