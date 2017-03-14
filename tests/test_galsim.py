# Test the GalSim interface to a PixelMapCollection
from __future__ import print_function
import pixmappy
import time
import numpy as np
import os
import galsim

def test_basic():
    """Test basic operation of the GalSimWCS class """

    # Check that building a GalSimWCS builds successfully and has some useful attributes.
    input_dir = 'input'
    file_name = 'test.astro'
    exp = 375294
    ccdnum = 14
    ccdname = 'S15'

    t0 = time.time()
    wcs = pixmappy.GalSimWCS(file_name, dir=input_dir, exp=exp, ccdnum=ccdnum)
    t1 = time.time() - t0
    print('wcs = ',wcs)
    print('wcs.exp = ',wcs.exp)
    print('wcs.ccdnum = ',wcs.ccdnum)
    print('wcs.ccdname = ',wcs.ccdname)
    print('wcs.wcs_name = ',wcs.wcs_name)
    print('time to load = ',t1)

    assert wcs.exp == exp
    assert wcs.ccdnum == ccdnum
    assert wcs.ccdname == ccdname
    assert wcs.wcs_name == 'D%s/%s'%(exp,ccdname)

    t0 = time.time()
    wcs2 = pixmappy.GalSimWCS(file_name, dir=input_dir, exp=252223, ccdnum=11)
    t2 = time.time() - t0
    print('wcs2 = ',wcs2)
    print('time to load = ',t2)
    assert t2 < 0.1  # This should be fast since already in cache
    assert wcs2.exp == 252223
    assert wcs2.ccdnum == 11

    # Check that invalid initializations raise the appropriate errors
    np.testing.assert_raises(TypeError, pixmappy.GalSimWCS, file_name=file_name, pmc=wcs.pmc,
                             wcs_name=wcs.wcs_name)
    np.testing.assert_raises(TypeError, pixmappy.GalSimWCS, wcs_name=wcs.wcs_name)
    np.testing.assert_raises(TypeError, pixmappy.GalSimWCS, pmc=wcs.pmc,
                             exp=exp, ccdnum=ccdnum, wcs_name=wcs.wcs_name)
    np.testing.assert_raises(TypeError, pixmappy.GalSimWCS, pmc=wcs.pmc)
    np.testing.assert_raises(IOError, pixmappy.GalSimWCS, file_name=file_name)
    np.testing.assert_raises(IOError, pixmappy.GalSimWCS, file_name=file_name,
                             wcs_name=wcs.wcs_name)
    np.testing.assert_raises(IOError, pixmappy.GalSimWCS, file_name=file_name, dir='foo',
                             wcs_name=wcs.wcs_name)
    np.testing.assert_raises(IOError, pixmappy.GalSimWCS, file_name='foo', dir=input_dir,
                             wcs_name=wcs.wcs_name)


def test_tpv():
    """Test that reading a tpv file is equivalent to a regular TPV FITS wcs"""

    yaml_name = os.path.join('input','tpv.yaml')
    wcs1 = pixmappy.GalSimWCS(yaml_name, wcs_name='testwcs')
    fits_name = os.path.join('input','tpv.fits')
    wcs2 = galsim.FitsWCS(fits_name)

    coords = [ (1322.1, 857.2), (1,1), (0,0), (943.234, 234.943), (2048, 2048) ]
    for coord in coords:
        print('coord = ',coord)
        pos = galsim.PositionD(*coord)
        sky1 = wcs1.toWorld(pos)
        sky2 = wcs2.toWorld(pos)
        print('  GalSimWCS: ',sky1)
        print('  FitsWCS: ',sky2)
        np.testing.assert_allclose(sky1.ra.rad(), sky2.ra.rad(), rtol=1.e-8)
        np.testing.assert_allclose(sky1.dec.rad(), sky2.dec.rad(), rtol=1.e-8)

    # Now do all the coords at once
    xy = np.array(coords).T
    all_sky1 = wcs1._radec(xy[0], xy[1])
    all_sky2 = wcs2._radec(xy[0], xy[1])
    for sky1, sky2 in zip(all_sky1, all_sky2):
        np.testing.assert_allclose(sky1[0], sky2[0], rtol=1.e-8)
        np.testing.assert_allclose(sky1[1], sky2[1], rtol=1.e-8)


def test_complex():
    """Test a complex PMC file against some reference values"""

    wcs = pixmappy.GalSimWCS(os.path.join('input', 'complex_wcs.astro'), wcs_name='TEST/N1')
    ref = np.genfromtxt(os.path.join('input', 'complex_wcs.results'), names=True)

    for row in ref:
        print(row)
        pos = galsim.PositionD(row['xpix'], row['ypix'])
        sky = wcs.toWorld(pos, color=row['color'])
        ra = sky.ra / galsim.degrees
        dec = sky.dec / galsim.degrees
        np.testing.assert_allclose(ra, row['RA'], rtol=1.e-6)
        np.testing.assert_allclose(dec, row['Dec'], rtol=1.e-6)

    # This WCS requires a color term.  Raises an exception if you don't provide it.
    np.testing.assert_raises(Exception, wcs.toWorld, pos)

def test_cache():
    """Test the caching features of GalSimWCS"""

    class MockGalSimWCS(pixmappy.GalSimWCS):
        # Everything is the same, but we have our own cache dict so we don't have to
        # worry about messing up other tests that may be using the cache.
        cache = dict()

    wcs = MockGalSimWCS('input/test.astro', wcs_name='D231890/N1')
    assert len(wcs.cache) == 1
    assert 'input/test.astro' in wcs.cache
    wcs2 = MockGalSimWCS('input/test.astro', wcs_name='D469524/S13')
    assert len(wcs2.cache) == 1
    wcs3 = MockGalSimWCS('input/tpv.yaml', wcs_name='testwcs')
    assert len(wcs3.cache) == 2
    assert 'input/test.astro' in wcs.cache
    assert 'input/tpv.yaml' in wcs.cache
    wcs4 = MockGalSimWCS('input/complex_wcs.astro', wcs_name='TEST/N1', cache=False)
    assert len(wcs3.cache) == 2
    assert 'input/complex_wcs.astro' not in wcs.cache
    wcs.clear_cache()
    assert len(MockGalSimWCS.cache) == 0
    # Can also call as a class method:
    MockGalSimWCS.clear_cache()
    assert len(MockGalSimWCS.cache) == 0

def test_sky():
    """Test using the GalSimWCS to fill an image with constant surface brightness from the sky
    """
    wcs = pixmappy.GalSimWCS(dir='input', file_name='test.astro', exp=375294, ccdnum=14)
    print('wcs = ',wcs)

    sky_level = 177
    im = galsim.Image(2048, 4096)
    print('start making sky')
    wcs.makeSkyImage(im, sky_level)
    print('made sky')
    im.write('sky.fits')

    for x,y in [ (im.bounds.xmin, im.bounds.ymin),
                 (im.bounds.xmax, im.bounds.ymin),
                 (im.bounds.xmin, im.bounds.ymax),
                 (im.bounds.xmax, im.bounds.ymax),
                 (im.center().x, im.center().y) ]:
        val = im(x,y)
        area = wcs.pixelArea(galsim.PositionD(x,y))
        np.testing.assert_almost_equal(val/(area*sky_level), 1., 6,
                                       "Sky image at %d,%d is wrong"%(x,y))

if __name__ == '__main__':
    test_basic()
    test_tpv()
    test_complex()
    test_cache()
    test_sky()
