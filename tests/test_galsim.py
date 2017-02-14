# Test the GalSim interface to a PixelMapCollection
from __future__ import print_function
import pixmappy
import time
import numpy as np

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
    print('wcs.pmc = ',wcs.pmc)
    print('wcs.exp = ',wcs.exp)
    print('wcs.ccdnum = ',wcs.ccdnum)
    print('wcs.ccdname = ',wcs.ccdname)
    print('wcs.wcs_name = ',wcs.wcs_name)
    print('wcs.wcs_name = ',wcs.wcs_name)

    assert wcs.exp == exp
    assert wcs.ccdnum == ccdnum
    assert wcs.ccdname == ccdname
    assert wcs.wcs_name == 'D%s/%s'%(exp,ccdname)

    t0 = time.time()

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

if __name__ == '__main__':
    test_basic()
