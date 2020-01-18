# Test that PixelMapCollection gives expected answers for some reference WCS files.

from __future__ import print_function
import pixmappy
import numpy as np
import astropy
import os
import warnings

def test_tpv():
    """Test that reading a tpv file is equivalent to a regular TPV FITS wcs"""

    pmc = pixmappy.PixelMapCollection(os.path.join('input', 'tpv.yaml'), use_json=False)
    wcs1 = pmc.getWCS('testwcs')
    with warnings.catch_warnings():
        # Ignore warnings about RADECSYS -> RADESYSa.
        warnings.simplefilter("ignore")
        wcs2 = astropy.wcs.WCS(os.path.join('input','tpv.fits'))


    coords = [ (1322.1, 857.2), (1,1), (0,0), (943.234, 234.943), (2048, 2048) ]
    for x,y in coords:
        print('coord = ',x,y)
        sky1 = wcs1.toSky(x,y)
        sky2 = wcs2.wcs_pix2world(x, y, 1, ra_dec_order=True)
        print('  PMC: ',sky1[0], sky1[1])
        print('  astropy: ',sky2[0], sky2[1])
        np.testing.assert_allclose(sky1, sky2, rtol=1.e-8)

        # And reverse
        x1,y1 = wcs1.toPix(*sky1)
        np.testing.assert_allclose([x1,y1], [x,y], rtol=1.e-6, atol=1.e-8)

    # Now do all the coords at once
    x = [ c[0] for c in coords ]
    y = [ c[1] for c in coords ]
    ra, dec = wcs1.toSky(x,y)
    ra2, dec2 = wcs2.wcs_pix2world(x, y, 1, ra_dec_order=True)
    np.testing.assert_allclose(ra, ra2, rtol=1.e-8)
    np.testing.assert_allclose(dec, dec2, rtol=1.e-8)

    # This time, a fresh clone will not have the .json file yet, so maybe no difference from pmc.
    pmc2 = pixmappy.PixelMapCollection(os.path.join('input', 'tpv.yaml'), use_json=True)
    # This time will definitely read from the json file.
    pmc3 = pixmappy.PixelMapCollection(os.path.join('input', 'tpv.yaml'), use_json=True)

    wcs3 = pmc3.getWCS('testwcs')

    for x,y in coords:
        sky1 = wcs1.toSky(x,y)
        sky3 = wcs3.toSky(x,y)
        np.testing.assert_allclose(sky1, sky3, rtol=1.e-15)

        # And reverse
        x1,y1 = wcs1.toPix(*sky1)
        x3,y3 = wcs3.toPix(*sky1)
        np.testing.assert_allclose([x1,y1], [x3,y3], rtol=1.e-15, atol=1.e-15)



def test_complex():
    """Test a complex PMC file against some reference values"""

    pmc = pixmappy.PixelMapCollection(os.path.join('input', 'complex_wcs.astro'))
    wcs = pmc.getWCS('TEST/N1')

    ref = np.genfromtxt(os.path.join('input', 'complex_wcs.results'), names=True)

    for row in ref:
        print(row)
        ra, dec = wcs.toSky( row['xpix'], row['ypix'], c=row['color'])
        np.testing.assert_allclose(ra, row['RA'], rtol=1.e-6)
        np.testing.assert_allclose(dec, row['Dec'], rtol=1.e-6)

        x1, y1 = wcs.toPix(ra, dec, c=row['color'])
        np.testing.assert_allclose([x1, y1], [row['xpix'], row['ypix']], rtol=1.e-6)

if __name__ == '__main__':
    test_tpv()
    test_complex()
