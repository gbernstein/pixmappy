# Test that PixelMapCollection gives expected answers for some reference WCS files.

from __future__ import print_function
import pixmappy
import numpy as np
import astropy
import os

def test_tpv():
    """Test that reading a tpv file is equivalent to a regular TPV FITS wcs"""

    pmc = pixmappy.PixelMapCollection(os.path.join('input', 'tpv.yaml'))
    wcs1 = pmc.getWCS('testwcs')
    wcs2 = astropy.wcs.WCS(os.path.join('input','tpv.fits'))


    coords = [ (1322.1, 857.2), (1,1), (0,0), (943.234, 234.943), (2048, 2048) ]
    for coord in coords:
        print('coord = ',coord)
        sky1 = wcs1.toSky(coord).icrs
        sky2 = wcs2.wcs_pix2world(coord[0], coord[1], 1, ra_dec_order=True)
        print('  PMC: ',sky1.ra.deg, sky1.dec.deg)
        print('  astropy: ',sky2[0], sky2[1])
        np.testing.assert_allclose(sky1.ra.deg, sky2[0], rtol=1.e-8)
        np.testing.assert_allclose(sky1.dec.deg, sky2[1], rtol=1.e-8)

        # And reverse
        pix1 = wcs1.toPix(sky1)
        np.testing.assert_allclose(pix1, coord, rtol=1.e-6, atol=1.e-8)

    # Now do all the coords at once
    all_sky = wcs1.toSky(coords)
    for sky1, coord in zip(all_sky, coords):
        sky2 = wcs2.wcs_pix2world(coord[0], coord[1], 1, ra_dec_order=True)
        np.testing.assert_allclose(sky1.icrs.ra.deg, sky2[0], rtol=1.e-8)
        np.testing.assert_allclose(sky1.icrs.dec.deg, sky2[1], rtol=1.e-8)


def test_complex():
    """Test a complex PMC file against some reference values"""

    pmc = pixmappy.PixelMapCollection(os.path.join('input', 'complex_wcs.astro'))
    wcs = pmc.getWCS('TEST/N1')

    ref = np.genfromtxt(os.path.join('input', 'complex_wcs.results'), names=True)

    for row in ref:
        print(row)
        coords = wcs.toSky( (row['xpix'], row['ypix']), c=row['color'])
        ra = coords.icrs.ra.deg
        dec = coords.icrs.dec.deg
        np.testing.assert_allclose(ra, row['RA'], rtol=1.e-6)
        np.testing.assert_allclose(dec, row['Dec'], rtol=1.e-6)

        pix1 = wcs.toPix(coords, c=row['color'])
        np.testing.assert_allclose(pix1, [row['xpix'], row['ypix']], rtol=1.e-6)

if __name__ == '__main__':
    test_tpv()
    test_complex()
