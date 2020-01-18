Y6A1 Astrometry Solutions

*Release History:*

Version 1.1 28 July 2019 GMB:
* Re-solve with g-band differential chromatic refraction coefficient changed 45.0->39.2 mas/mag
* Also lower the minimum turbulence error ellipse eigenvalues from 5->3 mas.

Version 1.0 26 July 2019 GMB:
* First complete release for Y6
* Includes per-CCD residual maps and time evolution of affine transforms per CCD
* solutions for wide survey and shallow SNe

*Summary*
The files in this directory encode the astrometric solutions derived for all of the "good" wide-survey exposures in the Y6A1 internal data release.  These were derived by matching the Y6A1 FINALCUT cataloged positions of high-S/N stars for all exposures, along with the Gaia DR2 catalog, using the `WCSFit` software described in Bernstein et al (PASP 129:074503 2017) and further documented with the code at https://github.com/gbernstein/gbdes.  The code is updated to allow for proper motion and parallax of all stars, jointly constrained by DECam and DES data.

These files are meant to be used with the `pixmappy` Python package available at https://github.com/gbernstein/pixmappy.  See the documentation there for instructions for use.  The `pixmappy` code makes the use of these files a turnkey operation, so you do not need to read this doc any further in order to use these solutions.

*File Contents:*

_y6a1.guts.astro_: This is the YAML-format specification of the instrumental portion of the astrometric solutions.  This includes the per-CCD polynomials describing the optical distortions; the lateral color terms (g and r bands only); the "tree ring" and edge distortions in the detectors; and the "epoch ccdshift" terms specifying the positional shifts of the individual CCDs between temperature cyclings of the camera.

_astrorings4.yaml, astroedges1.yaml_: These give the empirically derived profiles of the tree-ring and edge distortions, respectively.

_y6a1.exposureinfo.fits_: A FITS binary table giving exposure-specific components of the astrometric solution.  The columns of the table are as follows:

* _expnum_: exposure number
* _xpoly, ypoly_: 10-element arrays giving the coefficients of cubic polynomials in (x,y) that describe the x and y components of the large-scale distortion induced by atmospheric refraction, stellar aberration, or other time-dependent large-scale effects.
* _dcr_: 3-element array describing the differential chromatic refraction.  The first two elements are the x (East) and y (N) coefficients of a linear stretch that is proportional to (c-c0), where c is the object g-i color, and c0 is a reference color given as the 3rd element of the array.  Units of the first two elements are degrees of sky position per mag of color.
* _ra, dec_: the coordinates (in degrees, ICRS) of the optic axis of the exposure, also used as the center of a gnomonic projection coordinate system for this exposure.
* _zone_: the zone number for which this solution was derived.  While an exposure may overlap multiple zones, we retain only the solution from the zone with the most overlap.
* _band_: the filter used for the exposure
* _epoch_: the (string-valued) camera epoch in which the exposure was taken, which specifies which set of ccdshift coefficients were used.
* _cov_: the covariance matrix of astrometric errors induced by atmospheric turbulence during this exposure. The array has the 3 elements (covxx,covyy,covxy). Each is in (mas)^2. The x and y axes point east and north, respectively.
* _observatory_: The Cartesian equatorial ICRS coordinates of the telescope at the midpoint of the exposure (in AU).
* _mjd:_ MJD at the midpoint of the exposure
* _nstars:_ Number of stars used in the astrometric solution
* _npairs:_ Number of pairs of stars used to calculate the error covariance
* _chisqred:_ Reduced chisquared of astrometric fit to this exposure.
* _clippct:_ Percent of stars clipped during astrometric fit.


_y6a1.astroresids.fits_: A FITS file containing lookup tables of small astrometric corrections for each CCD.  Each extension of the FITS file is an image of numpy shape (2,256,128).  The 256x128 dimensions represent the CCD (y,x) coordinates when binned into regions of 16x16 pixels.  The [0,:,:] matrix gives the x displacement (in pixels) which should be _subtracted_ from the x positions of objects in each pixel bin, and the [1,:,:] image gives the amount to be subtracted from the y positions.  These are derived from the mean astrometric residuals of bright stars in all griz band survey exposures.  They should probably be scaled down somewhat for Y band (and very slightly in z) but this is not done because of the small effect and the relative unimportance of Y for astrometry.
Each extension of the file is named by the DETPOS of the detector it maps.

_y6a1.affine.fits_: A FITS table giving a time-dependent affine transformation that should be applied to each CCD, derived again by the residuals seen in bright-star fits.  The pixel coordinates of an object should be transformed as

x -> x - (x0 + (mag+e1)*x + (rot + e2)*y)
y -> y - (y0 + (e2-rot)*x + (mag-e1)*y

The columns of the table are
* x0, y0, mag, e1, e2, rot as used above
* sig0, sigmag: formal errors in the determinations of (x0,y0) and (mag,rot,e1,e2) respectively.
* detpos:  DETPOS of relevant detector
* mjd: the starting MJD in which to use this correction.  It is valid until the next MJD present for this DETPOS.

_README.txt:_:  This file
