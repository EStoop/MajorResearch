pro flux_compare
ra1 = 273.04000000
dec1 = 63.49508333

badpix = [-9., 9.] * 1.D8

get_centroids, '/Users/jkrick/irac_warm/calstars/r37993728/ch1/bcd/SPITZER_I1_37993728_0002_0000_1_cbcd.fits', t, dt, x0, y0, f, xs,ys,fs,b,/warm,/aper,ra = ra1, dec=dec1

print, 'get_centroids', f * 1000., x0, y0


fits_read,  '/Users/jkrick/irac_warm/calstars/r37993728/ch1/bcd/SPITZER_I1_37993728_0002_0000_1_cbcd.fits', bcddata, bcdheader
adxy, bcdheader, ra1, dec1, bcdx, bcdy
cntrd, bcddata, bcdx, bcdy, bcdxcen, bcdycen, 5,/silent 
aper, bcddata, x0,y0, bcdflux, errap, sky, skyerr, 1, [10], [12.,20.],badpix, /NAN,/flux,/silent,/exact


scale = 1.22D * 1.22D
; convert scale from arcsec^2 to sr and scale to Jy
scale = scale * !DPI * !DPI / (3600.D * 3600.D * 180.D * 180.D) * 1.0D+06

bcdflux = bcdflux * scale * 1000.

print, 'aper', bcdflux, bcdxcen, bcdycen


end
