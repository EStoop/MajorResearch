pro test_centroid
 
; staring mode ch1
ra_ref = 160.89641
dec_ref = 12.219967

fitsname = '/Users/jkrick/Desktop/SPITZER_I1_43336448_0001_0000_1_bcd.fits'
buncname = '/Users/jkrick/irac_warm/pcrs_planets/SD1043/r43336448/ch1/bcd/SPITZER_I1_43336448_0001_0000_1_bunc.fits'

;fitsname = '/Users/jkrick/irac_warm/pcrs_planets/SD1043/r43336448/ch1/bcd/SPITZER_I1_43336448_0001_0000_1_cbcd.fits'
;buncname = '/Users/jkrick/irac_warm/pcrs_planets/SD1043/r43336448/ch1/bcd/SPITZER_I1_43336448_0001_0000_1_cbunc.fits'

fits_read, fitsname, im, h
fits_read, buncname, unc, hunc

get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                     x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                     x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                     xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                     xfwhm, yfwhm, /WARM
      ;choose 3 pixel aperture, 3-7 pixel background
abcdflux = f[*,9]               ;put it back in same nomenclature
fs = fs[*,9]
x_center = x3
y_center = y3
back =  b[*,1]                  ;for 3, 11-15.5 ;b[*,2]  
backerr = bs[*,1]
        ;xfw = xfwhm
        ;yfw = yfwhm
;--------------------------------

print, x3, y3

end
