pro run_fixed_pos_aper

fits_read, '/Users/jkrick/external/dither_preaor/epic202083828/r54486528/ch2/bcd/SPITZER_I2_54486528_0010_0000_1_bcd.fits', im, h
fits_read, '/Users/jkrick/external/dither_preaor/epic202083828/r54486528/ch2/bcd/SPITZER_I2_54486528_0010_0000_1_bunc.fits', unc, hunc
xcen = 52.
ycen = 227.
fixed_pos_aper, im, h, unc, xcen, ycen, t, dt, hjd, xft, x3, y3, $
                       x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                       x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                       xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                       xfwhmarr, yfwhmarr, bb, /warm

;get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
;                                   x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
;                                   x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
;                                   xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
;                                   xfwhm, yfwhm,  /WARM

end
