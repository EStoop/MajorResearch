pro make_star

array = PSF_GAUSSIAN( Npixel=256, FWHM=[4.3,3.6], centroid = [124.0, 124.0], /NORMAL )
mkhdr, hd, array
sxaddpar, hd, 'CHNLNUM' ,     1 ,' 1 digit instrument channel number  '            
sxaddpar, hd, 'GAIN' ,     3.7 ,' 1 digit instrument channel number  '            
sxaddpar, hd, 'EXPTIME' ,     1 ,' 1 digit instrument channel number  '            
sxaddpar, hd, 'FLUXCONV' ,     1 ,' 1 digit instrument channel number  '            

fits_read, '/Users/jkrick/irac_warm/latent/r38763008/ch1/bcd/SPITZER_I1_38763008_0257_0000_1_cbcd.fits', junk, header

fits_write, '/Users/jkrick/test.fits', array, header

 get_centroids,'/Users/jkrick/test.fits', t, dt, bcdxcen, bcdycen, bcdflux, xs, ys, fs, b, /WARM, /APER;, RA=ra(r), DEC=dec(r),/silent

print, bcdxcen, bcdycen

end
