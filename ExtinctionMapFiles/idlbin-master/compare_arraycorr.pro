pro compare_arraycorr
fits_read, '/Users/jkrick/iwic/ch1_photcorr_ap_3.fits',ap3data, ap3head
fits_read, '/Users/jkrick/iwic/ch1_photcorr_ap_10.fits',ap10data, ap10head

diff = ap3data - ap10data
div = ap3data / ap10data

fits_write, '/Users/jkrick/iwic/testdiff.fits', diff, ap3head
fits_write, '/Users/jkrick/iwic/testdiv.fits', div, ap3head

end
