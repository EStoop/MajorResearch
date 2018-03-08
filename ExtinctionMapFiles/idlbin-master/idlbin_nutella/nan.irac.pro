fits_read, '/Users/jkrick/spitzer/irac/ch4/mosaic.fits', data, header
  data(where(~finite(data))) = 0
 fits_write, '/Users/jkrick/spitzer/irac/ch4/mosaic.nonan.fits', data, header
