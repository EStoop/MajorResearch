pro size_of_image

fits_read, '/Users/jkrick/Virgo/IRAC/ch1/ch1_Combine-mosaic/mosaic.fits', data, header
pixel_scale = abs(sxpar(header, 'PXSCAL2'))  ;need the absolute value, sometimes they are negative.

a = where(finite(data) gt 0)

print, n_elements(a)*(pixel_scale^2) / 60.^2 / 60.^2, 'square degrees'



end
