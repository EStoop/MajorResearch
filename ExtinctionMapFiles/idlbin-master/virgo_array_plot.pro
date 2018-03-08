pro virgo_array_plot
!P.multi = [0,1,2]
ps_open, filename='/Users/jkrick/Virgo/IRAC/array_location/array_hist.ps',/portrait,/square,/color



;now read in the combined array location dependence image and do some statistics
fits_read, '/Users/jkrick/Virgo/IRAC/array_location/ch1/ch1_Combine-mosaic/mosaic.fits', data, header
fits_read, '/Users/jkrick/Virgo/IRAC/array_location/ch1/ch1_Combine-mosaic/mosaic_cov.fits', covdata, covheader

;get rid of the nan's
a = where(finite(data) ne 1)
data(a) = 0

;get rid of low coverage
b = where(covdata lt 12)
data(b) = 0

plothist, data, xarr, yarr, bin = 0.01, xrange = [0.9, 1.0],/ylog, xtitle = 'Array dependent correction (%)', ytitle = 'Number of pixels', title = 'CH1'

;now read in the combined array location dependence image and do some statistics
fits_read, '/Users/jkrick/Virgo/IRAC/array_location/ch2/ch2_Combine-mosaic/mosaic.fits', data, header
fits_read, '/Users/jkrick/Virgo/IRAC/array_location/ch2/ch2_Combine-mosaic/mosaic_cov.fits', covdata, covheader

;get rid of the nan's
a = where(finite(data) ne 1)
data(a) = 0

;get rid of low coverage
b = where(covdata lt 12)
data(b) = 0

plothist, data, xarr, yarr, bin = 0.01, xrange = [0.9, 1.0],/ylog, xtitle = 'Array dependent correction (%)', ytitle = 'Number of pixels', title = 'CH2'
ps_close, /noprint,/noid

end
