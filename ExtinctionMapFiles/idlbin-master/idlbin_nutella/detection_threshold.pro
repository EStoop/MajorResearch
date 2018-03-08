pro detection_threshold

;fits_read, '/Users/jkrick/palomar/lfc/coadd_i.fits', data, ihead

;blank regions on the i-band image
xcenter = [3829.852,4462.3353,5130.2024,2958.5287,3431.7855,4608.2931,4847.1329,3515.8217,2839.1088,2666.6133]
ycenter = [4841.6314,4642.5982,4700.0967,4461.2568,4094.1511,4116.2659,3448.3988,3510.3202,3253.7885,4058.7674]
sig = xcenter

fits_read, '/Users/jkrick/hst/raw/smallacs.fits', data, acshead
;blank regions on the acs frame

xcenter=[400.56852,558.61918,132.91717,179.10406,279.19169,575.17303,823.07825,793.29751,109.47566,124.55716]

ycenter=[874.88436,831.06671,572.57533,410.62851,334.74395,351.08549,433.37159,243.3371,88.547955,884.56589]


;xcenter = [9581.5213,9892.7088,12349.012,11374.526,12726.839,9511.3768,11184.723,12255.744,11259.997,9579.3151,12660.107]
;ycenter = [10581.655,11286.015,11024.849,9711.7234,9441.2836,8973.573,7807.5031,8656.7395,11903.593,12493.654,11809.208]


;try a 6 pixel region
for i = 0, n_elements(xcenter) -1 do begin
   x = data[xcenter(i) - 10: xcenter(i) +10, ycenter(i)-10:ycenter(i) + 10]
   sig(i) = rms(x)
endfor

print, 'average sig', mean(sig)

end
