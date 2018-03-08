pro noise_exptime
!P.multi = [0,1,1]
;!P.charsize = 1.5
ps_open, filename='/Users/jkrick/irac_darks/warm_darks/noise_exptime_ch2.ps',/portrait,/color

; setting up input files & run mopex
cd, '/Users/jkrick/irac_darks/warm_darks'
command = './make_filelist'
;spawn, command

;read in the big image
fits_read,  '/Users/jkrick/irac_darks/warm_darks/ch2_mopex/Combine-mosaic/mosaic.fits', data, header
fits_read,  '/Users/jkrick/irac_darks/warm_darks/ch2_mopex/Combine-mosaic/mosaic_cov.fits', covdata, covheader

;mask sources
cd, '/Users/jkrick/irac_darks/warm_darks/ch2_mopex/Combine-mosaic'
command = 'sex mosaic.fits -c ~/irac_darks/warm_darks/default.sex'
;spawn, command

;first use the segmentation image as a mask for all of them
fits_read, '/Users/jkrick/irac_darks/warm_darks/ch2_mopex/Combine-mosaic/segmentation.fits', segdata, segheader
segmask = where(segdata gt 0)
data(segmask) = alog10(-1)

;fits_write, '/Users/jkrick/irac_darks/warm_darks/ch2_mopex/Combine-mosaic/mosaic_mask.fits', data, header


;-------------------------------
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
greencolor = FSC_COLOR("GREEN", !D.Table_Size-3)

covarr = fltarr(100)
stddevarr = fltarr(100)
countarr = fltarr(100)
i = 0
for c = 12, 300, 12 do begin
   a = where (finite(data) gt 0 and covdata gt c - 1. and covdata lt c + 1., count)
   covarr[i] = c
   
   plothist, data[a], xarr, yarr, bin = 0.0005, xrange = [-0.02,0.02],/noplot
   start = [0.00, 0.003, 1000]
   noise = fltarr(n_elements(yarr))
   noise[*] = 1
   result= MPFITFUN('mygauss',xarr,yarr, noise, start,/quiet) ;ICL
;   oplot, xarr,result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.), $
;          thick = 1, color = redcolor
   print, 'c, count', c, count, result(1)

   stddevarr[i] = result(1); robust_sigma(data_8[a])
   countarr[i] = count
   i = i + 1
endfor

covarr = covarr[0:i-1]
stddevarr = stddevarr[0:i-1]
countarr = countarr[0:i-1]

plot, covarr, stddevarr, psym = 2, charthick = 3, yrange = [0.00, 0.005],xtitle = 'Number of Frames',$
      ytitle = 'Noise (Mjy/sr)', xthick = 3, ythick = 3, thick = 3
;try fitting it with sqrtfunc
start = [0.015]
noise = fltarr(n_elements(stddevarr))
noise[*] = 1
;noise = countarr
result= MPFITFUN('sqrtfunc',covarr,stddevarr, noise, start)    ;ICL
oplot, covarr, result(0)/sqrt(covarr),  thick = 3
;ps_close, /noprint,/noid

;----------------------------------------------------------

;----------------------------------------------------------
;look at the big cryo dark field
print, 'cryo-------------'
;read in the mosaics & their coverage maps
fits_read, '/Users/jkrick/nutella/spitzer/irac/ch2/mosaic.fits', data, header
fits_read, '/Users/jkrick/nutella/spitzer/irac/ch2/mosaic_cov.fits', covdata, covheader

;mask sources
cd, '/Users/jkrick/nutella/spitzer/irac/ch2'
command = 'sex mosaic.fits -c ~/irac_darks/warm_darks/default.sex'
;spawn, command

;first use the segmentation image as a mask
fits_read, '/Users/jkrick/nutella/spitzer/irac/ch2/segmentation.fits', segdata, segheader
segmask = where(segdata gt 0)
data(segmask) = alog10(-1)

covarr = fltarr(100)
stddevarr = fltarr(100)
countarr = fltarr(100)
i = 0
for c = 12, 300, 12 do begin
   a = where (finite(data) gt 0 and covdata gt c - 1. and covdata lt c + 1., count)
   covarr[i] = c
   
   plothist, data[a], xarr, yarr, bin = 0.0005,/noplot
   start = [0.22, 0.003, 5000]
   ;start = [0.04, 0.003, 5000]

   noise = fltarr(n_elements(yarr))
   noise[*] = 1
   result= MPFITFUN('mygauss',xarr,yarr, noise, start,/quiet) ;ICL
;   oplot, xarr,result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.), $
;          thick = 1, color = redcolor
   print, 'c, count', c, count, result(1)

   stddevarr[i] = result(1); robust_sigma(data[a])
   countarr[i] = count
   i = i + 1
endfor

covarr = covarr[0:i-1]
stddevarr = stddevarr[0:i-1]
countarr = countarr[0:i-1]

oplot, covarr, stddevarr, psym = 6, thick = 3
;try fitting it with sqrtfunc
start = [0.015]
noise = fltarr(n_elements(stddevarr))
noise[*] = 1
;noise = countarr
result= MPFITFUN('sqrtfunc',covarr,stddevarr, noise, start)    ;ICL
oplot, covarr, result(0)/sqrt(covarr),  thick =3, linestyle = 2


ps_close, /noprint,/noid

end
