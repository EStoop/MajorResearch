pro noise_exptime
!P.multi = [0,1,1]
!P.charsize = 1.5

;run multiple mosaics of the warm dark field with increasing exposure time
;to see if noise really does decrease as expected with exptime.

;first run column pulldown code on the dark field data
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38008576'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38008832'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38009088'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38009344'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38009600'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38009856'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38010112'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38010368'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38010624'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38010880'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38054912'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38073600'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38094080'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38160128'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38168064'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38178304'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38392320'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38395904'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38408448'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38412032'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38689280'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38692864'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38708224'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38711808'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38810624'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38814208'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38927360'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r38930944'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39105280'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39108864'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39122944'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39126528'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39143680'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39147264'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39205632'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39209216'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39219968'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39223552'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39339776'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39343360'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39353600'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39357184'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39426816'
;irac_warm_correct_pulldown_pipe_v2, '/Users/jkrick/irac_darks/warm_darks/r39430400'

; setting up input files & run mopex
cd, '/Users/jkrick/irac_darks/warm_darks'
command = './make_filelist'
;spawn, command

;for the warm image then I need to know which pixels are unmasked in
;all images.  Since they are all on the same FIF table this should be
;straightforward. we can just use the mask from 5 tto apply to the
;rest.

;read in the mosaics & their coverage maps
fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_8/Combine-mosaic/mosaic.fits', data_8, header_8
fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_8/Combine-mosaic/mosaic_cov.fits', covdata_8, covheader_8

fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_7/Combine-mosaic/mosaic.fits', data_7, header_7
fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_7/Combine-mosaic/mosaic_cov.fits', covdata_7, covheader_7

fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_6/Combine-mosaic/mosaic.fits', data_6, header_6
fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_6/Combine-mosaic/mosaic_cov.fits', covdata_6, covheader_6

fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_5/Combine-mosaic/mosaic.fits', data_5, header_5
fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_5/Combine-mosaic/mosaic_cov.fits', covdata_5, covheader_5

fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_4/Combine-mosaic/mosaic.fits', data_4, header_4
fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_4/Combine-mosaic/mosaic_cov.fits', covdata_4, covheader_4

fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_3/Combine-mosaic/mosaic.fits', data_3, header_3
fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_3/Combine-mosaic/mosaic_cov.fits', covdata_3, covheader_3

fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_2/Combine-mosaic/mosaic.fits', data_2, header_2
fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_2/Combine-mosaic/mosaic_cov.fits', covdata_2, covheader_2

fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex/Combine-mosaic/mosaic.fits', data_1, header_1
fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex/Combine-mosaic/mosaic_cov.fits', covdata_1, covheader_1

;mask sources
cd, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_5/Combine-mosaic'
command = 'sex mosaic.fits -c ~/irac_darks/warm_darks/default.sex'
;spawn, command

;first use the segmentation image as a mask for all of them
fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_5/Combine-mosaic/segmentation.fits', segdata, segheader
segmask = where(segdata gt 0)
data_8(segmask) = alog10(-1)
data_7(segmask) = alog10(-1)
data_6(segmask) = alog10(-1)
data_5(segmask) = alog10(-1)
data_4(segmask) = alog10(-1)
data_3(segmask) = alog10(-1)
data_2(segmask) = alog10(-1)
data_1(segmask) = alog10(-1)

;fits_write, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_5/Combine-mosaic/mosaic_mask.fits', data_5, header_5
;fits_write, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_4/Combine-mosaic/mosaic_mask.fits', data_4, header_4
;fits_write, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_3/Combine-mosaic/mosaic_mask.fits', data_3, header_3
;fits_write, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_2/Combine-mosaic/mosaic_mask.fits', data_2, header_2
;fits_write, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex/Combine-mosaic/mosaic_mask.fits', data_1, header_1


;-------------------------
a = where(finite(data_5) gt 0 and covdata_1 gt 11 and covdata_1 lt 13 and covdata_2 gt 23 and covdata_2 lt 35 and covdata_3 gt 35 and covdata_3 lt 37 and covdata_4 gt 47 and covdata_4 lt 49 and covdata_5 gt 59 and covdata_5 lt 61 and covdata_6 gt 71 and covdata_6 lt 73 and covdata_7 gt 83 and covdata_7 lt 85 and covdata_8 gt 95 and covdata_8 lt 97,count)
;a = where(finite(data_5) gt 0 and covdata_1 gt 14 and covdata_1 lt 16 and covdata_2 gt 29 and covdata_2 lt 31 and covdata_3 gt 44 and covdata_3 lt 46 and covdata_4 gt 59 and covdata_4 lt 61 and covdata_5 gt 74 and covdata_5 lt 76 and covdata_6 gt 89 and covdata_6 lt 91 and covdata_7 gt 104 and covdata_7 lt 106 and covdata_8 gt 119 and covdata_8 lt 121,count)
print, 'number of pixels that are finite and have coverage on all images', count

;what is the coverage of those same pixels in the mosaic_2?
;plothist, covdata_6(a)

;what is the noise in those same 83 pixels across all 8 images?
covarr = fltarr(8)
stddevarr = fltarr(8)

;------
covarr(0) = median(covdata_1(a))
plothist, data_1(a), xarr, yarr, bin = 0.001, xrange = [-0.02,0.02]
start = [0.00, 0.003, 1000]
noise = fltarr(n_elements(yarr))
noise[*] = 1
result= MPFITFUN('mygauss',xarr,yarr, noise, start,/quiet)    ;ICL
oplot, xarr,result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.), $
   thick = 1, color = redcolor
stddevarr(0) =abs(result(1))
;------
covarr(1) = median(covdata_2(a))
plothist, data_2(a), xarr, yarr, bin = 0.001, xrange = [-0.02,0.02]
start = [0.00, 0.003, 1000]
noise = fltarr(n_elements(yarr))
noise[*] = 1
result= MPFITFUN('mygauss',xarr,yarr, noise, start,/quiet)    ;ICL
oplot, xarr,result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.), $
   thick = 1, color = redcolor
stddevarr(1) =abs(result(1))
;------
covarr(2) = median(covdata_3(a))
plothist, data_3(a), xarr, yarr, bin = 0.001, xrange = [-0.02,0.02]
start = [0.00, 0.003, 1000]
noise = fltarr(n_elements(yarr))
noise[*] = 1
result= MPFITFUN('mygauss',xarr,yarr, noise, start,/quiet)    ;ICL
oplot, xarr,result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.), $
   thick = 1, color = redcolor
stddevarr(2) =abs(result(1))
;------
covarr(3) = median(covdata_4(a))
plothist, data_4(a), xarr, yarr, bin = 0.001, xrange = [-0.02,0.02]
start = [0.00, 0.003, 1000]
noise = fltarr(n_elements(yarr))
noise[*] = 1
result= MPFITFUN('mygauss',xarr,yarr, noise, start,/quiet)    ;ICL
oplot, xarr,result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.), $
   thick = 1, color = redcolor
stddevarr(3) =abs(result(1))
;------
covarr(4) = median(covdata_5(a))
plothist, data_5(a), xarr, yarr, bin = 0.001, xrange = [-0.02,0.02]
start = [0.00, 0.003, 1000]
noise = fltarr(n_elements(yarr))
noise[*] = 1
result= MPFITFUN('mygauss',xarr,yarr, noise, start,/quiet)    ;ICL
oplot, xarr,result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.), $
   thick = 1, color = redcolor
stddevarr(4) =abs(result(1))
;------
covarr(5) = median(covdata_6(a))
plothist, data_6(a), xarr, yarr, bin = 0.001, xrange = [-0.02,0.02]
start = [0.00, 0.003, 1000]
noise = fltarr(n_elements(yarr))
noise[*] = 1
result= MPFITFUN('mygauss',xarr,yarr, noise, start,/quiet)    ;ICL
oplot, xarr,result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.), $
   thick = 1, color = redcolor
stddevarr(5) =abs(result(1))
;------
covarr(6) = median(covdata_7(a))
plothist, data_7(a), xarr, yarr, bin = 0.001, xrange = [-0.02,0.02]
start = [0.00, 0.003, 1000]
noise = fltarr(n_elements(yarr))
noise[*] = 1
result= MPFITFUN('mygauss',xarr,yarr, noise, start,/quiet)    ;ICL
oplot, xarr,result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.), $
   thick = 1, color = redcolor
stddevarr(6) =abs(result(1))
;------
covarr(7) = median(covdata_8(a))
plothist, data_8(a), xarr, yarr, bin = 0.001, xrange = [-0.02,0.02]
start = [0.00, 0.003, 1000]
noise = fltarr(n_elements(yarr))
noise[*] = 1
result= MPFITFUN('mygauss',xarr,yarr, noise, start,/quiet)    ;ICL
oplot, xarr,result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.), $
   thick = 1, color = redcolor
stddevarr(7) =abs(result(1))




;-------------------------------
ps_open, filename='/Users/jkrick/irac_darks/warm_darks/noise_exptime_ch1.ps',/portrait,/color
redcolor = FSC_COLOR("Red", !D.Table_Size-2)

plot, covarr, stddevarr, psym = 2, xtitle = 'number of frames',$
      ytitle = 'Noise (Mjy/sr)', ystyle = 1,  yrange = [0.00, 0.005]


;try fitting it with sqrtfunc
start = [0.15]
noise = fltarr(n_elements(stddevarr))
noise[*] = 1
result= MPFITFUN('sqrtfunc',covarr,stddevarr, noise, start,/quiet)    ;ICL
oplot, covarr, result(0)/sqrt(covarr),  thick = 1, color = redcolor

ps_close, /noprint,/noid

;----------------------------------------------------------
;----------------------------------------------------------
;look at just the 8 mosaic and see what the different pixels show.

covarr = fltarr(100)
stddevarr = fltarr(100)
i = 0
for c = 12, 144, 12 do begin
   a = where (finite(data_8) gt 0 and covdata_8 gt c - 1. and covdata_8 lt c + 1., count)
   print, 'c, count', c, count
   covarr[i] = c
   
   plothist, data_8[a], xarr, yarr, bin = 0.001, xrange = [-0.02,0.02]
   start = [0.00, 0.003, 1000]
   noise = fltarr(n_elements(yarr))
   noise[*] = 1
   result= MPFITFUN('mygauss',xarr,yarr, noise, start,/quiet) ;ICL
   oplot, xarr,result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.), $
          thick = 1, color = redcolor
   
   stddevarr[i] = result(1); robust_sigma(data_8[a])
   i = i + 1
endfor

covarr = covarr[0:i-1]
stddevarr = stddevarr[0:i-1]

ps_open, filename='/Users/jkrick/irac_darks/warm_darks/noise_exptime_single_ch1.ps',/portrait,/color
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
plot, covarr, stddevarr, psym = 2, charthick = 1, yrange = [0.00, 0.005],xtitle = 'number of frames',$
      ytitle = 'Noise (Mjy/sr)'
;try fitting it with sqrtfunc
start = [0.015]
noise = fltarr(n_elements(stddevarr))
noise[*] = 1
result= MPFITFUN('sqrtfunc',covarr,stddevarr, noise, start,/quiet)    ;ICL
oplot, covarr, result(0)/sqrt(covarr),  thick = 1, color = redcolor
ps_close, /noprint,/noid

;----------------------------------------------------------
;----------------------------------------------------------
;look at the big cryo dark field

;read in the mosaics & their coverage maps
fits_read, '/Users/jkrick/nutella/spitzer/irac/ch1/mosaic.fits', data, header
fits_read, '/Users/jkrick/nutella/spitzer/irac/ch1/mosaic_cov.fits', covdata, covheader

;mask sources
cd, '/Users/jkrick/nutella/spitzer/irac/ch1'
command = 'sex mosaic.fits -c ~/irac_darks/warm_darks/default.sex'
spawn, command

;first use the segmentation image as a mask
fits_read, '/Users/jkrick/nutella/spitzer/irac/ch1/segmentation.fits', segdata, segheader
segmask = where(segdata gt 0)
data(segmask) = alog10(-1)

covarr = fltarr(100)
stddevarr = fltarr(100)
i = 0
for c = 12, 500, 12 do begin
   a = where (finite(data) gt 0 and covdata gt c - 1. and covdata lt c + 1., count)
   print, 'c, count', c, count
   covarr[i] = c
   
   plothist, data[a], xarr, yarr, bin = 0.0005
   ;start = [0.22, 0.003, 5000]
   start = [0.04, 0.003, 5000]

   noise = fltarr(n_elements(yarr))
   noise[*] = 1
   result= MPFITFUN('mygauss',xarr,yarr, noise, start,/quiet) ;ICL
   oplot, xarr,result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.), $
          thick = 1, color = redcolor
   
   stddevarr[i] = result(1); robust_sigma(data[a])
   i = i + 1
endfor

covarr = covarr[0:i-1]
stddevarr = stddevarr[0:i-1]

ps_open, filename='/Users/jkrick/irac_darks/warm_darks/noise_exptime_cryo_ch1.ps',/portrait,/color
redcolor = FSC_COLOR("Red", !D.Table_Size-2)

plot, covarr, stddevarr, psym = 2, charthick = 1, yrange = [0.000, 0.005],xtitle = 'number of frames',$
      ytitle = 'Noise (Mjy/sr)'
;try fitting it with sqrtfunc
start = [0.015]
noise = fltarr(n_elements(stddevarr))
noise[*] = 1
result= MPFITFUN('sqrtfunc',covarr,stddevarr, noise, start,/quiet)    ;ICL
oplot, covarr, result(0)/sqrt(covarr),  thick = 1, color = redcolor
ps_close, /noprint,/noid

end


;ok, but those are all different pixels, aka different wcs pixels, is
;there a way to do this on the same pixels.

;first use the segmentation image as a mask
;fits_read, '/Users/jkrick/irac_darks/warm_darks/ch1_mopex_5/Combine-mosaic/segmentation.fits', segdata, segheader
;segmask = where(segdata gt 0)
;data(segmask) = alog10(-1)

;covarr = fltarr(100)
;stddevarr = fltarr(100)
;i = 0
;for c = 10, 90, 10 do begin
;   a = where (finite(data) gt 0 and covdata gt c - 2. and covdata lt c + 2., count)
;   print, 'c, count', c, count
;   covarr[i] = c
;   stddevarr[i] = stddev(data[a])
;   i = i + 1
;endfor

;covarr = covarr[0:i-1]
;stddevarr = stddevarr[0:i-1]
;plot, covarr, stddevarr, psym = 2, charthick = 1
;oplot, covarr, .018*1/(sqrt(covarr)) 



;now try to do this with the real dark field image = cold processing
;with a lot more frames

