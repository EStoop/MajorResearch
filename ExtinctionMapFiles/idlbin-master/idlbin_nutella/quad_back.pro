pro quad_back
;determines the background values in each quadrant of ACS data.
;input a list of *_flt images
;outputs a file which can be used as input to iraf to subtract the determined backgrounds.

device, true=24
device, decomposed=0
colors = GetColor(/load, Start=1)
close, /all
ps_open, filename='/Users/jkrick/hst/bkg.ps',/portrait,/square,/color


openw, outlun, "/Users/jkrick/hst/raw/quad_back.minmax24.cl", /get_lun


readcol, '/Users/jkrick/hst/raw/bkglist24', filename,format="A"

for num  = 0, n_elements(filename) -1 do begin


print, "working on ", filename(num)
fits_read, filename(num), data, header,exten_no=1

;first half
plothistminmax, data[1:1350,900:*], xhist, yhist, bin=0.05, min=-100, max = 150, /noprint,/noplot;,xrange=[0,480]
plot, xhist, yhist,psym = 3, title = filename(num),xrange=[-100,100]
start = [mean(data[1:2048,*]),8.,24000.0]
err = fltarr(n_elements(xhist)) +1.0 
result1= MPFITFUN('gauss',xhist,yhist, err, start,/quiet)    ;ICL
oplot, xhist,  result1(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result1(0))/result1(1))^2.), $
   thick = 3, color = colors.green

;error check
if result1(1) lt 5 then print, "not a good fit 1",result1(0),  result1(1)
if result1(2) lt 10000 then print, "not a good fit 1", result1(0), result1(2)

;second half
plothistminmax, data[3300:*,*], xhist, yhist, bin=0.05, min=-100, max = 150, /noprint,/noplot;,xrange=[0,480]
plot, xhist, yhist, xrange=[-100,100], psym = 3
start = [mean(data[2049:*,*]),8.,24000.0]
err = fltarr(n_elements(xhist)) +1.0 
result2= MPFITFUN('gauss',xhist,yhist, err, start,/quiet)    ;ICL
oplot, xhist,  result2(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result2(0))/result2(1))^2.), $
  thick = 3, color = colors.green
;error check
if result2(1) lt 5 then print, "not a good fit 2",result2(0), result2(1)
if result2(2) lt 10000 then print, "not a good fit 2",result2(0), result2(2)

printf, outlun, "imdel test.fits"
printf, outlun, strcompress('imcalc ' + strcompress(filename(num)+'[1]',/remove_all) + ' test.fits "if x .lt. 2049 then ' +string(result1(0)) + ' else' + string(result2(0)) + '"')
printf, outlun, strcompress('imarith ' + strcompress(filename(num)+'[1]',/remove_all) +' - test.fits[1] ' + strcompress(filename(num)+'[1,overwrite]',/remove_all))



fits_read, filename(num), data2, header2,exten_no=4


;third half
plothistminmax, data2[1:2048,1:2047], xhist, yhist, bin=0.05, min=-100, max = 150, /noprint,/noplot,/nan
plot, xhist, yhist, xrange=[-100,100], psym = 3, title = filename(num)
start = [mean(data2[1:2048,*]),8.,24000.0]
err = fltarr(n_elements(xhist)) +1.0 
;print, n_elements(err)
result1= MPFITFUN('gauss',xhist,yhist, err, start,/quiet)    ;ICL
oplot, xhist,  result1(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result1(0))/result1(1))^2.), $
   thick = 3, color = colors.green
;error check
if result1(1) lt 5 then print, "not a good fit 3",result1(0), result1(1)
if result1(2) lt 10000 then print, "not a good fit 3",result1(0), result1(2)

;fourth half
plothistminmax, data2[2049:*,*], xhist, yhist, bin=0.05, /noprint,/noplot,/nan, min=-100, max = 150 
plot, xhist, yhist, xrange=[-100,100], psym = 3, title = filename(num)
start = [mean(data[2049:*,*]),8.,24000.0]
err = fltarr(n_elements(xhist)) +1.0 
result2= MPFITFUN('gauss',xhist,yhist, err, start,/quiet)    ;ICL
oplot, xhist,  result2(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result2(0))/result2(1))^2.), $
  thick = 3, color = colors.green
;error check
if result1(1) lt 5 then print, "not a good fit 4",result2(0), result2(1)
if result1(2) lt 10000 then print, "not a good fit 4",result2(0), result2(2)

printf, outlun, "imdel test.fits"
printf, outlun, strcompress('imcalc ' + strcompress(filename(num)+'[4]',/remove_all) + ' test.fits "if x .lt. 2049 then ' +string(result1(0)) + ' else' + string(result2(0)) + '"')
printf, outlun, strcompress('imarith ' + strcompress(filename(num)+'[4]',/remove_all) +' - test.fits[1] ' + strcompress(filename(num)+'[4,overwrite]',/remove_all))

endfor


close, outlun
free_lun, outlun

ps_close, /noprint,/noid


end
