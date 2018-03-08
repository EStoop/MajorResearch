pro areavsexptime
!P.multi = [0,1,1]
ps_open, filename='/Users/jkrick/nep/datapaper/coverage.ps',/portrait,/square,/color

;make a cumulative distribution of area vs. exposure time

fits_read, '/Users/jkrick/spitzer/irac/ch1/mosaic_cov.fits', data_ch1_orig, head_ch1
fits_read, '/Users/jkrick/spitzer/irac/ch2/mosaic_cov.fits', data_ch2_orig, head_ch2

data_ch1_orig = data_ch1_orig*100. /60. /60.
data_ch2_orig = data_ch2_orig*100. /60. /60.


;get rid of zeros and nans
data_ch1 = data_ch1_orig(where(data_ch1_orig gt 0))
;data_ch1 = data_ch1(where( finite(data_ch1) gt 0))
data_ch2 = data_ch2_orig(where(data_ch2_orig gt 0))
;data_ch2 = data_ch2(where( finite(data_ch2) gt 0))


;sort data
sortarea_ch1 = data_ch1(sort(data_ch1))
sortarea_ch2 = data_ch2(sort(data_ch2))


N1 = n_elements(sortarea_ch1)
N2 = n_elements(sortarea_ch2)

f1 = (findgen(N1) + 1.)/ N1
f2 = (findgen(N2) + 1.)/ N2

;plot, sortarea_ch1, f1*100, xtitle = 'Exposure Time (hours)', ytitle = 'Fraction of pixels', charthick = 3, xthick=3, ythick=3, thick=3
;oplot, sortarea_ch2, f2*100, linestyle = 2, thick = 3

;---------------------------------------
;only include the central 15'

naxis1 = fxpar(head_ch1, 'NAXIS1')
naxis2 = fxpar(head_ch1, 'NAXIS2')

keeper = fltarr(naxis1*naxis2)
count = long(0)

for x = 0, naxis1 -1 do begin
   for y = 0, naxis2 - 1 do begin
      if sqrt((x-1457.)^2. + (y-1680.)^2.)  le 750 then begin
         keeper[count] = data_ch1_orig[x, y]
         count = count + 1
;         print, 'x,y,data', x, y, sqrt((x-1457.)^2. + (y-1680.)^2.) ,data_ch1_orig[x,y]
      endif
   endfor
endfor

keeper = keeper[0:count-1]

sortarea_ch1 = keeper(sort(keeper))
;print, 'sortarea_ch1', sortarea_ch1
N1 = n_elements(sortarea_ch1)
f1 = (findgen(N1) + 1.)/ N1
plot, sortarea_ch1, f1*100,  thick=3,  xtitle = 'Exposure Time (hours)', ytitle = 'Fraction of Pixels', charthick = 3, xthick=3, ythick=3
;-----
naxis1 = fxpar(head_ch2, 'NAXIS1')
naxis2 = fxpar(head_ch2, 'NAXIS2')

keeper = fltarr(naxis1*naxis2)
count = long(0)

for x = 0, naxis1 -1 do begin
   for y = 0, naxis2 - 1 do begin
      if sqrt((x-1457.)^2. + (y-1680.)^2.)  le 750 then begin
         keeper[count] = data_ch2_orig[x, y]
         count = count + 1
;         print, 'x,y,data', x, y, sqrt((x-1457.)^2. + (y-1680.)^2.) ,data_ch1_orig[x,y]
      endif
   endfor
endfor

keeper = keeper[0:count-1]

sortarea_ch2 = keeper(sort(keeper))
;print, 'sortarea_ch1', sortarea_ch1
N2 = n_elements(sortarea_ch2)
f2 = (findgen(N2) + 1.)/ N2
oplot, sortarea_ch2, f2*100,  thick=3, linestyle = 2

legend, ['3.6 & 5.8 ', '4.5 & 8.0'], linestyle = [0,2], charthick=3, thick = 3, position = [20,40]

ps_close, /noprint,/noid


end
