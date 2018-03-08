pro fit


device, true=24
device, decomposed=0
colors = GetColor(/Load, Start=1)
close, /all		;close all files = tidiness

;set up for plotting
mydevice = !D.NAME
!p.multi = [0, 0, 3]
SET_PLOT, 'ps'


device, filename = '/n/Sheriff3/jkrick/sep99/linearity.ps', /portrait, $
                BITS=8, scale_factor=0.9 , /color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;FIRST LINTEST
OPENR, lun,'/n/Sheriff3/jkrick/sep99/sep02/lintest1', /GET_LUN

;read in the radial profile
rows= 5
exptime = FLTARR(rows)
counts = FLTARR(rows)

e = 0.0				;radius
c = 0.0				;mean counts

FOR j=0,rows-1 DO BEGIN
      READF, lun, e,c
      exptime(j) = e
      counts(j) = c
ENDFOR

close, lun
free_lun, lun

err = dindgen(rows) - dindgen(rows) + 1

sortindex = Sort(exptime)
sortedexptime = exptime[sortindex]
sortedcounts = counts[sortindex]

yaxis = FLTARR(rows)
yaxis = sortedcounts/sortedexptime

normalized = yaxis(j-1)

FOR n = 0,j-1,1 DO BEGIN
    yaxis(n) = yaxis(n) / normalized
ENDFOR

;this is the curve fitting part
;give it starting values, and then ask it to find the best fit
start = [1E-4, 0.99]
result1 = MPFITFUN('linear',sortedcounts, yaxis, err, start)
start = [60,.0001,1E-6]
result2 = MPFITFUN('polynomial',sortedcounts, yaxis, err, start)
print, "linear result", result1
print, "polynomial result", result2

;SECOND TAKE THROUGH TIME SERIES
OPENR, lun2,'/n/Sheriff3/jkrick/sep99/sep02/lintest2', /GET_LUN

;read in the radial profile
rows= 5
exptime2 = FLTARR(rows)
counts2 = FLTARR(rows)

e2 = 0.0				;radius
c2 = 0.0				;mean counts

FOR j=0,rows-1 DO BEGIN
      READF, lun, e2,c2
      exptime2(j) = e2
      counts2(j) = c2
ENDFOR

close, lun2
free_lun, lun2

err2 = dindgen(rows) - dindgen(rows) + 1

sortindex2 = Sort(exptime2)
sortedexptime2 = exptime2[sortindex2]
sortedcounts2 = counts2[sortindex2]

yaxis2 = FLTARR(rows)
yaxis2 = sortedcounts2/sortedexptime2

normalized2 = yaxis2(j-1)

FOR n = 0,j-1,1 DO BEGIN
    yaxis2(n) = yaxis2(n) / normalized2
ENDFOR

;this is the curve fitting part
;give it starting values, and then ask it to find the best fit
start2 = [1E-4, 0.99]
result21 = MPFITFUN('linear',sortedcounts2, yaxis2, err2, start2)
start2 = [60,.0001,1E-6]
result22 = MPFITFUN('polynomial',sortedcounts2, yaxis2, err2, start2)
print, "linear result", result21
print, "polynomial result", result22




;plot the data = black, and the fit = blue
plot, sortedcounts , yaxis, psym = 2,thick = 3, $
	YRANGE = [0.96,1.04] ,xrange = [0,1.2E4],$
	ytitle = 'normalized counts/second', title = 'Sep02,1999 V', $
  xtitle = 'counts'
	


oplot, sortedcounts ,result1(0)*sortedcounts + result1(1),$
	thick = 3,color = colors.blue
;oplot, sortedcounts, result2(0) + result2(1)*sortedcounts + $
;  result2(2)*sortedcounts^2,thick = 3, color = colors.blue

;;;;
oplot, sortedcounts2 , yaxis2, psym = 2,thick = 3
;, $
;	YRANGE = [0.96,1.04] ,xrange = [0,1.2E4],$
;	ytitle = 'normalized counts/second', title = 'Sep02,1999 V', $
;  xtitle = 'counts'
	


;oplot, sortedcounts2 ,result21(0)*sortedcounts2 + result21(1),$
;  thick = 3,color = colors.red
;oplot, sortedcounts2, result22(0) + result22(1)*sortedcounts2 + $
;  result22(2)*sortedcounts2^2,thick = 3, color = colors.blue

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;SECOND LINTEST

OPENR, lun0,'/n/Sheriff3/jkrick/sep99/sep06/lintestV1', /GET_LUN

;read in the radial profile
rows= 12
exptime = FLTARR(rows)
counts = FLTARR(rows)

e = 0.0				;radius
c = 0.0				;mean counts

FOR j=0,rows-1 DO BEGIN
      READF, lun, e,c
      exptime(j) = e
      counts(j) = c
ENDFOR

close, lun0
free_lun, lun0

err = dindgen(rows) - dindgen(rows) + 1

sortindex = Sort(exptime)
sortedexptime = exptime[sortindex]
sortedcounts = counts[sortindex]

yaxis = FLTARR(rows)
yaxis = sortedcounts/sortedexptime

normalized = yaxis(j-1)

FOR n = 0,j-1,1 DO BEGIN
    yaxis(n) = yaxis(n) / normalized
ENDFOR

;this is the curve fitting part
;give it starting values, and then ask it to find the best fit
start = [1E-4, 0.99]
result1 = MPFITFUN('linear',sortedcounts, yaxis, err, start)
start = [60,.0001,1E-6]
result2 = MPFITFUN('polynomial',sortedcounts, yaxis, err, start)
print, "linear result", result1
print, "polynomial result", result2


;plot the data = black, and the fit = blue
plot, sortedcounts , yaxis, psym = 2,thick = 3, $
	yrange = [0.96,1.04],xrange = [0,1.2E4],$
	ytitle = 'normalized counts/second', title = 'Sep06,1999 V', $
  xtitle = 'counts'
	

;sortedcounts=[0,500,1000,1500,2000,4000,8000,12000]
oplot, sortedcounts ,result1(0)*sortedcounts + result1(1),$
	thick = 3,color = colors.red
;oplot, sortedcounts, result2(0) + result2(1)*sortedcounts + $
;  result2(2)*sortedcounts^2,thick = 3, color = colors.blue

oplot, sortedcounts2 ,result21(0)*sortedcounts2 + result21(1),$
  thick = 3,color = colors.blue


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;AFTER MAKING THE CORRECTION

OPENR, lun0,'/n/Sheriff3/jkrick/sep99/sep02/lintestafter', /GET_LUN

;read in the radial profile
rows= 10
exptime = FLTARR(rows)
counts = FLTARR(rows)

e = 0.0				;radius
c = 0.0				;mean counts

FOR j=0,rows-1 DO BEGIN
      READF, lun, e,c
      exptime(j) = e
      counts(j) = c
ENDFOR

close, lun0
free_lun, lun0

err = dindgen(rows) - dindgen(rows) + 1

sortindex = Sort(exptime)
sortedexptime = exptime[sortindex]
sortedcounts = counts[sortindex]

yaxis = FLTARR(rows)
yaxis = sortedcounts/sortedexptime

normalized = yaxis(j-1)

FOR n = 0,j-1,1 DO BEGIN
    yaxis(n) = yaxis(n) / normalized
ENDFOR

;this is the curve fitting part
;give it starting values, and then ask it to find the best fit
start = [1E-4, 0.99]
result1 = MPFITFUN('linear',sortedcounts, yaxis, err, start)
start = [60,.0001,1E-6]
result2 = MPFITFUN('polynomial',sortedcounts, yaxis, err, start)
print, "linear result", result1
print, "polynomial result", result2


;plot the data = black, and the fit = blue
plot, sortedcounts , yaxis, psym = 2,thick = 3, $
	yrange = [0.96,1.04],xrange = [0,1.2E4],$
	ytitle = 'normalized counts/second', title = 'Sep06,1999 V', $
  xtitle = 'counts'
	

;sortedcounts=[0,500,1000,1500,2000,4000,8000,12000]
oplot, sortedcounts ,result1(0)*sortedcounts + result1(1),$
	thick = 3,color = colors.red
;oplot, sortedcounts, result2(0) + result2(1)*sortedcounts + $
;  result2(2)*sortedcounts^2,thick = 3, color = colors.blue



device, /close
set_plot, mydevice

END

