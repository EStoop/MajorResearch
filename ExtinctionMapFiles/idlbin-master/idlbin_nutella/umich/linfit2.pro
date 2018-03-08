pro linfit2



colors = GetColor(/Load, Start=1)
close, /all		;close all files = tidiness


OPENR, lun,'/n/Sheriff1/jkrick/sep00/lintest/lintest', /GET_LUN

;read in the radial profile
rows= 13
exptime = FLTARR(rows)
counts = FLTARR(rows)
correct = FLTARR(rows)
average = FLTARR(rows)
midpt = FLTARR(rows)


e = 0.0				;radius
c = 0.0				;average counts
r = 0.0				;correct counts given 69.357783	
a = 0.0
m = 0.0	


FOR j=0,rows-1 DO BEGIN
      READF, lun, c,e,r,a,m

      exptime(j) = e
      counts(j) = c

ENDFOR

close, lun
free_lun, lun

err = dindgen(rows) - dindgen(rows) + 1

;this is the curve fitting part
;give it starting values, and then ask it to find the best fit
sortindex = Sort(exptime)
sortedexptime = exptime[sortindex]
sortedcounts = counts[sortindex]

print, "sortedcounts before normalization", sortedcounts

start = [60,.0001,1E-6]
result57 = MPFITFUN('polynomial',sortedcounts, sortedcounts/sortedexptime, err, start)
print, "result57", result57

yaxis = FLTARR(rows)
yaxis = sortedcounts/sortedexptime
normalized = yaxis(j-1)
FOR n = 0,j-1,1 DO BEGIN
	yaxis(n) = yaxis(n) / normalized
ENDFOR

start = [1E-4, 0.99]
result1 = MPFITFUN('linear',sortedcounts, yaxis, err, start)
start = [60,.0001,1E-6]
result2 = MPFITFUN('polynomial',sortedcounts, yaxis, err, start)
print, "result1", result1
print, "result2", result2


;set up for plotting
mydevice = !D.NAME
!p.multi = [0, 0, 3]
SET_PLOT, 'ps'


device, filename = '/n/Sheriff1/jkrick/sep00/lintest/linearity.ps', /portrait, $
                BITS=8, scale_factor=0.9 , /color


;plot the data = black, and the fit = blue
plot, sortedcounts , yaxis, psym = 2,thick = 3, $
	YRANGE = [0.96,1.02],XRANGE = [0,1.6E4],$
	ytitle = 'normalized counts/second', title = 'Sep00 '
	


oplot, sortedcounts ,result1(0)*sortedcounts + result1(1),$
	thick = 3,color = colors.red

oplot, sortedcounts, result2(0) + result2(1)*sortedcounts + result2(2)*sortedcounts^2,$
	thick = 3,color = colors.blue

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun0,'/n/Sheriff1/jkrick/sep00/lintest/lintestl', /GET_LUN
rows= 13
exptime = FLTARR(rows)
counts = FLTARR(rows)
correct = FLTARR(rows)
average = FLTARR(rows)
midpt = FLTARR(rows)


e = 0.0				;radius
c = 0.0				;average counts
r = 0.0				;correct counts given 69.357783	
a = 0.0
m = 0.0	


FOR j=0,rows-1 DO BEGIN
      READF, lun0, e,c

      exptime(j) = e
      counts(j) = c
;      correct(j) = r
;      average(j) = a
;      midpt(j) = m

ENDFOR


err = dindgen(rows) - dindgen(rows) + 1

;this is the curve fitting part
;give it starting values, and then ask it to find the best fit
sortindex = Sort(exptime)
sortedexptime = exptime[sortindex]
sortedcounts = counts[sortindex]

yaxis = FLTARR(rows)
yaxis = sortedcounts/sortedexptime
normalized = yaxis(j-1)
FOR n = 0,j-1,1 DO BEGIN
	yaxis(n) = yaxis(n) / normalized
ENDFOR

start = [1, 60]
result = MPFITFUN('linear',sortedcounts, yaxis, err, start)

print, result



plot, sortedcounts , yaxis, psym = 2,thick = 3, $
	YRANGE = [0.96,1.02],XRANGE = [0,1.6E4],$
	ytitle = 'normalized counts/second', $
	title = 'using the polynomial fit to correct for nonlinearity',$
	color =colors.blue


;oplot, sortedcounts ,result(1) + result(0)*sortedcounts ,$
;	thick = 3,color = colors.blue



;oplot, sortedcounts ,result1(0) + result1(1)*sortedcounts + result1(2)*(sortedcounts^2),$
;	thick = 2,color = colors.red, linestyle = 1
;oplot, sortedcounts ,result1(0)*sortedcounts + result1(1),$
;	thick = 3,color = colors.red, linestyle = 1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun1,'/n/Sheriff1/jkrick/sep00/lintest/lintestl2', /GET_LUN
rows= 13
exptime = FLTARR(rows)
counts = FLTARR(rows)

e = 0.0				;radius
c = 0.0				;average counts


FOR j=0,rows-1 DO BEGIN
      READF, lun1, e,c

      exptime(j) = e
      counts(j) = c

ENDFOR

print, "counts", counts

err = dindgen(rows) - dindgen(rows) + 1

;this is the curve fitting part
;give it starting values, and then ask it to find the best fit
sortindex = Sort(exptime)
sortedexptime = exptime[sortindex]
sortedcounts = counts[sortindex]

yaxis = FLTARR(rows)
yaxis = sortedcounts/sortedexptime
normalized = yaxis(j-1)
FOR n = 0,j-1,1 DO BEGIN
	yaxis(n) = yaxis(n) / normalized
ENDFOR

start = [60,.0001,1E-9]
result = MPFITFUN('polynomial',sortedcounts,yaxis, err, start)

print, result



plot, sortedcounts , yaxis, psym = 2,thick = 3,color = colors.red, $
	YRANGE = [0.96,1.02],XRANGE = [0,1.6E4],xtitle = 'counts',$
	ytitle = 'normalized counts/second',$
	 title = 'Using the linear fit to correct for nonlinearity'


;oplot, sortedcounts ,result(0) + result(1)*sortedcounts + result(2)*(sortedcounts^2), $
;	thick = 3,color = colors.blue

;oplot, sortedcounts ,result1(0) + result1(1)*sortedcounts + result1(2)*(sortedcounts^2),$
;	thick = 2,color = colors.red, linestyle = 1
;oplot, sortedcounts ,result1(0)*sortedcounts + result1(1),$
;	thick = 3,color = colors.red, linestyle = 1
close, lun0
free_lun, lun0

close, lun1
free_lun, lun1


device, /close
set_plot, mydevice
END
