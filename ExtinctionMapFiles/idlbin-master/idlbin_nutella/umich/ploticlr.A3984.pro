PRO ploticl
;device, true=24
;device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva4/jkrick/A3984/profiles.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;r  curve
OPENR, lun,'/n/Godiva4/jkrick/A3984/iclrtab', /GET_LUN

;read in the radial profile
rows= 17  - 3
rradius = FLTARR(rows)
rcounts = FLTARR(rows)
rerr = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      rradius(j) = r
      rcounts(j) = c 
      rerr(j) = ellerr
ENDFOR

close, lun
free_lun, lun


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun,'/n/Godiva4/jkrick/A3984/iclr1tab', /GET_LUN

;read in the radial profile
rows= 15
rradius1 = FLTARR(rows)
rcounts1 = FLTARR(rows)
rerr1 = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      rradius1(j) = r
      rcounts1(j) = c 
      rerr1(j) = ellerr
ENDFOR

close, lun
free_lun, lun


;print, rcounts, rradius,rerr
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun,'/n/Godiva4/jkrick/A3984/iclr3tab', /GET_LUN

;read in the radial profile
rows= 12
rradius3 = FLTARR(rows)
rcounts3 = FLTARR(rows)
rerr3 = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      rradius3(j) = r
      rcounts3(j) = c 
      rerr3(j) = ellerr
ENDFOR

close, lun
free_lun, lun


;print, rcounts, rradius,rerr
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva4/jkrick/A3984/galrtab', /GET_LUN

;read in the radial profile
rows= 17 - 4
gradius = FLTARR(rows)
gcounts = FLTARR(rows)
gerr = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      gradius(j) = r
      gcounts(j) = c 
      gerr(j) = ellerr
ENDFOR

close, lun
free_lun, lun


;print, rcounts, rradius,rerr
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


out = textoidl('Surface Brightness (mag/arcsec^{-2})')
plot, rradius* 0.259,24.6 - 2.5*alog10(rcounts/(0.259^2)),thick = 3, $
YRANGE = [30, 24], xrange = [0,250],ytitle = out, $
psym = 3, xstyle = 9,ystyle = 1,$
xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3

axis, 0, 24.0, xaxis=1, xrange=[0,0.766], xstyle = 1, xthick = 3,charthick = 3

;;;;;;;;;;


oplot,  gradius* 0.259,24.6 - 2.5*alog10(gcounts/(0.259^2)),thick = 3, color = colors.red

pxval = [rradius* 0.259,reverse(rradius1* 0.259)]
pyval = [24.6 - 2.5*alog10(rcounts/(0.259^2)),reverse(24.6 - 2.5*alog10(rcounts1/(0.259^2)))]
polyfill, pxval,pyval, color = colors.lightsalmon

pxval = [rradius* 0.259,reverse(rradius3* 0.259)]
pyval = [24.6 - 2.5*alog10(rcounts/(0.259^2)),reverse(24.6 - 2.5*alog10(rcounts3/(0.259^2)))]
polyfill, pxval,pyval, color = colors.indianred

oplot, rradius* 0.259,24.6 - 2.5*alog10(rcounts/(0.259^2)),thick = 3, color = colors.black

lower2 = 24.6 - 2.5*alog10((rcounts - 0.00123)/(0.259^2))
upper2 = 24.6 - 2.5*alog10((rcounts + 0.00123)/(0.259^2)) 
errplot, rradius* 0.259,lower2, upper2

xarrjunk = fltarr(10)+202
yarrjunk = findgen(10) + 28.2;29.8
oplot, xarrjunk, yarrjunk
;;;;;;;;
;Fitting OF the data 

err = dindgen(j) - dindgen(j) + 1
start = [0.03,200.0]

result = MPFITFUN('exponential',rradius,rcounts, err, start)   
;result2 = MPFITFUN('exponential',rradius[3:16],rcounts[3:16], err[3:16], start)    ;ICL


start=[200., 500.]
;result3 = mpfitfun('devauc', rradius, rcounts, err, start)

;exponential plotting
oplot, findgen(900)*0.259, 24.6 - 2.5*alog10(result(0)*exp(-findgen(900)/result(1))/(0.259^2) ),$
  thick = 3, color = colors.red
;oplot, findgen(900)*0.259, 24.6 - 2.5*alog10(result2(0)*exp(-findgen(900)/result2(1))/(0.259^2) ),$
;  thick = 3, color = colors.tomato

;exp2
;start  = [30, 20., 0.3, 200.]
;result3 = MPFITFUN('exp2',rradius,rcounts, err, start)
;numbers = result3(0)*exp(-findgen(900)/(result3(1))) + $
;  result3(2)*exp(-findgen(900)/(result3(3)))
;oplot, findgen(900)*0.259, 24.6 - 2.5*alog10(numbers/(0.259^2)), $
;  thick = 3, color = colors.red



;reynolds plotting
;numbers = result2(0) / ((1 + (radius2/result2(1)))^2)
;oplot, radius2[0:13]*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.green

;numbers = result3(0) / ((1 + (radius3/result3(1)))^2)
;oplot, radius3[1:19]*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.green


;sersic plotting
;numbers = (result2(0)) * (exp(-7.67*(((radius2[0:13]/result2(1))^(1.0/result2(2))) - 1.0)))
;oplot, radius2[0:13]*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.green

;numbers = (result3(0)) * (exp(-7.67*(((radius3[1:19]/result3(1))^(1.0/result3(2))) - 1.0)))
;oplot, radius3[1:19]*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.green

;devauc plotting
;numbers = (result2(0)) * (exp(-7.67*(((radius2[0:13]/result2(1))^(1.0/4.0)) - 1.0)))
;oplot, radius2[0:13]*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.green

;numbers = (result3(0)) * (exp(-7.67*(((rradius/result3(1))^(1.0/4.0)) - 1.0)))
;oplot, rradius*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.green


;;;;;;;;;;;;;;;;;;;;;;;;;;;


;oplot, findgen(400) +230, findgen(400)  - findgen(400) + 29.5, thick = 3, linestyle = 2;, color = colors.seagreen
;oplot, findgen(100) +50, findgen(400)  - findgen(400) + 28.3, thick = 3, linestyle = 2;, color = colors.seagreen
;oplot, findgen(400), findgen(400)  - findgen(400) + 30.9, color = colors.darkgreen, thick = 3

;oplot, findgen(400) + 165 , findgen(400)  - findgen(400) + 28.8, thick = 3, linestyle = 2;, color = colors.seagreen


;########################################################################3

OPENR, lun,'/n/Godiva4/jkrick/A3984/iclVtab', /GET_LUN

;read in the radial profile
rows= 18
radius = FLTARR(rows)
counts = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius(j) = r
      counts(j) = c
      Verr(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

counts = counts[0:j-1] ;- 2    
radius = radius[0:j-1]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva4/jkrick/A3984/iclVtab3', /GET_LUN

;read in the radial profile
rows= 16-2
radius3 = FLTARR(rows)
counts3 = FLTARR(rows)
Verr3 = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius3(j) = r
      counts3(j) = c
      Verr3(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun,'/n/Godiva4/jkrick/A3984/iclVtab1', /GET_LUN

;read in the radial profile
rows= 16 
radius1 = FLTARR(rows)
counts1 = FLTARR(rows)
Verr1 = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius1(j) = r
      counts1(j) = c
      Verr1(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva4/jkrick/A3984/galVtab', /GET_LUN

;read in the radial profile
rows= 18
radius4 = FLTARR(rows)
counts4 = FLTARR(rows)
Verr4 = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius4(j) = r
      counts4(j) = c
      Verr4(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



pxval = [radius* 0.259,reverse(radius1* 0.259)]
pyval = [24.3 - 2.5*alog10(counts/(0.259^2)),reverse(24.3 - 2.5*alog10(counts1/(0.259^2)))]
polyfill, pxval,pyval, color = colors.palegreen

pxval = [radius* 0.259,reverse(radius3* 0.259)]
pyval = [24.3 - 2.5*alog10(counts/(0.259^2)),reverse(24.3 - 2.5*alog10(counts3/(0.259^2)))]
polyfill, pxval,pyval, color = colors.forestgreen

oplot, radius* 0.259,24.3 - 2.5*alog10(counts/(0.259^2)),$
  thick = 3;, color = colors.green  

lower2 = 24.3 - 2.5*alog10((counts - 0.000555)/(0.259^2))
upper2 = 24.3 - 2.5*alog10((counts + 0.000555)/(0.259^2)) 
errplot, radius* 0.259,lower2, upper2

;re-plot r-error bars so they are on top of V shading
lower2 = 24.6 - 2.5*alog10((rcounts - 0.00123)/(0.259^2))
upper2 = 24.6 - 2.5*alog10((rcounts + 0.00123)/(0.259^2)) 
errplot, rradius* 0.259,lower2, upper2

oplot, radius4* 0.259,24.3 - 2.5*alog10(counts4/(0.259^2)),$
  thick = 3, color = colors.green  

;exp
start = [0.03,200.0]
result = MPFITFUN('exponential',radius,counts, Verr, start)   
oplot, findgen(900)*0.259, 24.3 - $
  2.5*alog10(result(0)*exp(-findgen(900)/result(1))/(0.259^2) ), $
  thick = 3, color = colors.seagreen

;exp2
;start  = [30, 20., 0.3, 200.]
;result3 = MPFITFUN('exp2',radius[1:17],counts[1:17], Verr[1:17], start)
;numbers = result3(0)*exp(-findgen(900)/(result3(1))) + $
;  result3(2)*exp(-findgen(900)/(result3(3)))
;oplot, findgen(900)*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), $
;  thick = 3, color = colors.lawngreen

;DeVauc2
start = [0.005, 200.0, 0.4, 40.0]
pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
pi(2).limited(0) = 1
pi(2).limits(0) = 0.0
;pi(1).limited(1) = 1
;pi(1).limits(1) = 800.0
;resultdevi = MPFITFUN('devauc2',radius[0:17],counts[0:17], Verr[0:17], $
;                      start,parinfo=pi)    
;numbers = (resultdevi(0)) * (exp(-7.67*(((findgen(900)/resultdevi(1))^$
;          (1.0/4.0)) - 1.0))) + (resultdevi(2)) * (exp(-7.67*$
;          (((findgen(900)/resultdevi(3))^(1.0/4.0)) - 1.0)))
;oplot, findgen(900)*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), $
;          thick = 3, color = colors.hotpink









xyouts, 260, 24, "A3984"

device, /close
set_plot, mydevice

END



