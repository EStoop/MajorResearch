PRO ploticlV
;device, true=24
;device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/final/twofitV.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva1/jkrick/A3888/final/iclVtab', /GET_LUN

;read in the radial profile
rows= 20-2
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

counts = counts[0:j-1] - 2    ;get rid of the points which have bad stop codes
radius = radius[0:j-1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;lower limit, 30% larger galaxy masks

OPENR, lun,'/n/Godiva1/jkrick/A3888/final/iclVtab2', /GET_LUN  ;iclVtab2

;read in the radial profile
rows=21
radius2 = FLTARR(rows)
counts2 = FLTARR(rows)
Verr2 = FLTARR(rows)
stop2 = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
i = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      IF c GT 2.0 THEN BEGIN
          radius2(i) = r
          counts2(i) = c
          Verr2(i) = ellerr 
          stop2(i) = s
          i = i +1
      ENDIF
ENDFOR

close, lun
free_lun, lun

counts2 = counts2[0:i - 2] -  2
radius2 = radius2[0:i-2]
;print, counts2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;extreme lower limit, 50% larger galaxy masks
OPENR, lun,'/n/Godiva1/jkrick/A3888/final/iclVtab3', /GET_LUN

;read in the radial profile
rows= 19-1
radius4 = FLTARR(rows)
counts4 = FLTARR(rows)
Verr4 = FLTARR(rows)
stop4 = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
l = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      IF c GT 2.0 THEN BEGIN
          radius4(l) = r
          counts4(l) = c
          Verr4(l) = ellerr 
          stop4(l) = s
          l = l +1
      ENDIF
ENDFOR

close, lun
free_lun, lun

counts4 = counts4[0:l - 2] -  2
radius4 = radius4[0:l-2]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; ellipse2 fits applied to
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;unmasked image

OPENR, lun,'/n/Godiva1/jkrick/A3888/final/galVtab2', /GET_LUN  ;galVtab2
;
rows= 21-1
radius3 = FLTARR(rows)
counts3 = FLTARR(rows)
Verr3 = FLTARR(rows)
stop3 = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius3(j) = r
      counts3(j) = c
      Verr3(j) = ellerr 
      stop3(j) = s
ENDFOR

close, lun
free_lun, lun

counts3 = counts3[0:j - 2] - 2
radius3 = radius3[0:j-2]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

out = textoidl('Surface Brightness (mag/arcsec^{-2})')
plot, radius2[2:i-2]* 0.259,24.29 - 2.5*alog10(counts2[2:i-2]/(0.259^2)),thick = 3, $
YRANGE = [30, 26], xrange = [0,250],ytitle = out, $
psym = 3,xstyle = 8,ystyle = 1,$
xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3

axis, 0, 26.0, xaxis=1, xrange=[0,1.034], xstyle = 1, xthick = 3,charthick = 3
oplot, radius2* 0.259,24.3 - 2.5*alog10(counts2/(0.259^2)),thick = 5;
;;;;;;
;pxval = [radius2[0:i - 2]* 0.259,reverse(radius4[0:l-2]* 0.259)]
;pyval = [24.29 - 2.5*alog10(counts2[0:i - 2]/(0.259^2)),reverse(24.29 - 2.5*alog10(counts4[0:l-2]/(0.259^2)))]
;;polyfill, pxval,pyval, color = colors.lightgray
;;;;polyfill, pxval,pyval, color = colors.aquamarine
;;polyfill, pxval,pyval, /line_fill, orientation = 45

;pxval = [radius2[0:i - 2]* 0.259,reverse(radius* 0.259)]
;pyval = [24.29 - 2.5*alog10(counts2[0:i - 2]/(0.259^2)),reverse(24.29 - 2.5*alog10(counts/(0.259^2)))]
;;polyfill, pxval,pyval, color = colors.darkgray
;;polyfill, pxval,pyval, color = colors.white, /line_fill, orientation = 45
;polyfill, pxval,pyval, color = colors.seagreen


;;;;;;;
;;oplot, radius2[0:i-2]* 0.259,24.29 - 2.5*alog10(counts2[0:i-2]/(0.259^2)),thick = 3;,linestyle = 2
;;oplot, radius3* 0.259,24.29 - 2.5*alog10(counts3/(0.259^2)),thick = 5;, linestyle = 2

;oplot, radius4* 0.259,24.29 - 2.5*alog10(counts4/(0.259^2)),thick = 3

;;;;;;;;
;plot the error bars
;;;;;;;;

;lower2 = 24.29 - 2.5*alog10((counts2[0:i - 2] - Verr2[0:i - 2])/(0.259^2))
;upper2 = 24.29 - 2.5*alog10((counts2[0:i - 2] + Verr2[0:i - 2])/(0.259^2)) 
lower2 = 24.29 - 2.5*alog10((counts2[0:i - 2] - 0.000555)/(0.259^2))
upper2 = 24.29 - 2.5*alog10((counts2[0:i - 2] + 0.000555)/(0.259^2)) 
errplot, radius2* 0.259,lower2, upper2
;print, lower2
;print, upper2

;newupper = (counts2[0:i - 2] + 0.000555)
;newlower = counts2[0:i - 2] - 0.000555
;lower3 = 24.29 - 2.5*alog10((counts3- Verr3)/(0.259^2))
;upper3 = 24.29 - 2.5*alog10((counts3 + Verr3)/(0.259^2)) 
;errplot, radius3* 0.259,lower3, upper3

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;profile fitting
;;;;;;;;;;;;;;;;;;;;;;;;;;;

err = dindgen(j) - dindgen(j) + .000555


;exp2
start = [0.5,105.0,0.004,545.0]  
;pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
;pi(0).limited(0) = 1
;pi(0).limits(0) = 0.01
;pi(2).limited(1) = 1
;pi(2).limits(1) = 0.01
print,radius2(16)
result2exp = MPFITFUN('exp2',radius2[1:16],counts2[1:16], err[1:16], start)      ;ICL

;exp
start = [0.03,300.0]
resultexp = MPFITFUN('exponential',radius2[1:16],counts2[1:16], err[1:16], start)    ;ICL

;result3 = MPFITFUN('exponential',radius3[1:18],counts3[1:18], err[0:18], start)    ;galaxies
;result2 = MPFITFUN('exponential',radius2[1:18],counts2[1:18], err[1:18], start)    ;ICL

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;profile plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
bigr = findgen(8300) / 10.

;2exponential
numbers = result2exp(0)*exp(-bigr/(result2exp(1))) + result2exp(2)*exp(-bigr/(result2exp(3)))
oplot, bigr*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), thick = 3, color = colors.gray


;exponential
oplot, bigr*0.259, 24.3 - 2.5*alog10(resultexp(0)*exp(-bigr/resultexp(1))/(0.259^2) ), $
   thick = 3, linestyle = 2;, color = colors.green


;exponential plotting
;;oplot, radius2[1:15]*0.259, 24.3 - 2.5*alog10(result2(0)*exp(-radius2[1:15]/result2(1))/(0.259^2) ), thick = 3;, color = colors.green
;;oplot, radius3[1:15]*0.259, 24.3 - 2.5*alog10(result3(0)*exp(-radius3[1:15]/result3(1))/(0.259^2) ), thick = 3;, color = colors.green

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

;numbers = (result3(0)) * (exp(-7.67*(((radius3[0:19]/result3(1))^(1.0/4.0)) - 1.0)))
;oplot, radius3[0:19]*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.green


;;;;;;;;;;;;;;;;;;;;;;;;;;;


;oplot, findgen(400) +230, findgen(400)  - findgen(400) + 29.5, thick = 3, linestyle = 2;, color = colors.seagreen
;oplot, findgen(100) +50, findgen(400)  - findgen(400) + 28.3, thick = 3, linestyle = 2;, color = colors.seagreen
;oplot, findgen(400), findgen(400)  - findgen(400) + 30.9, color = colors.darkgreen, thick = 3

;oplot, findgen(400) + 165 , findgen(400)  - findgen(400) + 28.8, thick = 3, linestyle = 2;, color = colors.seagreen


xarrjunk = fltarr(10)+215
yarrjunk = findgen(10) + 28.88;29.8
oplot, xarrjunk, yarrjunk
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


device, /close
set_plot, mydevice

END




