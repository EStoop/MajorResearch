PRO ploticlV
;device, true=24
;device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva2/jkrick/A141/icl.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva2/jkrick/A141/original/iclVtab1', /GET_LUN

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

counts = counts[0:j-1]     ;get rid of the points which have bad stop codes
radius = radius[0:j-1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva2/jkrick/A141/original/iclVtab2', /GET_LUN

;read in the radial profile
rows= 10
radius2 = FLTARR(rows)
counts2 = FLTARR(rows)
Verr2 = FLTARR(rows)
stop2 = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius2(j) = r
      counts2(j) = c
      Verr2(j) = ellerr 
      stop2(j) = s
ENDFOR

close, lun
free_lun, lun

counts2 = counts2[0:j-1]     ;get rid of the points which have bad stop codes
radius2 = radius2[0:j-1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva2/jkrick/A141/original/iclVtab3', /GET_LUN

;read in the radial profile
rows= 11
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

counts3 = counts3[0:j-1]     ;get rid of the points which have bad stop codes
radius3 = radius3[0:j-1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


out = textoidl('Surface Brightness (mag/arcsec^{-2})')
plot, radius* 0.259,24.29 - 2.5*alog10(counts/(0.259^2)),thick = 3, $
YRANGE = [30, 26], xrange = [0,200],ytitle = out, $
psym = 3,xstyle = 8,ystyle = 1,$
xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3

axis, 0, 26.0, xaxis=1, xrange=[0,1.120], xstyle = 1, xthick = 3,charthick = 3
oplot, radius* 0.259,24.3 - 2.5*alog10(counts/(0.259^2)),thick = 5;
oplot, radius2* 0.259,24.3 - 2.5*alog10(counts2/(0.259^2)),thick = 5, color = colors.orange;
oplot, radius3* 0.259,24.3 - 2.5*alog10(counts3/(0.259^2)),thick = 5, color = colors.orchid;

;;;;;;


;;;;;;;

;;;;;;;;
;plot the error bars
;;;;;;;;

lower = 24.29 - 2.5*alog10((counts - 0.000518)/(0.259^2))
upper = 24.29 - 2.5*alog10((counts + 0.000518)/(0.259^2)) 
errplot, radius* 0.259,lower, upper

lower2 = 24.29 - 2.5*alog10((counts2 - 0.000518)/(0.259^2))
upper2 = 24.29 - 2.5*alog10((counts2 + 0.000518)/(0.259^2)) 
errplot, radius2* 0.259,lower2, upper2

lower3 = 24.29 - 2.5*alog10((counts3 - 0.000518)/(0.259^2))
upper3 = 24.29 - 2.5*alog10((counts3 + 0.000518)/(0.259^2)) 
errplot, radius3* 0.259,lower3, upper3

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;profile fitting
;;;;;;;;;;;;;;;;;;;;;;;;;;;

err = dindgen(j) - dindgen(j) + 1


;exp2
start = [0.03,235.0,0.008,690.0]  
;pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},4)
;pi(0).limited(0) = 1
;pi(0).limits(0) = 0.01
;pi(2).limited(1) = 1
;pi(2).limits(1) = 0.01
;result2exp = MPFITFUN('exp2',radius2[1:15],counts2[1:15], Verr2[1:15], start)      ;ICL

;exp
start = [0.03,300.0]
;resultexp = MPFITFUN('exponential',radius2[1:15],counts2[1:15], Verr2[1:15], start)    ;ICL

;result3 = MPFITFUN('exponential',radius3[1:15],counts3[1:15], err[0:18], start)    ;galaxies
;result2 = MPFITFUN('exponential',radius2[1:15],counts2[1:15], err[1:15], start)    ;ICL

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;profile plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
bigr = findgen(8880) / 10.

;2exponential
;numbers = result2exp(0)*exp(-bigr/(result2exp(1))) + result2exp(2)*exp(-bigr/(result2exp(3)))
;oplot, bigr*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), thick = 3, color = colors.gray


;exponential
;oplot, bigr*0.259, 24.3 - 2.5*alog10(resultexp(0)*exp(-bigr/resultexp(1))/(0.259^2) ), $
;   thick = 3, linestyle = 2;, color = colors.green


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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


device, /close
set_plot, mydevice

END




