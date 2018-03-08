PRO ploticl
device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A4059/test.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;OPENR, lun,'/n/Godiva1/jkrick/A4059/original/icl2rtab', /GET_LUN
OPENR, lun,'/n/Godiva1/jkrick/A4059/original/cdtab', /GET_LUN

;read in the radial profile
rows= 34
radius = FLTARR(rows)
counts = FLTARR(rows)
Verr1 = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius(j) = r
      counts(j) = c
      Verr1(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

counts = counts[0:j-1]    ;get rid of the points which have bad stop codes
radius = radius[0:j-1]
Verr1 = Verr1[0:j-1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;lower limit, 30% larger galaxy masks


out = textoidl('Surface Brightness (mag/arcsec^{-2})')
plot, radius* 0.435,22.04 - 2.5*alog10(counts/(0.435^2)),thick = 3, $
YRANGE = [30., 20.0], xrange = [10,400],ytitle = out, $
xstyle = 9,ystyle = 1,/xlog, $
xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3
;
axis, 0, 20.0, xaxis=1, xrange=[0,0.424], xstyle = 1, xthick = 3,charthick = 3


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
err = dindgen(j) - dindgen(j) + 1
start = [0.03,200.0]

;exponential fitting
;result2 = MPFITFUN('exponential',radius[0:j-1],counts[0:j-1], Verr1[0:j-1], start)    ;ICL

;exponential plotting
;;oplot, findgen(900)*0.435, 22.04 - 2.5*alog10(result2(0)*exp(-findgen(900)/result2(1))/(0.435^2) ), thick = 3, color = colors.red


;sersic fitting and plotting
start = [0.03,200.0, 4.0]
;;result2 = MPFITFUN('sersic',radius[0:j-1],counts[0:j-1], Verr, start)    ;ICL

;;numbers = (result2(0)) * (exp(-7.67*(((radius[0:j-1]/result2(1))^(1.0/result2(2))) - 1.0)))
;;oplot, radius[0:j-1]*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)), thick = 5, color = colors.hotpink


;devauc fitting and plotting
start = [0.05,20.0]
result2 = MPFITFUN('devauc',radius[0:j-16],counts[0:j-16], Verr, start)    ;ICL

numbers = (result2(0)) * (exp(-7.67*(((radius[0:j-1]/result2(1))^(1.0/4.0)) - 1.0)))
oplot, radius[0:j-1]*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)), thick = 5, color = colors.orangered

;double devauc
start= [0.3, 40., 0.01, 200.]
result2 = MPFITFUN('devauc2',radius[0:j-1],counts[0:j-1], Verr, start)    ;ICL

numbers = (result2(0)) * (exp(-7.67*(((radius[0:j-1]/result2(1))^(1.0/4.0)) - 1.0))) + (result2(2)) * (exp(-7.67*(((radius[0:j-1]/(result2(3)))^(1.0/4.0)) - 1.0)))

oplot, radius[0:j-1]*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)), thick = 5, color = colors.orangered


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva1/jkrick/A4059/original/icl2rtab', /GET_LUN
;OPENR, lun,'/n/Godiva1/jkrick/A4059/original/cdtab', /GET_LUN

;read in the radial profile
rows= 20
radius = FLTARR(rows)
counts = FLTARR(rows)
Verr1 = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius(j) = r
      counts(j) = c
      Verr1(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

counts = counts[0:j-1]    ;get rid of the points which have bad stop codes
radius = radius[0:j-1]
Verr1 = Verr1[0:j-1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;lower limit, 30% larger galaxy masks

oplot, radius* 0.435,22.04 - 2.5*alog10(counts/(0.435^2)),thick = 3;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva1/jkrick/A4059/original/icl2r2tab', /GET_LUN

;read in the radial profile
rows= 26
radius2 = FLTARR(rows)
counts2 = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius2(j) = r
      counts2(j) = c
      Verr(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

counts2 = counts2[0:j-1]    ;get rid of the points which have bad stop codes
radius2 = radius2[0:j-1]

;;oplot, radius2* 0.435,22.04 - 2.5*alog10(counts2/(0.435^2)),thick = 3;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;lower limit, 30% larger galaxy masks
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva1/jkrick/A4059/original/icl2r3tab', /GET_LUN

;read in the radial profile
rows= 13
radius3 = FLTARR(rows)
counts3 = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius3(j) = r
      counts3(j) = c
      Verr(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

counts3 = counts3[0:j-1]    ;get rid of the points which have bad stop codes
radius3 = radius3[0:j-1]

;;oplot, radius3* 0.435,22.04 - 2.5*alog10(counts3/(0.435^2)),thick = 3;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pxval = [radius* 0.435,reverse(radius3* 0.435)]
pyval = [22.04 - 2.5*alog10(counts/(0.435^2)),reverse(22.04 - 2.5*alog10(counts3/(0.435^2)))]

;;polyfill, pxval,pyval, color = colors.salmon
;polyfill, pxval,pyval, /line_fill, orientation = 45

pxval = [radius* 0.435,reverse(radius2* 0.435)]
pyval = [22.04 - 2.5*alog10(counts/(0.435^2)),reverse(22.04 - 2.5*alog10(counts2/(0.435^2)))]
;;
;;polyfill, pxval,pyval, color = colors.firebrick
;polyfill, pxval,pyval, color = colors.white, /line_fill, orientation = 45

oplot, radius* 0.435,22.04 - 2.5*alog10(counts/(0.435^2)),thick = 3;
;;;;;;;;
;error bars
;;;;;;;;
lower2 = fltarr(n_elements(counts))
FOR counter = 0, n_elements(counts) - 1, 1 DO BEGIN
    IF (counts(counter) - .0009877 LE 0. )THEN BEGIN
        lower2(counter) = 35.
    ENDIF ELSE BEGIN
        lower2(counter) = 22.04 - 2.5*alog10((counts(counter) - .0009877)/(0.435^2))
    ENDELSE
ENDFOR

upper2 = 22.04 - 2.5*alog10((counts + .0009877)/(0.435^2)) 
errplot, radius* 0.435,lower2, upper2

;devauc fitting and plotting
start = [0.03,200.0]
;;result2 = MPFITFUN('devauc',radius[0:19],counts[0:19], Verr1[0:19], start)    ;ICL

;;numbers = (result2(0)) * (exp(-7.67*(((radius[0:19]/result2(1))^(1.0/4.0)) - 1.0)))
;oplot, radius[0:19]*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)), thick = 5, color = colors.black


;lower2 = 22.04 - 2.5*alog10((counts - .00104)/(0.435^2))
;upper2 = 22.04 - 2.5*alog10((counts + .00104)/(0.435^2)) 
;errplot, radius* 0.435,lower2, upper2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva1/jkrick/A4059/original/gal2rtab', /GET_LUN

;read in the radial profile
rows= 20
radius3 = FLTARR(rows)
counts3 = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius3(j) = r
      counts3(j) = c
      Verr(j) = ellerr 
      stop(j) = s
ENDFOR;

close, lun
free_lun, lun

counts3 = counts3[0:j-1]    ;get rid of the points which have bad stop codes
radius3 = radius3[0:j-1]

;oplot, radius3* 0.435,22.04 - 2.5*alog10(counts3/(0.435^2)),thick = 3, color = colors.red;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva1/jkrick/A4059/original/icl2Btab', /GET_LUN

;read in the radial profile
rows= 20 - 2
radiusB = FLTARR(rows)
countsB = FLTARR(rows)
VerrB = FLTARR(rows)
stopB = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radiusB(j) = r
      countsB(j) = c
      VerrB(j) = ellerr 
      stopB(j) = s
ENDFOR

close, lun
free_lun, lun

countsB = countsB[0:j-1]    ;get rid of the points which have bad stop codes
radiusB = radiusB[0:j-1]

;oplot, radiusB* 0.435,22.19 - 2.5*alog10(countsB/(0.435^2)),thick = 3;

err = dindgen(j) - dindgen(j) + 1
start = [0.03,200.0]

;exponential fitting
;;result2 = MPFITFUN('exponential',radiusB,countsB, VerrB, start)    ;ICL

;exponential plotting
;oplot, findgen(900)*0.435, 22.19 - 2.5*alog10(result2(0)*exp(-findgen(900)/result2(1))/(0.435^2) ), thick = 3, color = colors.blue


;sersic fitting and plotting
start = [0.03,200.0, 4.0]
;result2 = MPFITFUN('sersic',radiusB,countsB, Verr, start)    ;ICL

;numbers = (result2(0)) * (exp(-7.67*(((radiusB/result2(1))^(1.0/result2(2))) - 1.0)))
;oplot, radiusB*0.435, 22.19 - 2.5*alog10(numbers/(0.435^2)), thick = 5, color = colors.cyan


;devauc fitting and plotting
;start = [0.02,60.0]
;result2 = MPFITFUN('devauc',radiusB,countsB, VerrB, start)    ;ICL
;numbers = (result2(0)) * (exp(-7.67*(((radiusB/result2(1))^(1.0/4.0)) - 1.0)))
;oplot, radiusB*0.435, 22.19 - 2.5*alog10(numbers/(0.435^2)), thick = 5, color = colors.red



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva1/jkrick/A4059/original/icl2B2tab', /GET_LUN

;read in the radial profile
rows= 24 - 4
radiusB2 = FLTARR(rows)
countsB2 = FLTARR(rows)
VerrB = FLTARR(rows)
stopB = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radiusB2(j) = r
      countsB2(j) = c
      VerrB(j) = ellerr 
      stopB(j) = s
ENDFOR;

close, lun
free_lun, lun

countsB2 = countsB2[0:j-1]    ;get rid of the points which have bad stop codes
radiusB2 = radiusB2[0:j-1]

;oplot, radiusB2* 0.435,22.19 - 2.5*alog10(countsB2/(0.435^2)),thick = 3;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva1/jkrick/A4059/original/icl2B3tab', /GET_LUN

;read in the radial profile
rows= 14 - 1
radiusB3 = FLTARR(rows)
countsB3 = FLTARR(rows)
VerrB = FLTARR(rows)
stopB = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radiusB3(j) = r
      countsB3(j) = c
      VerrB(j) = ellerr 
      stopB(j) = s
ENDFOR

close, lun
free_lun, lun

countsB3 = countsB3[0:j-1]    ;get rid of the points which have bad stop codes
radiusB3 = radiusB3[0:j-1]

;oplot, radiusB3* 0.435,22.19 - 2.5*alog10(countsB3/(0.435^2)),thick = 3;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pxval = [radiusB* 0.435,reverse(radiusB3* 0.435)]
pyval = [22.19 - 2.5*alog10(countsB/(0.435^2)),reverse(22.19 - 2.5*alog10(countsB3/(0.435^2)))]
polyfill, pxval,pyval, color = colors.powderblue
;polyfill, pxval,pyval, /line_fill, orientation = 45

pxval = [radiusB* 0.435,reverse(radiusB2* 0.435)]
pyval = [22.19 - 2.5*alog10(countsB/(0.435^2)),reverse(22.19 - 2.5*alog10(countsB2/(0.435^2)))]
;;
polyfill, pxval,pyval, color = colors.steelblue
;polyfill, pxval,pyval, color = colors.white, /line_fill, orientation = 45

;oplot, radiusB* 0.435,22.19 - 2.5*alog10(countsB/(0.435^2)),thick = 3;
;;;;;;;;
;error bars
;;;;;;;;

lower2 = fltarr(n_elements(countsB))
FOR counter = 0, n_elements(countsB) - 1, 1 DO BEGIN
    IF (countsB(counter) - .000153 LE 0. )THEN BEGIN
        lower2(counter) = 35.
    ENDIF ELSE BEGIN
        lower2(counter) = 22.19 - 2.5*alog10((countsB(counter) - .000153)/(0.435^2))
    ENDELSE
ENDFOR

upper2 = 22.19 - 2.5*alog10((countsB + .000153)/(0.435^2)) 
errplot, radiusB* 0.435,lower2, upper2

;devauc fitting and plotting
start = [0.03,200.0]
;;result2 = MPFITFUN('devauc',radiusB,countsB, VerrB, start)    ;ICL

;;numbers = (result2(0)) * (exp(-7.67*(((radiusB/result2(1))^(1.0/4.0)) - 1.0)))
;oplot, radiusB*0.435, 22.19 - 2.5*alog10(numbers/(0.435^2)), thick = 5, color = colors.black;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun,'/n/Godiva1/jkrick/A4059/original/gal2Btab', /GET_LUN

;read in the radial profile
rows= 20
radiusB = FLTARR(rows)
countsB = FLTARR(rows)
VerrB = FLTARR(rows)
stopB = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radiusB(j) = r
      countsB(j) = c
      VerrB(j) = ellerr 
      stopB(j) = s
ENDFOR

close, lun
free_lun, lun

countsB = countsB[0:j-1]    ;get rid of the points which have bad stop codes
radiusB = radiusB[0:j-1]

;oplot, radiusB* 0.435,22.19 - 2.5*alog10(countsB/(0.435^2)),thick = 3, color = colors.blue;

;err = dindgen(j) - dindgen(j) + 1
;start = [0.03,200.0]

;exponential fitting
;result2 = MPFITFUN('exponential',radiusB,countsB, VerrB, start)    ;ICL

;exponential plotting
;oplot, findgen(460)*0.435, 22.19 - 2.5*alog10(result2(0)*exp(-findgen(460)/result2(1))/(0.435^2) ), thick = 3, color = colors.black


xyouts, 30,30,"A4059", charthick = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

device, /close
set_plot, mydevice

END






