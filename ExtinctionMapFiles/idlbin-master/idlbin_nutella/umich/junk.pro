PRO ploticl
device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/Users/jkrick/umich/icl/A2734/profiles.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/Users/jkrick/umich/icl/A2734/icl2rtab', /GET_LUN

;read in the radial profile
rows= 17 -2
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

counts = counts[0:j-1]    ;get rid of the points which have bad stop codes
radius = radius[0:j-1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;lower limit, 30% larger galaxy masks


out = textoidl('Surface Brightness (mag/arcsec^{-2})')
plot, radius* 0.435,22.04 - 2.5*alog10(counts/(0.435^2)),thick = 3, $
YRANGE = [29, 20.0], xrange = [0,250],ytitle = out, $
psym = 3,xstyle = 9,ystyle = 1,$
xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3
;
axis, 0, 20.0, xaxis=1, xrange=[0,0.361], xstyle = 1, xthick = 3,charthick = 3


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
err = dindgen(j) - dindgen(j) + 1
start = [0.03,200.0]

;exponential fitting
result2 = MPFITFUN('exponential',radius[0:j-1],counts[0:j-1], Verr, start)    ;ICL

;exponential plotting
;oplot, findgen(900)*0.435, 22.04 - 2.5*alog10(result2(0)*exp(-findgen(900)/result2(1))/(0.435^2) ), thick = 3, color = colors.red


;sersic fitting and plotting
start = [0.03,200.0, 4.0]
;result2 = MPFITFUN('sersic',radius[0:j-1],counts[0:j-1], Verr, start)    ;ICL

;numbers = (result2(0)) * (exp(-7.67*(((radius[0:j-1]/result2(1))^(1.0/result2(2))) - 1.0)))
;oplot, radius[0:j-1]*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)), thick = 5, color = colors.hotpink


;devauc fitting and plotting
;start = [0.03,200.0]
;result2 = MPFITFUN('devauc',radius[0:j-1],counts[0:j-1], Verr, start)    ;ICL

;numbers = (result2(0)) * (exp(-7.67*(((radius[0:j-1]/result2(1))^(1.0/4.0)) - 1.0)))
;oplot, radius[0:j-1]*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)), thick = 5, color = colors.orangered




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   add the cD galaxy
OPENR, lun,'/Users/jkrick/umich/icl/A2734/cdtab', /GET_LUN

;read in the radial profile
rows= 40
radiusc = FLTARR(rows)
countsc = FLTARR(rows)
Verrc = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radiusc(j) = r
      countsc(j) = c
      Verrc(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

countsc = countsc[0:j-1]    ;get rid of the points which have bad stop codes
radiusc = radiusc[0:j-1]
Verrc = Verrc[0:j-1]
oplot, radiusc* 0.435,22.04 - 2.5*alog10(countsc/(0.435^2)),thick = 3;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun,'/Users/jkrick/umich/icl/A2734/icl2r2tab', /GET_LUN

;read in the radial profile
rows= 19 -2
radius2 = FLTARR(rows)
counts2 = FLTARR(rows)
Verr2 = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      radius2(j) = r
      counts2(j) = c
      Verr2(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

counts2 = counts2[0:j-1]    ;get rid of the points which have bad stop codes
radius2 = radius2[0:j-1]

;oplot, radius2* 0.435,22.04 - 2.5*alog10(counts2/(0.435^2)),thick = 3;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;lower limit, 30% larger galaxy masks
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/Users/jkrick/umich/icl/A2734/icl2r3tab', /GET_LUN

;read in the radial profile
rows= 15 -2
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

counts3 = counts3[0:j-1]    ;get rid of the points which have bad stop codes
radius3 = radius3[0:j-1]

;oplot, radius3* 0.435,22.04 - 2.5*alog10(counts3/(0.435^2)),thick = 3;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;pxval = [radius* 0.435,reverse(radius3* 0.435)]
;pyval = [22.04 - 2.5*alog10(counts/(0.435^2)),reverse(22.04 - 2.5*alog10(counts3/(0.435^2)))]

;polyfill, pxval,pyval, color = colors.salmon

;pxval = [radius* 0.435,reverse(radius2* 0.435)]
;pyval = [22.04 - 2.5*alog10(counts/(0.435^2)),reverse(22.04 - 2.5*alog10(counts2/(0.435^2)))]
;;
;polyfill, pxval,pyval, color = colors.firebrick

;oplot, radius* 0.435,22.04 - 2.5*alog10(counts/(0.435^2)),thick = 3;
;;;;;;;;
;error bars
;;;;;;;;
lower2 = fltarr(n_elements(countsc))
FOR counter = 0, n_elements(countsc) - 1, 1 DO BEGIN
    IF (countsc(counter) - .000399 LE 0. )THEN BEGIN
        lower2(counter) = 35.
    ENDIF ELSE BEGIN
        lower2(counter) = 22.04 - 2.5*alog10((countsc(counter) - .000399)/(0.435^2))
    ENDELSE
ENDFOR

upper2 = 22.04 - 2.5*alog10((countsc + .000399)/(0.435^2)) 
errplot, radiusc* 0.435,lower2, upper2

;;;;;;;;
;error bars
;;;;;;;;
lower2 = fltarr(n_elements(counts))
FOR counter = 0, n_elements(counts) - 1, 1 DO BEGIN
    IF (counts(counter) - .000399 LE 0. )THEN BEGIN
        lower2(counter) = 35.
    ENDIF ELSE BEGIN
        lower2(counter) = 22.04 - 2.5*alog10((counts(counter) - .000399)/(0.435^2))
    ENDELSE
ENDFOR

upper2 = 22.04 - 2.5*alog10((counts + .000399)/(0.435^2)) 
errplot, radius* 0.435,lower2, upper2



;devauc fitting and plotting
print, "single devauc on ICL"
start = [0.03,200.0]
result2 = MPFITFUN('devauc',radius[0:j-1],counts[0:j-1], Verr, start)    ;ICL

;numbers = (result2(0)) * (exp(-7.67*(((radius[0:j-1]/result2(1))^(1.0/4.0)) - 1.0)))
;oplot, radius[0:j-1]*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)), thick = 5, color = colors.black

numbers = (result2(0)) * (exp(-7.67*(((findgen(700)/result2(1))^(1.0/4.0)) - 1.0)))
oplot, findgen(700)*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)) , thick = 5, color = colors.black

print, "single devauc on cD"
result27 =  MPFITFUN('devauc',radiusc,countsc, Verrc, start)    ;ICL
numbers = (result27(0)) * (exp(-7.67*(((findgen(700)/result27(1))^(1.0/4.0)) - 1.0)))
oplot, findgen(700)*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)) , thick = 5, color = colors.green

print, "double devauc on cD"
start=[0.05,83.0,0.43, 50.]
result37 =  MPFITFUN('devauc2',radiusc,countsc, Verrc, start)    ;ICL
numbers = (result37(0)) * (exp(-7.67*(((findgen(700)/result37(1))^(1.0/4.0)) - 1.0))) + (result37(2)) * (exp(-7.67*(((findgen(700)/result37(3))^(1.0/4.0)) - 1.0)))
;oplot, findgen(700)*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)) , thick = 5, color = colors.orange

print, "two fits on cD"
start=[0.05,83.0]
result27 =  MPFITFUN('devauc',radiusc[0:19],countsc[0:19], Verrc[0:19], start)    ;ICL
numbers = (result27(0)) * (exp(-7.67*(((findgen(700)/result27(1))^(1.0/4.0)) - 1.0)))
oplot, findgen(700)*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)) , thick = 5, color = colors.green

start=[0.05,83.0]
result27 =  MPFITFUN('devauc',radiusc[20:39],countsc[20:39], Verrc[20:39], start)    ;ICL
numbers = (result27(0)) * (exp(-7.67*(((findgen(700)/result27(1))^(1.0/4.0)) - 1.0)))
oplot, findgen(700)*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)) , thick = 5, color = colors.green


;numbers = (0.43) * (exp(-7.67*(((findgen(700)/50.94)^(1.0/4.0)) - 1.0))) + (0.05) * (exp(-7.67*(((findgen(700)/83)^(1.0/4.0)) - 1.0)))
;oplot, findgen(700)*0.435, 22.04 - 2.5*alog10(numbers/(0.435^2)) , thick = 5, color = colors.purple



;lower2 = 22.04 - 2.5*alog10((counts - .00104)/(0.435^2))
;upper2 = 22.04 - 2.5*alog10((counts + .00104)/(0.435^2)) 
;errplot, radius* 0.435,lower2, upper2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
xyouts, 30,30,"A2734", charthick=3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

device, /close
set_plot, mydevice

END






