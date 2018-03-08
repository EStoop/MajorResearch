PRO ploticl
device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva7/jkrick/A114/profiles2.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva7/jkrick/A114/original/icl2rtab', /GET_LUN

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

counts = counts[0:j-1]    ;get rid of the points which have bad stop codes
radius = radius[0:j-1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;lower limit, 30% larger galaxy masks


out = textoidl('Surface Brightness (mag/arcsec^{-2})')
plot, radius* 0.259,24.6 - 2.5*alog10(counts/(0.259^2)),thick = 3, $
YRANGE = [29, 26.0], xrange = [0,150],ytitle = out, $
psym = 3,xstyle = 9,ystyle = 1,$
xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3
;
axis, 0, 26.0, xaxis=1, xrange=[0,0.681], xstyle = 1, xthick = 3,charthick = 3

oplot, radius* 0.259,24.6 - 2.5*alog10(counts/(0.259^2)),thick = 3;
;;;;;;;;
;error bars
;;;;;;;;

lower2 = 24.6 - 2.5*alog10((counts - 0.000555)/(0.259^2))
upper2 = 24.6 - 2.5*alog10((counts + 0.000555)/(0.259^2)) 
errplot, radius* 0.259,lower2, upper2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva7/jkrick/A114/original/icl2r2tab', /GET_LUN

;read in the radial profile
rows= 20
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

;;oplot, radius2* 0.259,24.6 - 2.5*alog10(counts2/(0.259^2)),thick = 3, color = colors.red;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;OPENR, lun,'/n/Godiva7/jkrick/A114/original/iclr3tab', /GET_LUN

;;read in the radial profile
;rows= 17 - 9
;radius3 = FLTARR(rows)
;counts3 = FLTARR(rows)
;Verr = FLTARR(rows)
;stop = FLTARR(rows)
;r = 0.0				;radius
;c = 0.0				;mean counts
;ellerr = 0.0
;s = 0
;FOR j=0,rows-1 DO BEGIN
;      READF, lun, r,c, ellerr,s
;      radius3(j) = r
;      counts3(j) = c
;      Verr(j) = ellerr 
;      stop(j) = s
;ENDFOR;

;close, lun
;free_lun, lun

;counts3 = counts3[0:j-1]    ;get rid of the points which have bad stop codes
;radius3 = radius3[0:j-1]

;;oplot, radius3* 0.259,24.6 - 2.5*alog10(counts3/(0.259^2)),thick = 3, color = colors.salmon;


;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva7/jkrick/A114/original/gal2rtab', /GET_LUN

;read in the radial profile
rows= 21
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

oplot, radius3* 0.259,24.6 - 2.5*alog10(counts3/(0.259^2)),thick = 3, color = colors.red;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;pxval = [radius* 0.259,reverse(radius3* 0.259)]
;;pyval = [24.6 - 2.5*alog10(counts/(0.259^2)),reverse(24.6 - 2.5*alog10(counts3/(0.259^2)))]

;;polyfill, pxval,pyval, color = colors.salmon
;polyfill, pxval,pyval, /line_fill, orientation = 45

pxval = [radius* 0.259,reverse(radius2* 0.259)]
pyval = [24.6 - 2.5*alog10(counts/(0.259^2)),reverse(24.6 - 2.5*alog10(counts2/(0.259^2)))]
;;
polyfill, pxval,pyval, color = colors.firebrick
;polyfill, pxval,pyval, color = colors.white, /line_fill, orientation = 45

oplot, radius* 0.259,24.6 - 2.5*alog10(counts/(0.259^2)),thick = 3;

lower2 = 24.6 - 2.5*alog10((counts - 0.000555)/(0.259^2))
upper2 = 24.6 - 2.5*alog10((counts + 0.000555)/(0.259^2)) 
errplot, radius* 0.259,lower2, upper2


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
err = dindgen(j) - dindgen(j) + 1
start = [0.03,200.0]

;;result3 = MPFITFUN('exponential',radius3[1:15],counts3[1:15], err[0:18], start)    ;galaxies
result2 = MPFITFUN('exponential',radius,counts, Verr, start)    ;ICL

;exponential plotting
;oplot, radius*0.259, 24.6 - 2.5*alog10(result2(0)*exp(-radius/result2(1))/(0.259^2) ), thick = 3, color = colors.red
oplot, findgen(580)*0.259, 24.6 - 2.5*alog10(result2(0)*exp(-findgen(580)/result2(1))/(0.259^2) ), thick = 3, color = colors.red

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
;numbers = (result2(0)) * (exp(-7.67*(((radius2/result2(1))^(1.0/4.0)) - 1.0)))
;oplot, radius2*0.259, 24.6 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.red
;numbers = (result2(0)) * (exp(-7.67*(((radius/result2(1))^(1.0/4.0)) - 1.0)))
;oplot, radius*0.259, 24.6 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.red

;numbers = (result3(0)) * (exp(-7.67*(((radius3[0:19]/result3(1))^(1.0/4.0)) - 1.0)))
;oplot, radius3[0:19]*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.green
;devauc fitting and plotting
start = [0.03,200.0]
result2 = MPFITFUN('devauc',radius,counts, Verr, start)    ;ICL

numbers = (result2(0)) * (exp(-7.67*(((radius/result2(1))^(1.0/4.0)) - 1.0)))
oplot, radius*0.259, 24.6 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.black


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva7/jkrick/A114/original/icl2Vtab', /GET_LUN

;read in the radial profile
rows= 10
Vradius = FLTARR(rows)
Vcounts = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      Vradius(j) = r
      Vcounts(j) = c
      Verr(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

Vcounts = Vcounts[0:j-1]    ;get rid of the points which have bad stop codes
Vradius = Vradius[0:j-1]

;plot, Vradius* 0.259,24.3 - 2.5*alog10(Vcounts/(0.259^2)),thick = 3, $
;YRANGE = [29, 26.0], xrange = [0,150],ytitle = out, $
;psym = 3,xstyle = 8,ystyle = 1,$
;xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3

;axis, 0, 26.0, xaxis=1, xrange=[0,1.158], xstyle = 1, xthick = 3,charthick = 3

oplot, Vradius* 0.259,24.3 - 2.5*alog10(Vcounts/(0.259^2)),thick = 3;

lower2 = 24.3 - 2.5*alog10((Vcounts - 0.00043)/(0.259^2))
upper2 = 24.3 - 2.5*alog10((Vcounts + 0.00043)/(0.259^2)) 
errplot, Vradius* 0.259,lower2, upper2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
err = dindgen(j) - dindgen(j) + 1
start = [0.0003,200.0]

;;result3 = MPFITFUN('exponential',radius3[1:15],counts3[1:15], err[0:18], start)    ;galaxies
result2 = MPFITFUN('exponential',Vradius,Vcounts, Verr, start)    ;ICL

;numbers = (result2(0)) * (exp(-7.67*(((Vradius/result2(1))^(1.0/4.0)) - 1.0)))
;oplot, Vradius*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.green

;oplot, Vradius*0.259, 24.3 - 2.5*alog10(result2(0)*exp(-Vradius/result2(1))/(0.259^2) ), thick = 3, color = colors.green


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun,'/n/Godiva7/jkrick/A114/original/icl2V2tab', /GET_LUN

;read in the radial profile
rows= 20 - 3
Vradius2 = FLTARR(rows)
Vcounts2 = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      Vradius2(j) = r
      Vcounts2(j) = c
      Verr(j) = ellerr 
      stop(j) = s
ENDFOR

close, lun
free_lun, lun

Vcounts2 = Vcounts2[0:j-1]    ;get rid of the points which have bad stop codes
Vradius2 = Vradius2[0:j-1]

;oplot, Vradius2* 0.259,24.3 - 2.5*alog10(Vcounts2/(0.259^2)),thick = 3, color = colors.limegreen;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;OPENR, lun,'/n/Godiva7/jkrick/A114/original/iclV3tab', /GET_LUN

;;read in the radial profile
;rows= 9
;Vradius3 = FLTARR(rows)
;Vcounts3 = FLTARR(rows)
;Verr = FLTARR(rows)
;stop = FLTARR(rows)
;r = 0.0				;radius
;c = 0.0				;mean counts
;ellerr = 0.0
;s = 0
;FOR j=0,rows-1 DO BEGIN
;      READF, lun, r,c, ellerr,s
;      Vradius3(j) = r
;      Vcounts3(j) = c
;      Verr(j) = ellerr 
;      stop(j) = s
;ENDFOR;

;close, lun
;free_lun, lun

;Vcounts3 = Vcounts3[0:j-1]    ;get rid of the points which have bad stop codes
;Vradius3 = Vradius3[0:j-1]

;oplot, Vradius3* 0.259,24.3 - 2.5*alog10(Vcounts3/(0.259^2)),thick = 3, color = colors.green;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;pxval = [Vradius* 0.259,reverse(Vradius3* 0.259)]
;pyval = [24.29 - 2.5*alog10(Vcounts/(0.259^2)),reverse(24.29 - 2.5*alog10(Vcounts3/(0.259^2)))]

;;polyfill, pxval,pyval, color = colors.aquamarine
;;polyfill, pxval,pyval, /line_fill, orientation = 45

pxval = [Vradius* 0.259,reverse(Vradius2* 0.259)]
pyval = [24.29 - 2.5*alog10(Vcounts/(0.259^2)),reverse(24.29 - 2.5*alog10(Vcounts2/(0.259^2)))]

polyfill, pxval,pyval, color = colors.darkgreen
;;polyfill, pxval,pyval, color = colors.white, /line_fill, orientation = 45

oplot, Vradius* 0.259,24.3 - 2.5*alog10(Vcounts/(0.259^2)),thick = 3;

lower2 = 24.3 - 2.5*alog10((Vcounts - 0.00043)/(0.259^2))
upper2 = 24.3 - 2.5*alog10((Vcounts + 0.00043)/(0.259^2)) 
errplot, Vradius* 0.259,lower2, upper2

;pxval = [radius* 0.259,reverse(radius3* 0.259)]
;pyval = [24.6 - 2.5*alog10(counts/(0.259^2)),reverse(24.6 - 2.5*alog10(counts3/(0.259^2)))]

;;polyfill, pxval,pyval, color = colors.salmon
;polyfill, pxval,pyval, /line_fill, orientation = 45

pxval = [radius* 0.259,reverse(radius2* 0.259)]
pyval = [24.6 - 2.5*alog10(counts/(0.259^2)),reverse(24.6 - 2.5*alog10(counts2/(0.259^2)))]

;;polyfill, pxval,pyval, color = colors.firebrick
;polyfill, pxval,pyval, color = colors.white, /line_fill, orientation = 45

oplot, radius* 0.259,24.6 - 2.5*alog10(counts/(0.259^2)),thick = 3;

lower2 = 24.6 - 2.5*alog10((counts - 0.000555)/(0.259^2))
upper2 = 24.6 - 2.5*alog10((counts + 0.000555)/(0.259^2)) 
errplot, radius* 0.259,lower2, upper2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva7/jkrick/A114/original/gal2Vtab', /GET_LUN

;read in the radial profile
rows= 15
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

oplot, radius3* 0.259,24.3 - 2.5*alog10(counts3/(0.259^2)),thick = 3, color = colors.green;
oplot, findgen(420)*0.259, 24.3 - 2.5*alog10(result2(0)*exp(-findgen(420)/result2(1))/(0.259^2) ), thick = 3, color = colors.green

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;devauc fitting and plotting
start = [0.03,200.0]
result2 = MPFITFUN('devauc',Vradius,Vcounts, Verr, start)    ;ICL

numbers = (result2(0)) * (exp(-7.67*(((Vradius/result2(1))^(1.0/4.0)) - 1.0)))
oplot, Vradius*0.259, 24.3 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.black

xyouts, 10,26.3,"AC114", charthick = 3

device, /close
set_plot, mydevice

END






