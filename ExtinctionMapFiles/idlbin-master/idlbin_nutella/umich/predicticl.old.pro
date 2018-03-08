PRO predicticl
;device, true=24% ;device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva2/jkrick/A141/expect.3888.ps', /portrait, $
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
rows=21-1
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
plot, radius2[2:i-2]* 0.259,24.84 - 2.5*alog10(counts2[2:i-2]/(0.259^2)),thick = 3, $
YRANGE = [31, 23], xrange = [0,300],ytitle = out, $
psym = 3,xstyle = 4,ystyle = 1,$
xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3

axis, 0, 23.0, xaxis=1, xrange=[0,1.034], xstyle = 1, xthick =3,charthick = 3  ;A3888
axis, 0, 31.0, xaxis=0, xrange=[0,184.5], xstyle = 1, xthick =3,charthick = 3 ,xtitle = 'Semi-major Axis (arcseconds)' ;A3888
;axis, 0, 23.0, xaxis=1, xrange=[0,1.681], xstyle = 1, xthick =3,charthick = 3  ;A141

;;;;;;
pxval = [radius2[0:i - 2]* 0.259,reverse(radius4[0:l-2]* 0.259)]
pyval = [24.84 - 2.5*alog10(counts2[0:i - 2]/(0.259^2)),reverse(24.84 - 2.5*alog10(counts4[0:l-2]/(0.259^2)))]
;polyfill, pxval,pyval, color = colors.lightgray
polyfill, pxval,pyval, color = colors.aquamarine
polyfill, pxval,pyval, /line_fill, orientation = 45

pxval = [radius2[0:i - 2]* 0.259,reverse(radius* 0.259)]
pyval = [24.84 - 2.5*alog10(counts2[0:i - 2]/(0.259^2)),reverse(24.84 - 2.5*alog10(counts/(0.259^2)))]
;polyfill, pxval,pyval, color = colors.darkgray
polyfill, pxval,pyval, color = colors.white, /line_fill, orientation = 45
polyfill, pxval,pyval, color = colors.seagreen


;;;;;;;
oplot, radius2[0:i-2]* 0.259,24.84 - 2.5*alog10(counts2[0:i-2]/(0.259^2)),thick = 3;,linestyle = 2
oplot, radius3* 0.259,24.84 - 2.5*alog10(counts3/(0.259^2)),thick = 5, linestyle = 2, color = colors.seagreen

;oplot, radius4* 0.259,24.84 - 2.5*alog10(counts4/(0.259^2)),thick = 3

;;;;;;;;

;lower2 = 24.84 - 2.5*alog10((counts2[0:i - 2] - Verr2[0:i - 2])/(0.259^2))
;upper2 = 24.84 - 2.5*alog10((counts2[0:i - 2] + Verr2[0:i - 2])/(0.259^2)) 
lower2 = 24.84 - 2.5*alog10((counts2[0:i - 2] - 0.00051)/(0.259^2))
upper2 = 24.84 - 2.5*alog10((counts2[0:i - 2] + 0.00051)/(0.259^2)) 
errplot, radius2* 0.259,lower2, upper2
print, lower2
print, upper2
;lower3 = 24.84 - 2.5*alog10((counts3- Verr3)/(0.259^2))
;upper3 = 24.84 - 2.5*alog10((counts3 + Verr3)/(0.259^2)) 
;errplot, radius3* 0.259,lower3, upper3

;;;;;;;;
;Fitting OF the data 

err = dindgen(j) - dindgen(j) + 1
start = [0.03,200.0]

;;result3 = MPFITFUN('exponential',radius3[1:15],counts3[1:15], err[0:18], start)    ;galaxies
;;result2 = MPFITFUN('exponential',radius2[1:15],counts2[1:15], Verr2[0:13], start)    ;ICL

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


xarrjunk = fltarr(10)+228
yarrjunk = findgen(10) + 28.9;29.8
oplot, xarrjunk, yarrjunk
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;r  curve
OPENR, lun,'/n/Godiva1/jkrick/A3888/final/iclrtab', /GET_LUN

;read in the radial profile
rows= 19-2
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


rcounts = rcounts - 2;[0:j-5] - 2
rradius = rradius;[0:j-5]
;print, rcounts, rradius,rerr
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;r  curve
OPENR, lun,'/n/Godiva1/jkrick/A3888/final/iclrtab2', /GET_LUN 

;read in the radial profile
rows= 18-2
rradius2 = FLTARR(rows)
rcounts2 = FLTARR(rows)
rerr2 = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      rradius2(j) = r
      rcounts2(j) = c
      rerr2(j) = ellerr
ENDFOR

close, lun
free_lun, lun

rcounts2 = rcounts2 - 2;[0:j-3] - 2
rradous2 = rradius2;[0:j-3]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;r  curve
OPENR, lun,'/n/Godiva1/jkrick/A3888/final/galrtab2', /GET_LUN

;read in the radial profile
rows= 18-2
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

rcounts3 = rcounts3 -2;[0:j-3] - 2
rradius3 = rradius3;[0:j-3]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun,'/n/Godiva1/jkrick/A3888/final/iclrtab3', /GET_LUN  ;4

;read in the radial profile
rows= 19-2
rradius4 = FLTARR(rows)
rcounts4 = FLTARR(rows)
rerr4 = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      rradius4(j) = r
      rcounts4(j) = c
      rerr4(j) = ellerr
ENDFOR

close, lun
free_lun, lun

rcounts4 = rcounts4 -2;[0:j-6] - 2
rradius4 = rradius4;[0:j-6]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;device, filename = '/n/Godiva1/jkrick/A3888/iclrexpfit.ps', /portrait, $
;  BITS=8, scale_factor=0.9 , /color
out = textoidl('Surface Brightness (mag/arcsec^{-2})')

;;plot, rradius2* 0.259,25.03 - 2.5*alog10(rcounts2/(0.259^2)),thick = 3, $
;YRANGE = [31, 23], xrange = [0,300],ytitle = out, $
;psym = 3,xstyle = 8,ystyle = 1,$
;xtitle = 'Semi-major Axis (arcseconds, Mpc)', charthick = 3, xthick = 3, ythick = 3

;axis, 0, 23.0, xaxis=1, xrange=[0,1.034], xstyle = 1, xthick = 3,charthick = 3


;print, rradius2, rcounts2

oplot, rradius2* 0.259,25.03 - 2.5*alog10(rcounts2/(0.259^2)),thick = 3;, $

;YRANGE = [32, 24], xrange = [0,320],ytitle = 'Surface Brightness (mag/squarcsec)', $
;title = 'r-band ICL Radial Profile ',psym = 3,xstyle = 1,$
;xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3


;;;;;;
pxval = [rradius2* 0.259,reverse(rradius4* 0.259)]
pyval = [25.03 - 2.5*alog10(rcounts2/(0.259^2)),reverse(25.03 - 2.5*alog10(rcounts4/(0.259^2)))]
;polyfill, pxval,pyval, color = colors.gray
polyfill, pxval,pyval, color = colors.lightsalmon

pxval = [rradius2* 0.259,reverse(rradius* 0.259)]
pyval = [25.03 - 2.5*alog10(rcounts2/(0.259^2)),reverse(25.03 - 2.5*alog10(rcounts/(0.259^2)))]
;polyfill, pxval,pyval, color = colors.darkgray
polyfill, pxval,pyval, color = colors.orangered



;;;;;;;
;oplot, rradius* 0.259,25.03 - 2.5*alog10(rcounts/(0.259^2)),thick = 3
oplot, rradius2* 0.259,25.03 - 2.5*alog10(rcounts2/(0.259^2)),thick = 3
oplot, rradius3* 0.259,25.03 - 2.5*alog10(rcounts3/(0.259^2)),thick = 5,color = colors.red

;oplot, rradius4* 0.259,25.03 - 2.5*alog10(rcounts4/(0.259^2)),thick = 3

;out = textoidl('h_{70}^{-1}kpc')


;;;;;;;;

;lower2 = 25.03 - 2.5*alog10((rcounts2 - rerr2)/(0.259^2))
;upper2 = 25.03 - 2.5*alog10((rcounts2 + rerr2)/(0.259^2)) 
lower2 = 25.03 - 2.5*alog10((rcounts2 - 0.00158)/(0.259^2))
upper2 = 25.03 - 2.5*alog10((rcounts2 + 0.00158)/(0.259^2)) 
errplot, rradius2* 0.259,lower2, upper2

lower2 = 24.84 - 2.5*alog10((counts2[0:i - 2] - 0.00051)/(0.259^2))
upper2 = 24.84 - 2.5*alog10((counts2[0:i - 2] + 0.00051)/(0.259^2)) 
errplot, radius2* 0.259,lower2, upper2

;lower3 = 25.03 - 2.5*alog10((counts3- Verr3)/(0.259^2))
;upper3 = 25.03 - 2.5*alog10((counts3 + Verr3)/(0.259^2)) 
;errplot, radius3* 0.259,lower3, upper3

;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Fitting OF the data with exponentials
newcounts3 = interpolate(rcounts3, [0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
newradius3 = interpolate(rradius3, [0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])

err = dindgen(20) - dindgen(20) + 1
start = [0.02,300]
;start=[26.0,0.04]
;;result3 = MPFITFUN('exponential',rradius3[1:12],rcounts3[1:12],err[1:12], start,maxiter=300)        ;galaxies
;;result2 = MPFITFUN('exponential',rradius2[1:12],rcounts2[1:12], rerr2[1:12], start)      ;ICL

;exponential plotting
;;oplot, rradius2[1:12]*0.259, 25.03 - 2.5*alog10(result2(0)*exp(-rradius2[1:12]/result2(1))/(0.259^2) ), thick = 3;, color = colors.red
;;oplot, newradius3[1:12]*0.259, 25.03 - 2.5*alog10(result3(0)*exp(-newradius3[1:12]/result3(1))/(0.259^2) ), thick = 3;, color = colors.red


;reynolds plotting
;numbers = result2(0) / ((1 + (rradius3/result2(1)))^2)
;oplot, rradius2[0:11]*0.259, 25.03 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.red

;numbers =result3(0) / ((1 + (newradius3/result3(1)))^2)
;oplot, newradius3[0:19]*0.259, 25.03 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.red


;sersic plotting
;numbers = (result2(0)) * (exp(-7.67*(((rradius2[0:11]/result2(1))^(1.0/result2(2))) - 1.0)))
;oplot, rradius2[0:11]*0.259, 25.03 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.red

;numbers = (result3(0)) * (exp(-7.67*(((newradius3[0:19]/result3(1))^(1.0/result3(2))) - 1.0)))
;oplot, newradius3[0:19]*0.259, 25.03 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.red



;devauc plotting
;numbers = (result2(0)) * (exp(-7.67*(((rradius2[0:11]/result2(1))^(1.0/4.0)) - 1.0)))
;oplot, rradius2[0:11]*0.259, 25.03 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.red

;numbers = (result3(0)) * (exp(-7.67*(((newradius3[0:19]/result3(1))^(1.0/4.0)) - 1.0)))
;oplot, newradius3[0:19]*0.259, 25.03 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.red



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
;oplot, findgen(400), findgen(400)  - findgen(400) + 30.6, color = colors.firebrick, thick = 3
;oplot, findgen(400) + 230, findgen(400)  - findgen(400) + 28.8, thick = 3;, color = colors.orangered
;oplot, findgen(150) + 50, findgen(400)  - findgen(400) + 27.6, thick = 3;, color = colors.orangered

;oplot, findgen(400)+160, findgen(400)  - findgen(400) + 28.0, thick = 3;, color = colors.orangered



;xyouts, 30, 30.9, ' 1 sigma Accuracy Limit', color = colors.firebrick, charthick = 3
;xyouts, 30, 29.7, ' 3 sigma Accuracy Limit', color = colors.firebrick, charthick = 3

;xyouts, 145, 28,'Exponential Fit to the Data', color = colors.red, charthick = 3
device, /close
set_plot, mydevice

END
