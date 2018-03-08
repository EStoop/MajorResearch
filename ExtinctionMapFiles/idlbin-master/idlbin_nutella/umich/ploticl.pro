PRO ploticl
device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/iclnewr.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva1/jkrick/A3888/ellipseold/iclVtab', /GET_LUN

;read in the radial profile
rows= 22
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

OPENR, lun,'/n/Godiva1/jkrick/A3888/ellipseold//iclVtab2', /GET_LUN  ;iclVtab2

;read in the radial profile
rows=33
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
OPENR, lun,'/n/Godiva1/jkrick/A3888/ellipseold/iclVtab3', /GET_LUN

;read in the radial profile
rows= 24
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

OPENR, lun,'/n/Godiva1/jkrick/A3888/ellipseold//galVbkgtab', /GET_LUN  ;galVtab2

;read in the radial profile
rows= 33
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
;Fitting OF the data with exponentials

err = dindgen(j) - dindgen(j) + 1
start = [0.01,200.0]
;result = MPFITFUN('exponential',radius,counts, err, start)
;result2 = MPFITFUN('exponential',radius2[10:i - 2],counts2[10:i - 2], err, start)
radiusnew= radius3[10:*]
countsnew = counts3[10:*]
start = [0.1, 100];
;result3 = MPFITFUN('exponential',radiusnew,countsnew, err, start)
;result4 = MPFITFUN('exponential',radius4[0:l - 2],counts4[0:l - 2], err, start)

;print, "exponential result", result2, result3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


plot, radius2* 0.259,24.29 - 2.5*alog10(counts2/(0.259^2)),thick = 3, $
YRANGE = [32, 24], xrange = [0,320],ytitle = 'Surface Brightness (mag/squarcsec)', $
title = 'V-band ICL Radial Profile ',psym = 3,xstyle = 1,$
xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3


;;;;;;
pxval = [radius2[0:i - 2]* 0.259,reverse(radius4[0:l-2]* 0.259)]
pyval = [24.29 - 2.5*alog10(counts2[0:i - 2]/(0.259^2)),reverse(24.29 - 2.5*alog10(counts4[0:l-2]/(0.259^2)))]
polyfill, pxval,pyval, color = colors.aquamarine

pxval = [radius2[0:i - 2]* 0.259,reverse(radius* 0.259)]
pyval = [24.29 - 2.5*alog10(counts2[0:i - 2]/(0.259^2)),reverse(24.29 - 2.5*alog10(counts/(0.259^2)))]
polyfill, pxval,pyval, color = colors.seagreen


;;;;;;;
oplot, radius2* 0.259,24.29 - 2.5*alog10(counts2/(0.259^2)),thick = 3
oplot, radius3* 0.259,24.29 - 2.5*alog10(counts3/(0.259^2)),thick = 3, linestyle = 1
;oplot, radius4* 0.259,24.29 - 2.5*alog10(counts4/(0.259^2)),thick = 3

;;;;;;;;

lower2 = 24.29 - 2.5*alog10((counts2[0:i - 2] - Verr2[0:i - 2])/(0.259^2))
upper2 = 24.29 - 2.5*alog10((counts2[0:i - 2] + Verr2[0:i - 2])/(0.259^2)) 
errplot, radius2* 0.259,lower2, upper2


;lower3 = 24.29 - 2.5*alog10((counts3- Verr3)/(0.259^2))
upper3 = 24.29 - 2.5*alog10((counts3 + Verr3)/(0.259^2)) 
;errplot, radius3* 0.259,lower3, upper3

;;;;;;;;

;oplot, radius*0.259, 24.3 - 2.5*alog10(result(0)*exp(-radius/result(1))/(0.259^2) ), thick = 3
;;;oplot, radius2[10:i-2]*0.259, 24.29 - 2.5*alog10(result2(0)*exp(-radius2[10:i-2]/result2(1))/(0.259^2) ), thick = 3
;;;oplot, radiusnew*0.259, 24.29 - 2.5*alog10(result3(0)*exp(-radiusnew/result3(1))/(0.259^2) ), thick = 3
;oplot, radius4*0.259, 24.29 - 2.5*alog10(result4(0)*exp(-radius4/result4(1))/(0.259^2) ), thick = 3

oplot, findgen(400), findgen(400)  - findgen(400) + 29.7, color = colors.green, thick = 3, linestyle = 2
oplot, findgen(400), findgen(400)  - findgen(400) + 30.9, color = colors.green, thick = 3

xyouts, 30, 31.2, '1 sigma Accuracy Limit', color = colors.green, charthick = 3
xyouts, 30, 30.0, '3 sigma Accuracy Limit', color = colors.green, charthick = 3

xyouts, 43.48, 24.4, '100', charthick = 3, alignment = 0.5
xyouts, 87, 24.4, '200', charthick = 3, alignment = 0.5
xyouts, 130.4, 24.4, '300', charthick = 3, alignment = 0.5
xyouts, 173.9, 24.4, '400', charthick = 3, alignment = 0.5
xyouts, 217.4, 24.4, '500', charthick = 3, alignment = 0.5
xyouts, 260.9, 24.4, '600', charthick = 3, alignment = 0.5
xyouts, 278, 24.4, 'kpc', charthick = 3, alignment = 0.5

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;r  curve
OPENR, lun,'/n/Godiva1/jkrick/A3888/ellipse/iclrtab', /GET_LUN

;read in the radial profile
rows= 31-5
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;r  curve
OPENR, lun,'/n/Godiva1/jkrick/A3888/ellipse/iclrtab2', /GET_LUN 

;read in the radial profile
rows= 22-3
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
OPENR, lun,'/n/Godiva1/jkrick/A3888/ellipseold/galrbkgtab', /GET_LUN

;read in the radial profile
rows= 22-3
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
rradous3 = rradius3;[0:j-3]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun,'/n/Godiva1/jkrick/A3888/ellipse/iclrtab3', /GET_LUN;4

;read in the radial profile
rows= 27-6
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
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Fitting OF the data with exponentials

err = dindgen(j) - dindgen(j) + 1
start = [0.01,200.0]
;result = MPFITFUN('exponential',rradius,rcounts,rerr, start)
;result2 = MPFITFUN('exponential',rradius2,rcounts2, rerr2, start)
;radiusnew= rradius3[6:19]
;countsnew = rcounts3[6:19]
;print, radiusnew, countsnew
start = [0.01, 1000.]
;result3 = MPFITFUN('exponential',radiusnew,countsnew, rerr3, start)
;result4 = MPFITFUN('exponential',rradius4,rcounts4, rerr4, start)

;print, "exponential result", result,result2, result3, result4


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
plot, rradius2* 0.259,24.6 - 2.5*alog10(rcounts2/(0.259^2)),thick = 3, $
YRANGE = [32, 24], xrange = [0,320],ytitle = 'Surface Brightness (mag/squarcsec)', $
title = 'r-band ICL Radial Profile ',psym = 3,xstyle = 1,$
xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3


;;;;;;
pxval = [rradius2* 0.259,reverse(rradius4* 0.259)]
pyval = [24.6 - 2.5*alog10(rcounts2/(0.259^2)),reverse(24.6 - 2.5*alog10(rcounts4/(0.259^2)))]
polyfill, pxval,pyval, color = colors.bisque

pxval = [rradius2* 0.259,reverse(rradius* 0.259)]
pyval = [24.6 - 2.5*alog10(rcounts2/(0.259^2)),reverse(24.6 - 2.5*alog10(rcounts/(0.259^2)))]
polyfill, pxval,pyval, color = colors.salmon



;;;;;;;
;oplot, rradius* 0.259,24.6 - 2.5*alog10(rcounts/(0.259^2)),thick = 3
oplot, rradius2* 0.259,24.6 - 2.5*alog10(rcounts2/(0.259^2)),thick = 3
oplot, rradius3* 0.259,24.6 - 2.5*alog10(rcounts3/(0.259^2)),thick = 3, linestyle = 1
;oplot, rradius4* 0.259,24.6 - 2.5*alog10(rcounts4/(0.259^2)),thick = 3

;;;;;;;;

lower2 = 24.6 - 2.5*alog10((rcounts2 - rerr2)/(0.259^2))
upper2 = 24.6 - 2.5*alog10((rcounts2 + rerr2)/(0.259^2)) 
errplot, rradius2* 0.259,lower2, upper2

;lower3 = 24.6 - 2.5*alog10((counts3- Verr3)/(0.259^2))
;upper3 = 24.6 - 2.5*alog10((counts3 + Verr3)/(0.259^2)) 
;errplot, radius3* 0.259,lower3, upper3

;;;;;;;;

;oplot, radius*0.259, 24.3 - 2.5*alog10(result(0)*exp(-radius/result(1))/(0.259^2) ), thick = 3
;;;oplot, rradius2*0.259, 24.6 - 2.5*alog10(result2(0)*exp(-rradius2/result2(1))/(0.259^2) ), thick = 3
;;;;oplot, radiusnew*0.259, 24.6 - 2.5*alog10(result3(0)*exp(-radiusnew/result3(1))/(0.259^2) ), thick = 3
;oplot, radius4*0.259, 24.6 - 2.5*alog10(result4(0)*exp(-radius4/result4(1))/(0.259^2) ), thick = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
oplot, findgen(400), findgen(400)  - findgen(400) + 30.6, color = colors.firebrick, thick = 3
oplot, findgen(400), findgen(400)  - findgen(400) + 29.4, color = colors.firebrick, thick = 3, linestyle = 2
xyouts, 43.48, 24.4, '100', charthick = 3, alignment = 0.5
xyouts, 87, 24.4, '200', charthick = 3, alignment = 0.5
xyouts, 130.4, 24.4, '300', charthick = 3, alignment = 0.5
xyouts, 173.9, 24.4, '400', charthick = 3, alignment = 0.5
xyouts, 217.4, 24.4, '500', charthick = 3, alignment = 0.5
xyouts, 260.9, 24.4, '600', charthick = 3, alignment = 0.5
xyouts, 278, 24.4, 'kpc', charthick = 3, alignment = 0.5

xyouts, 30, 30.9, ' 1 sigma Accuracy Limit', color = colors.firebrick, charthick = 3
xyouts, 30, 29.7, ' 3 sigma Accuracy Limit', color = colors.firebrick, charthick = 3

;xyouts, 145, 28,'Exponential Fit to the Data', color = colors.red, charthick = 3
device, /close
set_plot, mydevice

END

;radius1 = radius(0:87)
;radius2 = radius(87:97)
;counts1 = counts(0:87)
;counts2 = counts(87:97)
;plot, radius* 0.259,24.3 - 2.5*alog10(counts/(0.259^2)),thick = 3, $
;YRANGE = [32, 24], xrange = [0,300],ytitle = 'Surface Brightness (mag/squarcsec)', title = 'V-band ICL Radial Profile ', $
; xtitle = 'Radius (arcseconds)', charthick = 3, xthick = 3, ythick = 3


;result = MPFITFUN('exponential',rad(4:l-1),averagearr(4:l-1), err, start);,PARINFO=pi);, weights=weight)


; somehow demarcate the points that Ellipse doesn't even believe
;FOR j = 0, rows - 1, 1 DO BEGIN
;    IF stop(j) GT 0 THEN BEGIN
;       xyouts, 0.259*radius(j), 24.29 - 2.5*alog10(counts(j)/(0.259^2)), ".",$
;color = colors.red, charthick = 5,alignment = 0.5
;    ENDIF
;ENDFOR


;plot, rradius* 0.259,24.3 - 2.5*alog10(rcounts/(0.259^2)),thick = 3, $
;YRANGE = [32, 24], xrange = [0,300],ytitle = 'Surface Brightness (mag/squarcsec)', $
;title = 'r-band ICL Radial Profile ', xtitle = 'Radius (arcseconds)', charthick = 3, xthick = 3, ythick = 3

;oplot, radius* 0.259,24.62 - 2.5*alog10(counts/(0.259^2)),thick = 3, linestyle = 2, color = colors.red



;oplot, radius2* 0.259,24.3 - 2.5*alog10(counts2/(0.259^2)),thick = 5, color = colors.blue

;oplot, radius*0.259, 24.3 - 2.5*alog10(result(0)*exp(-radius/result(1))/(0.259^2) ), color= colors.red, thick = 3

;oplot, radius*0.259, 24.62 - 2.5*alog10(result(0)*exp(-radius/result(1))/(0.259^2) ), $
;color= colors.red, thick = 3, linestyle = 2




;newcounts2= congrid(rcounts2, 34)
;newradius2= congrid(rradius2, 34)
;plotcolorfill, rradius* 0.259, 24.6 - 2.5*alog10(rcounts/(0.259^2)), color = colors.bisque,$
;bottom = 24.6 - 2.5*alog10(newcounts2/(0.259^2)), /noerase

;newcounts4 = congrid(rcounts4, 36)
;newradius4 = congrid(rradius4, 36)
;plotcolorfill, rradius2* 0.259, 24.6 - 2.5*alog10(rcounts2/(0.259^2)), color = colors.salmon,$
;bottom = 24.6 - 2.5*alog10(newcounts4/(0.259^2)), /noerase
;print, newcounts4 - rcounts4

;newcounts = congrid(rcounts2, 32)
;newradius = congrid(rradius2, 32)
;plotcolorfill, rradius4* 0.259, 24.6 - 2.5*alog10(rcounts4/(0.259^2)), color = colors.salmon,$
;bottom = 24.6 - 2.5*alog10(newcounts/(0.259^2)), /noerase
;print, newcounts4 - rcounts4



;plotcolorfill, radius* 0.259, 24.29 - 2.5*alog10(counts/(0.259^2)), color = colors.seagreen,$
;bottom = 24.29 - 2.5*alog10(counts2/(0.259^2)), /noerase, /notrace

;rebin counts4
;newcounts4 = congrid(counts4, 48)
;newradius4 = congrid(radius4, 48)
;plotcolorfill, radius2* 0.259, 24.29 - 2.5*alog10(counts2/(0.259^2)), color = colors.aquamarine,$
;bottom = 24.29 - 2.5*alog10(newcounts4/(0.259^2)), /noerase, /notrace
