PRO ploticlronefit
;device, true=24
;device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/final/onefit.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

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
galradius = FLTARR(rows)
galcounts = FLTARR(rows)
galerr = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      galradius(j) = r
      galcounts(j) = c
      galerr(j) = ellerr
ENDFOR

close, lun
free_lun, lun

galcounts = galcounts -2;[0:j-3] - 2
galradius = galradius ;[0:j-3]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun,'/n/Godiva1/jkrick/A3888/final/iclrtab3', /GET_LUN  ;4

;read in the radial profile
rows= 19-2
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

rcounts3 = rcounts3 -2;[0:j-6] - 2
rradius3 = rradius3;[0:j-6]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
out = textoidl('Surface Brightness (mag/arcsec^{-2})')

plot, rradius2* 0.259,24.6 - 2.5*alog10(rcounts2/(0.259^2)),thick = 3, $
YRANGE = [29, 26], xrange = [50,250],ytitle = out, $
psym = 3,xstyle = 8,ystyle = 1,$
xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3

axis, 0, 23.0, xaxis=1, xrange=[0,1.034], xstyle = 1, xthick = 3,charthick = 3

oplot, rradius2* 0.259,24.6 - 2.5*alog10(rcounts2/(0.259^2)),thick = 5;
oplot, galradius* 0.259,24.6 - 2.5*alog10(galcounts/(0.259^2)),thick = 5;,color = colors.black

;;;;;;;;;;;;;;;;;;;;;;;
;plot the shaded regions
;;;;;;;;;;;;;;;;;;;;;;;;
;pxval = [rradius2* 0.259,reverse(rradius4* 0.259)]
;pyval = [24.6 - 2.5*alog10(rcounts2/(0.259^2)),reverse(24.6 - 2.5*alog10(rcounts4/(0.259^2)))]
;;polyfill, pxval,pyval, color = colors.gray
;polyfill, pxval,pyval, color = colors.lightsalmon

;pxval = [rradius2* 0.259,reverse(rradius* 0.259)]
;pyval = [24.6 - 2.5*alog10(rcounts2/(0.259^2)),reverse(24.6 - 2.5*alog10(rcounts/(0.259^2)))]
;;polyfill, pxval,pyval, color = colors.darkgray
;polyfill, pxval,pyval, color = colors.orangered

;;;;;;;;;;;;;;;;;;;
;plot the error bars
;;;;;;;;;;;;;;;;;;;

lower2 = 24.6 - 2.5*alog10((rcounts2 - 0.001429)/(0.259^2))
upper2 = 24.6 - 2.5*alog10((rcounts2 + 0.001429)/(0.259^2)) 
errplot, rradius2* 0.259,lower2, upper2


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;profile fitting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

err = dindgen(20) - dindgen(20) + 1

;Sersic
start=[0.1,10.0,4.0]  
;resultserg= MPFITFUN('sersic',galradius[0:15],galcounts[0:15],err[0:15], start,maxiter=300)        ;galaxies
resultseri = MPFITFUN('sersic',rradius2[0:12],rcounts2[0:12], rerr2[0:12], start)      ;ICL

;Exponential
start = [0.02,300]  
;resultexpg = MPFITFUN('exponential',galradius[0:15],galcounts[0:15],err[0:15], start,maxiter=300)        ;galaxies
resultexpi = MPFITFUN('exponential',rradius2[0:12],rcounts2[0:12], rerr2[0:12], start)      ;ICL

;DeVauc
start = [.1, 10.0]
resultdevi = MPFITFUN('devauc',rradius2[0:12],rcounts2[0:12], rerr2[0:12], start)      ;ICL

;reynolds
start=[.1,100.0]
resultreyi = MPFITFUN('reynolds',rradius2[0:12],rcounts2[0:12], rerr2[0:12], start)      ;ICL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;profile plotting
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

newcounts3 = interpolate(galcounts, [0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])
newradius3 = interpolate(galradius, [0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19])

;sersic 
numbers = (resultseri(0)) * (exp(-7.67*(((rradius2[0:15]/resultseri(1))^(1.0/resultseri(2))) - 1.0)))
oplot, rradius2[0:15]*0.259, 24.6 - 2.5*alog10(numbers/(0.259^2)), thick = 3, color = colors.red

;numbers = (resultserg(0)) * (exp(-7.67*(((newradius3[0:19]/resultserg(1))^(1.0/resultserg(2))) - 1.0)))
;oplot, newradius3[0:19]*0.259, 24.6 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.red


;exponential 
oplot, rradius2[0:15]*0.259, 24.6 - 2.5*alog10(resultexpi(0)*exp(-rradius2[0:15]/resultexpi(1))/(0.259^2) ), $
     thick = 3, color = colors.blue
;oplot, newradius3[0:15]*0.259, 24.6 - 2.5*alog10(resultexpg(0)*exp(-newradius3[0:15]/resultexpg(1))/(0.259^2) ), thick = 3

;devauc 
numbers = (resultdevi(0)) * (exp(-7.67*(((rradius2[0:15]/resultdevi(1))^(1.0/4.0)) - 1.0)))
oplot, rradius2[0:15]*0.259, 24.6 - 2.5*alog10(numbers/(0.259^2)), thick = 3, color = colors.hotpink

numbers = (resultdevi(0)) * (exp(-7.67*(((findgen(1000)/resultdevi(1))^(1.0/4.0)) - 1.0)))
oplot, findgen(1000)*0.259, 24.6 - 2.5*alog10(numbers/(0.259^2)), thick = 3, color = colors.hotpink

;numbers = (result3(0)) * (exp(-7.67*(((newradius3[0:19]/result3(1))^(1.0/4.0)) - 1.0)))
;oplot, newradius3[0:19]*0.259, 24.6 - 2.5*alog10(numbers/(0.259^2)), thick = 5;, color = colors.red

;reynolds 
numbers = resultreyi(0) / ((1 + (rradius2[0:15]/resultreyi(1)))^2)
oplot, rradius2[0:15]*0.259, 24.6 - 2.5*alog10(numbers/(0.259^2)), thick = 3, color = colors.green

;numbers =result3(0) / ((1 + (newradius3/result3(1)))^2)
;oplot, newradius3[0:19]*0.259, 24.6 - 2.5*alog10(numbers/(0.259^2)), thick = 5, color = colors.red

xyouts, 60, 28, "exponential", color = colors.blue
xyouts, 60, 28.2, "sersic", color = colors.red
xyouts, 60, 28.4, "devauc", color = colors.hotpink
xyouts, 60, 28.6, "reynolds", color = colors.green


device, /close
set_plot, mydevice

END



