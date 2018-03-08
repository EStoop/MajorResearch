PRO lumfunc
device, decomposed=0
colors = GetColor(/load, Start=1)

close,/all
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva4/jkrick/A114/lumfunc.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color


mag = fltarr(12)
n = fltarr(12)
err = fltarr(12)

mag = [-27,-26,-25.5,-25,-24.5,-24,-23.5,-23,-22.2,-22,-21.5,-21]
n = [4.3,6.5,20,21,30,43,40.5,51,41,50,68,80]
err=[3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.]

;start = [-22.5, -1.0,51.]
start = [50.]
result1 = MPFITFUN('schechter1',mag, n, err, start, perr = testerr)


plot, mag, n, psym = 2, thick = 3,/ylog, xtitle = "Absolute Magnitude(K)", ytitle = "Number per 5 x 5 arcmin",$
title = "andreon 2001 Luminosity Distribution", charthick = 3, xthick = 3, ythick = 3
errplot, mag,n-err,n+err

oplot, mag,(result1(0)*0.4*alog(10))*(exp(-10^(0.4*((-25.3)-mag))))*(10^(0.4*((-25.3)-mag)*((-1.03)+1))), thick = 3

;oplot, mag,(31.21*0.4*alog(10))*(exp(-10^(0.4*((-25.3)-mag))))*(10^(0.4*((-25.3)-mag)*((-1.1)+1))), thick = 3
;oplot, mag,(102*0.4*alog(10))*(exp(-10^(0.4*((-25.3)-mag))))*(10^(0.4*((-25.3)-mag)*((-1.1)+1))), thick = 3

s = "M* -25.3"
xyouts, -22.9,10, s, charthick = 3
s = "alpha -1.03"
xyouts, -22.9, 5, s, charthick = 3
s = "phistar"+ string(result1(0))
xyouts, -22.9, 3, s, charthick = 3


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



device, /close
set_plot, mydevice


END
