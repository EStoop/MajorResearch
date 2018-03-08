PRO lumfunc
device, decomposed=0
colors = GetColor(/load, Start=1)

close,/all
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva4/jkrick/A118/lumfunc.ps', /portrait, $
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

mag2 = fltarr(9)
n2 = fltarr(9)
err2 = fltarr(9)

mag2 = [18.5,19,19.5,20,20.5,21,21.5,22,22.5]
n2 = [6.4,5.4,15.,31.,50.,49.,89.,87.,80.]
err2 = [2.0,2.,2.,2.,2.,2.,2.,2.,2.]

start = [70.]
result2 = MPFITFUN('schechter2',mag2, n2, err2, start, perr = testerr)


plot, mag2, n2, psym = 2, thick = 3,/ylog, xtitle = "magnitude(V)", ytitle = "Number per 8.6 x 8.6 arcmin",$
title = "Busarello 2002 Luminosity Distribution", charthick = 3, xthick = 3, ythick = 3
errplot, mag2,n2-err2,n2+err2

oplot, mag2,(result2(0)*0.4*alog(10))*(exp(-10^(0.4*(20.1-mag2))))*(10^(0.4*(20.1-mag2)*((-1.02)+1))), thick = 3


s = "M* 20.1"
xyouts, 21.2,10, s, charthick = 3
s = "alpha -1.02"
xyouts, 21.2, 5, s, charthick = 3
s = "phistar"+ string(result2(0))
xyouts, 21.2, 3, s, charthick = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

mag3 = fltarr(10)
n3 = fltarr(10)
err3 = fltarr(10)

mag3 = [17.5,18.0,18.5,19,19.5,20,20.5,21,21.5,22]
n3 =   [1,5.8,10,25,48,50,87,69,90,80]
err3 = [3.,3.,3.,3.,3.,3.,3.,3.,3.,3.]

start = [70.]
result3 = MPFITFUN('schechter3',mag3, n3, err3, start, perr = testerr)


plot, mag3, n3, psym = 2, thick = 3,/ylog, xtitle = "magnitude(R)", ytitle = "Number per 8.6 x 8.6 arcmin",$
title = "Busarello 2002 Luminosity Distribution", charthick = 3, xthick = 3, ythick = 3
errplot, mag3,n3-err3,n3+err3

oplot, mag3,(result3(0)*0.4*alog(10))*(exp(-10^(0.4*(19.4-mag3))))*(10^(0.4*(19.4-mag3)*((-0.94)+1))), thick = 3


s = "M* 19.4"
xyouts, 20.2,10, s, charthick = 3
s = "alpha -0.94"
xyouts, 20.2, 5, s, charthick = 3
s = "phistar"+ string(result3(0))
xyouts, 20.2, 3, s, charthick = 3


device, /close
set_plot, mydevice


END
