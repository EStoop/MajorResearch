PRO lumfunc
device, decomposed=0
;colors = GetColor(/load, Start=1)

close,/all

mag = fltarr(12)
n = fltarr(12)
err = fltarr(12)

;for A3888 in Driver etal. 1998
;mag = [-24.25,-23.75,-23.25,-22.75,-22.25,-21.75,-21.25,-20.75,-20.25,-19.75,-19.25,-18.75,-18.25,-17.75,-17.25]
;n = [4.,3.,8.,13.,22.,24.,26.,38.,35.,29.,27.,24.,43.,51.,11.]
;err = [2.,2.,3.,4.,5.,5.,6.,7.,8.,8.,10.,11.,14.,16.,19.]

mag = [-27,-26,-25.5,-25,-24.5,-24,-23.5,-23,-22.2,-22,-21.5,-21]
n = [4.3,6.5,20,21,30,43,40.5,51,41,50,68,80]
err=[3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.]

mag2 = fltarr(9)
n2 = fltarr(9)
err2 = fltarr(9)

mag2 = [18.5,19,19.5,20,20.5,21,21.5,22,22.5]
n2 = [6.4,5.4,15.,31.,50.,49.,89.,87.,80.]
err2 = [2.0,2.,2.,2.,2.,2.,2.,2.,2.]

;start = [-22.5, -1.0,51.]
start = [50.]
result1 = MPFITFUN('schechter1',mag, n, err, start, perr = testerr)

start = [20.1,-1.02,70.]
result2 = MPFITFUN('schechter',mag2, n2, err2, start, perr = testerr)





mydevice = !D.NAME
!p.multi = [0, 0, 2]
SET_PLOT, 'ps'

device, filename = '/Users/jkrick/Virgo/testlumfunc.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot, mag, n, psym = 2, thick = 3,/ylog, xtitle = "Absolute Magnitude(K)", ytitle = "Number per fov",$
title = "andreon 2001 Luminosity Distribution", charthick = 3, xthick = 3, ythick = 3
errplot, mag,n-err,n+err

oplot, mag,(result1(0)*0.4*alog(10))*(exp(-10^(0.4*((-25.3)-mag))))*(10^(0.4*((-25.3)-mag)*((-1.1)+1))), thick = 3

oplot, mag,(31.21*0.4*alog(10))*(exp(-10^(0.4*((-25.3)-mag))))*(10^(0.4*((-25.3)-mag)*((-1.1)+1))), thick = 3
oplot, mag,(102*0.4*alog(10))*(exp(-10^(0.4*((-25.3)-mag))))*(10^(0.4*((-25.3)-mag)*((-1.1)+1))), thick = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

plot, mag2, n2, psym = 2, thick = 3,/ylog, xtitle = "magnitude(V)", ytitle = "Number per fov",$
title = "Busarello 2002 Luminosity Distribution", charthick = 3, xthick = 3, ythick = 3
errplot, mag2,n2-err2,n2+err2

oplot, mag2,(result2(2)*0.4*alog(10))*(exp(-10^(0.4*(result2(0)-mag))))*(10^(0.4*(result2(0)-mag)*(result2(1)+1))), thick = 3

;oplot, mag,(50.*0.4*alog(10))*(exp(-10^(0.4*(result1(0)-mag))))*(10^(0.4*(result1(0)-mag)*(result1(1)+1))), thick = 3, color = colors.blue

;s = "M*"+ string(result1(0))
;xyouts, -21,2., s, charthick = 3
;s = "alpha"+ string(result1(1))
;xyouts, -21.2, 1.6, s, charthick = 3
;s = "phistar"+ string(result1(2))
;xyouts, -21.4, 1.3, s, charthick = 3

oplot, mag, fltarr(15)+result1(0)

device, /close
set_plot, mydevice

print, int_tabulated(mag,n)


print, schechter(mag,result1)
print, qsimp(wvar,-24.25,-17.25,/double)

print, schechter(mag, result1)
print, schechter(findgen(27), result1)


print, qsimp('testfunc2',1,12)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



END

FUNCTION wvar, X
     return, (10^(-0.4*(result(1)+1.)*X))*exp(-10^(0.4*(result1(0)-X)))
END
FUNCTION testfunc2,X
  return, X^2
END
