PRO beyondlf

close,/all
colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva7/jkrick/lf.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

l=[-2,0,2,0,-2]
m=[0,2,0,-2,0]
usersym,l,m,/fill

x = findgen(100)/ 10. + 18
;x = x / 10.0 - 30
plot,  x, (100.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-x))))*(10^(0.4*(21.0-x)*(-1.0+1))), thick = 3, /ylog ,$
 yrange =[1E-2,1E3],ystyle=1, ytitle = "mag^-1 Mpc^-1", $
charthick = 3, xthick = 3, ythick = 3
oplot,  x, (100.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-x))))*(10^(0.4*(21.0-x)*(-0.6+1))), thick = 3
oplot,  x, (100.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-x))))*(10^(0.4*(21.0-x)*(-0.9+1))), thick = 3
oplot,  x, (100.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-x))))*(10^(0.4*(21.0-x)*(-1.25+1))), thick = 3
xyouts, 24, 1, "Vary alpha", charthick = 3
p = fltarr(2) + 21.0
p1 = fltarr(2) + (100.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-21.0))))*(10^(0.4*(21.0-21.0)*(-1.0+1)))
oplot, p, p1, psym = 8

plot,  x, (100.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-x))))*(10^(0.4*(21.0-x)*(-1.0+1))), thick = 3, $
color = colors.blue, /ylog , yrange =[1E-2,1E3],ystyle=1, $
ytitle = "mag^-1 Mpc^-1", charthick = 3, xthick = 3, ythick = 3
oplot,  x, (100.0*0.4*alog(10))*(exp(-10^(0.4*(20.0-x))))*(10^(0.4*(20.0-x)*(-1.0+1))), thick = 3, color = colors.blue
oplot,  x, (100.0*0.4*alog(10))*(exp(-10^(0.4*(22.0-x))))*(10^(0.4*(22.0-x)*(-1.0+1))), thick = 3, color = colors.blue
p = fltarr(2) + 21.0
p1 = fltarr(2) + (100.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-21.0))))*(10^(0.4*(21.0-21.0)*(-1.0+1)))
oplot, p, p1, psym = 8
p = fltarr(2) + 20.0
p1 = fltarr(2) + (100.0*0.4*alog(10))*(exp(-10^(0.4*(20.0-20.0))))*(10^(0.4*(20.0-20.0)*(-1.0+1)))
oplot, p, p1, psym = 8
p = fltarr(2) + 22.0
p1 = fltarr(2) + (100.0*0.4*alog(10))*(exp(-10^(0.4*(22.0-22.0))))*(10^(0.4*(22.0-22.0)*(-1.0+1)))
oplot, p, p1, psym = 8
xyouts, 24, 1, "Vary M*", charthick = 3, color = colors.blue

plot,  x, (100.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-x))))*(10^(0.4*(21.0-x)*(-1.0+1))), thick = 3, $
color = colors.red, /ylog , yrange =[1E-2,1E3],ystyle=1, xtitle = "apparent magnitude", $
ytitle = "mag^-1 Mpc^-1", charthick = 3, xthick = 3, ythick = 3
oplot,  x, (400.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-x))))*(10^(0.4*(21.0-x)*(-1.0+1))), thick = 3, color = colors.red
oplot,  x, (200.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-x))))*(10^(0.4*(21.0-x)*(-1.0+1))), thick = 3, color = colors.red
p = fltarr(2) + 21.0
p1 = fltarr(2) + (100.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-21.0))))*(10^(0.4*(21.0-21.0)*(-1.0+1)))
oplot, p, p1, psym = 8
p = fltarr(2) + 21.0
p1 = fltarr(2) + (200.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-21.0))))*(10^(0.4*(21.0-21.0)*(-1.0+1)))
oplot, p, p1, psym = 8
p = fltarr(2) + 21.0
p1 = fltarr(2) + (400.0*0.4*alog(10))*(exp(-10^(0.4*(21.0-21.0))))*(10^(0.4*(21.0-21.0)*(-1.0+1)))
oplot, p, p1, psym = 8
xyouts, 24, 1.0, "Vary normalization", charthick = 3, color = colors.red
device, /close
set_plot, mydevice
END
