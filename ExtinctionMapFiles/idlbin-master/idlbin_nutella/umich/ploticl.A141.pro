PRO ploticl

;device, true=24
;device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva2/jkrick/A141/icl2.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

area = fltarr(6)
area = [123502.3,128893.3,57726.8,53146.3,1318966,123659]
aver = [.00653,.003295,.005026,.003162,0.004802,.002062]
aveV = [.002723,.00091,.001986,.0008308,.001909,.0003458]
xarr = [273,390,147,204,328, 457]


newarr = [xarr[0]/2., (xarr[1] - xarr[0])/2. + xarr[0]]
plot, newarr*0.259, 24.6 -2.5*(alog10(aver[0:1]/(0.259^2))), thick = 3, $
YRANGE = [30, 26], xrange = [0,200],ytitle = out, $
xstyle = 8,ystyle = 1,$
xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3
axis, 0, 26.0, xaxis=1, xrange=[0,1.120], xstyle = 1, xthick = 3,charthick = 3
xerr=[(xarr[0]/2.)*0.259,((xarr[1] - xarr[0])/2.)*0.259 ]
yerr = [0,0]

oploterror, newarr*0.259, 24.6 -2.5*(alog10(aver[0:1]/(0.259^2))),$
  xerr, yerr, psym = 2



newarr = [xarr[2]/2., (xarr[3] - xarr[2])/2. + xarr[2]]
oplot,newarr*0.259, 24.6 -2.5*(alog10(aver[2:3]/(0.259^2))), thick = 3,$
  color = colors.gray
xerr=[(xarr[2]/2.)*0.259,((xarr[3] - xarr[2])/2.)*0.259 ]
yerr = [0,0]

oploterror, newarr*0.259, 24.6 -2.5*(alog10(aver[2:3]/(0.259^2))),$
  xerr, yerr, psym = 2



newarr = [xarr[4]/2., (xarr[5] - xarr[4])/2. + xarr[4]]
oplot,newarr*0.259, 24.6 -2.5*(alog10(aver[4:5]/(0.259^2))), thick = 3,$
  color = colors.darkgray
xerr=[(xarr[4]/2.)*0.259,((xarr[5] - xarr[4])/2.)*0.259 ]
yerr = [0,0]

oploterror, newarr*0.259, 24.6 -2.5*(alog10(aver[4:5]/(0.259^2))),$
  xerr, yerr, psym = 4

;-------------------------------------------

newarr = [xarr[0]/2., (xarr[1] - xarr[0])/2. + xarr[0]]
oplot, newarr*0.259, 24.3 -2.5*(alog10(aveV[0:1]/(0.259^2))), thick = 3,$
  color = colors.black, linestyle = 2

xerr=[(xarr[0]/2.)*0.259,((xarr[1] - xarr[0])/2.)*0.259 ]
yerr = [0,0]

oploterror, newarr*0.259, 24.3 -2.5*(alog10(aveV[0:1]/(0.259^2))),$
  xerr, yerr, psym = 2



newarr = [xarr[2]/2., (xarr[3] - xarr[2])/2. + xarr[2]]
oplot,newarr*0.259, 24.3 -2.5*(alog10(aveV[2:3]/(0.259^2))), thick = 3,$
  color = colors.gray, linestyle = 2
xerr=[(xarr[2]/2.)*0.259,((xarr[3] - xarr[2])/2.)*0.259 ]
yerr = [0,0]

oploterror, newarr*0.259, 24.3 -2.5*(alog10(aveV[2:3]/(0.259^2))),$
  xerr, yerr, psym = 2



newarr = [xarr[4]/2., (xarr[5] - xarr[4])/2. + xarr[4]]
oplot,newarr*0.259, 24.3 -2.5*(alog10(aveV[4:5]/(0.259^2))), thick = 3,$
  color = colors.darkgray,linestyle = 2
xerr=[(xarr[4]/2.)*0.259,((xarr[5] - xarr[4])/2.)*0.259 ]
yerr = [0,0]

oploterror, newarr*0.259, 24.3 -2.5*(alog10(aveV[4:5]/(0.259^2))),$
  xerr, yerr, psym = 4

;----------------------------------------
;make a color plot as well
color = [0.65,1.1,1.0,1.1,0.7,1.6]
newarr = [xarr[0]/2., (xarr[1] - xarr[0])/2. + xarr[0]]
plot, newarr*0.259, color[0:1],  thick = 3,$
  YRANGE = [0, 2], xrange = [0,200],ytitle = out, $
  xstyle = 8,ystyle = 1,charthick = 3,$
  xtitle = 'Semi-major Axis (arcseconds)', xthick = 3, ythick = 3
axis, 0, 2.0, xaxis=1, xrange=[0,1.120], xstyle = 1, xthick = 3,charthick = 3

xerr=[(xarr[0]/2.)*0.259,((xarr[1] - xarr[0])/2.)*0.259 ]
yerr = [0,0]

oploterror,  newarr*0.259, color[0:1], xerr, yerr, psym = 2



newarr = [xarr[2]/2., (xarr[3] - xarr[2])/2. + xarr[2]]
oplot,newarr*0.259, color[2:3], thick = 3,$
  color = colors.gray
xerr=[(xarr[2]/2.)*0.259,((xarr[3] - xarr[2])/2.)*0.259 ]
yerr = [0,0];

oploterror, newarr*0.259, color[2:3],$
  xerr, yerr, psym = 2



newarr = [xarr[4]/2., (xarr[5] - xarr[4])/2. + xarr[4]]
oplot,newarr*0.259, color[2:3], thick = 3,$
  color = colors.darkgray
xerr=[(xarr[4]/2.)*0.259,((xarr[5] - xarr[4])/2.)*0.259 ]
yerr = [0,0]

oploterror, newarr*0.259, color[2:3],$
  xerr, yerr, psym = 4



device, /close
set_plot, mydevice

END
