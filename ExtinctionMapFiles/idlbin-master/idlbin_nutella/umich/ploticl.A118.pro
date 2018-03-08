PRO ploticl

;device, true=24
;device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva7/jkrick/A118/iclplot.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

area = fltarr(4)
area = [21871.8,58358.2,27401.,93032.]
aver = [.01461,.0059,.01075,.004679]
aveV = [.00631,.0022,.004604,.001794]
xarr = [118.,226.,98.,205.]


newarr = [xarr[0]/2., (xarr[1] - xarr[0])/2. + xarr[0]]
plot, newarr*0.259, 24.6 -2.5*(alog10(aver[0:1]/(0.259^2))), thick = 3, $
YRANGE = [30, 26], xrange = [0,200],ytitle = out, $
xstyle = 8,ystyle = 1,$
xtitle = 'Semi-major Axis (arcseconds)', charthick = 3, xthick = 3, ythick = 3
axis, 0, 26.0, xaxis=1, xrange=[0,1.544], xstyle = 1, xthick = 3,charthick = 3
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




;----------------------------------------
;make a color plot as well
color = [0.6,0.8,0.6,0.7]
newarr = [xarr[0]/2., (xarr[1] - xarr[0])/2. + xarr[0]]
plot, newarr*0.259, color[0:1],  thick = 3,$
  YRANGE = [0, 1], xrange = [0,200],ytitle = out, $
  xstyle = 8,ystyle = 1,charthick = 3,$
  xtitle = 'Semi-major Axis (arcseconds)', xthick = 3, ythick = 3
axis, 0, 1.0, xaxis=1, xrange=[0,1.544], xstyle = 1, xthick = 3,charthick = 3

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





device, /close
set_plot, mydevice

END
