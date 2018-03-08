PRO colorgradient
;device, true=24
;device, decomposed=0

;colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 2]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/colorgradient.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color, /inches, ysize = 8, yoffset = 0.5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva1/jkrick/A3888/ellipse/iclVtab2', /GET_LUN
rows= 33
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

OPENR, lun,'/n/Godiva1/jkrick/A3888/ellipse/iclrtab2', /GET_LUN
rows = 22
rradius = FLTARR(rows)
rcounts = FLTARR(rows)
rerr = FLTARR(rows)
rstop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      rradius(j) = r
      rcounts(j) = c
      rerr(j) = ellerr 
      rstop(j) = s
ENDFOR

close, lun
free_lun, lun

;interpolate so that I am sampling r and V at the same radii
;using a spline interpolation
rx = rradius 
ry = rcounts - 2
rt = [110.,170.,230.,290., 350.,410.,620., 650., 680.,710.,740.,800.,830.,860.,890.,920.,950.,980.,1010.,1040.,1070.]

rz = spline(rx,ry,rt)

rx = [rx ,rt]
ry = [ry ,rz]
sortindex  = sort(rx)
sortrx = rx[sortindex]
sortry = ry[sortindex]

;----------------------
Vx = Vradius 
Vy = Vcounts - 2
Vt =  [630., 670.,  720., 820.,  880. , 940.  , 1000. ,1060. ,1120.,  1180. ]

Vz = spline(Vx,Vy,Vt)

Vx = [Vx ,Vt]
Vy = [Vy ,Vz]
sortindex  = sort(Vx)
sortVx = Vx[sortindex]
sortVy = Vy[sortindex]


plot, sortry[5:43-5] , sortVy[5:43-5], psym = 2, thick = 3, ytitle = 'V Intensity ', $
title = 'Color gradient in ICL ',charthick = 3, xthick = 3, ythick = 3, xstyle = 12

print, sortVx[5:43-5]

;fit a line to the color color plot
err = dindgen(12) - dindgen(12) + 1
start = [0.5,0.001]
result = MPFITFUN('linear',sortry[5:43-5],sortVy[5:43-5], err, start)
print, result
oplot,sortry[5:43-5], result(0)*(sortry[5:43-5]) + result(1), thick = 3
 
;make a plot of the residuals
plot, sortry[5:43-5], (result(0)*(sortry[5:43-5]) + result(1)) - (sortVy[5:43-5]), psym = 2, thick = 3,$
ytitle = 'residual', xtitle = 'r intensity', charthick = 3, xthick = 3, ythick = 3, POS=[0.125, 0.35, 0.964, 0.57]

oplot, findgen(2) - findgen(2), thick = 3
device, /close
set_plot, mydevice


END
