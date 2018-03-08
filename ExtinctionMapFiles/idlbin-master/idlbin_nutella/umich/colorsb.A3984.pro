PRO colorsb
close,/all

device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva2/jkrick/A3984/iclcolor2.ps', /portrait, $
  BITS=8, scale_factor=0.9 ;, /color


fits_read, "/n/Godiva2/jkrick/A3984/centersmo.r.fits", data, header

data = data[400:2200,600:3000]


FOR x =0,1799, 1 DO BEGIN
    FOR y=0, 2400, 1 DO BEGIN
        IF data[x,y] GT 0.011657 THEN data[x,y] = 7
        IF data[x,y] GT 0.007355 AND data[x,y] LE 0.011657 THEN data[x,y] = 6
        IF data[x,y] GT 0.004641 AND data[x,y] LE 0.007355 THEN data[x,y] = 5
        IF data[x,y] GT 0.002982 AND data[x,y] LE 0.004641 THEN data[x,y] = 4
        IF data[x,y] GT 0.0014 AND data[x,y] LE 0.002982 THEN data[x,y] = 3
        IF data[x,y] GT -0.05 AND data[x,y] LE 0.0014 THEN data[x,y] = 2
        IF data[x,y] LT-0.05 THEN data[x,y] = 1
    ENDFOR
ENDFOR

plotimage, bytscl(data,min = 1, max = 7) ,$; bytscl(data, min=1,max=3),$
 /preserve_aspect, /noaxes, ncolors=7


oplot, findgen(273)+140 , findgen(273) -findgen(273)+ 190, thick = 3,color = 245

xyouts,110,100,'300 kpc', color = 245, charthick = 3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun,'/n/Godiva2/jkrick/A3984/rtab', /GET_LUN 
rows= 17
sma = FLTARR(rows)
ellip = FLTARR(rows)
r = 0.0				;radius
el = 0.0				;mean counts
pa = 0.0
x0=0.
y0=0.

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,el,pa,x0,y0
      sma(j) = r
      ellip(j) = el
ENDFOR
x0 = x0 -400.
y0 = y0 - 600.
smb = sma*(1-ellip)
pa = (pa+90.)*(!PI/180.)   ;in radians


npoints = 120
phi = 2 * !PI * (Findgen(npoints) / (npoints-1))

FOR k =0, rows-1, 1 DO BEGIN
    x = sma(k)*cos(phi)
    y = smb(k)*sin(phi)

    xprime = x0 + (x*cos(pa)) - (y*sin(pa))
    yprime = y0 + (x*sin(pa)) + (y*cos(pa))
 
    oplot, xprime, yprime, thick = 3, color = 245
    print, sma(k)
ENDFOR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;







device, /close
set_plot, mydevice


END
