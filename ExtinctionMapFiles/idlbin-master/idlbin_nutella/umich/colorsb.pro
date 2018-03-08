PRO colorsb
close,/all

device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/final/samemask/iclcolor.ps', /portrait, $
  BITS=8, scale_factor=0.9 ;, /color


;fits_read, "/n/Godiva1/jkrick/A3888/final/samemask/centersmoV2.fit", data, header
fits_read, "/n/Godiva1/jkrick/A3888/final/centersmoV2.fit", data, header

;data = data[700:2400,800:2500]
data = data[200:2400,400:2500]

;FOR x =1000, 1020, 1 DO BEGIN
;    FOR y = 1000, 1020, 1 DO BEGIN
;        print, data[x,y]
;    ENDFOR
;ENDFOR

FOR x =1,2199, 1 DO BEGIN
    FOR y=1, 2099, 1 DO BEGIN
        IF data[x,y] GT 2.011657 THEN data[x,y] = 7
        IF data[x,y] GT 2.007355 AND data[x,y] LE 2.011657 THEN data[x,y] = 6
        IF data[x,y] GT 2.004641 AND data[x,y] LE 2.007355 THEN data[x,y] = 5
        IF data[x,y] GT 2.002982 AND data[x,y] LE 2.004641 THEN data[x,y] = 4
        IF data[x,y] GT 2.0014 AND data[x,y] LE 2.002982 THEN data[x,y] = 3
 ;       IF data[x,y] GT 2.0 AND data[x,y] LE 2.001848 THEN data[x,y] = 1
        IF data[x,y] GT 1.95 AND data[x,y] LE 2.0014 THEN data[x,y] = 2
        IF data[x,y] LT 1.95 THEN data[x,y] = 1
    ENDFOR
ENDFOR

plotimage, bytscl(data) ,$; bytscl(data, min=1,max=3),$
 /preserve_aspect, /noaxes, ncolors=10

OPENR, lun,'/n/Godiva1/jkrick/A3888/final/samemask/vtab2', /GET_LUN 
rows= 25
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
x0 = x0 -700.
y0 = y0 - 800.
smb = sma*(1-ellip)
pa = (pa+90.)*(!PI/180.)   ;in radians


npoints = 120
phi = 2 * !PI * (Findgen(npoints) / (npoints-1))

FOR k =1, rows-7, 2 DO BEGIN
    x = sma(k)*cos(phi)
    y = smb(k)*sin(phi)

    xprime = x0 + (x*cos(pa)) - (y*sin(pa))
    yprime = y0 + (x*sin(pa)) + (y*cos(pa))
 
;    oplot, xprime, yprime, thick = 3, color = colors.white
    print, sma(k)
ENDFOR

oplot, findgen(441)+1000 , findgen(441) -findgen(441)+ 150, thick = 3,color = 245
xyouts, 1060,80,'300 kpc', color = 245, charthick = 3

device, /close
set_plot, mydevice


END
