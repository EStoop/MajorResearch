PRO colorsb
close,/all

device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva4/jkrick/A118/iclcolor.ps', /portrait, $
  BITS=8, scale_factor=0.9 ;, /color


fits_read, "/n/Godiva4/jkrick/A118/original/centersmo.fits", data, header

data = data[400:1900,200:1800]


FOR x =0,1499, 1 DO BEGIN
    FOR y=0, 1600, 1 DO BEGIN
        IF data[x,y] GT 0.011657 THEN data[x,y] = 7
        IF data[x,y] GT 0.007355 AND data[x,y] LE 0.011657 THEN data[x,y] = 6
        IF data[x,y] GT 0.004641 AND data[x,y] LE 0.007355 THEN data[x,y] = 5
        IF data[x,y] GT 0.002982 AND data[x,y] LE 0.004641 THEN data[x,y] = 4
        IF data[x,y] GT 0.0014 AND data[x,y] LE 0.002982 THEN data[x,y] = 3
 ;       IF data[x,y] GT 2.0 AND data[x,y] LE 2.001848 THEN data[x,y] = 1
        IF data[x,y] GT -0.05 AND data[x,y] LE 0.0014 THEN data[x,y] = 2
        IF data[x,y] LT-0.05 THEN data[x,y] = 1
    ENDFOR
ENDFOR

plotimage, xrange=[1,1500],yrange=[1600,200], bytscl(data, min = 0, max = 7) ,$; bytscl(data, min=1,max=3),$
 /preserve_aspect, /noaxes, ncolors=8

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

OPENR, lun,'/n/Godiva4/jkrick/A118/original/iclVtab', /GET_LUN 
rows= 7
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
y0 = y0 - 200.
smb = sma*(1-ellip)
pa = (pa+90.)*(!PI/180.)   ;in radians


npoints = 120
phi = 2 * !PI * (Findgen(npoints) / (npoints-1))

FOR k =0, rows-1, 1 DO BEGIN
    x = sma(k)*cos(phi)
    y = smb(k)*sin(phi)

    xprime = x0 + (x*cos(pa)) - (y*sin(pa))
    yprime = y0 + (x*sin(pa)) + (y*cos(pa))
 
;    oplot, xprime, yprime, thick = 3, color = 245
    print, sma(k)
ENDFOR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
OPENR, lun,'/n/Godiva4/jkrick/A118/original/iclV2tab', /GET_LUN 
rows= 6
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
y0 = y0 - 200.
smb = sma*(1-ellip)
pa = (pa+90.)*(!PI/180.)   ;in radians


npoints = 120
phi = 2 * !PI * (Findgen(npoints) / (npoints-1))

FOR k =0, rows-1, 1 DO BEGIN
    x = sma(k)*cos(phi)
    y = smb(k)*sin(phi)

    xprime = x0 + (x*cos(pa)) - (y*sin(pa))
    yprime = y0 + (x*sin(pa)) + (y*cos(pa))
 
;    oplot, xprime, yprime, thick = 3, color = 245
    print, sma(k)
ENDFOR
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
oplot, findgen(200)+140 , findgen(200) -findgen(200)+ 230, thick = 3,color = 245
xyouts,150,280,'400 kpc', color = 245, charthick = 3
xyouts, 1200, 1570, 'AC118', color = 245, charthick = 3
device, /close
set_plot, mydevice


END
