PRO colorsb
close,/all

device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva2/jkrick/A141/original/shell.ps', /portrait, $
  BITS=8, scale_factor=0.9 ;, /color


fits_read, "/n/Godiva2/jkrick/A141/original/larger.fits", data, header

;data = data[700:2400,800:2500]

;FOR x =1000, 1020, 1 DO BEGIN
;    FOR y = 1000, 1020, 1 DO BEGIN
;        print, data[x,y]
;    ENDFOR
;ENDFOR

FOR x =1400,1620, 1 DO BEGIN
    FOR y=615,900, 1 DO BEGIN
;        IF data[x,y] GT 1.0 THEN data[x,y] = 4
;        IF data[x,y] GT 0.6 AND data[x,y] LE 1.0 THEN data[x,y] = 3
;        IF data[x,y] GT 0.2 AND data[x,y] LE 0.6 THEN data[x,y] = 2
;        IF data[x,y] GT 0.0 AND data[x,y] LE 0.2 THEN data[x,y] = 1
;        IF data[x,y] LE 0.0 THEN data[x,y] = 0

        IF data[x,y] GT 0.5 THEN data[x,y] = 8
;        IF data[x,y] GT 0.2 AND data[x,y] LE 0.5 THEN data[x,y] = 9
;        IF data[x,y] GT 0.15 AND data[x,y] LE 0.2 THEN data[x,y] = 8
        IF data[x,y] GT 0.1 AND data[x,y] LE 0.05 THEN data[x,y] = 7
        IF data[x,y] GT 0.07 AND data[x,y] LE 0.1 THEN data[x,y] = 6
;        IF data[x,y] GT 0.1 AND data[x,y] LE 0.12 THEN data[x,y] = 5
;
        IF data[x,y] GT 0.05 AND data[x,y] LE 0.07 THEN data[x,y] = 5
;        IF data[x,y] GT 0.03 AND data[x,y] LE 0.05 THEN data[x,y] = 6
        IF data[x,y] GT 0.011657 AND data[x,y] LE 0.05 THEN data[x,y] = 4
;        IF data[x,y] GT 0.007355 AND data[x,y] LE 0.011657 THEN data[x,y] = 4
        IF data[x,y] GT 0.004641 AND data[x,y] LE 0.011657 THEN data[x,y] = 3
        IF data[x,y] GT 0.002982 AND data[x,y] LE 0.004641 THEN data[x,y] = 2
;        IF data[x,y] GT -0.05 AND data[x,y] LE 0.00464 THEN data[x,y] = 1
 ;       IF data[x,y] GT 2.0 AND data[x,y] LE 2.001848 THEN data[x,y] = 1
       IF data[x,y] GT -0.05 AND data[x,y] LE 0.002982 THEN data[x,y] = 1
;        IF data[x,y] GT -0.05 AND data[x,y] LE 0.0009 THEN data[x,y] = 1
;        IF data[x,y] LT-0.05 THEN data[x,y] = 0

    ENDFOR
ENDFOR

plotimage, bytscl(data[1400:1620,615:900],min=0,max = 7) ,$; bytscl(data, min=1,max=3),$
 /preserve_aspect, /noaxes, ncolors=8

;OPENR, lun,'/n/Godiva1/jkrick/A3888/final/vtab2', /GET_LUN 
;rows= 25
;sma = FLTARR(rows)
;ellip = FLTARR(rows)
;r = 0.0				;radius
;el = 0.0				;mean counts
;pa = 0.0
;x0=0.
;y0=0.

;FOR j=0,rows-1 DO BEGIN
;      READF, lun, r,el,pa,x0,y0
;      sma(j) = r
;      ellip(j) = el
;ENDFOR
;x0 = x0 -700.
;y0 = y0 - 800.
;smb = sma*(1-ellip)
;pa = (pa+90.)*(!PI/180.)   ;in radians;


;npoints = 120
;phi = 2 * !PI * (Findgen(npoints) / (npoints-1))

;FOR k = 3, rows-5, 2 DO BEGIN
;    x = sma(k)*cos(phi)
;    y = smb(k)*sin(phi)

;    xprime = x0 + (x*cos(pa)) - (y*sin(pa))
;    yprime = y0 + (x*sin(pa)) + (y*cos(pa))
; 
;    oplot, xprime, yprime, thick = 3, color = colors.white
;    print, sma(k)
;ENDFOR

oplot, findgen(207)+1200 , findgen(207) -findgen(207)+ 150, thick = 3,color = 245
xyouts, 1200,100,'300 kpc', color = 245, charthick = 3

device, /close
set_plot, mydevice


END
