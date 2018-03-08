PRO areavsaxis

close,/all

aarr = fltarr(8000)
area = fltarr(8000)
i = 0

openr, lun,"/n/Godiva1/jkrick/A3888/final/fullr2.cat", /get_lun
WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, o, xcenter, ycenter, a, b, e, f, m, isoarea, fwhm, pa, bkgd, apflux,isocorflux,isoflux
 
    IF (f GT 0) AND (apflux GE 0.14 )AND (isoarea GE 14.0) AND (fwhm GT 5.4) AND (m LT 21) THEN BEGIN

        IF sqrt((xcenter -1534)^2 + (ycenter - 1800)^2) GT 900 THEN BEGIN
            area[i] = sqrt(isoarea /(!PI*(1-e)))
            aarr[i] = a
            i = i + 1
            print, xcenter, ycenter
        ENDIF
    ENDIF

ENDWHILE
area = area[0:i-1]
aarr = aarr[0:i-1]

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/final/isoarea.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot, area, aarr, psym = 2, thick = 3, xrange=[0,40], xtitle = "area to get a", ytitle = "a from 2nd moment"

x = findgen(80)
y = 0.5*x
oplot, x, y, thick = 3
oplot, x, (1./3.)*x, thick =3, linestyle =2
device, /close
set_plot, mydevice


END
