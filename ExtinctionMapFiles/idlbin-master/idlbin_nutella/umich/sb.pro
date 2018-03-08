PRO sb
device, true=24
device, decomposed=0

close, /all
magarr = fltarr(11000)
centsbarr = fltarr(11000)
i = 0

openr, lun, "/n/Godiva2/jkrick/A141/SExtractor2.cat", /get_lun
WHILE (NOT EOF(lun)) DO BEGIN
    readf, lun, num, x,y,a,b,e, flux, mag, magerr, isoarea, fwhm,theta, back,fluxaper
    IF ( e GT 0.1) AND (fwhm GT 3.7) AND (mag LT 25.)  AND (fluxaper GT 0.)THEN BEGIN
        magarr[i] = mag
        var = fluxaper /(!Pi*(7.7^2))   ;counts per second per pixel
        centsbarr[i] = 24.6-2.5*(alog10(var /(0.259^2)))
        i = i + 1
    ENDIF

ENDWHILE

magarr = magarr(0:i-1)
centsbarr = centsbarr(0:i-1)
colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva2/jkrick/A141/sb.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot, centsbarr, magarr, psym = 2, thick = 3, xrange=[15,30], yrange=[15,30],$
 xtitle = "central 2 arcsec diamter SB", ytitle = "isophotal mag", title = "A141 galaxies"

oplot, findgen(30), findgen(30)-2, color = colors.red

device, /close
set_plot, mydevice


END
