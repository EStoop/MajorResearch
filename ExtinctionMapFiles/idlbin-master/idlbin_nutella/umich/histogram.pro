PRO histogram
close, /all
colors = GetColor(/load, Start=1)

magarr = fltarr(32000)
i = 0

openr, sexlun, "/Users/jkrick/palomar/lfc/catalog/SExtractor.r.cat", /get_lun
WHILE (NOT EOF(sexlun)) DO BEGIN
    READF, sexlun, o, xwcs, ywcs, xcenter, ycenter, fluxauto, magauto, magerr, fwhm, isoarea
   

    IF (magauto GT 0) and (magauto LT 50) THEN BEGIN
        magarr(i) = magauto
        i = i + 1
    ENDIF

ENDWHILE

magarr = magarr(0:i-1)
close,sexlun
free_lun, sexlun

device, true=24
device, decomposed=0

!p.multi = [0, 0, 1]
ps_open, file = "/Users/jkrick/palomar/LFC/catalog/rhist.ps", /portrait, xsize = 4, ysize = 4,/color


plothist, magarr, xhist, yhist, bin=0.10, /noprint;,xrange=[0,480]

plot, xhist, yhist,title = "rmag distribution ",$
thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "r mag",$
xstyle = 1


err = fltarr(i) + 1.
start = [25.,5., 5000.]
result = MPFITFUN('gauss',xhist,yhist, err, start) 
xarr = (findgen(1000.))

oplot, xarr, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.),$
  thick = 3, color = colors.red
;========================================
magarr = fltarr(42000)
i = 0.

openr, sexlun, "/Users/jkrick/palomar/lfc/catalog/SExtractor.g.cat", /get_lun
WHILE (NOT EOF(sexlun)) DO BEGIN
    READF, sexlun, o, xwcs, ywcs, xcenter, ycenter, fluxauto, magauto, magerr, fwhm, isoarea
   

    IF (magauto GT 0) and (magauto LT 50) THEN BEGIN
        magarr(i) = magauto
        i = i + 1.
    ENDIF

ENDWHILE

magarr = magarr(0:i-1)
close,sexlun
free_lun, sexlun


plothist, magarr, xhist, yhist, bin=0.10, /noprint,/noplot;,xrange=[0,480]

plot, xhist, yhist,title = "gmag distribution ",$
thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "g mag",$
xstyle = 1


err = fltarr(i) + 1.
start = [25.,5., 5000.]
result = MPFITFUN('gauss',xhist,yhist, err, start) 
xarr = (findgen(1000.))

oplot, xarr, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.),$
  thick = 3, color = colors.red
;========================================
magarr = fltarr(42000)
i = 0.

openr, sexlun, "/Users/jkrick/palomar/lfc/catalog/SExtractor.i.cat", /get_lun
WHILE (NOT EOF(sexlun)) DO BEGIN
    READF, sexlun, o, xwcs, ywcs, xcenter, ycenter, fluxauto, magauto, magerr, fwhm, isoarea
   

    IF (magauto GT 0) and (magauto LT 50) THEN BEGIN
        magarr(i) = magauto
        i = i + 1.
    ENDIF

ENDWHILE

magarr = magarr(0:i-1)
close,sexlun
free_lun, sexlun


plothist, magarr, xhist, yhist, bin=0.10, /noprint,/noplot;,xrange=[0,480]

plot, xhist, yhist,title = "imag distribution ",$
thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "i mag",$
xstyle = 1


err = fltarr(i) + 1.
start = [25.,5., 5000.]
result = MPFITFUN('gauss',xhist,yhist, err, start) 
xarr = (findgen(1000.))

oplot, xarr, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.),$
  thick = 3, color = colors.red
;========================================
magarr = fltarr(162000)
i = 0.

openr, sexlun, "/Users/jkrick/palomar/lfc/catalog/SExtractor.u.cat", /get_lun
WHILE (NOT EOF(sexlun)) DO BEGIN
    READF, sexlun, o, xwcs, ywcs, xcenter, ycenter, fluxauto, magauto, magerr, fwhm, isoarea
   

    IF (magauto GT 0) and (magauto LT 50) THEN BEGIN
        magarr(i) = magauto
        i = i + 1.
    ENDIF

ENDWHILE

magarr = magarr(0:i-1)
close,sexlun
free_lun, sexlun


plothist, magarr, xhist, yhist, bin=0.10, /noprint,/noplot;,xrange=[0,480]

plot, xhist, yhist,title = "umag distribution ",$
thick = 3, charthick = 3, xthick = 3, ythick = 3, xtitle = "u mag",$
xstyle = 1


err = fltarr(i) + 1.
start = [25.,5., 5000.]
result = MPFITFUN('gauss',xhist,yhist, err, start) 
xarr = (findgen(1000.))

oplot, xarr, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - result(0))/result(1))^2.),$
  thick = 3, color = colors.red

ps_close, /noprint, /noid

end
