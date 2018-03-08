PRO fwhmcutoff
close, /all

device, true=24
device, decomposed=0
colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva4/jkrick/A4010/fwhmcutoff.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

magarr = fltarr(30000)
fwhmarr = fltarr(30000)
i = 0

commandline = 'sex ' + '/n/Godiva6/jkrick/A3880/original/fullB.fits'+ ' -c  /n/Godiva4/jkrick/A3984/default.sex'
spawn, commandline

openr, lun, "/n/Godiva4/jkrick/A3984/SExtractor.cat", /GET_LUN
WHILE (NOT eof(lun)) DO BEGIN
    READF,lun,junk, x, y, junk,junk,junk,junk,m,junk,junk,junk,f,junk,junk
    magarr(i) = m
    fwhmarr(i) = f
    i = i + 1
ENDWHILE
magarr = magarr[0:i - 1]
fwhmarr = fwhmarr[0:i - 1]
close, lun
free_lun, lun


plot, magarr, fwhmarr,  thick = 1, psym = 2, xthick = 3, ythick = 3,charthick = 3,$
          xrange = [15,26], yrange = [0,8], ystyle = 1,$
	  ytitle = "fwhm", title = 'SExtractored objects from the combined V image', $
          xtitle = "r magnitude"

k = 0
starfarr = fltarr(12000)
starmarr = fltarr(12000)
FOR j = 0, i - 1, 1 DO BEGIN
    IF (fwhmarr(j) LT 7.0) AND (fwhmarr(j) GE 3.2 ) THEN BEGIN
       starfarr(k) = fwhmarr(j)
       starmarr(k) = magarr(j)
       k = k + 1
    ENDIF
ENDFOR

starfarr = starfarr[0:k-1]
starmarr = starmarr[0:k-1]

oplot, starmarr, starfarr, psym = 2, color = colors.red
device, /close
set_plot, mydevice

close, /all

END
