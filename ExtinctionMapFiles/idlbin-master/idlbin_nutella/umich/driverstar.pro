PRO driverstar
;device, true=24
;device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/final/junk.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color


openr, lun, "/n/Godiva1/jkrick/A3888/final/junk.prof", /get_lun
radarr = fltarr(201)
countarr = fltarr(201)
i = 0
WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, rad,counts
    radarr(i) = rad
    countarr(i) = counts
    i = i + 1
ENDWHILE
print, radarr
print, countarr
close, lun
free_lun, lun
;print, radarr


plot, alog10(radarr*0.259), alog10(countarr), yrange=[-3,0.5]

start = [230060.0]
err = dindgen(i) - dindgen(i) + 1.
result = mpfitfun('rcubed',radarr[102,*],countarr[102,*], err, start)

oplot, alog10(radarr*0.259),  alog10(result(0)*(radarr^(-3.0)))

;------------------------------------


;the homemade star
datafile = "/n/Godiva1/jkrick/A3888/final/satstarprof.prof"
openr, lun1,datafile, /get_lun
radarr1 = fltarr(205)
countarr1 = fltarr(205)
i = 0
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, rad,counts
    radarr1(i) = rad
    countarr1(i) = counts
    i = i + 1
ENDWHILE


oplot, alog10(radarr1*0.259), alog10(countarr1), color = colors.red

start = [153000.0]
err = dindgen(i) - dindgen(i) + 1.
result1 = mpfitfun('rcubed',radarr1,countarr1, err, start)

oplot, alog10(radarr1*0.259),  alog10(result1(0)*((radarr1)^(-3.))),color = colors.blue



device, /close
set_plot, mydevice

END
