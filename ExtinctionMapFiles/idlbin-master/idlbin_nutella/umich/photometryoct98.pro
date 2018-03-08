PRO photometry

close, /all
colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva4/jkrick/oct98/phot.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

openr, lun, "/n/Godiva4/jkrick/oct98/oct20/phot.log", /get_lun

airmassarr = fltarr(20)
mag1arr = fltarr(20)
mag2arr = fltarr(20)
mag3arr = fltarr(20)
i = 0

WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, airmass, mag1, mag2, mag3
    airmassarr[i] = airmass
    mag1arr[i] = mag1
    mag2arr[i] = mag2
    mag3arr[i] = mag3
    i = i + 1
ENDWHILE

airmassarr = airmassarr[0:i-1]
mag1arr = mag1arr[0:i-1] - mag1arr[fix(i/2)]
mag2arr = mag2arr[0:i-1]- mag2arr[fix(i/2)]
mag3arr = mag3arr[0:i-1]- mag3arr[fix(i/2)]

close, lun
free_lun, lun


plot, airmassarr, mag1arr , thick = 3,charthick = 3, xthick = 3, ythick = 3,$
title = "r band  photometric?", yrange=[-0.1,0.1], xrange = [1.0,1.6],psym = 2
oplot, airmassarr, mag1arr, thick = 3

oplot, airmassarr, mag2arr, thick = 3, psym = 2
oplot, airmassarr, mag2arr, thick = 3
oplot, airmassarr, mag3arr, thick = 3,psym = 2
oplot, airmassarr, mag3arr, thick = 3

airmassarr20 = [airmassarr,airmassarr,airmassarr]
magarr20 = [mag1arr, mag2arr, mag3arr]
err = dindgen(i) - dindgen(i) + 1
start = [0.06,-0.1]
result20 = MPFITFUN('linear',airmassarr20,magarr20,err, start)


;#######################################
openr, lun, "/n/Godiva4/jkrick/oct98/oct21/phot.log", /get_lun

airmassarr = fltarr(20)
mag1arr = fltarr(20)
mag2arr = fltarr(20)
mag3arr = fltarr(20)
i = 0

WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, airmass, mag1, mag2, mag3
    airmassarr[i] = airmass
    mag1arr[i] = mag1
    mag2arr[i] = mag2
    mag3arr[i] = mag3
    i = i + 1
ENDWHILE

airmassarr = airmassarr[0:i-1]
mag1arr = mag1arr[0:i-1] - mag1arr[fix(i-2/2)]
mag2arr = mag2arr[0:i-1]- mag2arr[fix(i-2/2)]
mag3arr = mag3arr[0:i-1]- mag3arr[fix(i-2/2)]

close, lun
free_lun, lun


oplot, airmassarr, mag1arr , psym = 2, thick = 3, color = colors.blue
oplot, airmassarr, mag1arr, thick = 3, color = colors.blue

oplot, airmassarr, mag2arr, thick = 3, color = colors.blue, psym = 2
oplot, airmassarr, mag2arr, thick = 3, color = colors.blue
oplot, airmassarr, mag3arr, thick = 3, color = colors.blue,psym = 2
oplot, airmassarr, mag3arr, thick = 3, color = colors.blue

airmassarr21 = [airmassarr,airmassarr,airmassarr]
magarr21 = [mag1arr, mag2arr, mag3arr]
err = dindgen(i) - dindgen(i) + 1
start = [0.06,-0.1]
result21 = MPFITFUN('linear',airmassarr21,magarr21,err, start)

;##########################################

openr, lun, "/n/Godiva4/jkrick/oct98/oct24/phot.log", /get_lun

airmassarr = fltarr(20)
mag1arr = fltarr(20)
mag2arr = fltarr(20)
mag3arr = fltarr(20)
i = 0

WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, airmass, mag1, mag2, mag3
    airmassarr[i] = airmass
    mag1arr[i] = mag1
    mag2arr[i] = mag2
    mag3arr[i] = mag3
    i = i + 1
ENDWHILE

airmassarr = airmassarr[0:i-1]
mag1arr = mag1arr[0:i-1] - mag1arr[fix(i/2)]
mag2arr = mag2arr[0:i-1]- mag2arr[fix(i/2)]
mag3arr = mag3arr[0:i-1]- mag3arr[fix(i/2)]

close, lun
free_lun, lun


oplot, airmassarr, mag1arr , psym = 2, thick = 3, color = colors.red
oplot, airmassarr, mag1arr, thick = 3, color = colors.red

oplot, airmassarr, mag2arr, thick = 3, color = colors.red, psym = 2
oplot, airmassarr, mag2arr, thick = 3, color = colors.red
oplot, airmassarr, mag3arr, thick = 3, color = colors.red,psym = 2
oplot, airmassarr, mag3arr, thick = 3, color = colors.red

airmassarr24 = [airmassarr,airmassarr,airmassarr]
magarr24 = [mag1arr, mag2arr, mag3arr]
err = dindgen(i) - dindgen(i) + 1
start = [0.06,-0.1]
result24= MPFITFUN('linear',airmassarr24,magarr24,err, start)
;##########################################
openr, lun, "/n/Godiva4/jkrick/oct98/oct25/phot.log", /get_lun

airmassarr = fltarr(20)
mag1arr = fltarr(20)
mag2arr = fltarr(20)
mag3arr = fltarr(20)
i = 0

WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, airmass, mag1, mag2, mag3
    airmassarr[i] = airmass
    mag1arr[i] = mag1
    mag2arr[i] = mag2
    mag3arr[i] = mag3
    i = i + 1
ENDWHILE

airmassarr = airmassarr[0:i-1]
mag1arr = mag1arr[0:i-1] - mag1arr[fix(i/2)]
mag2arr = mag2arr[0:i-1]- mag2arr[fix(i/2)]
mag3arr = mag3arr[0:i-1]- mag3arr[fix(i/2)]

close, lun
free_lun, lun


oplot, airmassarr, mag1arr , psym = 2, thick = 3, color = colors.green
oplot, airmassarr, mag1arr, thick = 3, color = colors.green

oplot, airmassarr, mag2arr, thick = 3, color = colors.green, psym = 2
oplot, airmassarr, mag2arr, thick = 3, color = colors.green
;oplot, airmassarr, mag3arr, thick = 3, color = colors.green,psym = 2
;oplot, airmassarr, mag3arr, thick = 3, color = colors.green
airmassarr25 = [airmassarr,airmassarr,airmassarr]
magarr25 = [mag1arr, mag2arr, mag3arr]
err = dindgen(i) - dindgen(i) + 1
start = [0.06,-0.1]
result25 = MPFITFUN('linear',airmassarr25,magarr25,err, start)
;##########################################

xyouts, 1.05, 0.09, "oct20", color = colors.black, charthick = 3
xyouts, 1.05, 0.08, "oct21", color = colors.blue, charthick = 3
xyouts, 1.05, 0.07, "oct24", color = colors.red, charthick = 3
xyouts, 1.05, 0.06, "oct25", color = colors.green, charthick = 3

plot, airmassarr20, result20(0)*airmassarr20 + result20(1), color = colors.black, thick = 3,$
charthick = 3, xthick = 3, ythick = 3,$
title = "r band  photometric?", yrange=[-0.1,0.1], xrange = [1.0,1.6]

oplot, airmassarr21, result21(0)*airmassarr21 + result21(1), color = colors.blue, thick = 3
oplot, airmassarr24, result24(0)*airmassarr24 + result24(1), color = colors.red, thick = 3
oplot, airmassarr25, result25(0)*airmassarr25 + result25(1), color = colors.green, thick = 3

xyouts, 1.05, 0.09, "oct20"+ string(result20(0)), color = colors.black, charthick = 3
xyouts, 1.05, 0.08, "oct21"+ string(result21(0)), color = colors.blue, charthick = 3
xyouts, 1.05, 0.07, "oct24"+ string(result24(0)), color = colors.red, charthick = 3
xyouts, 1.05, 0.06, "oct25"+ string(result25(0)), color = colors.green, charthick = 3


;##########################################
openr, lun, "/n/Godiva4/jkrick/oct98/oct22/phot.log", /get_lun
;B
airmassarr = fltarr(20)
mag1arr = fltarr(20)
mag2arr = fltarr(20)
mag3arr = fltarr(20)
i = 0

WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, airmass, mag1, mag2, mag3
    airmassarr[i] = airmass
    mag1arr[i] = mag1
    mag2arr[i] = mag2
    mag3arr[i] = mag3
    i = i + 1
ENDWHILE

airmassarr = airmassarr[0:i-1]
mag1arr = mag1arr[0:i-1] - mag1arr[fix(i/2)]
mag2arr = mag2arr[0:i-1]- mag2arr[fix(i/2)]
mag3arr = mag3arr[0:i-1]- mag3arr[fix(i/2)]

close, lun
free_lun, lun


plot, airmassarr, mag1arr , psym = 2, thick = 3, color = colors.magenta,$
yrange=[-0.1, 0.1], charthick=3, xthick=3,ythick=3, title = "B photometric?"
oplot, airmassarr, mag1arr, thick = 3, color = colors.magenta

oplot, airmassarr, mag2arr, thick = 3, color = colors.magenta, psym = 2
oplot, airmassarr, mag2arr, thick = 3, color = colors.magenta
;oplot, airmassarr, mag3arr, thick = 3, color = colors.magenta,psym = 2
;oplot, airmassarr, mag3arr, thick = 3, color = colors.magenta
airmassarr22 = [airmassarr,airmassarr,airmassarr]
magarr22 = [mag1arr, mag2arr, mag3arr]
err = dindgen(i) - dindgen(i) + 1
start = [0.06,-0.1]
result22 = MPFITFUN('linear',airmassarr22,magarr22,err, start)
;##########################################
device, /close
set_plot, mydevice


END
