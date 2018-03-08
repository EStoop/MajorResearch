PRO photometry

close, /all
;colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/Users/jkrick/palomar/LFC/phot/phot.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

openr, lun, "/Users/jkrick/palomar/LFC/phot/09aug04.phot", /get_lun

airmassarr = fltarr(20)
mag1arr = fltarr(20)
mag2arr = fltarr(20)
mag3arr = fltarr(20)
i = 0

WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, airmass, mag1, mag2, mag3
;    print, airmass, mag1, mag2, mag3
    airmassarr[i] = airmass
    mag1arr[i] = mag1
    mag2arr[i] = mag2
    mag3arr[i] = mag3
    i = i + 1
ENDWHILE

airmassarr = airmassarr[0:i-1]
mag1arr = mag1arr[0:i-1] ;- mag1arr[fix(i/2)]
mag2arr = mag2arr[0:i-1];- mag2arr[fix(i/2)]
mag3arr = mag3arr[0:i-1];- mag3arr[fix(i/2)]

;print, mag1arr

close, lun
free_lun, lun


plot, airmassarr, mag1arr , thick = 3,charthick = 3, xthick = 3, ythick = 3,$
      yrange=[21.4,22], title="aug 09 i band"

oplot, airmassarr, mag1arr, thick = 3, psym = 2

oplot, airmassarr, mag2arr, thick = 3, psym = 2
oplot, airmassarr, mag2arr, thick = 3
oplot, airmassarr, mag3arr, thick = 3,psym = 2
oplot, airmassarr, mag3arr, thick = 3

airmassarr20 = [airmassarr,airmassarr,airmassarr]
magarr20 = [mag1arr, mag2arr, mag3arr]
err = dindgen(i) - dindgen(i) + 1
start = [0.06,-0.1]
result20 = MPFITFUN('linear',airmassarr20,magarr20,err, start)
result = MPFITFUN('linear', airmassarr, mag3arr, err, start)

xyouts, 1.21, 21.95, string(result20(0))

oplot, airmassarr20, result20(0)*airmassarr20 + result20(1), thick = 3
oplot, airmassarr, result(0)*airmassarr + result(1), thick = 3, linestyle = 2


;##########################################
openr, lun, "/Users/jkrick/palomar/LFC/phot/10aug04.phot", /get_lun

airmassarr = fltarr(20)
mag1arr = fltarr(20)
mag2arr = fltarr(20)
mag3arr = fltarr(20)
i = 0

WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, airmass, mag1, mag2, mag3
;    print, airmass, mag1, mag2, mag3
    airmassarr[i] = airmass
    mag1arr[i] = mag1
    mag2arr[i] = mag2
    mag3arr[i] = mag3
    i = i + 1
ENDWHILE

airmassarr = airmassarr[0:i-1]
mag1arr = mag1arr[0:i-1] ;- mag1arr[fix(i/2)]
mag2arr = mag2arr[0:i-1];- mag2arr[fix(i/2)]
mag3arr = mag3arr[0:i-1];- mag3arr[fix(i/2)]

;print, mag1arr

close, lun
free_lun, lun


plot, airmassarr, mag1arr , thick = 3,charthick = 3, xthick = 3, ythick = 3,$
      yrange=[21.4,22], title="aug 10 i band"


oplot, airmassarr, mag1arr, thick = 3, psym = 2

oplot, airmassarr, mag2arr, thick = 3, psym = 2
oplot, airmassarr, mag2arr, thick = 3
oplot, airmassarr, mag3arr, thick = 3,psym = 2
oplot, airmassarr, mag3arr, thick = 3

airmassarr20 = [airmassarr,airmassarr,airmassarr]
magarr20 = [mag1arr, mag2arr, mag3arr]
err = dindgen(i) - dindgen(i) + 1
start = [0.06,-0.1]
result20 = MPFITFUN('linear',airmassarr20,magarr20,err, start)
result = MPFITFUN('linear', airmassarr, mag3arr, err, start)

xyouts, 1.5, 21.45, string(result20(0))

oplot, airmassarr20, result20(0)*airmassarr20 + result20(1), thick = 3
oplot, airmassarr, result(0)*airmassarr + result(1), thick = 3, linestyle = 2


;##########################################
openr, lun, "/Users/jkrick/palomar/LFC/phot/23jul04.phot", /get_lun

airmassarr = fltarr(20)
mag1arr = fltarr(20)
mag2arr = fltarr(20)
mag3arr = fltarr(20)
i = 0

WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, airmass, mag1, mag2, mag3
;    print, airmass, mag1, mag2, mag3
    airmassarr[i] = airmass
    mag1arr[i] = mag1
    mag2arr[i] = mag2
    mag3arr[i] = mag3
    i = i + 1
ENDWHILE

airmassarr = airmassarr[0:i-1]
mag1arr = mag1arr[0:i-1] ;- mag1arr[fix(i/2)]
mag2arr = mag2arr[0:i-1];- mag2arr[fix(i/2)]
mag3arr = mag3arr[0:i-1];- mag3arr[fix(i/2)]

;print, mag1arr

close, lun
free_lun, lun


plot, airmassarr, mag1arr , thick = 3,charthick = 3, xthick = 3, ythick = 3,$
      yrange=[22.4,23.4], title="jul 23 r band"


oplot, airmassarr, mag1arr, thick = 3, psym = 2

oplot, airmassarr, mag2arr, thick = 3, psym = 2
oplot, airmassarr, mag2arr, thick = 3
oplot, airmassarr, mag3arr, thick = 3,psym = 2
oplot, airmassarr, mag3arr, thick = 3

airmassarr20 = [airmassarr,airmassarr,airmassarr]
magarr20 = [mag1arr, mag2arr, mag3arr]
err = dindgen(i) - dindgen(i) + 1
start = [0.06,-0.1]
result20 = MPFITFUN('linear',airmassarr20,magarr20,err, start)
result = MPFITFUN('linear', airmassarr, mag3arr, err, start)

xyouts, 1.52, 23.35, string(result20(0))

oplot, airmassarr20, result20(0)*airmassarr20 + result20(1), thick = 3
oplot, airmassarr, result(0)*airmassarr + result(1), thick = 3, linestyle = 2

;##########################################
openr, lun, "/Users/jkrick/palomar/LFC/phot/24jul04.phot", /get_lun

airmassarr = fltarr(20)
mag1arr = fltarr(20)
mag2arr = fltarr(20)
mag3arr = fltarr(20)
i = 0

WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, airmass, mag1, mag2, mag3
;    print, airmass, mag1, mag2, mag3
    airmassarr[i] = airmass
    mag1arr[i] = mag1
    mag2arr[i] = mag2
    mag3arr[i] = mag3
    i = i + 1
ENDWHILE

airmassarr = airmassarr[0:i-1]
mag1arr = mag1arr[0:i-1] ;- mag1arr[fix(i/2)]
mag2arr = mag2arr[0:i-1];- mag2arr[fix(i/2)]
mag3arr = mag3arr[0:i-1];- mag3arr[fix(i/2)]

;print, mag1arr

close, lun
free_lun, lun


plot, airmassarr, mag1arr , thick = 3,charthick = 3, xthick = 3, ythick = 3,$
      yrange=[21.8,23.0], title="jul 24 g band"


oplot, airmassarr, mag1arr, thick = 3, psym = 2

oplot, airmassarr, mag2arr, thick = 3, psym = 2
oplot, airmassarr, mag2arr, thick = 3
oplot, airmassarr, mag3arr, thick = 3,psym = 2
oplot, airmassarr, mag3arr, thick = 3

airmassarr20 = [airmassarr,airmassarr,airmassarr]
magarr20 = [mag1arr, mag2arr, mag3arr]
err = dindgen(i) - dindgen(i) + 1
start = [0.06,-0.1]
result20 = MPFITFUN('linear',airmassarr20,magarr20,err, start)
result = MPFITFUN('linear', airmassarr, mag3arr, err, start)

xyouts, 1.25, 22.9, string(result20(0))

oplot, airmassarr20, result20(0)*airmassarr20 + result20(1), thick = 3
oplot, airmassarr, result(0)*airmassarr + result(1), thick = 3, linestyle = 2

;##########################################
openr, lun, "/Users/jkrick/palomar/LFC/phot/08aug04.phot", /get_lun

airmassarr = fltarr(20)
mag1arr = fltarr(20)
mag2arr = fltarr(20)
mag3arr = fltarr(20)
i = 0

WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, airmass, mag1, mag2, mag3
;    print, airmass, mag1, mag2, mag3
    airmassarr[i] = airmass
    mag1arr[i] = mag1
    mag2arr[i] = mag2
    mag3arr[i] = mag3
    i = i + 1
ENDWHILE

airmassarr = airmassarr[0:i-1]
mag1arr = mag1arr[0:i-1] ;- mag1arr[fix(i/2)]
mag2arr = mag2arr[0:i-1];- mag2arr[fix(i/2)]
mag3arr = mag3arr[0:i-1];- mag3arr[fix(i/2)]

;print, mag1arr

close, lun
free_lun, lun


plot, airmassarr, mag1arr , thick = 3,charthick = 3, xthick = 3, ythick = 3,$
      yrange=[19.8,21.4], title="aug 08 g band"


oplot, airmassarr, mag1arr, thick = 3, psym = 2

oplot, airmassarr, mag2arr, thick = 3, psym = 2
oplot, airmassarr, mag2arr, thick = 3
oplot, airmassarr, mag3arr, thick = 3,psym = 2
oplot, airmassarr, mag3arr, thick = 3

airmassarr20 = [airmassarr,airmassarr,airmassarr]
magarr20 = [mag1arr, mag2arr, mag3arr]
err = dindgen(i) - dindgen(i) + 1
start = [0.06,-0.1]
result20 = MPFITFUN('linear',airmassarr20,magarr20,err, start)
result = MPFITFUN('linear', airmassarr, mag3arr, err, start)

xyouts, 1.46, 21.4, string(result20(0))

oplot, airmassarr20, result20(0)*airmassarr20 + result20(1), thick = 3
oplot, airmassarr, result(0)*airmassarr + result(1), thick = 3, linestyle = 2


;##########################################
openr, lun, "/Users/jkrick/palomar/LFC/phot/22jul04.phot", /get_lun

airmassarr = fltarr(20)
mag1arr = fltarr(20)
mag2arr = fltarr(20)
mag3arr = fltarr(20)
i = 0

WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, airmass, mag1, mag2, mag3
;    print, airmass, mag1, mag2, mag3
    airmassarr[i] = airmass
    mag1arr[i] = mag1
    mag2arr[i] = mag2
    mag3arr[i] = mag3
    i = i + 1
ENDWHILE

airmassarr = airmassarr[0:i-1]
mag1arr = mag1arr[0:i-1] ;- mag1arr[fix(i/2)]
mag2arr = mag2arr[0:i-1];- mag2arr[fix(i/2)]
mag3arr = mag3arr[0:i-1];- mag3arr[fix(i/2)]

;print, mag1arr

close, lun
free_lun, lun


plot, airmassarr, mag1arr , thick = 3,charthick = 3, xthick = 3, ythick = 3,$
      yrange=[20.,21.5], title="jul 22 u band"


oplot, airmassarr, mag1arr, thick = 3, psym = 2

oplot, airmassarr, mag2arr, thick = 3, psym = 2
oplot, airmassarr, mag2arr, thick = 3
oplot, airmassarr, mag3arr, thick = 3,psym = 2
oplot, airmassarr, mag3arr, thick = 3

airmassarr20 = [airmassarr,airmassarr,airmassarr]
magarr20 = [mag1arr, mag2arr, mag3arr]
err = dindgen(i) - dindgen(i) + 1
start = [0.06,-0.1]
result20 = MPFITFUN('linear',airmassarr20,magarr20,err, start)
result = MPFITFUN('linear', airmassarr, mag3arr, err, start)

xyouts, 1.25, 21.2, string(result20(0))

oplot, airmassarr20, result20(0)*airmassarr20 + result20(1), thick = 3
oplot, airmassarr, result(0)*airmassarr + result(1), thick = 3, linestyle = 2

;##############

device, /close
set_plot, mydevice


END
