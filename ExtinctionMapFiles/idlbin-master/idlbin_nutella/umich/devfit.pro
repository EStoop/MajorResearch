
pro devfit


device, true=24
device, decomposed=0
colors = GetColor(/Load, Start=1)
close, /all		;close all files = tidiness

OPENR, lun,'/n/Sheriff1/jkrick/sep00/A3888/ellipse/centerfit.txt', /GET_LUN

rows = 90
rad = FLTARR(rows)
intens = FLTARR(rows)
intens_err = FLTARR(rows)

r = 0.0
i = 0.0
e = 0.0

FOR j=0,rows-1 DO BEGIN
    READF, lun, r,i,e
    rad(j) = r
    intens(j) = i
    intens_err(j) = e
ENDFOR

close, lun
free_lun, lun

start = [0.008,1500]
err = dindgen(rows) - dindgen(rows) + 1
result = MPFITFUN('devauc',rad, intens, intens_err, start, maxiter = 300)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Sheriff1/jkrick/sep00/A3888/ellipse/devauc.ps', /portrait, $
                BITS=8, scale_factor=0.9 , /color


plot, rad,intens, psym = 4, thick=3
oploterr, rad, intens,intens_err

a = rad/result(1)
b = (result(0) *exp((-7.67)*((rad/result(1))^(1./4.)-1)))
print, rad,b
oplot, rad,b, color = colors.blue
;oplot, rad, result(0) +result(1)*rad + result(2)*rad^2, color = colors.blue
device, /close
set_plot, mydevice
END
