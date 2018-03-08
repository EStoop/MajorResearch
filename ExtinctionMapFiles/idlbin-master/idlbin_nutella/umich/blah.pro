PRO blah
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/final/junk.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color


radarr = findgen(100) + 1.
newarr = 25000.0*(radarr^(-3.))
plot, alog10(radarr), alog10(newarr), thick = 3

start = [9600.0]
err = dindgen(100) - dindgen(100) + 1.
result1 = mpfitfun('rcubed',radarr,newarr, err, start)

oplot, alog10(radarr), alog10(result1(0)*(radarr^(-3.0)))
device, /close
set_plot, mydevice

END
