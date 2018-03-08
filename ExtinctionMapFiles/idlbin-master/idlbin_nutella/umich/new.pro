PRO new
close, /all
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/junk.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;r  curve
OPENR, lun,'/n/Godiva1/jkrick/A3888/final/iclrtab2', /GET_LUN 

;read in the radial profile
rows= 18-2
rradius2 = FLTARR(rows)
rcounts2 = FLTARR(rows)
rerr2 = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      rradius2(j) = r
      rcounts2(j) = c
      rerr2(j) = ellerr
ENDFOR

close, lun
free_lun, lun

rcounts2 = rcounts2 - 2;[0:j-3] - 2

;plot, rradius2* 0.259,24.6 - 2.5*alog10(rcounts2/(0.259^2)),thick = 3, yrange = [31,26]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
newradius = rradius2*0.259
newcounts = 24.6 - 2.5*alog10(rcounts2 / (0.259^2))
print, newradius, newcounts

start = [26.0, -400.]
err2 = dindgen(20) - dindgen(20) + 1

result2 = MPFITFUN('exponential',newradius,newcounts, rerr2, start)      ;ICL

plot, newradius, newcounts,thick = 3, yrange = [31,26]

oplot, newradius, result2(0)*exp(-newradius/result2(1)), thick = 3, linestyle = 2

device, /close
set_plot, mydevice

END
