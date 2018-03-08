PRO colorplot
close,/all



;read in the V radial profile
OPENR, lun,'/Users/jkrick/umich/icl/A3984/iclVtab', /GET_LUN  ;iclVtab2
rows=18
radius = FLTARR(rows)
counts = FLTARR(rows)
Verr = FLTARR(rows)
stop = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0
s = 0
i = 0
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      
          radius(i) = r
          counts(i) = c
          Verr(i) = ellerr 
          stop(i) = s
          i = i +1
     
ENDFOR

close, lun
free_lun, lun

counts = counts[0:i - 1] 
radius = radius[0:i-1]



;iclrtab2 vs samemask/iclVtab2.2
interpcounts = fltarr(29)
interpradius = fltarr(29)
interpradius = interpolate(radius, [0,1,2,3,3.75,4,4.66,5,6,7,7.66,8,9,9.66,10,10.66,11,11.66,12,12.66,13,13.66,14,14.66,15,15.66,16,16.66,17])
interpcounts = interpolate(counts, [0,1,2,3,3.75,4,4.66,5,6,7,7.66,8,9,9.66,10,10.66,11,11.66,12,12.66,13,13.66,14,14.66,15,15.66,16,16.66,17])
print, interpradius





;read in the r sb profile
OPENR, lun,'/Users/jkrick/umich/icl/A3984/iclrtab', /GET_LUN 
rows= 17 
rradius = FLTARR(rows)
rcounts = FLTARR(rows)
rerr = FLTARR(rows)
r = 0.0				;radius
c = 0.0				;mean counts
ellerr = 0.0

FOR j=0,rows-1 DO BEGIN
      READF, lun, r,c, ellerr,s
      rradius(j) = r
      rcounts(j) = c
      rerr(j) = ellerr
ENDFOR

close, lun
free_lun, lun

rcounts = rcounts 
rradius = rradius


;iclrtab2 vs samemask/iclVtab2.2
interprcounts = fltarr(15)
interprradius = fltarr(15)
interprradius = interpolate(rradius, [0,0.33,0.66,1,2,2.33,3,3.11,3.44,3.77,4,4.16,4.66,5,5.33,6,6.33,7,7.33,8,8.33,9,9.33,10,10.33,11,11.33,12,12.33])
interprcounts = interpolate(rcounts, [0,0.33,0.66,1,2,2.33,3,3.11,3.44,3.77,4,4.16,4.66,5,5.33,6,6.33,7,7.33,8,8.33,9,9.33,10,10.33,11,11.33,12,12.33])

print, interprradius

Vyarr = fltarr(7)
ryarr = fltarr(7)
radarr = fltarr(7)

Vyarr(0) = mean(interpcounts(0:3))
Vyarr(1) = mean(interpcounts(4:7))
Vyarr(2) = mean(interpcounts(8:11))
Vyarr(3) = mean(interpcounts(12:15))
Vyarr(4) = mean(interpcounts(16:19))
Vyarr(5) = mean(interpcounts(20:23))
Vyarr(6) = mean(interpcounts(24:28))




radarr(0) = mean(interpradius(0:3))
radarr(1) = mean(interpradius(4:7))
radarr(2) = mean(interpradius(8:10))
radarr(3) = mean(interpradius(11:14))
radarr(4) = mean(interpradius(16:19))
radarr(5) = mean(interpradius(20:23))
radarr(6) = mean(interpradius(24:28))

ryarr(0) = mean(interprcounts(0:3))
ryarr(1) = mean(interprcounts(4:7))
ryarr(2) = mean(interprcounts(8:10))
ryarr(3) = mean(interprcounts(11:14))
ryarr(4) = mean(interprcounts(16:19))
ryarr(5) = mean(interprcounts(20:23))
ryarr(6) = mean(interprcounts(24:28))

;print, Vyarr, ryarr

Vyarr = 24.3 -2.5*(alog10(Vyarr/(0.259^2)))
ryarr = 24.6 -2.5*(alog10(ryarr/(0.259^2)))
;print, Vyarr, ryarr


f0v = 3.75E-9
f0r = 1.75E-9
noiser = fltarr(3)
noiser= [29.1,29.1,29.1, 29.1,29.1,29.1,29.1]
noisev = noiser -0.1


fv = f0v*10^(Vyarr/(-2.5))
sigv = f0v*10^(noisev/(-2.5))
fr = f0r*10^(ryarr/(-2.5))
sigr = f0r*10^(noiser/(-2.5))


;because I have averaged together 4(sometimes 3 values, the sigma is
;better than the individual sigma)

sigv = sqrt(1./4.)*sigv
;sigv(4) = sqrt(1./3.)*sigv(4)

sigr = sqrt(1./4.)*sigr
;sigr(1) = sqrt(1./3.)*sigr(1)
;sigr(2) = sqrt(1./3.)*sigr(2)
;sigr(3) = (1./2.)*sigr(3)

color = Vyarr - ryarr
;print, "color", color

sig = sqrt(((fv/fr)^2)*(((sigv/fv)^2)+((sigr/fr)^2)))
percent = sig/(fv/fr)
plus = -2.5*alog10(1+percent)
minus = -2.5*alog10(1-percent)

;print, color;, plus, minus
;print, plus
;print, minus

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/Users/jkrick/umich/icl/A3984/color.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot, radarr*0.259, color, psym = 2, xrange = [0,250], ytitle ="V-r",thick = 3,$
charthick = 3, xthick = 3, ythick = 3,xtitle = 'Semi-major Axis (arcseconds)',$
yrange=[-1.2,1.0],ystyle = 1, xstyle = 8

axis, 0, 1.0, xaxis=1, xrange=[0,1.06], xstyle = 1, xthick = 3,charthick = 3

;upper = [0.2,0.2,0.2,0.2]
errplot, radarr*0.259, color-minus, color-plus

;fit this with a linear function
err = minus
start = [-0.0025,0.6]
result = MPFITFUN('linear',radarr,color,err, start,perror =test)
print, color, plus, minus
oplot, findgen(1200)*0.259, result(0)*findgen(1200) + result(1), thick = 3, linestyle = 2
;oplot, findgen(1200)*0.259, (result(0) -2*test(0))*findgen(1200) + (result(1) - 2*test(1)), linestyle = 2
;oplot, findgen(1200)*0.259, (result(0)+2*test(0) )*findgen(1200) + (result(1) + 2*test(1)), linestyle = 2

out = textoidl('h_{70}^{-1}kpc')


;the red cluster sequence
x= findgen(100) + 205
y = x - (findgen(100)+205) + 0.3
y2 = x - (findgen(100)+205) + 0.0

oplot, x, y, thick = 3
;oplot, x, y2, thick = 3
xyouts,   205, 0.22, 'RCS', charthick = 3
;xyouts, 230, 0.02, 'Tidal Features', charthick = 3


xyouts, 20,0.8,"A3984", charthick  = 3


sum = total(color)
av = sum / n_elements(color)

sig2 = sig^2
sumsig = sqrt(total(sig2))

sigav = av*(sqrt(sumsig^2/sum^2))


print, "sum, av, sumsig, sigav", sum, av, sumsig, sigav



device, /close
set_plot, mydevice


END
