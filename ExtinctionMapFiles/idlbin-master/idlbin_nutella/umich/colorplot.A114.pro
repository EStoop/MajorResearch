PRO colorplot
close,/all



;read in the V radial profile
OPENR, lun,'/Users/jkrick/umich/icl/A114/icl2Vtab', /GET_LUN  ;iclVtab2
rows=10
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

counts = counts[0:i -1 ] 
radius = radius[0:i-1]

;print, radius

;iclrtab2 vs samemask/iclVtab2.2
interpcounts = fltarr(11)
interpradius = fltarr(11)
interpradius = interpolate(radius, [0,1,2,3,4,5,6,7,8,8.3,9])
interpcounts = interpolate(counts, [0,1,2,3,4,5,6,7,8,8.3,9])

print, "interpradius", interpradius





;read in the r sb profile
OPENR, lun,'/Users/jkrick/umich/icl/A114/icl2rtab', /GET_LUN 
rows= 18
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

rcounts = rcounts[0:j-4]; - 2
rradius = rradius[0:j-4]

;print, rradius
;iclrtab2 vs samemask/iclVtab2.2
interprcounts = fltarr(11)
interprradius = fltarr(11)
interprradius = interpolate(rradius, [1,2,3,4,5,6,7,8,8.67,9,10])
interprcounts = interpolate(rcounts, [1,2,3,4,5,6,7,8,8.67,9,10])

print, "interprradius", interprradius

Vyarr = fltarr(3)
ryarr = fltarr(3)
radarr = fltarr(3)

Vyarr(0) = mean(interpcounts(0:2))
Vyarr(1) = mean(interpcounts(3:6))
Vyarr(2) = mean(interpcounts(7:10))


radarr(0) = mean(interpradius(0:2))
radarr(1) = mean(interpradius(3:6))
radarr(2) = mean(interpradius(7:10))

ryarr(0) = mean(interprcounts(0:2))
ryarr(1) = mean(interprcounts(3:6))
ryarr(2) = mean(interprcounts(7:10))

;print, Vyarr, ryarr

Vyarr = 24.3 -2.5*(alog10(Vyarr/(0.259^2)))
ryarr = 24.6 -2.5*(alog10(ryarr/(0.259^2)))
;print, Vyarr, ryarr


f0v = 3.75E-9
;f0v = 6.40E-9 ;= f0B
f0r = 1.75E-9
noiser = fltarr(3)
noiser= [27.7,27.7,27.7]
noisev = noiser + 1.7


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
print, "color", color

sig = sqrt(((fv/fr)^2)*(((sigv/fv)^2)+((sigr/fr)^2)))
percent = sig/(fv/fr)
FOR gah =0, 2 DO BEGIN
    IF percent[gah] GT 1 THEN percent[gah] =0.99
ENDFOR
plus = -2.5*alog10(1+percent)
minus = -2.5*alog10(1-percent)

;print, color;, plus, minus
;print, plus
;print, minus

sum = total(color)
av = sum / n_elements(color)

sig2 = sig^2
sumsig = sqrt(total(sig2))

sigav = av*(sqrt(sumsig^2/sum^2))


print, "sum, av, sumsig, sigav", sum, av, sumsig, sigav

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/Users/jkrick/umich/icl/A114/color.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot, radarr*0.259, color, psym = 2, xrange = [0,150], ytitle ="V-r",thick = 3,$
charthick = 3, xthick = 3, ythick = 3,xtitle = 'Semi-major Axis (arcseconds)',$
yrange=[-0.5,2.0],ystyle = 1, xstyle = 8

axis, 0, 2.0, xaxis=1, xrange=[0,1.158], xstyle = 1, xthick = 3,charthick = 3

;upper = [0.2,0.2,0.2,0.2]
errplot, radarr*0.259, color-minus, color-plus

;fit this with a linear function
err = plus
start = [-0.0025,0.6]
result = MPFITFUN('linear',radarr,color,err, start,perror =test)

oplot, findgen(1200)*0.259, result(0)*findgen(1200) + result(1), thick = 3, linestyle = 2
oplot, findgen(1200)*0.259, (result(0) -2*test(0))*findgen(1200) + (result(1) - 2*test(1)), linestyle = 2
oplot, findgen(1200)*0.259, (result(0)+2*test(0) )*findgen(1200) + (result(1) + 2*test(1)), linestyle = 2

out = textoidl('h_{70}^{-1}kpc')


;the red cluster sequence
x= findgen(100) + 130
y = x - (findgen(100)+130) + 0.5
;y2 = x - (findgen(100)+135) + 0.0

oplot, x, y, thick = 3
;oplot, x, y2, thick = 3
xyouts, 130, 0.4, 'RCS', charthick = 3
;xyouts, 230, 0.02, 'Tidal Features', charthick = 3
xyouts, 10, 1.8, "AC114", charthick=3
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;read in the V radial profile
OPENR, lun,'/Users/jkrick/umich/icl/A2556/galBtab', /GET_LUN  ;iclVtab2
rows=14
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

counts = counts[3:i -1 ] 
radius = radius[3:i-1]

print, radius

;read in the r sb profile
OPENR, lun,'/Users/jkrick/umich/icl/A2556/galrtab', /GET_LUN 
rows= 14
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

rcounts = rcounts[0:j-4]; - 2
rradius = rradius[0:j-4]

print, rradius

Vyarr = fltarr(3)
ryarr = fltarr(3)
radarr = fltarr(3)

Vyarr(0) = mean(counts(0:2))
Vyarr(1) = mean(counts(3:6))
Vyarr(2) = mean(counts(7:10))


radarr(0) = mean(radius(0:2))
radarr(1) = mean(radius(3:6))
radarr(2) = mean(radius(7:10))

ryarr(0) = mean(rcounts(0:2))
ryarr(1) = mean(rcounts(3:6))
ryarr(2) = mean(rcounts(7:10))

;print, Vyarr, ryarr

Vyarr = 22.19 -2.5*(alog10(Vyarr/(0.435^2)))
ryarr = 22.04 -2.5*(alog10(ryarr/(0.435^2)))
;print, Vyarr, ryarr


;f0v = 3.75E-9
f0v = 6.40E-9 ;= f0B
f0r = 1.75E-9
noiser = fltarr(3)
noiser= [27.7,27.7,27.7]
noisev = noiser + 1.7


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
print, "color", color

sig = sqrt(((fv/fr)^2)*(((sigv/fv)^2)+((sigr/fr)^2)))
percent = sig/(fv/fr)
plus = -2.5*alog10(1+percent)
minus = -2.5*alog10(1-percent)

;oplot, radarr*0.435, color, psym = 4,thick = 3



device, /close
set_plot, mydevice


END
