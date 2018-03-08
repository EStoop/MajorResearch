PRO colorplot
close,/all



;read in the V radial profile
OPENR, lun,'/Users/jkrick/umich/icl/A4010/icl2Btab', /GET_LUN  ;iclVtab2
rows=32
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
interpcounts = fltarr(25)
interpradius = fltarr(25)
interpradius = interpolate(radius, [0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24])
interpcounts = interpolate(counts, [0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24])
print, "B radius", interpradius





;read in the r sb profile
OPENR, lun,'/Users/jkrick/umich/icl/A4010/icl2rtab', /GET_LUN 
rows= 23
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

rcounts = rcounts[0:j-1]; - 2
rradius = rradius[0:j-1]

;print, "rradius", rradius
;iclrtab2 vs samemask/iclVtab2.2
interprcounts = fltarr(25)
interprradius = fltarr(25)
interprradius = interpolate(rradius, [0,1,2,2.67,3,4,4.5,5,6,7,8.5,9.5,10.5,12.5,13.5,14.5,15.25,15.75,16.5,17.5,18.5,19.5,20.5,21.5,22.5])
interprcounts = interpolate(rcounts, [0,1,2,2.67,3,4,4.5,5,6,7,8.5,9.5,10.5,12.5,13.5,14.5,15.25,15.75,16.5,17.5,18.5,19.5,20.5,21.5,22.5])

print, "r radius" , interprradius

Vyarr = fltarr(5)
ryarr = fltarr(5)
radarr = fltarr(5)

Vyarr(0) = mean(interpcounts(0:4))
Vyarr(1) = mean(interpcounts(5:9))
Vyarr(2) = mean(interpcounts(10:14))
Vyarr(3) = mean(interpcounts(15:19))
Vyarr(4) = mean(interpcounts(20:24))

radarr(0) = mean(interpradius(0:4))
radarr(1) = mean(interpradius(5:9))
radarr(2) = mean(interpradius(10:14))
radarr(3) = mean(interpradius(15:19))
radarr(4) = mean(interpradius(20:24))

ryarr(0) = mean(interprcounts(0:4))
ryarr(1) = mean(interprcounts(5:9))
ryarr(2) = mean(interprcounts(10:14))
ryarr(3) = mean(interprcounts(15:19))
ryarr(4) = mean(interprcounts(20:24))


;print, Vyarr, ryarr

Vyarr = 22.19 -2.5*(alog10(Vyarr/(0.435^2)))
ryarr = 22.04 -2.5*(alog10(ryarr/(0.435^2)))
;print, Vyarr, ryarr


;f0v = 3.75E-9
f0v = 6.40E-9 ;= f0B
f0r = 1.75E-9
noiser = fltarr(5)
noiser= [28.4,28.4,28.4,28.4,28.4]
noisev = noiser + 1.5


fv = f0v*10^(Vyarr/(-2.5))
sigv = f0v*10^(noisev/(-2.5))
fr = f0r*10^(ryarr/(-2.5))
sigr = f0r*10^(noiser/(-2.5))


;because I have averaged together 4(sometimes 3 values, the sigma is
;better than the individual sigma)

sigv = sqrt(1./5.)*sigv
;sigv(4) = sqrt(1./3.)*sigv(4)

sigr = sqrt(1./5.)*sigr
;sigr(1) = sqrt(1./3.)*sigr(1)
;sigr(2) = sqrt(1./3.)*sigr(2)
;sigr(3) = (1./2.)*sigr(3)

color = Vyarr - ryarr
print, "color", color

sig = sqrt(((fv/fr)^2)*(((sigv/fv)^2)+((sigr/fr)^2)))
percent = sig/(fv/fr)
plus = -2.5*alog10(1+percent)
minus = -2.5*alog10(1-percent)

print, "color, plus, minus",color, plus, minus
;print, plus
;print, minus
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/Users/jkrick/umich/icl/A4010/color.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot, radarr*0.435, color, psym = 2, xrange = [0,400], ytitle ="B-r",thick = 3,$
charthick = 3, xthick = 3, ythick = 3,xtitle = 'Semi-major Axis (arcseconds)',$
yrange=[1.0,3.0],ystyle = 1, xstyle = 8

axis, 0, 3.0, xaxis=1, xrange=[0,0.85], xstyle = 1, xthick = 3,charthick = 3

;upper = [0.2,0.2,0.2,0.2]
errplot, radarr*0.435, color-minus, color-plus

;fit this with a linear function
err = plus
start = [-0.0025,0.6]
result = MPFITFUN('linear',radarr,color,err, start,perror =test)

oplot, findgen(1200)*0.435, result(0)*findgen(1200) + result(1), thick = 3, $
  linestyle = 2
oplot, findgen(1200)*0.435, (result(0) -2*test(0))*findgen(1200) + (result(1) - 2*test(1)), linestyle = 2
oplot, findgen(1200)*0.435, (result(0)+2*test(0) )*findgen(1200) + (result(1) + 2*test(1)), linestyle = 2

out = textoidl('h_{70}^{-1}kpc')


;the red cluster sequence
x= findgen(100) + 360
y = x - (findgen(100)+360) + 1.8
;y2 = x - (findgen(100)+135) + 0.0

oplot, x, y, thick = 3
;oplot, x, y2, thick = 3
xyouts, 360, 1.70, 'RCS', charthick = 3
;xyouts, 230, 0.02, 'Tidal Features', charthick = 3
xyouts, 30,2.7,"A4010", charthick = 3
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
noiser= [28.4,28.4,28.4]
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


sum = total(color)
av = sum / n_elements(color)

sig2 = sig^2
sumsig = sqrt(total(sig2))

sigav = av*(sqrt(sumsig^2/sum^2))


print, "sum, av, sumsig, sigav", sum, av, sumsig, sigav

device, /close
set_plot, mydevice


END