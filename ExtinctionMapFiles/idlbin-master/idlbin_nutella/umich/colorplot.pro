PRO colorplot
close,/all



;read in the V radial profile
OPENR, lun,'/Users/jkrick/umich/icl/A3888/iclVtab2.2', /GET_LUN  ;iclVtab2
rows=26
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
      IF c GT 2.0 THEN BEGIN
          radius(i) = r
          counts(i) = c
          Verr(i) = ellerr 
          stop(i) = s
          i = i +1
      ENDIF
ENDFOR

close, lun
free_lun, lun

counts = counts[0:i - 1] -  2
radius = radius[0:i-1]

;iclrtab2 vs. iclVtab2
;interpradius = fltarr(35)
;interpcounts = fltarr(35)
;interpradius = interpolate(radius, [0,0.5,1,1.33,2,3,4,5,5.33,6,6.66,7,8,8.33,9,9.33,10,10.33,11,11.33,12,12.66,13,14,15,15.2,16,16.2,17,17.2,18,18.2,19,19.2,20])
;interpradius = interpradius[1:34]
;interpcounts = interpolate(counts, [0,0.5,1,1.33,2,3,4,5,5.33,6,6.66,7,8,8.33,9,9.33,10,10.33,11,11.33,12,12.66,13,14,15,15.2,16,16.2,17,17.2,18,18.2,19,19.2,20])
;interpcounts = interpcounts[1:34]

;iclVtab3 vs iclrtab2
;interpradius = fltarr(33)
;interpcounts = fltarr(33)
;interpradius = interpolate(radius, [0,0.33,1,2,3,4,4.33,5,5.66,6,7,7.33,8,8.33,9,9.33,10,10.33,11,11.5,12,12.6,13,13.6,14,14.6,15,15.5,16,16.4,17,17.4,18])
;interpcounts = interpolate(counts, [0,0.33,1,2,3,4,4.33,5,5.66,6,7,7.33,8,8.33,9,9.33,10,10.33,11,11.5,12,12.6,13,13.6,14,14.6,15,15.5,16,16.4,17,17.4,18])

;iclrtab vs iclVtab2
;interpcounts = fltarr(35)
;interpradius = fltarr(35)
;interpradius = interpolate(radius, [0,0.25,1,1.66,2,3,4,5,5.33,6,6.66,7,8,8.33,9,9.33,10,10.33,11,11.33,12,12.66,13,14,15,15.2,16,16.2,17,17.2,18,18.2,19,19.2,20])
;interpcounts = interpolate(counts, [0,0.25,1,1.66,2,3,4,5,5.33,6,6.66,7,8,8.33,9,9.33,10,10.33,11,11.33,12,12.66,13,14,15,15.2,16,16.2,17,17.2,18,18.2,19,19.2,20])
;interpradius = interpradius[0:33]
;interpcounts = interpcounts[0:33]



;iclrtab2 vs samemask/iclVtab2.2
interpcounts = fltarr(37)
interpradius = fltarr(37)
interpradius = interpolate(radius, [0,1,1.33,2,3,4,5,5.33,6,6.66,7,8,8.33,9,9.33,10,10.33,11,11.33,12,12.66,13,14,14.33,15,16,17,17.66,18,19,19.33,20,20.75,21,22,23,23.25])
interpcounts = interpolate(counts, [0,1,1.23,2,3,4,5,5.33,6,6.66,7,8,8.33,9,9.33,10,10.33,11,11.33,12,12.66,13,14,14.33,15,16,17,17.66,18,19,19.33,20,20.75,21,22,23,23.25])
print, interpradius





;read in the r sb profile
OPENR, lun,'/Users/jkrick/umich/icl/A3888/iclrtab2', /GET_LUN 
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

rcounts = rcounts - 2;[0:j-3] - 2
rradius = rradius;[0:j-3]

;iclrtab2 vs iclVtab2
;interprcounts = fltarr(35)
;interprradius = fltarr(35)
;interprradius = interpolate(rradius, [0,0.66,1,1.4,2,3,3.75,4,4.5,5,5.2,5.8,6,6.66,7,7.66,8,8.66,9,9.5,10.,10.2,11,11.8,12,12.8,13,13.8,14,14.8,15,15.8,16,16.8,17])
;interprradius = interprradius[0:33]
;interprcounts = interpolate(rcounts, [0,0.66,1,1.4,2,3,3.75,4,4.5,5,5.2,5.8,6,6.66,7,7.66,8,8.66,9,9.5,10.,10.2,11,11.8,12,12.8,13,13.8,14,14.8,15,15.8,16,16.8,17])
;interprcounts = interprcounts[0:33]

;iclVtab3 vs iclrtab2
;interprcounts = fltarr(35)
;interprradius = fltarr(35)
;interprradius = interpolate(rradius, [0,0.66,1,1.4,2,3,3.75,4,4.5,5,5.2,5.8,6,6.66,7,7.66,8,8.66,9,9.5,10,10.4,11,11.4,12,12.4,13,13.4,14,14.6,15,15.6,16,16.6,17])
;interprcounts = interpolate(rcounts, [0,0.66,1,1.4,2,3,3.75,4,4.5,5,5.2,5.8,6,6.66,7,7.66,8,8.66,9,9.5,10,10.4,11,11.4,12,12.4,13,13.4,14,14.6,15,15.6,16,16.6,17])
;interprradius= interprradius[1:34]
;interprcounts= interprcounts[1:34]

;iclrtab vs iclVtab2
;interprcounts = fltarr(35)
;interprradius = fltarr(35)
;interprradius = interpolate(rradius, [0,1,2,2.6,3,3.25,4,5,5.75,6,6.5,7,7.2,7.8,8,8.66,9,9.66,10,10.66,11,11.5,12,12.2,13,13.8,14,14.8,15,15.8,16,16.8,17,17.8,18])
;interprcounts = interpolate(rcounts, [0,1,2,2.6,3,3.25,4,5,5.75,6,6.5,7,7.2,7.8,8,8.66,9,9.66,10,10.66,11,11.5,12,12.2,13,13.8,14,14.8,15,15.8,16,16.8,17,17.8,18])
;interprradius = interprradius[1:34]
;interprcounts= interprcounts[1:34]

;iclrtab2 vs samemask/iclVtab2.2
interprcounts = fltarr(37)
interprradius = fltarr(37)
interprradius = interpolate(rradius, [0,0.66,1,1.4,2,3,3.75,4,4.5,5,5.2,5.8,6,6.66,7,7.66,8,8.66,9,9.5,10,10.2,10.8,11,11.4,12,12.6,13,13.2,13.8,14,14.4,15,15.2,16,16.8,17])
interprcounts = interpolate(rcounts, [0,0.66,1,1.4,2,3,3.75,4,4.5,5,5.2,5.8,6,6.66,7,7.66,8,8.66,9,9.5,10,10.2,10.8,11,11.4,12,12.6,13,13.2,13.8,14,14.4,15,15.2,16,16.8,17])

print, interprradius

Vyarr = fltarr(5)
ryarr = fltarr(5)
radarr = fltarr(5)

;iclVtab2 original
;Vyarr(0) = mean(counts(0:3))
;Vyarr(1) = mean(counts(4:7))
;Vyarr(2) = mean(counts(8:11))
;Vyarr(3) = mean(counts(12:15))
;Vyarr(4) = mean(counts(15:17))

;iclVtab2 to match iclrtab
;Vyarr(0) = mean(counts(0:3))
;Vyarr(1) = mean(counts(4:7))
;Vyarr(2) = mean(counts(8:11))
;Vyarr(3) = mean(counts(12:15))
;Vyarr(4) = mean(counts(16:18))

;iclVtab3
;Vyarr(0) = mean(counts(0:2))
;Vyarr(1) = mean(counts(3:6))
;Vyarr(2) = mean(counts(7:10))
;Vyarr(3) = mean(counts(11:13))
;Vyarr(4) = mean(counts(14:15))

Vyarr(0) = mean(interpcounts(0:4))
Vyarr(1) = mean(interpcounts(5:10))
Vyarr(2) = mean(interpcounts(11:15))
Vyarr(3) = mean(interpcounts(16:21))
Vyarr(4) = mean(interpcounts(22:27))



radarr(0) = mean(interpradius(0:4))
radarr(1) = mean(interpradius(5:10))
radarr(2) = mean(interpradius(11:15))
radarr(3) = mean(interpradius(16:21))
radarr(4) = mean(interpradius(22:27))

ryarr(0) = mean(interprcounts(0:4))
ryarr(1) = mean(interprcounts(5:10))
ryarr(2) = mean(interprcounts(11:15))
ryarr(3) = mean(interprcounts(16:21))
ryarr(4) = mean(interprcounts(22:27))

;iclrtab2 original
;ryarr(0) = mean(rcounts(0:2))
;ryarr(1) = mean(rcounts(3:5))
;ryarr(2) = mean(rcounts(6:8))
;ryarr(3) = mean(rcounts(9:12))
;ryarr(4) = mean(rcounts(12:14))

;iclrtab2 to match iclVtab
;ryarr(0) = mean(rcounts(0:2))
;ryarr(1) = mean(rcounts(3:5))
;ryarr(2) = mean(rcounts(6:8))
;ryarr(3) = mean(rcounts(9:11))
;ryarr(4) = mean(rcounts(12:14))

;iclrtab to match iclVtab2
;ryarr(0) = mean(rcounts(0:4))
;ryarr(1) = mean(rcounts(5:7))
;ryarr(2) = mean(rcounts(8:11))
;ryarr(3) = mean(rcounts(12:14))
;ryarr(4) = mean(rcounts(15:17))

Vyarr = 24.3 -2.5*(alog10(Vyarr/(0.259^2)))
ryarr = 24.6 -2.5*(alog10(ryarr/(0.259^2)))
;print, Vyarr, ryarr


f0v = 3.75E-9
f0r = 1.75E-9
noiser = fltarr(5)
noiser= [28.8,28.8,28.8,28.8,28.8]
noisev = noiser + 0.7


fv = f0v*10^(Vyarr/(-2.5))
sigv = f0v*10^(noisev/(-2.5))
fr = f0r*10^(ryarr/(-2.5))
sigr = f0r*10^(noiser/(-2.5))


;because I have averaged together 4(sometimes 3 values, the sigma is
;better than the individual sigma)

sigv(0:4) = sqrt(1./6.)*sigv(0:4)
;sigv(4) = sqrt(1./3.)*sigv(4)

sigr(0:4) = sqrt(1./6.)*sigr(0:4)
;sigr(1) = sqrt(1./3.)*sigr(1)
;sigr(2) = sqrt(1./3.)*sigr(2)
;sigr(3) = (1./2.)*sigr(3)

color = Vyarr - ryarr


sig = sqrt(((fv/fr)^2)*(((sigv/fv)^2)+((sigr/fr)^2)))
percent = sig/(fv/fr)
plus = -2.5*alog10(1+percent)
minus = -2.5*alog10(1-percent)

print, color;, plus, minus
print, plus
print, minus
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/Users/jkrick/umich/icl/A3888/color.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot, radarr*0.259, color, psym = 2, xrange = [0,300], ytitle ="V-r",thick = 3,$
charthick = 3, xthick = 3, ythick = 3,xtitle = 'Semi-major Axis (arcseconds)',$
yrange=[-0.2,1.4],ystyle = 1, xstyle = 8
axis, 0, 1.4, xaxis=1, xrange=[0,1.034], xstyle = 1, xthick = 3,charthick = 3

upper = [0.2,0.2,0.2,0.2]
errplot, radarr*0.259, color-minus, color-plus

;fit this with a linear function
err = color+plus
start = [0.0025,0.]
result = MPFITFUN('linear',radarr,color,err, start,perror =test)

oplot, findgen(1200)*0.259, result(0)*findgen(1200) + result(1), thick = 3, linestyle = 2
;oplot, findgen(1200)*0.259, (result(0) -1*test(0))*findgen(1200) + (result(1) - 1*test(1)), linestyle = 2
;oplot, findgen(1200)*0.259, (result(0)+1*test(0) )*findgen(1200) + (result(1) + 1*test(1)), linestyle = 2

out = textoidl('h_{70}^{-1}kpc')

;xyouts, 28.9, 1.3, '100', charthick = 3, alignment = 0.5
;xyouts, 57.8, 1.3, '200', charthick = 3, alignment = 0.5
;xyouts, 86.7, 1.3, '300', charthick = 3, alignment = 0.5
;xyouts, 115.6, 1.3, '400', charthick = 3, alignment = 0.5
;xyouts, 144.5, 1.3, '500', charthick = 3, alignment = 0.5
;xyouts, 173.4, 1.3, '600', charthick = 3, alignment = 0.5
;xyouts, 202.3, 1.3, '700', charthick = 3, alignment = 0.5
;xyouts, 231.2, 1.3, '800', charthick = 3, alignment = 0.5
;xyouts, 270, 1.3, out, charthick = 3, alignment = 0.5

;the red cluster sequence
x= findgen(100) + 230
y = x - (findgen(100)+230) + 0.3
y2 = x - (findgen(100)+230) + 0.0

oplot, x, y, thick = 3
oplot, x, y2, thick = 3
xyouts, 230, 0.32, 'RCS', charthick = 3
xyouts, 230, 0.02, 'Tidal Features', charthick = 3

xyouts,30,1.2, 'A3888', charthick = 3
sum = total(color)
av = sum / n_elements(color)

sig2 = sig^2
sumsig = sqrt(total(sig2))

sigav = av*(sqrt(sumsig^2/sum^2))


print, "sum, av, sumsig, sigav", sum, av, sumsig, sigav



device, /close
set_plot, mydevice


END
;xarr = fltarr(4)
;xarr = [80., 110., 140., 170.]  ;arcseconds
;xarr = xarr/0.259    ;pixels

;Vyarr = 0.0189*exp(-xarr/230.9)
;Vyarr = 24.3 - (2.5*alog10(Vyarr / (0.259^2)))

;ryarr = 0.02099*exp(-xarr/317.634)
;ryarr = 24.6 - (2.5*alog10(ryarr / (0.259^2)))
