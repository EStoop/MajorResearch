pro grasil_convolv


COMMON FILTERCURVE,filters
;get_filts
filename = '/Users/jkrick/Virgo/grasil/elli1_5.txt'
readcol, filename, lambda, cirrus, star1, star2, mc, total, star3, sub1, sub2,format="A"


lam_um = lambda/1.e4 ;microns

;worry about units of flux, and extinction
;what about normalizing the filter curves

;total is in 1e30 erg/s/A
;need total in Lsun
; 1 Lsun = 3.839 × 10^33 erg/s

lum = total *1D30* lambda ; now in erg/s
lum = lum / 3.839D33 ; now in Lsun

ch1_flux= obsv_filter(lam_um,lum,/box,lrange=[3.0, 4.0])

print, 'test', ch1_flux  ; in Lsun/Hz?

central_lambda = 3.6 ; microns
c = 3E8 ; m/s
c = c * 1E6 ; um/s
central_nu = c / central_lambda

ch1_flux = ch1_flux * central_nu ; now in Lsun
ch1_flux = ch1_flux *3.839D33 ; now in erg/s
END
 


