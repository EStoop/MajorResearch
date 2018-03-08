pro plotrcsratio
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)

!p.multi = [0, 1, 1]
ps_open, file = "/Users/jkrick/nep/clusters/rcsratio.ps", /portrait, /color,/square, xsize=4, ysize=4

;----------------------------------------------------------------
;plot summary of all 3 cluster measurements
;----------------------------------------------------------------
x = [0.1, 1.0]
;old
;nbright=[134.4,13.0,31.2]
;nbrighterr = [37.0,4.0,24.8]
;ndim = [469.0,24.6,1]
;ndimerr=[124.7,7.8,7]
;errors taken from stddev. of background counts.

;375kpc radius
nbright=[39,31.8]
nbrighterr = [sqrt(4.)*1.8,sqrt(3.)*1.4]
ndim = [41.,1.4]
ndimerr=[sqrt(4.)*3.5,sqrt(3.)*2.1]
y375 = ndim/(ndim+nbright)


;750kpc radius
;nbright750=[71.2,97.5]
;nbrighterr750=[sqrt(4.)*4.1,sqrt(3.)*6.7]
;ndim750=[66.7,4.5]
;ndimerr750=[sqrt(4.)*11.7,sqrt(3.)*13.9]
;y750 = ndim750/(ndim750+nbright750)

nbright750=[33.3,88.5]
nbrighterr750=[sqrt(4.)*5.8,sqrt(3.)*14.0]
ndim750=[25.0,0.1]
ndimerr750=[sqrt(4.)*13.4,sqrt(3.)*14.6]
y750 = ndim750/(ndim750+nbright750)

;500kpc radius
;nbright=[52.3,62.6]
;nbrighterr = [sqrt(4.)*2.7,sqrt(3.)*5.1]
;ndim = [52.2,18.2]
;ndimerr=[sqrt(4.)*5.8,sqrt(3.)*7.2]
;y500 = ndim/(ndim+nbright)

nbright=[37.8,56.8]
nbrighterr = [sqrt(4.)*3.7,sqrt(3.)*8.0]
;for an average, don't need the sqrt3 out front
nbrighterr = [3.7,8.0]
ndim = [27.0,0.2]
ndimerr=[sqrt(4.)*6.9,sqrt(3.)*7.3]
ndimerr=[6.9,7.3]
y500 = ndim/(ndim+nbright)

;;plot, x, y500, psym=5, thick=3, xthick=3, ythick=3, charthick=3,$
;;      xtitle = 'Redshift', ytitle = 'Nfaint/Ntotal', xrange=[0.0,1.2], xstyle=1, yrange=[0,1]
;oplot, x, y500, psym=6, thick=3
;;oplot, x, y500, thick=3
;;oplot, x, y750, psym = 2, thick=3
;error propagation
;sigy^2/y^2 = siga^2/a^2 + (siga^2 + sigb^2)/(a + b)^2
;sigy = sqrt((y^2)*( siga^2/a^2 + (siga^2 + sigb^2)/(a + b)^2))

sigy500 = sqrt((y500^2)*( ndimerr^2/ndim^2 + (ndimerr^2 + nbrighterr^2)/(ndim + nbright)^2))

;;errplot, x, y500 - sigy500,y500 + sigy500, thick=3, width=.04

;-----------------------------------------------------------------------
;consider faint/bright ratio.
faintbright = ndim/nbright

plot, x, faintbright, psym = 6, thick=3, xthick=3, ythick=3, charthick=3,$
      xtitle = 'Redshift', ytitle = 'Nfaint/Nbright', xrange=[0.0,1.2], xstyle=1, yrange=[0,1.3], ystyle=1

brighterr = (nbrighterr^2)/(nbright^2)
dimerr = (ndimerr^2) / (ndim^2)
;sigfaintbright = sqrt(faintbright^2* ((brighterr^2)/(nbright^2) + (dimerr^2)/(ndim^2)) )
sigfaintbright = sqrt(faintbright^2 * ((ndimerr^2)/(ndim^2) + (nbrighterr^2)/(nbright^2)))
print, sigfaintbright
errplot, x, ndim/nbright - sigfaintbright, ndim/nbright + sigfaintbright,  thick=3, width=.04

oplot, x, ndim750/nbright750, thick=3, psym = 2
ps_close, /noprint, /noid

end
