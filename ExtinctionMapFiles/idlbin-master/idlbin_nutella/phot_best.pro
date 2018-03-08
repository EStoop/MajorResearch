pro photbest
close, /all
!P.THICK = 3
readcol, '/Users/jkrick/palomar/lfc/catalog/SExtractor.r.cat', NUMBER,X_WORLD      ,Y_WORLD      ,X_IMAGE,Y_IMAGE,FLUX_AUTO,MAG_AUTO,MAGERR_AUTO,MAG_APER1,magerr1,FWHM_IMAGE,ISOAREA_IMAGE,FLUX_MAX,ELLIPTICITY, format="A"

autoarr = fltarr(n_elements(number))
bestarr = fltarr(n_elements(number))
aperarr1 = autoarr
aperarr2 = autoarr
aperarr3 = autoarr
isoarr = autoarr

k = 0
for i = 0l, n_elements(number) - 1 do begin
   if y_image(i) lt 6000 and y_image(i) gt 2200 and x_image(i) gt 3000 and x_image(i) lt 6000 and fwhm_image(i) gt 5.0 then begin
      autoarr(k) = mag_auto(i)
;      bestarr(k) = mag_best(i)
;      aperarr3(k) = magaper3(i)
;      aperarr2(k) = magaper2(i)
      aperarr1(k) = mag_aper1(i)
      isoarr(k) = isoarea_image(i)
      k = k + 1

      if mag_aper1(i) gt 30 and mag_aper1(i) lt 90 then print, mag_aper1(i), x_image(i), y_image(i), fwhm_image(i)
   endif
endfor

autoarr = autoarr[0:k-1]
bestarr = bestarr[0:k-1]
aperarr1 = aperarr1[0:k-1]
;aperarr2 = aperarr2[0:k-1]
;aperarr3 = aperarr3[0:k-1]
isoarr = isoarr[0:k-1]
ps_open, filename='/Users/jkrick/palomar/lfc/catalog/phot.ps', /portrait,/color

;plot, autoarr, bestarr, psym = 2, thick=3, xrange=[10,28],yrange=[10,28],xtitle='mag_auto', ytitle="mag_best"
;plothist, autoarr, xhist, yhist, bin=0.1, /noprint,yrange=[0,100], xtitle='mag_auto'
;plothist, aperarr, xhist, yhist, bin=0.1, /nopront, yrange=[0,100], xtitle = 'mag_aper'

;plot, autoarr, autoarr / aperarr, psym = 2, thick=3, xrange=[10,28], xtitle='mag_auto', ytitle='mag_auto / mag_aper3', charthick = 3, xthick=3, ythick=3
;plot, aperarr, autoarr / aperarr, psym = 2, thick=3, xrange=[10,28], xtitle='mag_aper', ytitle='mag_auto / mag_aper3', charthick = 3, xthick=3, ythick=3
;plot, autoarr, autoarr - aperarr2, psym = 2, thick=3, xrange=[10,28], xtitle='mag_auto', ytitle='mag_auto - mag_aper2', charthick = 3, xthick=3, ythick=3,yrange=[-10,10]
;plot, aperarr2, autoarr - aperarr2, psym = 2, thick=3, xrange=[10,28], xtitle='mag_aper', ytitle='mag_auto - mag_aper2', charthick = 3, xthick=3, ythick=3,yrange=[-10,10]

;plot, aperarr3, aperarr1, xtitle = 'magaper 3', ytitle='mag aper 1', psym = 2, charthick = 3, xthick=3, ythick=3,yrange=[10,30], xrange=[10,30]

plot, autoarr, autoarr - aperarr1, psym = 2, thick=3, xrange=[10,28], xtitle='mag_auto', ytitle='mag_auto - mag_aper1', charthick = 3, xthick=3, ythick=3,yrange=[-10,10]
plot, aperarr1, autoarr - aperarr1, psym = 2, thick=3, xrange=[10,28], xtitle='mag_aper', ytitle='mag_auto - mag_aper1', charthick = 3, xthick=3, ythick=3,yrange=[-10,10]

plot, isoarr,  autoarr - aperarr1, psym = 2, thick=3,  xtitle='isoarea', ytitle='mag_auto - magaper', charthick = 3, xthick=3, ythick=3,yrange=[-5,2], /xlog

ps_close, /noprint, /noid
end
