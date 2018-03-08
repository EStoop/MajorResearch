pro irac_phot

restore, '/Users/jkrick/idlbin/object.sav'

ps_open, filename='/Users/jkrick/spitzer/irac/phot.ps',/portrait,/square,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

vsym, /polygon, /fill

;readcol, "/Users/jkrick/palomar/LFC/catalog/combined_catalog.prt",  wcsra, wcsdec, xcenter, ycenter, tmassJa,tmassha,tmasska,irac1a,irac2a,irac3a, irac4a,mips24a


readlargecol, '/Users/jkrick/spitzer/irac/SExtractor.1.cat', NUMBER,X_WORLD      ,Y_WORLD      ,X_IMAGE,Y_IMAGE,FLUX_AUTO,MAG_AUTO,MAGERR_AUTO,FLUX_APER,MAG_APER,MAGERR_APER,FLUX_BEST,MAG_BEST,MAGERR_BEST,FWHM_IMAGE,ISOAREA_IMAGE,FLUX_MAX,ELLIPTICITY,CLASS_STAR,FLAGS

;match object.irac1flux with new sextractor photometry
;-----------------------------------------------
; create initial arrays
ir=n_elements(number)

irmatch=fltarr(ir)
irmatch[*]=-999

print,"Starting at "+systime()
dist=irmatch
dist[*]=0

;for q=0,ir-1 do begin

;   dist=sphdist(X_WORLD[q],Y_WORLD[q], wcsra,wcsdec,/degrees)
;   sep=min(dist,ind)
;   if (sep LT 0.0008) then begin
;;      irmatch[q]=ind
;   endif 
;endfor

;print,"Finished at "+systime()

;matched=where(irmatch GT 0)
;nonmatched = where(irmatch lt 0)

;print, n_elements(matched),"matched"
;print, n_elements(nonmatched),"nonmatched"

!P.multi=[0,1,1]

;openw, outlun, '/Users/jkrick/spitzer/irac/phot.txt', /get_lun
;for i = 0, n_elements(matched) - 1 do printf, matched(i), irmatch[matched(i)]

;plot, flux_best(matched)*8.46, irac1a(irmatch[matched]), psym = 2, xrange=[0,100], yrange=[0,100], xtitle ='new Sextractor flux_best', ytitle ='original 2" aperture photometry'
;oplot, findgen(200), findgen(200)

plot, flux_aper*8.46/0.87, flux_best*8.46,psym = 8, xtitle='flux aper', ytitle = ' flux_auto', /xlog, /ylog,yrange=[0.001,10000], xrange=[0.001, 10000], xstyle = 1, ystyle = 1, thick = 3, charthick = 3, xthick = 3, ythick = 3, symsize = 0.3


x = [1E-4, 1E5]
y = [1E-4,1E5]
oplot, x,y, color = redcolor, thick = 3


plot, flux_aper*8.46/0.87, (flux_aper/0.87) / flux_auto, psym =8, xtitle='flux aper', ytitle = 'flux_aper/ flux_auto', /xlog, yrange=[0,2], xrange=[0.001, 100000], xstyle = 1, ystyle = 1,thick = 3, charthick = 3, xthick = 3, ythick = 3, symsize = 0.3

oplot,x, [1,1], color = redcolor, thick = 3

plot, flux_aper*8.46/0.87, isoarea_image, psym = 3, xtitle ='flux_aper', ytitle = 'isoarea_image', /xlog, /ylog, xrange=[0.001, 10000],/nodata, xstyle = 1, ystyle = 1, thick = 3, charthick = 3, xthick = 3, ythick = 3
a = where(class_star gt 0.7)
b = where(class_star le 0.7)
oplot, flux_aper(a), isoarea_image(a), psym = 8, color = bluecolor, symsize = 0.3
oplot, flux_aper(b), isoarea_image(b), psym = 8, color = redcolor, symsize = 0.3

plot, flux_auto, class_star, psym = 8, yrange=[0,1], /xlog, xrange=[0.001,10000], symsize = 0.3

plot, fwhm_image, ellipticity, psym = 8, xrange=[0,50], symsize = 0.3
ps_close, /noprint,/noid

end
