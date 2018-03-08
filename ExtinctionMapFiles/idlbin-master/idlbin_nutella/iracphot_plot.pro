pro iracphot_plot
!p.multi=[0,0,1]

restore, '/Users/jkrick/idlbin/objectnew.sav'
restore, '/Users/jkrick/idlbin/object.sav'


;acs detected
a = where(object.acsmag gt 0 and object.acsmag lt 90)

;acs detected and extended
ea = where(objectnew[a].extended gt 0)   

;acs detected and point source
pa = where(objectnew[a].extended lt 1)

;non acs detected but irac detected
j = where(object.acsmag lt 0 or object.acsmag gt 90)
i = where(object[j].irac1mag gt 0 and object[j].irac1mag lt 90)

;irac extended
ei1 = where(objectnew[j[i]].extended gt 0)

;irac point source
pi1 = where(objectnew[j[i]].extended lt 1)



ps_open, filename='/Users/jkrick/spitzer/irac/iracphot.ps',/portrait,/square,/color
!p.multi=[0,2,2]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
vsym,  /polygon, /fill

plot, object.irac1mag, object.irac1mag - objectnew.irac1mag, psym = 8, xrange=[10,30], symsize=0.2, xtitle='original aperture magnitudes', ytitle='original - new',  charthick=3, xthick=3, ythick=3, title='acs detected, extended objects'
oplot, object[a[ea]].irac1mag, object[a[ea]].irac1mag - objectnew[a[ea]].irac1mag, psym = 8, color=redcolor, symsize=0.2

plot, object.irac1mag, object.irac1mag - objectnew.irac1mag, psym = 8, xrange=[10,30], symsize=0.2, xtitle='original aperture magnitudes', ytitle='original - new',  charthick=3, xthick=3, ythick=3, title='acs detected, point source objects'
oplot, object[a[pa]].irac1mag, object[a[pa]].irac1mag - objectnew[a[pa]].irac1mag, psym = 8, color=orangecolor, symsize=0.2

plot, object.irac1mag, object.irac1mag - objectnew.irac1mag, psym = 8, xrange=[10,30], symsize=0.2, xtitle='original aperture magnitudes', ytitle='original - new',  charthick=3, xthick=3, ythick=3, title='non acs detected, extended objects'
oplot, object[j[i[ei1]]].irac1mag, object[j[i[ei1]]].irac1mag - objectnew[j[i[ei1]]].irac1mag, psym = 8, color=bluecolor, symsize=0.2

plot, object.irac1mag, object.irac1mag - objectnew.irac1mag, psym = 8, xrange=[10,30], symsize=0.2, xtitle='original aperture magnitudes', ytitle='original - new',  charthick=3, xthick=3, ythick=3, title='non acs detected,point source objects'
oplot, object[j[i[pi1]]].irac1mag, object[j[i[pi1]]].irac1mag - objectnew[j[i[pi1]]].irac1mag, psym = 8, color=greencolor, symsize=0.2


ps_close, /noprint,/noid




end
;junk=where(object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac1mag - objectnew.irac1mag gt 4)
;openw, outlunred, '/Users/jkrick/spitzer/irac/phot.reg', /get_lun
;printf, outlunred, 'fk5'
;for rc=0, n_elements(junk) -1 do  printf, outlunred, 'circle( ', object[junk[rc]].ra, object[junk[rc]].dec, ' 3")'
;close, outlunred
;free_lun, outlunred


;readcol, '/Users/jkrick/spitzer/irac/SExtractor.1.cat', ch1NUMBER,ch1X_WORLD ,ch1Y_WORLD ,ch1X_IMAGE,ch1Y_IMAGE,ch1FLUX_AUTO,ch1MAG_AUTO,ch1MAGERR_AUTO,ch1FLUX_APER,ch1MAG_APER,ch1MAGERR_APER,ch1FLUX_BEST,ch1MAG_BEST,ch1MAGERR_BEST,ch1FWHM_IMAGE,ch1ISOAREA_IMAGE,ch1FLUX_MAX,ch1ELLIPTICITY,ch1CLASS_STAR,ch1FLAGS, ch1theta, ch1a, ch1b, format="A"

;ch1ratio = 1./(1.-ch1ellipticity)

;ap=sqrt(ch1isoarea_image/!PI)
;semi_minor=sqrt(ch1isoarea_image/(!PI*ch1ratio))
;semi_major=ch1ratio*semi_minor


;plot, semi_major, ch1a, psym = 8, symsize=0.2, xrange=[0,20], yrange=[0,20]
;oplot, findgen(30), findgen(30)
;plot, semi_minor, ch1b, psym = 8, symsize=0.2, xrange=[0,20], yrange=[0,20]
;oplot, findgen(30), findgen(30)
