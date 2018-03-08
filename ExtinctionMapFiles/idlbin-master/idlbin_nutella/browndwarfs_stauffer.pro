


pro browndwarfs
close, /all
help, /memory
;device, true_color=24
;device, decomposed=0
!P.charsize = 1

;colors = GetColor(/load, Start=1)

;AB - Vega
;ch 1 = 2.7
;ch2 = 3.25
;ch3 = 3.73
;ch4 = 4.37

ps_open, filename='/Users/jkrick/nep/bdwarf/bdwarfs_stauffer.ps',/portrait,/square;,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

!P.thick = 5
!P.charthick = 5

restore, '/Users/jkrick/idlbin/idlbin_nutella/objectnew.sav'
acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');

object = objectnew

notstars = [898,922,1181,1658,5990,6082,6239,6881,7519,7636,8164]
object[notstars].acsclassstar = 0.
object[notstars].uclassstar = 0.
object[notstars].gclassstar = 0.
object[notstars].rclassstar = 0.
object[notstars].iclassstar = 0.

;-----------------------------------------------------------------
; first try to recreate plots of Patten et al with only the stars.  
;be wary of flagged items
;point8 = where(object.acsclassstar eq 1 or object.uclassstar eq 1 or object.gclassstar eq 1 or object.rclassstar eq 1 or object.iclassstar eq 1 and object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 and object.irac3mag gt 0 and object.irac3mag lt 90 and object.irac4mag gt 0 and object.irac4mag lt 90  and object.acsflag eq 0 and object.uflags eq 0 and object.gflags eq 0 and object.rflags eq 0 and object.iflags eq 0)

;print, "point8", n_elements(point8)

;-----------------------------------------------------------------
;then loosen the criterion to be anything with some ellip and fwhm (let go of classstar)

;relax = where(object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 and object.irac3mag gt 0 and object.irac3mag lt 90 and object.acsflag eq 0 and object.uflags eq 0 and object.gflags eq 0 and object.rflags eq 0 and object.iflags eq 0)
;print, "relax", n_elements(relax)

;-----------------------------------------------------------------

;bring in near-IR K band

withk = where(object.acsmag gt 0 and object.acsmag lt 90 and object.acsclassstar eq 1 and object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 and object.wirckmag gt 0. and object.wirckmag lt 90. )

print, "withk", n_elements(withk)

good = where(  (object[withk].wirckmag) - (object[withk].irac2mag-3.25) gt 1.0  )
print, "good", n_elements(good)

;plothist, object[withk[good]].rmaga, xrange=[10,30], bin=0.5, xtitle='rmag'
;-----------------------------------------------------------------

plot, (object[withk].irac1mag-2.7) - (object[withk].irac2mag-3.25), (object[withk].wirckmag) - (object[withk].irac1mag - 2.7), psym = 2, xrange=[-0.5, 2.0], yrange=[0, 4.0], xtitle = '[3.6]-[4.5]', ytitle='K-[3.6]', xstyle = 1, ystyle = 1, xthick = 5, ythick = 5, ticklen = 1

;oplot, findgen(20)/ 10, 1.05*(findgen(20)/10) + 1.2

;oplot, (object[withk[good]].irac1mag-2.7) - (object[withk[good]].irac2mag-3.25), (object[withk[good]].wirckmag) - (object[withk[good]].irac2mag - 3.25), psym = 2,color = redcolor


;plot, (object[withk].irac1mag-2.7) - (object[withk].irac2mag-3.25), (object[withk].wirckmag) - (object[withk].irac1mag - 2.7), psym = 2, xrange=[-0.5, 2.0], yrange=[0, 1.8], xtitle = '[3.6]-[4.5]', ytitle='K-[3.6]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

;oplot, (object[withk[good]].irac1mag-2.7) - (object[withk[good]].irac2mag-3.25), (object[withk[good]].wirckmag) - (object[withk[good]].irac1mag - 2.7), psym = 2, color = redcolor


;-----------------------------------------------------------------

;now if these are our candidates, where do they fall in the other color color plots?
;fortunately they are all detected in ch3 and ch 4

;plot, (object[withk[good]].irac2mag-3.25) - (object[withk[good]].irac3mag-3.73), (object[withk[good]].irac1mag - 2.7) - (object[withk[good]].irac2mag - 3.25), psym = 2, xrange=[-0.8, 0.6], yrange=[-0.5, 2.5], xtitle = '[4.5]-[5.8]', ytitle='[3.6]-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

;plot, (object[withk[good]].irac3mag-3.73) - (object[withk[good]].irac4mag-4.37), (object[withk[good]].irac1mag - 2.7) - (object[withk[good]].irac2mag - 3.25), psym = 2, xrange=[-0.2, 1.0], yrange=[-0.5, 2.5], xtitle = '[5.8]-[8.0]', ytitle='[3.6]-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

;plot, (object[withk[good]].irac1mag-2.7) - (object[withk[good]].irac4mag-4.37), (object[withk[good]].irac2mag - 3.25) - (object[withk[good]].irac3mag - 3.73), psym = 2, xrange=[-0.2, 2.3], yrange=[-0.7, 0.6], xtitle = '[3.6]-[8.0]', ytitle='[4.5]-[5.8]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1










ps_close, /noprint,/noid

end


