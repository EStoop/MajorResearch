


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

ps_open, filename='/Users/jkrick/nep/bdwarf/bdwarfs3.ps',/portrait,/square;,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

!P.thick = 5
!P.charthick = 5

restore, '/Users/jkrick/idlbin/object.sav'
acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');


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

plothist, object[withk[good]].rmaga, xrange=[10,30], bin=0.5, xtitle='rmag'
;-----------------------------------------------------------------

plot, (object[withk].irac1mag-2.7) - (object[withk].irac2mag-3.25), (object[withk].wirckmag) - (object[withk].irac2mag - 3.25), psym = 2, xrange=[-0.5, 2.0], yrange=[0, 4.0], xtitle = '[3.6]-[4.5]', ytitle='K-[4.5]', xstyle = 1, ystyle = 1, xthick = 5, ythick = 5, ticklen = 1

oplot, findgen(20)/ 10, 1.05*(findgen(20)/10) + 1.2

;oplot, (object[withk[good]].irac1mag-2.7) - (object[withk[good]].irac2mag-3.25), (object[withk[good]].wirckmag) - (object[withk[good]].irac2mag - 3.25), psym = 2,color = redcolor


plot, (object[withk].irac1mag-2.7) - (object[withk].irac2mag-3.25), (object[withk].wirckmag) - (object[withk].irac1mag - 2.7), psym = 2, xrange=[-0.5, 2.0], yrange=[0, 1.8], xtitle = '[3.6]-[4.5]', ytitle='K-[3.6]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

oplot, (object[withk[good]].irac1mag-2.7) - (object[withk[good]].irac2mag-3.25), (object[withk[good]].wirckmag) - (object[withk[good]].irac1mag - 2.7), psym = 2, color = redcolor


;-----------------------------------------------------------------

;now if these are our candidates, where do they fall in the other color color plots?
;fortunately they are all detected in ch3 and ch 4

plot, (object[withk[good]].irac2mag-3.25) - (object[withk[good]].irac3mag-3.73), (object[withk[good]].irac1mag - 2.7) - (object[withk[good]].irac2mag - 3.25), psym = 2, xrange=[-0.8, 0.6], yrange=[-0.5, 2.5], xtitle = '[4.5]-[5.8]', ytitle='[3.6]-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

plot, (object[withk[good]].irac3mag-3.73) - (object[withk[good]].irac4mag-4.37), (object[withk[good]].irac1mag - 2.7) - (object[withk[good]].irac2mag - 3.25), psym = 2, xrange=[-0.2, 1.0], yrange=[-0.5, 2.5], xtitle = '[5.8]-[8.0]', ytitle='[3.6]-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

plot, (object[withk[good]].irac1mag-2.7) - (object[withk[good]].irac4mag-4.37), (object[withk[good]].irac2mag - 3.25) - (object[withk[good]].irac3mag - 3.73), psym = 2, xrange=[-0.2, 2.3], yrange=[-0.7, 0.6], xtitle = '[3.6]-[8.0]', ytitle='[4.5]-[5.8]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1









;-----------------------------------------------------------------

;try selecting only based on acs detection, acs starlike, and ch1 and ch2 detection

withoutk = where( object.acsmag gt 0 and object.acsmag lt 90 and object.acsfwhm lt 3.2 and object.acsellip lt 0.19 and object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 )

print, "withoutk", n_elements(withoutk)

good2 = where(  (object[withoutk].irac1mag - 2.7 ) - (object[withoutk].irac2mag-3.25) gt 0.2  )
print, "good2", n_elements(good2)

;-----------------------------------------------------------------
;how do they fit on those plots again
plot, (object[withoutk].irac2mag-3.25) - (object[withoutk].irac3mag-3.73), (object[withoutk].irac1mag - 2.7) - (object[withoutk].irac2mag - 3.25), psym = 2, xrange=[-0.8, 0.6], yrange=[-0.5, 2.5], xtitle = '[4.5]-[5.8]', ytitle='[3.6]-[4.5]', xstyle = 1, ystyle = 1, xthick = 3, ythick = 3, ticklen = 1

oplot, (object[withoutk[good2]].irac2mag-3.25) - (object[withoutk[good2]].irac3mag-3.73), (object[withoutk[good2]].irac1mag - 2.7) - (object[withoutk[good2]].irac2mag - 3.25), psym = 2, color = redcolor

;where are they on the acs image

adxy, acshead, object[withoutk[good2]].ra, object[withoutk[good2]].dec, xcenter3, ycenter3

openw, outlun, '/users/jkrick/nep/bdwarfs3.reg', /get_lun
for count = 0, n_elements(good2) - 1 do printf, outlun, 'circle( ', xcenter3(count), ycenter3(count), ' 100)'
close, outlun
free_lun, outlun 


;----------------------------------------------------------------
; try selecting objects detected in ch2, acs, and are acs starlike but not detected in ch1

good3 = where(object.irac2mag gt 0 and object.irac2mag lt 90 and object.irac1mag lt 0 and object.acsmag gt 0 and object.acsmag lt 90 and object.acsfwhm lt 3.2 and object.acsellip lt 0.19)

good4 = where(object.irac2mag gt 0 and object.irac2mag lt 90 and object.irac1mag gt 90 and object.acsmag gt 0 and object.acsmag lt 90 and object.acsfwhm lt 3.2 and object.acsellip lt 0.19)
print, good3, good4
print, 'good3', n_elements(good3) + n_elements(good4)

;adxy, acshead, object[good3].ra, object[good3].dec, xcenter3, ycenter3

;openw, outlun, '/users/jkrick/nep/bdwarfs4.reg', /get_lun
;for count = 0, n_elements(good3) - 1 do printf, outlun, 'circle( ', xcenter3(count), ycenter3(count), ' 100)'
;close, outlun
;free_lun, outlun 

;----------------------------------------------------------------
; how about selecting objects brighter in ch2 than ch1 and which have k - ch2 gt 1

a = where(object.irac2mag lt object.irac1mag and object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 and object.irac3mag gt 0 and object.irac3mag lt 90   )
print, "n a", n_elements(a)

b = where(object[a].wirckmag - object[a].irac2mag + 3.25 gt 1.0)
print, n_elements(b)
;----------------------------------------------------------------

a = where(object.irac1mag gt 0 and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 and object.wirckmag gt 0 and object.wirckmag lt 90 and object.irac2flux gt object.irac1flux and object.wirckmag - object.irac2mag + 3.25 gt 1.0)

 print, "now a", n_elements(a)


 c = where(object[a].acsmag gt 0 and object[a].acsmag lt 90 and object[a].acsfwhm lt 3.2 )
 print, n_elements(c)
d = where(object[a].acsmag lt 0);
e = where(object[a].acsmag gt 90)


;----------------------------------------------------------------
; select based on all patten criterion.  all four irac plus k detections
a = where(object.irac1mag gt 0  and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 and object.irac3mag gt 0 and object.irac3mag lt 90  and object.irac4mag gt 0 and object.irac4mag lt 90 and object.wirckmag gt 0 and object.wirckmag lt 90)

b =  where(  (object[a].wirckmag) - (object[a].irac2mag-3.25) ge 1.0 );and  (object[a].wirckmag) - (object[a].irac1mag-2.7) ge 0.9 and (object[a].irac1mag - 2.7) - (object[a].irac2mag - 3.25) ge 0.1 and (object[a].irac2mag - 3.25) - (object[a].irac3mag - 3.73) le 0.5 and (object[a].irac3mag - 3.73) - (object[a].irac4mag - 4.37) ge 0.15)


c = where(object[a[b]].acsmag gt 0 and object[a[b]].acsmag lt 90 and object[a[b]].acsfwhm lt 3.2 and object[a[b]].acsellip lt 0.19)



print, 'a,b,c patten', n_elements(a), n_elements(b), n_elements(c)

for cbd = 0, n_elements(c) -1 do begin
   print, object[a[b[c[cbd]]]].ra, object[a[b[c[cbd]]]].dec
endfor

;----------------------------------------------------------------
; select based on all patten criterion.  3 irac plus k detections
d = where(object.irac1mag gt 0  and object.irac1mag lt 90 and object.irac2mag gt 0 and object.irac2mag lt 90 and object.irac3mag gt 0 and object.irac3mag lt 90 and object.wirckmag gt 0 and object.wirckmag lt 90)

e =  where(  (object[d].wirckmag) - (object[d].irac2mag-3.25) ge 1.0 and  (object[d].wirckmag) - (object[d].irac1mag-2.7) ge 0.9 and (object[d].irac1mag - 2.7) - (object[d].irac2mag - 3.25) ge 0.1 and (object[d].irac2mag - 3.25) - (object[d].irac3mag - 3.73) le 0.5 )


f = where(object[d[e]].acsmag gt 0 and object[d[e]].acsmag lt 90 and object[d[e]].acsfwhm lt 3.2 and object[d[e]].acsellip lt 0.19)

print, 'd,e,f patten', n_elements(d), n_elements(e), n_elements(f)


ps_close, /noprint,/noid

;plotsed, withk[good],'/Users/jkrick/nep/bdwarfs.sed.1.ps'
plotsed, a[b[c]],'/Users/jkrick/nep/bdwarf/bdwarfs.sed.3.ps';
plotsed, d[e[f]],'/Users/jkrick/nep/bdwarf/bdwarfs.sed.4.ps'
;plotsed, a[e],'/Users/jkrick/nep/bdwarfs.sed.5.ps'

;!P.multi=[0,3,2]
;!P.thick = 1
;!P.charthick = 1
;!P.charsize = 3
;plothist, (object[a[b]].irac1mag - 2.7) - (object[a[b]].irac2mag - 3.25), title =' [1-2]', xrange=[-0.5,2], bin=0.05
;plothist, (object[a[b]].irac2mag - 3.25) - (object[a[b]].irac3mag - 3.73), title ='[2-3]', xrange=[-1,1], bin=0.05
;plothist, (object[a[b]].irac3mag - 3.73) - (object[a[b]].irac4mag - 4.37) , title='[3-4]', xrange=[0,1], bin=0.05
;plothist, (object[a[b]].wirckmag) - (object[a[b]].irac2mag-3.25) , title='K-[2]', xrange=[0,4], bin=0.05
;plothist, (object[a[b]].wirckmag) - (object[a[b]].irac1mag-2.7) , title='K-[1]', xrange=[0,2], bin=0.05

end


