


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

restore, '/Users/jkrick/idlbin/objectnew.sav'
acshead = headfits('/Users/jkrick/hst/raw/wholeacs.fits');



;----------------------------------------------------------------
; select based on all patten criterion.  all four irac plus k detections
a = where(objectnew.irac1mag gt 0  and objectnew.irac1mag lt 90 and objectnew.irac2mag gt 0 and objectnew.irac2mag lt 90 and objectnew.irac3mag gt 0 and objectnew.irac3mag lt 90  and objectnew.irac4mag gt 0 and objectnew.irac4mag lt 90 and objectnew.wirckmag gt 0 and objectnew.wirckmag lt 90)

b =  where(  (objectnew[a].wirckmag) - (objectnew[a].irac2mag-3.25) ge 1.0 );and  (objectnew[a].wirckmag) - (objectnew[a].irac1mag-2.7) ge 0.9 and (objectnew[a].irac1mag - 2.7) - (objectnew[a].irac2mag - 3.25) ge 0.1 and (objectnew[a].irac2mag - 3.25) - (objectnew[a].irac3mag - 3.73) le 0.5 and (objectnew[a].irac3mag - 3.73) - (objectnew[a].irac4mag - 4.37) ge 0.15)


c = where(objectnew[a[b]].acsmag gt 0 and objectnew[a[b]].acsmag lt 90 and objectnew[a[b]].acsfwhm lt 3.2 and objectnew[a[b]].acsellip lt 0.19)



print, 'a,b,c patten', n_elements(a), n_elements(b), n_elements(c)

for cbd = 0, n_elements(c) -1 do begin
   print, objectnew[a[b[c[cbd]]]].ra, objectnew[a[b[c[cbd]]]].dec
endfor

;----------------------------------------------------------------
; select based on all patten criterion.  3 irac plus k detections
d = where(objectnew.irac1mag gt 0  and objectnew.irac1mag lt 90 and objectnew.irac2mag gt 0 and objectnew.irac2mag lt 90 and objectnew.irac3mag gt 0 and objectnew.irac3mag lt 90 and objectnew.wirckmag gt 0 and objectnew.wirckmag lt 90)

e =  where(  (objectnew[d].wirckmag) - (objectnew[d].irac2mag-3.25) ge 1.0 and  (objectnew[d].wirckmag) - (objectnew[d].irac1mag-2.7) ge 0.9 and (objectnew[d].irac1mag - 2.7) - (objectnew[d].irac2mag - 3.25) ge 0.1 and (objectnew[d].irac2mag - 3.25) - (objectnew[d].irac3mag - 3.73) le 0.5 )


f = where(objectnew[d[e]].acsmag gt 0 and objectnew[d[e]].acsmag lt 90 and objectnew[d[e]].acsfwhm lt 3.2 and objectnew[d[e]].acsellip lt 0.19)

print, 'd,e,f patten', n_elements(d), n_elements(e), n_elements(f)

print, 'r', objectnew[d[e[f]]].rmag
print, 'k', objectnew[d[e[f]]].wirckmag

ps_close, /noprint,/noid

;plotsed, withk[good],'/Users/jkrick/nep/bdwarfs.sed.1.ps'
plotsed, a[b[c]],'/Users/jkrick/nep/bdwarf/bdwarfs.sed.3.ps';
plotsed, d[e[f]],'/Users/jkrick/nep/bdwarf/bdwarfs.sed.4.ps'
;plotsed, a[e],'/Users/jkrick/nep/bdwarfs.sed.5.ps'

;!P.multi=[0,3,2]
;!P.thick = 1
;!P.charthick = 1
;!P.charsize = 3
;plothist, (objectnew[a[b]].irac1mag - 2.7) - (objectnew[a[b]].irac2mag - 3.25), title =' [1-2]', xrange=[-0.5,2], bin=0.05
;plothist, (objectnew[a[b]].irac2mag - 3.25) - (objectnew[a[b]].irac3mag - 3.73), title ='[2-3]', xrange=[-1,1], bin=0.05
;plothist, (objectnew[a[b]].irac3mag - 3.73) - (objectnew[a[b]].irac4mag - 4.37) , title='[3-4]', xrange=[0,1], bin=0.05
;plothist, (objectnew[a[b]].wirckmag) - (objectnew[a[b]].irac2mag-3.25) , title='K-[2]', xrange=[0,4], bin=0.05
;plothist, (objectnew[a[b]].wirckmag) - (objectnew[a[b]].irac1mag-2.7) , title='K-[1]', xrange=[0,2], bin=0.05

end


