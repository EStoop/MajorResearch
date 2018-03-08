pro find_highz

!p.multi=[0,0,1]
ps_open, filename='/Users/jkrick/nep/highz/iraccolor.ps',/portrait,/square;,/color
vsym, /polygon, /fill

redcolor = FSC_COLOR("Red", !D.Table_Size-2);
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

restore, '/Users/jkrick/idlbin/objectnew_swirc.sav'
;restore, '/Users/jkrick/idlbin/object.sav'


good = where(objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.irac2mag gt 0 and objectnew.irac2mag lt 90 and objectnew.irac3mag gt 0 and objectnew.irac3mag lt 90 and objectnew.irac4mag gt 0 and objectnew.irac4mag lt 90  and objectnew.acsmag gt 90 and objectnew.zmagbest gt 24.2 and objectnew.mips24mag gt 90) ;

;good=[4016,10410,14147,14434,15638,8468,14341,15700,6722,8640,10894,14374,14739,15127,6192,3723,4950,6686,6828,8100,8602,11357,13180,13228,13269,13306,13317,13665,13739,14225,14268,14301,14368,14378,14497,14645,14810,14823,14850,15127,15409,15743,16154,11987,12011,12073,12114,12385]

print, "n_elements(good) ", n_elements(good)

openw, outlun, '/Users/jkrick/nep/highz/newlist.txt', /get_lun
;for i = 0, n_elements(good) - 1 do printf, outlun, format='(I10, F10.5, F10.5, F10.2, F10.2, F10.2,F10.2, F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2, F10.2, F10.2 )', good[i],  objectnew[good[i]].ra, objectnew[good[i]].dec, objectnew[good[i]].gmaga, objectnew[good[i]].gmagerra, objectnew[good[i]].zmagbest, objectnew[good[i]].zmagerrbest,objectnew[good[i]].irac1flux, objectnew[good[i]].irac1fluxerr, objectnew[good[i]].irac2flux, objectnew[good[i]].irac2fluxerr, objectnew[good[i]].irac3flux, objectnew[good[i]].irac3fluxerr, objectnew[good[i]].irac4flux, objectnew[good[i]].irac4fluxerr,objectnew[good[i]].mips24flux, objectnew[good[i]].mips24fluxerr
close, outlun
free_lun, outlun

;plothyperz, good, '/Users/jkrick/nep/highz/allcand.ps'

good2 = where(alog10(objectnew[good].irac1flux / objectnew[good].irac4flux) lt 0. and alog10(objectnew[good].irac2flux / objectnew[good].irac3flux) lt -1.)



;plot, objectnew[good].irac1mag - objectnew[good].irac4mag, objectnew[good].irac2mag - objectnew[good].irac3mag, psym = 8, xtitle = '3.6-8.0 (AB)', ytitle='4.5-5.8 (AB)', symsize=0.3, charthick=3, xthick=3,ythick=3

plot, alog10(objectnew[good].irac1flux / objectnew[good].irac4flux), alog10(objectnew[good].irac2flux / objectnew[good].irac3flux), psym = 8, xtitle = 'log 3.6/8.0', ytitle='log 4.5/5.8 ', symsize=0.3, charthick=3, xthick=3,ythick=3
;oplot, findgen(50) - findgen(50), findgen(50)
;oplot, findgen(50), findgen(50)- findgen(50)
;oplot, alog10(objectnew[good[good2]].irac1flux / objectnew[good[good2]].irac4flux), alog10(objectnew[good[good2]].irac2flux / objectnew[good[good2]].irac3flux), psym = 8,symsize=0.3, color=redcolor

;plot, objectnew[good].irac1mag - objectnew[good].irac4mag, objectnew[good].irac3mag - objectnew[good].irac4mag, psym = 8, xtitle = '3.6-8.0 (AB)', ytitle='5.8-8.0 (AB)', symsize=0.3, charthick=3, xthick=3,ythick=3
plot, alog10(objectnew[good].irac1flux / objectnew[good].irac4flux), alog10(objectnew[good].irac3flux / objectnew[good].irac4flux), psym = 8, xtitle = 'log 3.6/8.0 ', ytitle='log 5.8/8.0 ', symsize=0.3, charthick=3, xthick=3,ythick=3
;oplot, findgen(50) - findgen(50), findgen(50)
;oplot, findgen(50), findgen(50)- findgen(50)

plot, alog10(objectnew[good].irac1flux / objectnew[good].irac4flux), alog10(objectnew[good].irac3flux / objectnew[good].irac2flux), psym = 8, xtitle = 'log 3.6/8.0 ', ytitle='log 3.6/4.5 ', symsize=0.3, charthick=3, xthick=3,ythick=3, xrange=[-1,0.5], yrange=[-0.5,0.5]

;openw, outlunred, '/Users/jkrick/nep/highz/iraccolor.reg', /get_lun
;printf, outlunred, 'fk5'
;for rc=0, n_elements(good2) -1 do  printf, outlunred, 'circle( ', objectnew[good[good2[rc]]].ra, objectnew[good[good2[rc]]].dec, ' 3")'
;close, outlunred
;free_lun, outlunred

;-----------------------------------------------------------

onefour = where(alog10(objectnew[good].irac1flux / objectnew[good].irac4flux) lt 0 and alog10(objectnew[good].irac2flux / objectnew[good].irac3flux) lt 0 and alog10(objectnew[good].irac3flux / objectnew[good].irac4flux) gt 0 and alog10(objectnew[good].irac1flux / objectnew[good].irac2flux) ge 0)
print, '1/4 < 0 ', n_elements(onefour)

;plothyperz, good[onefour], '/Users/jkrick/nep/highz/onetwo.ps'
;-----------------------------------------------------------

a = where(alog10(objectnew[good].irac1flux / objectnew[good].irac4flux) le 0 and alog10(objectnew[good].irac1flux / objectnew[good].irac2flux) le 0 and alog10(objectnew[good].irac2flux / objectnew[good].irac3flux) le 0 and alog10(objectnew[good].irac3flux / objectnew[good].irac4flux) ge 0)
print, 'a ', n_elements(a)
;plothyperz, good[a], '/Users/jkrick/nep/highz/iraccolor.thumb.ps'


openw, outluna, '/Users/jkrick/nep/highz/iraccolorold.txt', /get_lun
for i = 0, n_elements(onefour) -1 do begin
   printf, outluna, format='(I10, F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2 )', good[onefour[i]],  objectnew[good[onefour[i]]].irac1flux, objectnew[good[onefour[i]]].irac1fluxerr, objectnew[good[onefour[i]]].irac2flux, objectnew[good[onefour[i]]].irac2fluxerr, objectnew[good[onefour[i]]].irac3flux, objectnew[good[onefour[i]]].irac3fluxerr, objectnew[good[onefour[i]]].irac4flux, objectnew[good[onefour[i]]].irac4fluxerr

endfor
close, outluna
free_lun, outluna
;-----------------------------------------------------------


;no optical detections

vred = where(objectnew.rmag gt 90 and   objectnew.acsmag gt 90 and  objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.irac2mag gt 0 and objectnew.irac2mag lt 90 )

print,'vred',  n_elements(vred)
;plothyperz, vred, '/Users/jkrick/nep/highz/vred.ps'

;--------------------------------------------------------------
;acs drops

acsdrop = where(objectnew.rmag gt 90 and objectnew.acsmag gt 90 and objectnew.irac1flux gt 0 and objectnew.irac2flux gt 0 and objectnew.irac3flux gt 0 and objectnew.irac4flux gt 0 and objectnew.swircjmagauto gt 0 and objectnew.swircjmagauto lt 90)
 
print, n_elements(acsdrop), 'acsdrop'

for newj = 0, n_elements(acsdrop) - 1 do begin
   print, newj, objectnew[acsdrop(newj)].swircjmagaper, objectnew[acsdrop(newj)].swircjmagauto
endfor

openw, outlunred, '/Users/jkrick/nep/highz/acsdrop.reg', /get_lun
printf, outlunred, 'fk5'
for c=0, n_elements(acsdrop) -1 do  printf, outlunred, 'circle( ', objectnew[acsdrop[c]].ra, objectnew[acsdrop[c]].dec, ' 3")'
close, outlunred
free_lun, outlunred

;plothyperz_highz, acsdrop, '/Users/jkrick/nep/highz/acsdrop.ps'



ps_close, /noprint,/noid

end

;ero = where(objectnew.rmaga gt 0 and objectnew.rmaga lt 90 and objectnew.wirckmag gt 0 and objectnew.wirckmag lt 90 and objectnew.rmaga - objectnew.wirckmag gt 5)

;print, 'ero', n_elements(ero)


;--------------------------------------------------------------
