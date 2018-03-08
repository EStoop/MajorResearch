pro mdwarfs
restore, '/Users/jkrick/idlbin/objectnew.sav'

;want things in vega for this work
;wirc already in vega
;need to convert irac magnitudes into vega

objectnew.irac1mag = objectnew.irac1mag-2.7
objectnew.irac2mag = objectnew.irac2mag-3.25
objectnew.irac3mag = objectnew.irac3mag-3.73
objectnew.irac4mag = objectnew.irac4mag-4.37


;all objects with k detections and 0.1<Ks-IRAC1<0.5 and
;0<IRAC(i)-IRAC(i+1)<0.1 and 0.1<IRAC1-IRAC4<0.2

;assuming mag errors are ~0.02
magerr=0.02
candidate1 = where(objectnew.wirckmag gt 0 and objectnew.wirckmag lt 90 and objectnew.wirckmag - objectnew.irac1mag lt 0.5+magerr and objectnew.wirckmag - objectnew.irac1mag gt 0.1-magerr  and objectnew.irac1mag - objectnew.irac2mag gt 0.-magerr and objectnew.irac1mag - objectnew.irac2mag lt 0.1+magerr and objectnew.irac2mag - objectnew.irac3mag gt 0.-magerr and objectnew.irac2mag - objectnew.irac3mag lt 0.1+magerr and objectnew.irac3mag - objectnew.irac4mag gt 0.-magerr and objectnew.irac3mag - objectnew.irac4mag lt 0.1+magerr  and objectnew.irac1mag - objectnew.irac4mag gt 0.1-magerr and objectnew.irac1mag - objectnew.irac2mag lt 0.2+magerr )

print, 'n candidate1', n_elements(candidate1)

openw, outlunred, '/Users/jkrick/nep/mdwarf/candidate1.reg', /get_lun
printf, outlunred, 'fk5'
for rc=0, n_elements(candidate1) -1 do  begin
   printf, outlunred, 'circle( ', objectnew[candidate1[rc]].ra, objectnew[candidate1[rc]].dec, ' 3")'
   print, format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', candidate1[rc], objectnew[candidate1[rc]].wircjmag, objectnew[candidate1[rc]].wirchmag, objectnew[candidate1[rc]].wirckmag, objectnew[candidate1[rc]].irac1mag, objectnew[candidate1[rc]].irac2mag, objectnew[candidate1[rc]].irac3mag,  objectnew[candidate1[rc]].irac4mag
endfor

close, outlunred
free_lun, outlunred

;plothyperz, candidate1, '/Users/jkrick/nep/mdwarf/candidate1.ps'


;---------------------------------------------------------------------------------------
;no k requirement

candidate2 = where(objectnew.irac1mag - objectnew.irac2mag gt 0.-magerr and objectnew.irac1mag - objectnew.irac2mag lt 0.1+magerr and objectnew.irac2mag - objectnew.irac3mag gt 0.-magerr and objectnew.irac2mag - objectnew.irac3mag lt 0.1+magerr and objectnew.irac3mag - objectnew.irac4mag gt 0.-magerr and objectnew.irac3mag - objectnew.irac4mag lt 0.1+magerr  and objectnew.irac1mag - objectnew.irac4mag gt 0.1-magerr and objectnew.irac1mag - objectnew.irac2mag lt 0.2+magerr );and objectnew.wirckmag gt 90)

print, 'n candidate2', n_elements(candidate2)

openw, outlunred, '/Users/jkrick/nep/mdwarf/candidate2.reg', /get_lun
printf, outlunred, 'fk5'
for rc=0, n_elements(candidate2) -1 do  begin
   printf, outlunred, 'circle( ', objectnew[candidate2[rc]].ra, objectnew[candidate2[rc]].dec, ' 3")'
   print, format='(I10,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)', candidate2[rc], objectnew[candidate2[rc]].wircjmag, objectnew[candidate2[rc]].wirchmag, objectnew[candidate2[rc]].wirckmag, objectnew[candidate2[rc]].irac1mag,  objectnew[candidate2[rc]].irac2mag,  objectnew[candidate2[rc]].irac3mag, objectnew[candidate2[rc]].irac4mag

endfor
close, outlunred
free_lun, outlunred

;plothyperz, candidate2, '/Users/jkrick/nep/mdwarf/candidate2.ps'

;---------------------------------------------------------------------------------------
;only K requirement
candidate3 = where(objectnew.wirckmag gt 0 and objectnew.wirckmag lt 90 and objectnew.wirckmag - objectnew.irac1mag lt 0.5+magerr and objectnew.wirckmag - objectnew.irac1mag gt 0.1-magerr and  objectnew.acsmag gt 0 and objectnew.acsmag lt 90 and objectnew.acsfwhm lt 3.2 and objectnew.acsellip lt 0.19 )

print, 'n3', n_elements(candidate3)

plothyperz, candidate3, '/Users/jkrick/nep/mdwarf/candidate3.ps'
end


;and  objectnew.acsmag gt 0 and objectnew.acsmag lt 90 and objectnew.acsfwhm lt 3.2 and objectnew.acsellip lt 0.19
