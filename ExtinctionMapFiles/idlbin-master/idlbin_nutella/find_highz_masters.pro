pro find_highz_masters


restore, '/Users/jkrick/idlbin/idlbin_nutella/objectnew.sav'
;restore, '/Users/jkrick/idlbin/object.sav'

;detected at ch1, ch2, ch3, but not ACS F814
good = where(objectnew.irac1mag gt 24. and objectnew.irac1mag lt 90 and objectnew.irac2mag gt 0 and objectnew.irac2mag lt 90 and objectnew.irac3mag gt 0 and objectnew.irac3mag lt 90  and objectnew.acsmag gt 90 and objectnew.rmag gt 90 and objectnew.gmag gt 90 and objectnew.imag gt 90 and objectnew.umag gt 90 and objectnew.irac3mag - objectnew.irac2mag gt 0. and objectnew.irac2mag - objectnew.irac3mag gt objectnew.irac1mag - objectnew.irac2mag + 0.2 , goodcount); , goodcount) ;and objectnew.swircjmagauto gt 90

print, "n_elements(good) ",goodcount


;good2 = where(objectnew[good].irac3flux - objectnew[good].irac2flux gt 0. and objectnew[good].irac2flux - objectnew[good].irac3flux gt objectnew[good].irac1flux - objectnew[good].irac2flux + 0.2, good2count )
;print, "n_elements(good2) ", good2count

openw, outlun, '/Users/jkrick/external/nep/highz/masters_catalog.txt', /get_lun
printf, outlun, 'number,  ra, dec, gmag, gmagerr, zmagbest, zmagerrbest,irac1mag, irac1magerr, irac2mag, irac2magerr, irac3mag, irac3magerr, irac4mag, irac4magerr,mips24mag, mips24magerr'

for i = 0, n_elements(good) - 1 do printf, outlun, format='(I10, F10.5, F10.5, F10.2, F10.2, F10.2,F10.2, F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2, F10.2, F10.2 )', good[i],  objectnew[good[i]].ra, objectnew[good[i]].dec, objectnew[good[i]].gmag, objectnew[good[i]].gmagerr, objectnew[good[i]].zmagbest, objectnew[good[i]].zmagerrbest,objectnew[good[i]].irac1mag, objectnew[good[i]].irac1magerr, objectnew[good[i]].irac2mag, objectnew[good[i]].irac2magerr, objectnew[good[i]].irac3mag, objectnew[good[i]].irac3magerr, objectnew[good[i]].irac4mag, objectnew[good[i]].irac4magerr,objectnew[good[i]].mips24mag, objectnew[good[i]].mips24magerr
close, outlun
free_lun, outlun

plothyperz, good, '/Users/jkrick/external/nep/highz/candidates_masters.ps'


end

