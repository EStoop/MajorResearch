pro find_highz_ranga


restore, '/Users/jkrick/idlbin/idlbin_nutella/objectnew.sav'
;restore, '/Users/jkrick/idlbin/object.sav'

;detected at ch3, ch4, but not R
good = where( objectnew.irac3mag gt 0 and objectnew.irac3mag lt 90  and objectnew.irac4mag gt 0 and objectnew.irac4mag lt 90 and objectnew.rmag gt 90 and objectnew.gmag gt 90 and objectnew.umag gt 90, goodcount); , goodcount) ;and objectnew.swircjmagauto gt 90

print, "n_elements(good) ",goodcount


;good2 = where(objectnew[good].irac3flux - objectnew[good].irac2flux gt 0. and objectnew[good].irac2flux - objectnew[good].irac3flux gt objectnew[good].irac1flux - objectnew[good].irac2flux + 0.2, good2count )
;print, "n_elements(good2) ", good2count

openw, outlun, '/Users/jkrick/external/nep/highz/ranga_catalog_16_01_22.txt', /get_lun
printf, outlun, 'number,  ra, dec, rmag, rmagerr, zmagbest, zmagerrbest,acsmag, acsmagerr, irac1mag, irac1magerr, irac2mag, irac2magerr, irac3mag, irac3magerr, irac4mag, irac4magerr,mips24mag, mips24magerr'

for i = 0, n_elements(good) - 1 do printf, outlun, format='(I10, F10.5, F10.5, F10.2, F10.2, F10.2,F10.2, F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2, F10.2, F10.2,F10.2, F10.2 )', good[i],  objectnew[good[i]].ra, objectnew[good[i]].dec, objectnew[good[i]].rmag, objectnew[good[i]].rmagerr, objectnew[good[i]].zmagbest, objectnew[good[i]].zmagerrbest,objectnew[good[i]].acsmag, objectnew[good[i]].acsmagerr, objectnew[good[i]].irac1mag, objectnew[good[i]].irac1magerr, objectnew[good[i]].irac2mag, objectnew[good[i]].irac2magerr, objectnew[good[i]].irac3mag, objectnew[good[i]].irac3magerr, objectnew[good[i]].irac4mag, objectnew[good[i]].irac4magerr,objectnew[good[i]].mips24mag, objectnew[good[i]].mips24magerr
close, outlun
free_lun, outlun

;;plothyperz, good, '/Users/jkrick/external/nep/highz/candidates_ranga.ps'


end

