pro agn

!P.multi = [0,0,1]
!P.charthick = 1
restore, '/Users/jkrick/idlbin/objectnew.sav'
;--------------------------------------------------
;mips detected, chandra undetected = 3200 sources...
print, 'mips not chandra', n_elements(where(objectnew.mips24mag gt 0 and objectnew.mips24mag lt 90 and objectnew.xnetrate le 0))

;--------------------------------------------------
;chandra detected, mips undetected

xnomips = where(objectnew.xnetrate gt 0 and objectnew.mips24mag gt 90)
;plothyperz, xnomips, '/Users/jkrick/nep/agn/x_no_mips.ps'

openw, outlunred, '/Users/jkrick/nep/agn/x_no_mips.reg', /get_lun
printf, outlunred, 'fk5'
for rc=0, n_elements(xnomips) -1 do  printf, outlunred, 'circle( ', objectnew[xnomips[rc]].ra, objectnew[xnomips[rc]].dec, ' 3")'
close, outlunred
free_lun, outlunred
;--------------------------------------------------
;extreme 24 - r
extreme24_r = where(objectnew.mips24mag gt 0 and objectnew.mips24mag lt 90 and objectnew.rmag gt 0 and objectnew.rmag lt 90 and objectnew.mips24mag - objectnew.rmag lt -6.0)
;plothyperz, extreme24_r, '/Users/jkrick/nep/agn/mips_r_extreme.ps'

;--------------------------------------------------
;blue u - r

blueu_r = where(objectnew.umag gt 0 and objectnew.umag lt 90 and objectnew.rmag gt 0 and objectnew.rmag lt 90 and objectnew.umag - objectnew.rmag lt -4.0)
;plothyperz, blueu_r, '/Users/jkrick/nep/agn/blueu_r.ps'
;--------------------------------------------------
;plot 24 / 8 vs 24

good = where(objectnew.mips24mag gt 0 and objectnew.mips24mag lt 90 and objectnew.irac4mag gt 0 and objectnew.irac4mag lt 90)

plot, objectnew[good].mips24mag, objectnew[good].mips24mag - objectnew[good].irac4mag, psym = 3
end
