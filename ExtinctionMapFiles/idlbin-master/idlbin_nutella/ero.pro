pro ero

restore, '/Users/jkrick/idlbin/objectnew.sav'

mips = where(objectnew.mips24flux gt 30 and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag gt 90 )

print, 'n_elements(mips)', n_elements(mips)


;color = objectnew[mips].irac1mag - objectnew[mips].irac2mag

;sort color
;sortcolor = color[sort(color)]
;sortmips = mips[sort(color)]
bright = objectnew[mips].irac1mag
sortbright = bright[sort(bright)]
sortmips = mips[sort(bright)]

print, sortbright
print, 'max', sortbright(0), sortmips(0);sortcolor(n_elements(mips) - 1), sortmips(n_elements(mips) - 1)

reddest = sortbright[10:20];sortcolor[n_elements(mips) - 10 : n_elements(mips) - 1 ]
reddestcat = sortmips[10:20];sortmips[n_elements(mips) - 10 : n_elements(mips) - 1 ]
print, 'top ten', reddest, reddestcat



plothyperz, reddestcat, '/Users/jkrick/nep/ero.ps'

openw, outlun, '/Users/jkrick/nep/ero.reg', /get_lun;
printf, outlun, 'fk5'
for i = 0, n_elements(reddestcat) - 1 do printf, outlun, 'circle(', objectnew[reddestcat(i)].ra, objectnew[reddestcat(i)].dec,' 5")'
close, outlun
free_lun, outlun

;plothyperz,red, '/Users/jkrick/nep/ero.ps'

end

;candidate = where(object.gmaga gt 90 and object.rmaga gt 90 and object.imaga gt 90 and object.umaga gt 90 and object.irac1mag lt 90 and object.irac1mag gt 0 and object.irac2mag gt 0 and object.irac2mag lt 99)
;candidate = where(object.acsmag gt 0 and object.acsmag gt 90 and object.irac1mag lt 90 and object.irac1mag gt 0 and object.irac2mag gt 0 and object.irac2mag lt 99 and object.gmaga gt 90)

;candidate =  where( object.acsmag gt 90 and object.irac1mag lt 90 and object.irac1mag gt 0 and object.irac2mag gt 0 and object.irac2mag lt 99 and object.gmaga gt 90)



;print, n_elements(candidate) , "candidate"
;plothist, 28.5 - object[candidate].irac1mag, bin = 0.1, thick=1, charthick=1

;

;b = where(28.5 - object[candidate].irac1mag gt 7 and 28.5 - object[candidate].irac1mag lt 9)



;plothyperz, candidate[b], '/Users/jkrick/nep/ero.ps'


