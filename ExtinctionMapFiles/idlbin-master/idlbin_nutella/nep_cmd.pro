pro nepcmd

restore, '/Users/jkrick/idlbin/idlbin_nutella/objectnew.sav'

goodall = where( objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.rmag gt 0 and objectnew.rmag lt 28 , goodcount)

print,goodcount

;a = plot (objectnew[goodall].rmag - objectnew[goodall].irac1mag ,objectnew[goodall].rmag,  '6r1.', xtitle ='r - ch1', ytitle='r mag ', title = 'NEP', yrange = [12, 28], xrange = [-4, 8])

;--------------------
goodall = where( objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.irac2mag gt 0 and objectnew.irac2mag lt 90 , goodcount)

print,goodcount

b = plot (objectnew[goodall].irac1mag - objectnew[goodall].irac2mag ,objectnew[goodall].irac1mag,  '6r1.', xtitle ='ch1 - ch2 (AB mag)', ytitle='ch1 mag AB', xrange = [-8, 5], yrange = [15, 27], title = 'NEP')

;b2 = plot (objectnew[goodall].irac1flux - objectnew[goodall].irac2flux ,objectnew[goodall].irac1flux,  '6r1.', xtitle ='ch1 - ch2 (flux)', ytitle='ch1 flux', /xlog,/ylog, title = 'NEP')

b3 = plot ( objectnew[goodall].irac2flux ,objectnew[goodall].irac1flux,  '6r1.', xtitle ='ch2 (flux)', ytitle='ch1 flux', /xlog,/ylog, title = 'NEP', xrange = [1E-4, 1E4], yrange = [1E-4, 1E4])
;--------------------
goodall = where( objectnew.rmag gt 0 and objectnew.rmag lt 90 and objectnew.irac2mag gt 0 and objectnew.irac2mag lt 90 , goodcount)
print,goodcount

c = plot (objectnew[goodall].rmag - objectnew[goodall].irac2mag ,objectnew[goodall].rmag,  '6r1.', xtitle ='r - ch2', ytitle='r mag ', title = 'NEP', yrange = [12, 28], xrange = [-4, 8])


end
