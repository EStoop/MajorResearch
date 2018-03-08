pro dog
restore, '/Users/jkrick/idlbin/objectnew.sav'

a = where(objectnew.rmaga gt 0 and objectnew.rmaga lt 90 and objectnew.mips24mag gt 0 and objectnew.mips24mag lt 90)

mips24vega = - 2.5*alog10(objectnew[a].mips24flux/7.29/1E6)
print, objectnew[a[0]].mips24flux, mips24vega(0)

good = where(objectnew[a].rmaga - mips24vega gt 14)

print, n_elements(good)
!p.multi=[0,0,1]
plothist, objectnew[a[good]].mips24flux

plothyperz, a[good], '/Users/jkrick/nep/dog.ps'
end
