pro junk

restore, '/Users/jkrick/nutella/idlbin/objectnew.sav'

a = where (objectnew.acsmag gt 0 and objectnew.acsmag lt 90)
sortindex = sort(objectnew[a].acsmag)
sortmag = objectnew[a(sortindex)].acsmag
print, sortmag[0:10]
print, min(objectnew[a].acsmag)
end
