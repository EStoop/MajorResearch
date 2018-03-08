pro irac2mass

restore, '/Users/jkrick/idlbin/objectnew.sav'

a = where(objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.tmasskmag gt 0 and objectnew.tmasskmag lt 90)

plot, objectnew[a].irac1mag, objectnew[a].tmasskmag, psym = 2, xrange=[12,20], yrange = [12,20]


end
