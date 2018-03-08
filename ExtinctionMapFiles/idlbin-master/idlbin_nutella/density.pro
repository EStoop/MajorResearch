pro density

restore, '/Users/jkrick/idlbin/objectnew.sav'
racenter = 264.88783
deccenter = 68.993018

d20 = where(objectnew.rmag lt 20.0 and objectnew.rmag gt 0 and sphdist(objectnew.ra, objectnew.dec,racenter,deccenter, /degrees) lt 0.122746)
d21 = where(objectnew.rmag lt 21.0 and objectnew.rmag gt 0 and sphdist(objectnew.ra, objectnew.dec, racenter,deccenter, /degrees) lt 0.122746)

print,n_elements(d20) / ( !PI*(7.3647)^2), n_elements(d21) /  (!PI*(7.3647)^2), !PI*(7.3647)^2
end
