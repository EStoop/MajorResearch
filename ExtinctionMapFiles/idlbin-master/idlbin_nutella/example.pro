pro example

restore, '/Users/jkrick/idlbin/object.sav'

;only work with those objects which are actually detected by acs and are decently bright
x = where(object.acsmag gt 0 and object.acsmag lt 25)

plot, object[x].acsellip, object[x].acsfwhm, psym = 3, yrange=[0,30]

;if you wanted to print out a list of ra and dec of all objects which appear starlike
star = where(object.acsmag gt 0 and object.acsmag lt 90 and object.acsellip lt 0.19 and object.acsfwhm lt 3.2)

for i = 0, n_elements(star) - 1 do print, object[star(i)].ra, object[star(i)].dec

;or if you want to go the easy route I have already classified things which appear starlike into the variable acsclassstar
;this does include some galaxies and is not complete, but if you just need a few test objects...
;object.acsclassstar = 1 for stars, lt 1 for everything else

print, "number of objects which are starlike",n_elements(where(object.acsclassstar eq 1.))
end
