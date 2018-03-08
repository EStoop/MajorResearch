pro cosmic_makecrossref

;readcol, '/Users/jkrick/palomar/cosmic/slitmasks/june08/masklist.txt', slitnum, hr, min, sec, dd, dm, ds, rmag, format="A"
readcol, '/Users/jkrick/palomar/cosmic/slitmasks/june07/bmasklist.txt', slitnum, hr, min, sec, dd, dm, ds, rmag, format="A"

;convert hr, min, sec into degrees.
ra = ((hr/24.)*360.) + ((min/60.)*15.) + ((sec/60/60)*15.)
dec = (((ds/60.) + dm ) /60.) + dd

;need to figure out which catalog number this corresponds to.
catnum = lonarr(n_elements(slitnum))

;write out the results to a file
;openw, outlunw, '/Users/jkrick/palomar/cosmic/slitmasks/june07/crossref.txt', /GET_LUN
;printf, outlun, 'printing to file'

for i = 0, n_elements(slitnum) - 1 do begin

   catnum(i) = findobject_return(ra(i), dec(i))
;   print, 'catnum', catnum(i)
   print, slitnum(i), fix(catnum(i)), ra(i), dec(i)
endfor

;close, outlunw
;free_lun, outlunw
end
