pro find_akari


;read in akari data  
readcol, '/Users/jkrick/akari/cat_jan2909.txt', hh, mm, ss, dd, dm, ds, flux11a, format="A"

;convert h,m,s into degrees
ra = 15.* tenv(hh, mm, ss)
dec = tenv(dd, dm, ds)



openw, outlun, '/Users/jkrick/akari/egami_24match.txt', /get_lun

for i = 0,n_elements(ra) do begin
   
   find_object_akari, ra(i), dec(i)

endfor





close, outlun
free_lun, outlun 


end
