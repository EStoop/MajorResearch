function prepare_fits, f
  readcol,f, unsortfitsname, format = 'A', /silent
  
  a = sort(unsortfitsname)
  fitsname = unsortfitsname(a)
  
     ;remove duplicates in processing number
  u = uniq(strmid(fitsname, 0, 57))
  fitsname = fitsname(u)
  
  openw, outlun, '/Users/jkrick/irac_warm/calstars/allch2bcdlist_sort.txt',/get_lun
  for fn = 0, n_elements(fitsname) - 1 do begin
     printf, outlun, fitsname(fn)
     
  endfor
  free_lun, outlun
  
return, fitsname
end
