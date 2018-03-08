pro make_hyperz_filter
readcol,'/Users/jkrick/hst/f814acs.trace',lambda,flux,format="A"
openw, outlun, '/Users/jkrick/hst/f814acs.txt',/get_lun

for i = 0, n_elements(flux) -1 do begin
   printf,outlun,  i+1, "  ", lambda[i]*10., flux[i]
endfor

close, outlun
free_lun, outlun
end
