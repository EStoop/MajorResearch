pro tweakmerge
close, /all
name = "hello"
readcol, '/Users/jkrick/hst/raw/shiftscomb.txt', name, x, y, rot, mag,format="A"
readcol, '/Users/jkrick/hst/raw/shifts4.txt', name, x2, y2, rot2, mag2, format = "A"
openw, outlun, '/Users/jkrick/hst/raw/shiftscomb.txt', /get_lun

for i = 0, n_elements(x) - 1 do begin
   x[i] = x[i] + x2[i]
   y[i] = y[i] + y2[i]
   rot[i] = rot[i] + rot2[i]
endfor

printf, outlun, "#frame: output"
printf, outlun, "#reference: tweak_wcs_11_6.fits[wcs]"
printf, outlun, "#form: delta"
printf, outlun, "#units: pixels"

for j = 0, n_elements(x) - 1 do begin
   printf, outlun, name[j], x[j], y[j], rot[j], mag[j]
endfor

close, outlun
free_lun, outlun
end

