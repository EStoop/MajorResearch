pro wircjcat

readcol, '/Users/jkrick/palomar/wirc/SExtractor.j.cat', junk, ra, dec, junk, junk, junk, mag, junk, junk,junk,junk,junk,junk,junk,junk,junk

plothist, mag,bin = 0.1, xrange=[10,30]

a = where(mag gt 21.)
openw, outlun, '/Users/jkrick/palomar/wirc/j.reg', /get_lun
printf, outlun, 'fk5'
for i = 0, n_elements(a) - 1 do begin
   printf, outlun, 'circle( ', ra(a(i)), dec(a(i)), ' 40)'
endfor


close, outlun
free_lun, outlun

end
