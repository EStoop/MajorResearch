pro readlexiobj

readcol, '/Users/jkrick/nep/lexiobj.txt', catnum, ra, dec, x, y

for i = 0, n_elements(ra) - 1 do begin
   if ra(i) lt 264 then begin
      ra(i) = ra(i) * 2
      dec(i) = dec(i)*2
   endif
endfor


findobject, ra, dec

end
