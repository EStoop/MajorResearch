pro cosmic_akari

;where do we have cosmic spectra and akari spectra?

;akari spectra
readlargecol,'/Users/jkrick/akari/nine_match.txt', akarinum,catalognumakari,akarira,akaridec,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk,junk

;palomar spectra
readcol, '/Users/jkrick/palomar/cosmic/slitmasks/observed.reg',palomarra, palomardec, format ='(X, F10.5, F10.5, X)'

for q = 0, n_elements(akarira) - 1 do begin
   dist = sphdist(akarira[q], akaridec[q], palomarra, palomardec,/degrees)
   sep = min(dist, ind)

   if sep lt 0.0001 then begin
      ;it did match

      print, format='(I10, I10, F10.5, F10.5)',akarinum[q], catalognumakari[q], akarira[q], akaridec[q]
   endif
endfor

end

