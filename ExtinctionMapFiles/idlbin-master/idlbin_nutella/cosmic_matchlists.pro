pro cosmic_matchlists

;readcol, '/Users/jkrick/palomar/cosmic/01june2008/listmask1.master', slitname, top, bottom, junk, format="A"
readcol, '/Users/jkrick/palomar/cosmic/slitmasks/june07/mask6.txt', slitname, xlocation,junk,junk,junk,junk, junk,junk,junk,junk,junk, format="A"
catname = lonarr(n_elements(slitname))
;readcol, '/Users/jkrick/palomar/cosmic/slitmasks/june08/crossref_good2.txt', slitnumber, catnumber, ra, dec, format="A"
readcol, '/Users/jkrick/palomar/cosmic/slitmasks/june07/crossref.txt', slitnumber, catnumber, ra, dec, format="A"
;openw, outlun, '/Users/jkrick/palomar/cosmic/01june2008/mask1.master', /append, /get_lun
openw, outlun, '/Users/jkrick/palomar/cosmic/june2007/mask6.master',  /get_lun

;put them in order
;aka sort on xlocation
slitname = slitname[sort(xlocation)]
xlocation = xlocation[sort(xlocation)]
for i = 0, n_elements(slitname) -1 do begin
   for j = 0, n_elements(slitnumber) - 1 do begin
      if slitnumber(j) eq slitname(i) then begin
         catname(i) = catnumber(j)
         printf, outlun, format='(I5, I5)', catname(i), slitname(i) ;format='(I5, 3X, I4, 3X,  I4, 3X,  I4, 3X, I1, 3X, F10.4, F10.4)',
      endif

   endfor
endfor

close, outlun
free_lun, outlun
;---------
;convert coord files to be the real catalog numbers

;readcol, '/Users/jkrick/palomar/cosmic/01june2008/mask1mosaic_tvmark.coo', x, y, slitname, format='(I2, I4, A)'
;slitname = strmid(slitname, 3)
;catname = lonarr(n_elements(slitname))
;readcol, '/Users/jkrick/palomar/cosmic/slitmasks/june08/crossref_good2.txt', slitnumber, catnumber, ra, dec, format="A"
;openw, outlun, '/Users/jkrick/palomar/cosmic/01june2008/mask1.reg', /get_lun

;printf, outlun, 'global color=green font="helvetica 10 normal" '
;printf, outlun, 'image'
;
;for i = 0, n_elements(slitname) -1 do begin
;   for j = 0, n_elements(slitnumber) - 1 do begin
;      if slitnumber(j) eq slitname(i) then begin
;         catname(i) = catnumber(j)
;         printf, outlun, '#text', x(i), ',', y(i),') text={', catname(i), '}'
;      endif

;   endfor
;endfor

;close, outlun
;free_lun, outlun






end
