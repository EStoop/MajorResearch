pro find_simbad_stars

readlargecol, '/Users/jkrick/irac_pc/Nov_opz.txt', obj, ra, dec, jmag, hmag, kmag, ref, notes, format='L, F,F,F,F,F,L,L'

;print, 'SEP K brightest', min(kmag)


;want to know the brightness of the other stars within +-2.5'
;of the bright stars.

a = sort(kmag)
kmag_sort = kmag[a]
jmag_sort = jmag[a]
hmag_sort = hmag[a]
ra_sort = ra[a]
dec_sort = dec[a]

;print, 'the brightest stars', kmag_sort[0:40]

;for each of the 50 brightest stars
count = 0
for i = 0, n_elements(kmag_sort) - 1 do begin

  ; what is the mag of the other
  ; stars within 2.5' of the
  ; center of this star

   dist=sphdist(ra_sort(i),dec_sort(i),ra,dec, /degrees)

   good = where(dist lt .08)
   
   if n_elements(good) gt 2 then begin
      count = count + 1
;   if n_elements(good) gt 1 then begin
      print, 'j,h,kmag', jmag_sort(i), hmag_sort(i), kmag_sort(i), ra_sort(i), dec_sort(i)
      
      for j = 0, n_elements(good) - 1 do begin
         if kmag(good(j)) lt 8.0 then   print, 'stars inside of FOV',  kmag(good(j)), dist(good(j))
      endfor
   endif

endfor

print, 'copunt', count

end
