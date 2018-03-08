pro find_2mass_stars

;readlargecol, '/Users/jkrick/irac_pc/SEP_2mass.txt', ra,       dec,err_maj,err_min,err_ang,      designation,   jmag,j_cmsig,j_msigcom,     j_snr,   hmag,h_cmsig,h_msigcom,     h_snr,   kmag,k_cmsig,k_msigcom,     k_snr,ph_qual,rd_flg,bl_flg,cc_flg,  ndet,  prox,    pxcntr,gal_contam,mp_flg, format = '    F10.6,    F10.6, F10.2, F10.2,    I5,             A,F10.3, F10.3,   F10.3,    F10.1,F10.3, F10.3,   F10.3,    F10.1,F10.3, F10.3,   F10.3,    F10.1,   A,  A,  A,  A,  A,F10.1,       I10,       I1,   I1'
readlargecol, '/Users/jkrick/irac_pc/all_2mass.txt', ra,       dec,  jmag,j_cmsig,j_msigcom,    hmag,h_cmsig,h_msigcom,     kmag,k_cmsig,k_msigcom, gal_contam,mp_flg, format = '    F10.6,    F10.6,F10.3, F10.3,   F10.3,F10.3, F10.3,   F10.3,   F10.3, F10.3,   F10.3,       I1,   I1'

 ;print, 'SEP K brightest', min(kmag)


;want to know the brightness of the other stars within +-2.5'
;of the bright stars.

a = sort(kmag)
kmag_sort = kmag[a]
jmag_sort = jmag[a]
hmag_sort = hmag[a]
ra_sort = ra[a]
dec_sort = dec[a]

print, 'the brightest stars', kmag_sort[0:40]

;for each of the 50 brightest stars
for i = 0, 500 do begin

  ; what is the mag of the other
  ; stars within 2.5' of the
  ; center of this star

   dist=sphdist(ra_sort(i),dec_sort(i),ra,dec, /degrees)

   good = where(dist lt .08)
   
   if n_elements(good) gt 1 then begin
      print, 'j,h,kmag', jmag_sort(i), hmag_sort(i), kmag_sort(i), ra_sort(i), dec_sort(i)
      
      for j = 0, n_elements(good) - 1 do begin
         if kmag(good(j)) lt 5.0 then   print, 'stars inside of FOV',  kmag(good(j)), dist(good(j))
      endfor
   endif

endfor

readlargecol, '/Users/jkrick/irac_pc/all_2mass_2.txt',ra,       dec,  jmag,j_cmsig,j_msigcom,    hmag,h_cmsig,h_msigcom,     kmag,k_cmsig,k_msigcom, gal_contam,mp_flg, rd_flag, jh,hk,jk,format = '    F10.6,    F10.6,F10.3, F10.3,   F10.3,F10.3, F10.3,   F10.3,   F10.3, F10.3,   F10.3,       I1,   I1,A,F10.2,F10.3,F10.3'



;want to know the brightness of the other stars within +-2.5'
;of the bright stars.

a = sort(kmag)
kmag_sort = kmag[a]
jmag_sort = jmag[a]
hmag_sort = hmag[a]
ra_sort = ra[a]
dec_sort = dec[a]

print, 'the brightest NEP stars', kmag_sort[0:10]

;for each of the 50 brightest stars
for i = 0, 5000 do begin

  ; what is the mag of the other
  ; stars within 2.5' of the
  ; center of this star

   dist=sphdist(ra_sort(i),dec_sort(i),ra,dec, /degrees)

   good = where(dist lt .08)
   
   if n_elements(good) gt 1 then begin
      print, 'j,h,kmag', jmag_sort(i), hmag_sort(i), kmag_sort(i), ra_sort(i), dec_sort(i)
      
      for j = 0, n_elements(good) - 1 do begin
         if kmag(good(j)) lt 5.0 then   print, 'stars inside of FOV',  kmag(good(j)), dist(good(j))
      endfor
   endif

endfor

end
