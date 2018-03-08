pro cosmicmatch

;restore, '/Users/jkrick/idlbin/objectnew.sav'

;read in 07 list
;readcol, '/Users/jkrick/palomar/cosmic/slitmasks/june07/observed.07.txt', junk, slitnum, rh, rm, rs, dd, mm, ss, mag, format="(A, I4, F10.2, F10.2, F10.2, F10.2, F10.2, F10.2, F10.2)"

;ra = 15.* tenv(rh, rm, rs)
;dec = tenv(dd, mm, ss)

;read in 08 list
;readcol, '/Users/jkrick/palomar/cosmic/01june2008/spectra/specz.txt', masknum, catnum08, photz, xccz, xcczerr,	zxc,	xcr,	xcgrade,	emcz,	emczerr,	zem,	emgrade,	whichmethod, format="(I1,L5,F10.3, F10.2,F10.2, F10.3, F10.2, I1, F10.2,F10.2,F10.3,I1,A)"


;what is on the 07 list that is not on the 08 list

;if it does not exist on 08 add it to an array

;first match the 07 ra and dec to catalog numbers.
;all of these should match.

;catnum07 = fltarr(n_elements(slitnum))
;for i = 0, n_elements(slitnum) - 1 do begin
;   catnum07(i) = findobject_return(ra(i), dec(i))
;endfor

;print, 'catnum07', catnum07

;save this.
;openw, outlunw, '/Users/jkrick/palomar/cosmic/slitmasks/june07/observed.07.cat.txt', /get_lun
;for ca = 0, n_elements(catnum07) - 1 do begin
;   printf, outlunw, catnum07(ca), slitnum(ca), ra(ca), dec(ca), mag(ca)
;endfor
;close, outlunw
;free_lun, outlunw

;restore, '/Users/jkrick/palomar/cosmic/slitmasks/june07/observed.07.cat.sav'

;;then match catalog numbers to catalog numbers
;nomatch = lonarr(n_elements(catnum07))
;q = 0
;bool07 = fltarr(n_elements(catnum07))
;for i = 0 , n_elements(catnum07) - 1 do begin
;   for j = 0, n_elements(catnum08) - 1 do begin
;;      print, 'start', catnum07(i), catnum08(j)
;      if catnum07(i) eq catnum08(j) then begin
;         bool07(i) = 1
;      endif  
;   endfor
;endfor



;nomatch = catnum07(where(bool07 lt 1))

;print, 'n_elements(nomatch', n_elements(nomatch)
;print, nomatch
;save, nomatch, filename='/Users/jkrick/palomar/cosmic/slitmasks/june07/nomatch.sav'

;make thumbnails and SED's, 
;are those likely to have emission lines.
;what are their phot-z's and mags.

;object 5915 is creating problems.  it has Nan's
;just skip it for now.
;nomatch = [nomatch[0:79],nomatch[81:85], nomatch[87:*]]

;plothyperz, nomatch, '/Users/jkrick/palomar/cosmic/slitmasks/june07/nomatchjune08.ps'

;get rid of stars
;nomatch = nomatch(where( objectnew[nomatch].acsclassstar lt 0.95))
;print, 'n_elements(nomatch', n_elements(nomatch)
;print, nomatch




;------------------------------------
mask3=[2489,7376,2071,2062,1893,1400,5096,1902,6445,5306,7410,2291,5845,5791,5077,5056,5070,6944,6694,7038,7075,27503,27839,27733]


mask4=[28888,8872,8299,9684,2637,731,800,11389,11058,10383,478,64320,453,563,482,761,719,358,680,1945,8726,8509,2358,350]

mask6=[2442,1985,12772,13342,1784,1061,13956,2436,2666,1102,2007,12940,1767,1498,2982,14709,13844,900,2142,14545,1788,13054]
allmasks07=[mask3,mask4,mask6]



readcol, '/Users/jkrick/palomar/cosmic/01june2008/spectra/specz.txt', masknum, catnum08, photz, xccz, xcczerr,	zxc,	xcr,	xcgrade,	emcz,	emczerr,	zem,	emgrade,	whichmethod, format="(I1,L5,F10.3, F10.2,F10.2, F10.3, F10.2, I1, F10.2,F10.2,F10.3,I1,A)"

;then match catalog numbers to catalog numbers
nomatch = lonarr(n_elements(allmasks07))
q = 0
bool07 = fltarr(n_elements(allmasks07))
for i = 0 , n_elements(allmasks07) - 1 do begin
   for j = 0, n_elements(catnum08) - 1 do begin
;      print, 'start', catnum07(i), catnum08(j)
      if allmasks07(i) eq catnum08(j) then begin
         bool07(i) = 1
      endif  
   endfor
endfor



nomatch = allmasks07(where(bool07 lt 1))
match = allmasks07(where(bool07 gt 0))
print, 'n_elements(nomatch', n_elements(nomatch)

for c = 0, n_elements(nomatch) - 1 do print, nomatch(c)
print, 'in between'
for c = 0, n_elements(match) - 1 do print, match(c)

end
