pro cosmic_done

restore, '/Users/jkrick/idlbin/object.sav'

filenamelist = ['/Users/jkrick/palomar/cosmic/slitmasks/mask1.txt', '/Users/jkrick/palomar/cosmic/slitmasks/mask2.txt','/Users/jkrick/palomar/cosmic/slitmasks/mask3.txt','/Users/jkrick/palomar/cosmic/slitmasks/mask4.txt','/Users/jkrick/palomar/cosmic/slitmasks/mask6.txt' ]

openw, outlun, '/Users/jkrick/palomar/cosmic/slitmasks/observed.reg', /get_lun
printf, outlun, 'fk5'
numlist = fltarr(500)
j = 0
for i = 0, n_elements(filenamelist) -1 do begin
   readcol, filenamelist(i), badnum, junk, junk, junk, junk, junk, junk, junk,junk,junk,junk
   for k = 0, n_elements(badnum) - 1 do begin
      if badnum(k) gt 0 then begin
         numlist(j) = badnum(k)
         j = j + 1
      endif

   endfor
 
endfor
numlist = numlist[0:j-1]


readcol, '/Users/jkrick/palomar/cosmic/slitmasks/bmasklist.txt', cosmicnum, hh, mm, ss, dd, dm, ds, mag
ra = ((ss /60. + mm )/60. + hh ) * 15.
dec = (ds /60. + dm )/60. + dd 

print, ra(0), dec(0)

;now match those ra & dec values with the catalog.  output a list of slitmask numbers, catalog numbers, ra, dec

; create initial arrays
m=long(n_elements(numlist))
ir=n_elements(object.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999

print,'Matching  to object'
dist=irmatch
dist[*]=0
for q=0,m-1 do begin
   a = where( cosmicnum eq  numlist[q] )
   print, format='(I10,I10,I10,F10.5,F10.5)',numlist[q], a, cosmicnum[a], ra[a], dec[a]
   printf, outlun, format='(A10,F10.5,F10.5,A5)', 'circle( ', ra[a], dec[a], '10")'

;   findobject, ra[a], dec[a]
   dist=sphdist( ra[a], dec[a], object.ra, object.dec,/degrees)
   sep=min(dist,ind)
   print, "sep", sep
   if (sep LT 0.01)  then begin
      mmatch[q]=ind
      print, numlist[q], ind, ra[cosmicnum[a]], dec[cosmicnum[a]], object[ind].ra, object[ind].dec


     
   endif 
   
endfor
;print,"mmatch",  n_elements(mmatch), mmatch
close, outlun
free_lun, outlun

end


