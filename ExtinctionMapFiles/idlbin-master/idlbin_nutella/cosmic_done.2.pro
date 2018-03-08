pro cosmic_done

restore, '/Users/jkrick/idlbin/object.sav'

openw, outlun, '/Users/jkrick/palomar/cosmic/objectmatch.txt', /get_lun


readcol, '/Users/jkrick/palomar/cosmic/slitmasks/bmasklist.txt', cosmicnum, hh, mm, ss, dd, dm, ds, mag
ra = ((ss /60. + mm )/60. + hh ) * 15.
dec = (ds /60. + dm )/60. + dd 


;now match those ra & dec values with the catalog.  output a list of slitmask numbers, catalog numbers, ra, dec

; create initial arrays
m=long(n_elements(cosmicnum))
ir=n_elements(object.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999

print,'Matching  to object'
dist=irmatch
dist[*]=0
for q=0,m-1 do begin
;   a = where( cosmicnum eq  numlist[q] )
;   print, format='(I10,I10,I10,F10.5,F10.5)',numlist[q], a, cosmicnum[a], ra[a], dec[a]
;   printf, outlun, format='(A10,F10.5,F10.5,A5)', 'circle( ', ra[a], dec[a], '10")'

;   findobject, ra[a], dec[a]
   dist=sphdist( ra[q], dec[q], object.ra, object.dec,/degrees)
   sep=min(dist,ind)
;   print, "sep", sep
   if (sep LT 0.01)  then begin
      mmatch[q]=ind
;      print, numlist[q], ind, ra[cosmicnum[a]], dec[cosmicnum[a]], object[ind].ra, object[ind].dec
      printf, outlun,  format='(I10,I10,F10.5,F10.5,F10.5,F10.5)' ,cosmicnum[q], ind, ra[q], dec[q], object[ind].ra, object[ind].dec

     
   endif 
   
endfor

close, outlun
free_lun, outlun


end


