  

function findobject_return, rawant, decwant
;device, true=24
;device, decomposed=0

close, /all

best = 1; do find only the best match, 0 = find all nearest matches


restore, '/Users/jkrick/idlbin/objectnew.sav'


; create initial arrays
m=n_elements(rawant)
ir=n_elements(objectnew.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


;print,'Matching  to objectnew'
;print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin
;   print, q, rawant[q]
;   print, 'working on new obj'
   dist=sphdist( rawant[q], decwant[q], objectnew.ra, objectnew.dec, /degrees )
   sep=min(dist,ind)
;   print, 'matched ', ind
   return, ind


endfor



END

