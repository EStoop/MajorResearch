pro star_acs
restore, '/Users/jkrick/idlbin/objectnew.sav'

readcol, '/Users/jkrick/nep/pop3/starlist.txt', ra, dec, format="A"


; create initial arrays
m=n_elements(ra)
ir=n_elements(objectnew.ra)

irmatch=fltarr(ir)
irmatch[*]=-999

;print,'Matching  to objectnew'

dist=irmatch
dist[*]=0

indexlist = fltarr(m)
for q=0,m-1 do begin
   dist=sphdist( ra[q], dec[q], objectnew.ra, objectnew.dec, /degrees )
   sep=min(dist,ind)
   indexlist(q) = ind
endfor




plothist, objectnew[indexlist].acsmag
print, 'max', max(objectnew[indexlist].acsmag)


;only want a certain brightness level
a = where(objectnew[indexlist].acsmag lt 20 and objectnew[indexlist].acsmag gt 15)
openw, outlunred, '/Users/jkrick/IRAC/iwic210/stars.txt', /get_lun
for i=0, n_elements(a) -1 do  printf, outlunred, objectnew[indexlist[a[i]]].ra, objectnew[indexlist[a[i]]].dec
close, outlunred
free_lun, outlunred


end
