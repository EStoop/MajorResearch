pro astrometry

close, /all

;
imagefilelist= ['/Users/jkrick/hst/raw/j9aj04010_drz.fits', '/Users/jkrick/hst/raw/j9aj01010_drz.fits']

;for c = 0, n_elements(imagefilelist) - 1 do begin
c = 0
;run Sextractor
commandline = '/SciApps/bin/sex ' + imagefilelist(c) + " -c /Users/jkrick/hst/raw/default.sex"
spawn, commandline


newname = strcompress("/Users/jkrick/hst/raw/match0104.3.txt",/remove_all)
print, newname

readcol, '/Users/jkrick/hst/raw/test.asc', x,y,mag,flag,fwhm,ellip, class_star,format="A"

openw, outlun, newname, /get_lun

fits_read, imagefilelist(c), data, header

j = 0
ra = fltarr(n_elements(x))
dec = ra
;clean sextractor catalog
for i = 0, n_elements(x) - 1 do begin
   if mag(i) lt 24 and flag(i) lt 1  then begin
      xyad, header, x(i), y(i), a, d
      ra(j) = a
      dec(j) = d
      j = j + 1
   endif

endfor
ra = ra[0:j-1]
dec = dec[0:j-1]


;---------
c = c + 1
;run Sextractor
commandline = '/SciApps/bin/sex ' + imagefilelist(c) + " -c /Users/jkrick/hst/raw/default.sex"
spawn, commandline


readcol, '/Users/jkrick/hst/raw/test.asc', x,y,mag,flag,fwhm,ellip, class_star,format="A"


fits_read, imagefilelist(c), data2, header2

j = 0
ra2 = dblarr(n_elements(x))
dec2 = ra2
;clean sextractor catalog
for i = 0, n_elements(x) - 1 do begin
   if mag(i) lt 24 and flag(i) lt 1  then begin
      xyad, header2, x(i), y(i), a, d
 ;     print, x(i), y(i), a, d
      ra2(j) = a
      dec2(j) = d
      j = j + 1
   endif

endfor
ra2 = ra2[0:j-1]
dec2 = dec2[0:j-1]



; create initial arrays
m=n_elements(ra)
ir=n_elements(ra2)

irmatch=fltarr(ir);
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


print,"Starting at "+systime()
dist=irmatch
dist[*]=0


for q=0,m-1 do begin
  
   dist=sphdist(ra[q],dec[q], ra2, dec2,/degrees)
   sep=min(dist,ind)

   if (sep LT 0.0005) then begin
      mmatch[q]=ind
      adxy, header2, ra2[ind], dec2[ind], refx, refy
      adxy, header,  ra[q], dec[q], x, y
;      printf, outlun, q, "  ", refra[ind],"  ", refdec[ind]
   printf, outlun, ra2[ind], dec2[ind], x, y
   print, ra2[ind], dec2[ind], x, y
   endif 
endfor

print,"Finished at "+systime()




close, outlun
free_lun, outlun


end
