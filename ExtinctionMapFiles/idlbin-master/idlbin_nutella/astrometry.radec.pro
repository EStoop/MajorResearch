pro astrometry

close, /all

;
imagefilelist= ['/Users/jkrick/hst/drz/j9aj01010_drz.fits', '/Users/jkrick/hst/drz/j9aj02010_drz.fits']

;for c = 0, n_elements(imagefilelist) - 1 do begin
for c = 0, 1 do begin

;run Sextractor
commandline = '/SciApps/bin/sex ' + imagefilelist(c) + " -c /Users/jkrick/hst/raw/default.sex"
spawn, commandline


newname = strcompress("/Users/jkrick/hst/drz/match" + strmid(imagefilelist(c), 24,5) +".cat",/remove_all)
print, newname

readcol, '/Users/jkrick/hst/raw/rband.cat', refra, refdec, refmag,format="A"
;readcol, '/Users/jkrick/hst/raw/tmass.cat', refra, refdec, refmag,junk,junk,junk,junk,junk,junk,junk,junk,format="A"
readcol, '/Users/jkrick/hst/raw/test.asc', x,y,mag,flag,fwhm,ellip, class_star,format="A"

openw, outlun, newname, /get_lun

fits_read, imagefilelist(c), data, header

j = 0
ra = fltarr(n_elements(x))
dec = ra
;clean sextractor catalog
for i = 0, n_elements(x) - 1 do begin
   if mag(i) lt 25 and flag(i) lt 1 and fwhm(i) lt 3.0 and ellip(i) lt 0.3 then begin
      xyad, header, x(i), y(i), a, d
      ra(j) = a
      dec(j) = d
      j = j + 1
   endif

endfor

ra = ra[0:j-1]
dec = dec[0:j-1]

; create initial arrays
m=n_elements(ra)
ir=n_elements(refra)

irmatch=fltarr(ir);
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin

   dist=sphdist(ra[q],dec[q], refra, refdec,/degrees)
   sep=min(dist,ind)

   if (sep LT 0.0005) then begin
      mmatch[q]=ind
      adxy, header, refra[ind], refdec[ind], refx, refy
      adxy, header,  ra[q], dec[q], x, y
      printf, outlun, q, "  ", refra[ind],"  ", refdec[ind]
   endif 
endfor

print,"Finished at "+systime()

close, outlun
free_lun, outlun
endfor


readcol, '/Users/jkrick/hst/drz/newmatch01.cat',refx, refy, x, y,format="A"
fits_read, '/Users/jkrick/hst/drz/j9aj01010_drz.fits', data, header
openw, outlunnew, '/Users/jkrick/hst/drz/match01.gaia.txt', /get_lun
xyad, header, x, y, a, d

for i = 0, n_elements(a) - 1 do printf,outlunnew, a(i), d(i)

close, outlunnew
free_lun, outlunnew

end
