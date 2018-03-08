pro astrometry

close, /all

;imagefilelist  = ['/Users/jkrick/hst/raw/j9aj25yeq_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj25yfq_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj25yhq_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj25yjq_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj25ylq_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj25ynq_single_sci.fits','/Users/jkrick/hst/raw/j9aj25ypq_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj25yrq_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj24lsq_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj24ltq_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj24lvq_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj24lxq_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj24lzq_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj24m1q_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj24m3q_single_sci.fits' ,'/Users/jkrick/hst/raw/j9aj24m5q_single_sci.fits']

;imagefilelist= ['/Users/jkrick/hst/drz/j9aj01010_drz.fits', '/Users/jkrick/hst/drz/j9aj02010_drz.fits']

imagefilelist = ['/Users/jkrick/hst/raw/d.fits']

for c = 0, n_elements(imagefilelist) - 1 do begin
;for c = 0, 0 do begin

;run Sextractor
commandline = '/SciApps/bin/sex ' + imagefilelist(c) + " -c /Users/jkrick/hst/raw/default.sex"
spawn, commandline

newname = strcompress("/Users/jkrick/hst/raw/matchd.txt",/remove_all)
;newname = strcompress("/Users/jkrick/hst/raw/match" + strmid(imagefilelist(c), 23,5) +".cat",/remove_all)
print, newname

readcol, '/Users/jkrick/hst/rband.cat', refra, refdec, refmag,format="A"
;readcol, '/Users/jkrick/hst/raw/tmass.cat', refra, refdec, refmag,junk,junk,junk,junk,junk,junk,junk,junk,format="A"
readcol, '/Users/jkrick/hst/raw/test.asc', x,y,mag,flag,fwhm,ellip, class_star,format="A"

openw, outlun, newname, /get_lun

fits_read, imagefilelist(c), data, header

j = 0
ra = fltarr(n_elements(x))
dec = ra
;clean sextractor catalog
for i = 0, n_elements(x) - 1 do begin
   if mag(i) lt 22 and flag(i) lt 1 and fwhm(i) lt 2.5 and ellip(i) lt 0.2 then begin
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
;      adxy, header, refra[ind], refdec[ind], refx, refy
      adxy, header,  ra[q], dec[q], x, y
      printf, outlun , refra[ind], refdec[ind], x,y 
;      printf, outlun, refra[ind], refdec[ind], "   ", refx, refy
;      printf, outlun, q, refra[ind], refra[dec]
   endif 
endfor

print,"Finished at "+systime()

close, outlun
free_lun, outlun
endfor



;fits_read,'/Users/jkrick/hst/drz/coadd.fits', data, coaddhead
;adxy, coaddhead, refra, refdec, x, y
;openw, outlunnew, '/Users/jkrick/hst/drz/testmatch.txt', /get_lun
;for i = 0, n_elements(refra) - 1 do printf, outlunnew, x(i), y(i)

;readcol, '/Users/jkrick/hst/drz/newmatch01.cat',refx, refy, x, y,format="A";
;fits_read, '/Users/jkrick/hst/drz/j9aj01010_drz.fits', data, header
;openw, outlunnew, '/Users/jkrick/hst/drz/match01.gaia.txt', /get_lun
;xyad, header, x, y, a, d

;for i = 0, n_elements(a) - 1 do printf,outlunnew, a(i), d(i)

;close, outlunnew;
;free_lun, outlunnew

;openw, out2, '/Users/jkrick/hst/drz/rband.gaia.cat', /get_lun
;for i = 0, n_elements(refra) - 1 do printf, out2, i, "  ", refra(i), refdec(i)
;close, out2
;free_lun, out2

end
