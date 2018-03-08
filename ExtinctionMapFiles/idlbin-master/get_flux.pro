 pro get_flux, filelist, pos1, pos2


;++++++++++++++++++++++++++++++++++++++
; Give it a file list of fits files and 
;it will compute the flux, centroids, backg, s/n
;and print them out.
;
INPUT:
;    file: string containing filename of fits data
;    ra, and dec of object
;
;PJL 09/09
;++++++++++++++++++++++++++++++++++++


spawn, 'wc '+filelist, wc
foo=strsplit(wc, /extract)
np=foo(0)

flux=dblarr(np)
fluxerr=dblarr(np)
xcent=dblarr(np)
ycent=dblarr(np)
bckg=dblarr(np)
sclk=dblarr(np)
sign=dblarr(np)
ra=pos1
dec=pos2

readcol, filelist, fitsname , format="A", /silent
for t=0, np-1 do begin 
   print, fitsname(t)
endfor

for i = 0, np-1 do begin
  get_centroids, fitsname(i), t, dt, x0, y0, f, xs, ys, fs, b, /WARM, /APER, RA=pos1, DEC=pos2
  flux[i] = f
  fluxerr[i] = fs
  xcent[i] = x0
  ycent[i] = y0 
  bckg[i] = b
  sclk[i] = t
  sign[i] = f/fs
endfor

print, 'file', ' RA       DEC       xcent  ', '   ycent   ', '   flux   ', '    fluxerr   ', '   S/N'
for j=0, np-1 do begin
    print, fitsname(j), ra, dec, xcent(j), ycent(j), flux(j), fluxerr(j), sign(j)
endfor
print, 'you are finished'
end
