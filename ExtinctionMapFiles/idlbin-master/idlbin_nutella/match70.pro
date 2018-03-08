pro match70

;read in the apex catalog
readcol, '/Users/jkrick/spitzer/mips/mips70/mosaic/test/mosaic_extract.txt', srcid, detid, apexra, dra, apexdec, ddec, drad, apexx, dx, apexy, dy, dxy, apexflux, dflux, ch2dof, psch2dof, status,snr, format="A"

match = fltarr(n_elements(srcid))
matchdist = fltarr(n_elements(srcid))
match[*] = -1
;read in the SExtractor catalog

readcol,  "/Users/jkrick/spitzer/mips/mips70/mosaic/Combine_old/mosaic.asc", num,ra, dec, x, y, fluxiso, fluxerriso,fluxauto, fluxerrauto,kron,back,isoarea,a, b, theta, flag, fwhm,elon, ellip, format="A"

;do some matching
m = n_elements(apexra)
apexmatchflux = fltarr(m)
sexmatchflux = fltarr(m)
i = 0

print, 'starting matching'
for q=0,n_elements(apexra) -1 do begin
   dist=sphdist( apexra(q), apexdec(q),ra,dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0011)  then begin  ;3"
      print, 'match', sep, apexx(q), x(ind), apexy(q), y(ind)
      apexmatchflux(i) = apexflux(q)
      sexmatchflux(i) = fluxiso(ind)
      i = i + 1
   endif else begin
      print, 'nomatch', sep, apexx(q), x(ind), apexy(q), y(ind)
   endelse

endfor

apexmatchflux = apexmatchflux[0:i-1]
sexmatchflux = sexmatchflux[0:i-1]

print, "matched", n_elements(apexmatchflux)
print, "nonmatched", n_elements(apexra) - n_elements(apexmatchflux)

;plot out the fluxes

plot, apexmatchflux, sexmatchflux * 376., psym = 2,/xlog,/ylog
oplot, findgen(1E6), findgen(1E6)


end
