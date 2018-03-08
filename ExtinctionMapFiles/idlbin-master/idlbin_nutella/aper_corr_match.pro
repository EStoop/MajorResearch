pro aper_corr_match

readcol, '/Users/jkrick/nutella/palomar/lfc/catalog/SExtractor.g.ap14.cat',  o, xwcs, ywcs, xcenter, ycenter, fluxauto, magauto, magerr, magaper, magapererr, fwhm, isoarea, fluxmax, ellip, class_star, flags, format="A"

readcol, '/Users/jkrick/nutella/palomar/lfc/catalog/SExtractor.g.ap33.cat',  o2, xwcs2, ywcs2, xcenter2, ycenter2, fluxauto2, magauto2, magerr2, magaper2, magapererr2, fwhm2, isoarea2, fluxmax2, ellip2, class_star2, flags2, format="A"



m = n_elements(xwcs)
mmatch = fltarr(m)
mmatch[*] = -999
ap14 = fltarr(m)
ap33  = fltarr(m)
t = 0
for q=0,m -1 do begin
   dist=sphdist( xwcs(q), ywcs(q),xwcs2,ywcs2,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin 
      ap14(t) = magaper(q)
      ap33(t) = magaper2(ind)
      t = t +1
   endif
endfor
ap14 = ap14[0:t-1]
ap33 = ap33[0:t-1]
plot, ap14, ap33, psym = 3

save, ap14, ap33, filename= 'junk.save'
end
