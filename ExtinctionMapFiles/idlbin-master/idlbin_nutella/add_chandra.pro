pro addchandra

      
ftab_ext, '/Users/jkrick/chandra/evt2/image_all_src_wav.fits', [1,2,3,4,5,6,7,8,9], ra, dec, ra_err, dec_err, xpix, ypix, x_errpix, y_errpix, npixsou
ftab_ext, '/Users/jkrick/chandra/evt2/image_all_src_wav.fits', [10,11,12,13,14,15,16,17,18], netcounts, netcountserr, bkgcounts, bkgcountserr,netrate, netrateerr, bkgrate, bkgrateerr, exptime
ftab_ext, '/Users/jkrick/chandra/evt2/image_all_src_wav.fits', [19,20,21,22,23,25,26], exptimeerr, src_significance, psfsize, multi_correl_max, shape,  rotang, psfratio

print, 'ra', ra

print, 'min', min(netrate), min(netcounts)

chandra  = { xra:0D, xdec:0D, xra_err:0D, xdec_err:0D, xxpix:0D, xypix:0D, xx_errpix:0D, xy_errpix:0D, xnpixsou:0D, xnetcounts:0D, xnetcountserr:0D, xbkgcounts:0D, xbkgcountserr:0D, xnetrate:0D, xnetrateerr:0D, xbkgrate:0D, xbkgrateerr:0D, xexptime:0D, xexptimeerr:0D, xsrc_significance:0D, xpsfsize:0D, xmulti_correl_max:0D, xshape:' ',  xrotang:0D, xpsfratio:0D, xmatch:0D, xmatchdist:0D, xecf:0D, xflux:0D}


restore, '/Users/jkrick/idlbin/objectnew.sav'

b = replicate(chandra, n_elements(objectnew.ra))
objectnew = struct_addtags(objectnew, b)

objectnew.xmatch = -999
objectnew.xmatchdist = 999

;do some matching
m = n_elements(ra)
mmatch = fltarr(m)
mmatch[*] = -999

for q=0,m -1 do begin
   dist=sphdist( ra(q), dec(q),objectnew.ra,objectnew.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.001)  then begin  ;4"
      if objectnew[ind].xmatch ge 0 then begin   ;if there was a previous match
         if sep lt objectnew[ind].xmatchdist then begin   ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[objectnew[ind].xmatch] = -999
            ;put the new match into the matched pile
            mmatch[q]=ind
            objectnew[ind].xra = ra(q)
            objectnew[ind].xdec = dec(q)
            objectnew[ind].xra_err = ra_err(q)
            objectnew[ind].xdec_err = dec_err(q)
            objectnew[ind].xxpix = xpix(q)
            objectnew[ind].xypix = ypix(q)
            objectnew[ind].xx_errpix = x_errpix(q)
            objectnew[ind].xy_errpix = y_errpix(q)
            objectnew[ind].xnpixsou = npixsou(q)
            objectnew[ind].xnetcounts = netcounts(q)
            objectnew[ind].xnetcountserr = netcountserr(q)
            objectnew[ind].xbkgcounts = bkgcounts(q)
            objectnew[ind].xbkgcountserr = bkgcountserr(q)
            objectnew[ind].xnetrate = netrate(q)
            objectnew[ind].xnetrateerr = netrateerr(q)
            objectnew[ind].xbkgrate = bkgrate(q)
            objectnew[ind].xbkgrateerr = bkgrateerr(q)
            objectnew[ind].xexptime = exptime(q)
            objectnew[ind].xexptimeerr = exptimeerr(q)
            objectnew[ind].xsrc_significance = src_significance(q)
            objectnew[ind].xpsfsize = psfsize(q)
            objectnew[ind].xmulti_correl_max = multi_correl_max(q)
            objectnew[ind].xshape = shape(q)
            objectnew[ind].xrotang = rotang(q)
            objectnew[ind].xpsfratio = psfratio(q)
            objectnew[ind].xmatch = q
            objectnew[ind].xmatchdist = sep
            objectnew[ind].xecf = 1.076E-11
            objectnew[ind].xflux = netrate(q) * 1.076E-11

         endif
      endif else begin          ;if there was no previous match
         mmatch[q]=ind
         objectnew[ind].xra = ra(q)
         objectnew[ind].xdec = dec(q)
         objectnew[ind].xra_err = ra_err(q)
         objectnew[ind].xdec_err = dec_err(q)
         objectnew[ind].xxpix = xpix(q)
         objectnew[ind].xypix = ypix(q)
         objectnew[ind].xx_errpix = x_errpix(q)
         objectnew[ind].xy_errpix = y_errpix(q)
         objectnew[ind].xnpixsou = npixsou(q)
         objectnew[ind].xnetcounts = netcounts(q)
         objectnew[ind].xnetcountserr = netcountserr(q)
         objectnew[ind].xbkgcounts = bkgcounts(q)
         objectnew[ind].xbkgcountserr = bkgcountserr(q)
         objectnew[ind].xnetrate = netrate(q)
         objectnew[ind].xnetrateerr = netrateerr(q)
         objectnew[ind].xbkgrate = bkgrate(q)
         objectnew[ind].xbkgrateerr = bkgrateerr(q)
         objectnew[ind].xexptime = exptime(q)
         objectnew[ind].xexptimeerr = exptimeerr(q)
         objectnew[ind].xsrc_significance = src_significance(q)
         objectnew[ind].xpsfsize = psfsize(q)
         objectnew[ind].xmulti_correl_max = multi_correl_max(q)
         objectnew[ind].xshape = shape(q)
         objectnew[ind].xrotang = rotang(q)
         objectnew[ind].xpsfratio = psfratio(q)
         objectnew[ind].xmatch = q
         objectnew[ind].xmatchdist = sep
         objectnew[ind].xecf = 1.076E-11
         objectnew[ind].xflux = netrate(q) * 1.076E-11

      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))

;print out any that do not have a match

for num= 0,  n_elements(nonmatched) -1. do begin

   print, ra(nonmatched[num]), dec(nonmatched[num]), x_errpix(nonmatched[num]), y_errpix(nonmatched[num])
endfor

;these are nonmatches becuase they are near a bright star, or off the
;edge of most of my images.




save, objectnew, filename='/Users/jkrick/idlbin/objectnew.sav'

end
