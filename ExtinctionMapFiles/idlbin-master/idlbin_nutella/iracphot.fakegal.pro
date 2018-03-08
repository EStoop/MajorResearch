pro iracphot_fake
 


;read in file list of fake galaxy images with images of all separation distances
readcol, '/Users/jkrick/spitzer/irac/fakegal/fakestarlist.txt', filename, format='A,F'

fluxa=fltarr(n_elements(filename))
fluxb=fltarr(n_elements(filename))
flluxa20 = 0.
fluxb20= 0.
separation = fltarr(n_elements(filename))
nr=0
magkeep = fltarr(2)
sepkeep = fltarr(2)
;for each image
for nim=0, n_elements(filename) - 1 do begin
   print,"working on ",  filename(nim)

; run sextractor on the image
   commandline = '/SciApps/bin/sex ' + filename(nim) + " -c /Users/jkrick/spitzer/irac/fakegal/defaultfake.sex"
   spawn, commandline

;read in the SExtractor file 
   readcol, '/Users/jkrick/spitzer/irac/SExtractor.f.cat', ch1NUMBER,ch1X_WORLD ,ch1Y_WORLD ,ch1X_IMAGE,ch1Y_IMAGE,ch1FLUX_AUTO,ch1MAG_AUTO,ch1MAGERR_AUTO,ch1FLUX_APER,ch1MAG_APER,ch1MAGERR_APER,ch1FLUX_BEST,ch1MAG_BEST,ch1MAGERR_BEST,ch1FWHM_IMAGE,ch1ISOAREA_IMAGE,ch1FLUX_MAX,ch1ELLIPTICITY,ch1CLASS_STAR,ch1FLAGS, ch1theta, ch1a, ch1b, format="A"

   ch1ratio = 1./(1.-ch1ellipticity)
   
;need to check to make sure there are 2 detections, otherwise, don't use this point.
   if n_elements(ch1number) gt 1 then begin
   a = ch1y_image[0] - ch1y_image[1]
   b = ch1x_image[0] - ch1x_image[1]
   sep = sqrt(a^2 + b^2)
   print, sep

;read in the image
   fits_read, filename(nim), ch1data, ch1head  
   naxis1=sxpar(ch1head,'NAXIS1')
   naxis2=sxpar(ch1head,'NAXIS2')

   psf=readfits('/Users/jason/psf/massimo/I1_psf_x4_best.fits',psf_hd)



; create arrays to hold output
   npoints=(n_elements(ch1number))
   print, 'there are ', n_elements(ch1number), ' detections'
;npoints=10
   iracflux_acs=fltarr(npoints)
   iracap_acs=fltarr(npoints)
   iracfluxerr = fltarr(npoints)
   finalflux = fltarr(npoints)
   finalfluxerr = fltarr(npoints)
   gain = [3.3]                 ; ch1, ch2-4=, 3.7, 3.8, 3.8]
   kron_factor=2.5
   min_radius=2.0               ; minimum radius for computing an ellipse (3.2 default)
   mask=ch1data                 ; initialize the mask to be the same as the image
   mask[*]=0
   mask_total=mask
   mask_flux=fltarr(npoints)
   mask_flux[*]=sqrt(-1)
   resid_flux=mask_flux
 
; compute the fiducial aperture correction ;
   aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,500,20,[40,50],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
   fiducial=pflux
   aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,500,6.4,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
   fiducial_cor=pflux/fiducial

;loop through all objects in the catalog
   for i = 0,  npoints- 1 do begin
   
  ; set the aperture size                                          ;

 ;           ap=sqrt(ch1isoareanow/3.141)
 ;           semi_minor=sqrt(ch1isoareanow/(3.141*ch1rationow))
;            semi_major=ch1rationow*semi_minor

;instead try using SExtractor's a & b and not isoarea with ellipticity
      semi_major = ch1a(i)
      semi_minor = ch1b(i)
      ap = (semi_major + semi_minor )/2.

; now set the minimum aperture
      if semi_major lt 3.2 then ap=3.2
;ap(where(finite(ap) NE 1))=3.2

; this fixes the ratio value to give a minimum semi_minor axis
; must be done _after_ semi-major axis is computed above
      semi_major_orig=semi_major
      if semi_major lt 3.2 then semi_major=3.2
      if semi_minor lt 2.0 then semi_minor = 2.0
            
                                ;if object is near to a bright object then give it a small aperture
                                ;tricky with irac detected things since sextractor knows that this is a confused image
      if ch1flags(i)  ge 32  then ap = 3.2


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; for all objects, measure a circular aperture
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      aper,ch1data,ch1x_image[i],ch1y_image[i],tflux,tfluxerr,tsky,tskyerr,3,ap,[ap+10,ap+20],[-100,1000],/nan,/exact,/flux,/silent

      tflux=tflux*8.46          ; convert to uJy


; derive the appropriate aperture correction
      apcor=fiducial_cor        ; default aperture
      if ((ap GT 3.2) and (ap LT 100))then begin
         aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,500,ap*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
         apcor=pflux/fiducial
      endif

      if (ap GT 100) then apcor=total(psf)/fiducial

      tflux=tflux/apcor         ; apply aperture correction

; register the fluxes and apertures
      iracflux_acs[i]=tflux
      iracap_acs[i]=ap
      iracfluxerr[i] = tfluxerr

      finalflux[i] = tflux
      finalfluxerr[i]=tfluxerr

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; do photometery on extended sources
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      if ap gt 3.2 then begin
;      if ((semi_major_orig GT min_radius)) then begin
         print, 'extended'
   
; first mask the aperture itself

         dist_ellipse,mask_ellipse,[naxis1,naxis2],ch1x_image[i],ch1y_image[i],ch1ratio[i],90+ch1theta[i]
         mask[*]=0
         mask(where(mask_ellipse LE kron_factor*semi_major))=1
         mask_total=mask_total OR mask
        
         mask_object=where(mask EQ 1)
         mask_object_area=total(mask)
         mask_object_flux=total(ch1data(mask_object))

; now process the sky
         mask[*]=0
         mask(where((mask_ellipse LE kron_factor*semi_major+15)and(mask_ellipse GE kron_factor*semi_major+5)))=1
         mask_sky=where(mask EQ 1)
         mmm,ch1data(mask_sky),background,background_sigma,background_skew
         mask_sky_flux=background*mask_object_area

         mask_flux[i]=8.46*(mask_object_flux-mask_sky_flux)


; derive the appropriate aperture correction

; first, recreate the effective aperture
         mask_radius=sqrt(mask_object_area / !PI)

         if (mask_radius LT 100) then begin
            aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,gain, mask_radius*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
            apcor=pflux/fiducial

         endif else begin
            apcor=total(psf)/fiducial
         endelse

         mask_flux[i]=mask_flux[i]/apcor
             
;            print, i, iracb(0), object[i].ra, object[i].dec, ch1isoareanow, ch1rationow, semi_major_orig, semi_major, semi_minor , iracap_acs(i), iracflux_acs(i), iracfluxerr(i), mask_flux[i], nearestdist, nearestflux, nearestisoarea
         if finite(mask_flux[i]) gt 0 then finalflux[i]=mask_flux[i] else finalflux[i] =iracflux_acs(i)


      endif

;the end
      finalmag = 8.926 - 2.5*alog10(finalflux[i]*1E-6)

      print, 'final', i, finalflux[i], finalfluxerr[i], finalmag, finalmag - (8.926 - 2.5*alog10((finalflux[i] + finalfluxerr[i])*1E-6))
      if sep gt 19 then begin
         fluxa20=0.1*finalflux[0] + finalflux[0]
         fluxb20=0.1*finalflux[1] + finalflux[1]
      endif

   endfor

   fluxa[nim] = finalflux[0]
   fluxb[nim] = finalflux[1]
   separation[nim] = sep

endif else begin
   fluxa[nim] = 0
   fluxb[nim] = 0
   separation[nim] = 0
endelse
endfor
;------
sortindex = sort(fluxa)
sortfluxa0 = fluxa(sortindex)
sortfluxa = sortfluxa0(where(sortfluxa0 gt 0))

sortsepa0 = separation(sortindex)
sortsepa = sortsepa0(where(sortfluxa0 gt 0))

limit = 0.1*sortfluxa(0) + sortfluxa(0)
print, 'limit ', limit
a = where(sortfluxa ge limit)
print, "sep val a", sortsepa(a(0))

sepkeep(nr) = sortsepa(a(0))
magkeep(nr) = 8.926 - 2.5*alog10(sortfluxa(0)*1E-6)
nr = nr + 1
;------
sortindex = sort(fluxb)
sortfluxb0 = fluxb(sortindex)
sortfluxb = sortfluxb0(where(sortfluxb0 gt 0))

sortsepb0 = separation(sortindex)
sortsepb = sortsepb0(where(sortfluxb0 gt 0))

limit = 0.05*sortfluxb(0) + sortfluxb(0)
print, 'limit ', limit
b = where(sortfluxb ge limit)
print, 'b', b
print, "sep val b", sortsepb(b(0))

sepkeep(nr) = sortsepb(b(0))
magkeep(nr) = 8.926 - 2.5*alog10(sortfluxb(0)*1E-6)
nr = nr + 1
;------



; save some results
save,/variables,filename='fake.sav'
ps_open, filename='/Users/jkrick/spitzer/irac/fakegal/separation.ps',/portrait,/color

!P.multi=[0,0,1]
plot, separation, fluxa, psym = 2,  symsize=0.5, yrange=[300,1000], ystyle = 1, xtitle='separation (irac pixels)', ytitle='flux', xthick=3, ythick=3, charthick=3
oplot, separation, fluxb, psym = 6,  symsize=0.5
oplot, findgen(100), findgen(100) - findgen(100) + fluxa20
oplot, findgen(100), findgen(100) - findgen(100) + fluxb20

plot, magkeep, sepkeep, psym = 2, xtitle='magnitude of star',ytitle= 'minimum separation value'

ps_close, /noprint,/noid

end
