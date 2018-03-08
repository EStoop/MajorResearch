pro iracphot


readcol, '/Users/jkrick/spitzer/irac/SExtractor.1.cat', ch1NUMBER,ch1X_WORLD ,ch1Y_WORLD ,ch1X_IMAGE,ch1Y_IMAGE,ch1FLUX_AUTO,ch1MAG_AUTO,ch1MAGERR_AUTO,ch1FLUX_APER,ch1MAG_APER,ch1MAGERR_APER,ch1FLUX_BEST,ch1MAG_BEST,ch1MAGERR_BEST,ch1FWHM_IMAGE,ch1ISOAREA_IMAGE,ch1FLUX_MAX,ch1ELLIPTICITY,ch1CLASS_STAR,ch1FLAGS, ch1theta, ch1a, ch1b, format="A"

ch1ratio = 1./(1.-ch1ellipticity)

fits_read, '/Users/jkrick/spitzer/irac/ch1/mosaic.fits', ch1data, ch1head  
fits_read, '/Users/jkrick/spitzer/irac/ch1/mosaic_cov.fits', ch1cov, ch1covhead  
naxis1=sxpar(ch1head,'NAXIS1')
naxis2=sxpar(ch1head,'NAXIS2')

fits_read, '/Users/jkrick/spitzer/irac/ch2/mosaic.fits', ch2data, ch2head  
fits_read, '/Users/jkrick/spitzer/irac/ch2/mosaic_cov.fits', ch2cov, ch2covhead  
fits_read, '/Users/jkrick/spitzer/irac/ch3/mosaic.fits', ch3data, ch3head  
fits_read, '/Users/jkrick/spitzer/irac/ch3/mosaic_cov.fits', ch3cov, ch3covhead  
fits_read, '/Users/jkrick/spitzer/irac/ch4/mosaic.fits', ch4data, ch4head  
fits_read, '/Users/jkrick/spitzer/irac/ch4/mosaic_cov.fits', ch4cov, ch4covhead  

psf1=readfits('/Users/jason/psf/massimo/I1_psf_x4_best.fits',psf_hd)
psf2=readfits('/Users/jason/psf/massimo/I2_psf_x4_best.fits',psf_hd)
psf3=readfits('/Users/jason/psf/massimo/I3_psf_x4_best.fits',psf_hd)
psf4=readfits('/Users/jason/psf/massimo/I4_psf_x4_best.fits',psf_hd)

restore, '/Users/jkrick/idlbin/object.sav'

;convert acs ra and dec to irac x, y
adxy, ch1head, object.ra, object.dec, ch1x, ch1y
adxy, ch2head, object.ra, object.dec, ch2x, ch2y
adxy, ch3head, object.ra, object.dec, ch3x, ch3y
adxy, ch4head, object.ra, object.dec, ch4x, ch4y


;onimage=where(((ch1x +15. LT naxis1) AND (ch1x- 15. GT 2) AND (ch1y + 15. LT naxis2) AND (ch1y - 15. GT 2)),count)

;want to keep both old and new
;would be better if move current irac1flux to irac1fluxold and put new flux into irac1flux

a  = { irac1fluxerr:0D, irac2fluxerr:0D, irac3fluxerr:0D, irac4fluxerr:0D, irac1magerr:0D,  irac2magerr:0D, irac3magerr:0D, irac4magerr:0D, nearestdist:0D, nearestnum:0, extended:0D, irac1fluxold:0D, irac2fluxold:0D, irac3fluxold:0D, irac4fluxold:0D, irac1magold:0D, irac2magold:0D, irac3magold:0D, irac4magold:0D }

b = replicate(a, n_elements(object.ra))
objectnew = struct_addtags(object, b)

objectnew.irac1fluxold = object.irac1flux
objectnew.irac2fluxold = object.irac2flux
objectnew.irac3fluxold = object.irac3flux
objectnew.irac4fluxold = object.irac4flux
objectnew.irac1magold = object.irac1mag
objectnew.irac2magold = object.irac2mag
objectnew.irac3magold = object.irac3mag
objectnew.irac4magold = object.irac4mag
;--------------------------------------------------------------------------
; create arrays to hold output
npoints=long(n_elements(object.ra))
npoints=100
iracap_acs=fltarr(npoints)
iracflux1_acs=fltarr(npoints)
irac1fluxerr = fltarr(npoints)
finalflux1 = fltarr(npoints)
finalfluxerr1 = fltarr(npoints)
iracflux2_acs=fltarr(npoints)
irac2fluxerr = fltarr(npoints)
finalflux2 = fltarr(npoints)
finalfluxerr2 = fltarr(npoints)
iracflux3_acs=fltarr(npoints)
irac3fluxerr = fltarr(npoints)
finalflux3 = fltarr(npoints)
finalfluxerr3 = fltarr(npoints)
iracflux4_acs=fltarr(npoints)
irac4fluxerr = fltarr(npoints)
finalflux4 = fltarr(npoints)
finalfluxerr4 = fltarr(npoints)
gain = [3.3, 3.7, 3.8, 3.8]
kron_factor=2.5
min_radius=2.0                  ; minimum radius for computing an ellipse (3.2 default)

mask=ch1data                    ; initialize the mask to be the same as the image
mask[*]=0
mask_total=mask
mask_flux=fltarr(npoints)
mask_flux[*]=sqrt(-1)
mask_flux1 = mask_flux
mask_flux2=mask_flux
mask_flux3=mask_flux
mask_flux4=mask_flux
;resid_flux=mask_flux


;testdeltamag = fltarr(7000)
;testskyerr = fltarr(7000)
g=0
acsratio=1./(1.-object.acsellip) ; convert sextrator ellipticity to ratio of smaj/smin

; compute the fiducial aperture correction ;
aper,psf1,521,521,pflux,pfluxerr,psky,pskyerr,500,20,[40,50],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
fiducial1=pflux
aper,psf1,521,521,pflux,pfluxerr,psky,pskyerr,500,6.4,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
fiducial_cor1=pflux/fiducial1

aper,psf2,521,521,pflux,pfluxerr,psky,pskyerr,500,20,[40,50],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
fiducial2=pflux
aper,psf2,521,521,pflux,pfluxerr,psky,pskyerr,500,6.4,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
fiducial_cor2=pflux/fiducial2

aper,psf3,521,521,pflux,pfluxerr,psky,pskyerr,500,20,[40,50],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
fiducial3=pflux
aper,psf3,521,521,pflux,pfluxerr,psky,pskyerr,500,6.4,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
fiducial_cor3=pflux/fiducial3

aper,psf4,521,521,pflux,pfluxerr,psky,pskyerr,500,20,[40,50],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
fiducial4=pflux
aper,psf4,521,521,pflux,pfluxerr,psky,pskyerr,500,6.4,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
fiducial_cor4=pflux/fiducial4

; define the bounding polygonal area for stars 
px=[12.5,25, 25, 10, 12.5]       ; acsmag 
py=[6e5 ,60, 35, 3e5, 6e5]       ; ascisoarea 
py=alog10(py)                    ; switch to log space 
roi = Obj_New('IDLanROI', px, py) 
star=roi->ContainsPoints(object.acsmag,alog10(object.acsisoarea)) 

;loop through all objects in the catalog
;for i = long(10000),  17000 do begin
for i = long(0),  npoints- 1 do begin
   
 
   if object[i].irac1mag gt 0.  then begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; figure out distance to nearest neighbor
;and the brightness and isoarea of that neighbor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;   print, 'working on ', i
   dist=sphdist( object[i].ra, object[i].dec,object.ra,object.dec,/degrees)
   sorteddist = dist(sort(dist))
   b = sort(dist)
   nearestdist = sorteddist(1)
   nearestflux  = object[b(1)].acsflux
   nearestisoarea = object[b(1)].acsisoarea
   objectnew[i].nearestdist = nearestdist
   objectnew[i].nearestnum = b(1)
   ;why isn't there an irac isoarea in the catalog?

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; calculate effective gain based on coverage map
;effective gain = N * gain
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   N1= ch1cov[ch1x[i], ch1y[i]]
   phpadu1 = N1*gain[0]
   N2= ch2cov[ch2x[i], ch2y[i]]
   phpadu2 = N2*gain[1]
   N3= ch3cov[ch3x[i], ch3y[i]]
   phpadu3 = N3*gain[2]
   N4= ch4cov[ch4x[i], ch4y[i]]
   phpadu4 = N4*gain[3]
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if the object is detected in acs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   if object[i].acsmag gt 0 and object[i].acsmag lt 90 then begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set the aperture size                                          ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ap=sqrt(object[i].acsisoarea/3.141)
      semi_minor=sqrt(object[i].acsisoarea/(3.141*acsratio[i]))
      semi_major=acsratio[i]*semi_minor
; convert that to IRAC mosaic units
      ap=ap*0.05/0.6
      semi_major=semi_major*0.05/0.6
      semi_minor=semi_minor*0.05/0.6
; now set the minimum aperture
      if ap lt 3.2 then ap =3.2
;ap(where(finite(ap) NE 1))=3.2

; this fixes the ratio value to give a minimum semi_minor axis
; must be done _after_ semi-major axis is computed above
      semi_major_orig=semi_major
      if semi_major lt 3.2 then semi_major=3.2
      if semi_minor lt 2.0 then semi_minor = 2.0

;make sure the stars get the circular apertures
      if star(i) gt 0 then ap = 3.2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; for all objects, measure a circular aperture
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      aper,ch1data,ch1x[i],ch1y[i],tflux,tfluxerr1,tsky,tskyerr,phpadu1,ap,[ap+10,ap+20],[-100,1000],/nan,/exact,/flux,/silent
      tflux1=tflux*8.46          ; convert to uJy
      aper,ch2data,ch2x[i],ch2y[i],tflux,tfluxerr2,tsky,tskyerr,phpadu2,ap,[ap+10,ap+20],[-100,1000],/nan,/exact,/flux,/silent
      tflux2=tflux*8.46          ; convert to uJy
      aper,ch3data,ch3x[i],ch3y[i],tflux,tfluxerr3,tsky,tskyerr,phpadu3,ap,[ap+10,ap+20],[-100,1000],/nan,/exact,/flux,/silent
      tflux3=tflux*8.46          ; convert to uJy
      aper,ch4data,ch4x[i],ch4y[i],tflux,tfluxerr4,tsky,tskyerr,phpadu4,ap,[ap+10,ap+20],[-100,1000],/nan,/exact,/flux,/silent
      tflux4=tflux*8.46          ; convert to uJy


; derive the appropriate aperture correction
      apcor1=fiducial_cor1        ; default aperture
      apcor2=fiducial_cor2        ; default aperture
      apcor3=fiducial_cor3        ; default aperture
      apcor4=fiducial_cor4        ; default aperture
      if ((ap GT 3.2) and (ap LT 100))then begin
         aper,psf1,521,521,pflux,pfluxerr,psky,pskyerr,gain[0],ap*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
         apcor1=pflux/fiducial1
         aper,psf2,521,521,pflux,pfluxerr,psky,pskyerr,gain[1],ap*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
         apcor2=pflux/fiducial2
         aper,psf3,521,521,pflux,pfluxerr,psky,pskyerr,gain[2],ap*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
         apcor3=pflux/fiducial3
         aper,psf4,521,521,pflux,pfluxerr,psky,pskyerr,gain[3],ap*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
         apcor4=pflux/fiducial4
      endif
      if (ap GT 100) then apcor1=total(psf1)/fiducial1
      if (ap GT 100) then apcor2=total(psf2)/fiducial2
      if (ap GT 100) then apcor3=total(psf3)/fiducial3
      if (ap GT 100) then apcor4=total(psf4)/fiducial4


; apply aperture correction
      tflux1=tflux1/apcor1       
      tflux2=tflux2/apcor2       
      tflux3=tflux3/apcor3       
      tflux4=tflux4/apcor4       

; register the fluxes and apertures
      iracap_acs[i]=ap
      iracflux1_acs[i]=tflux1
      irac1fluxerr[i] = tfluxerr1*8.46      
      finalflux1[i] =  iracflux1_acs[i]
      finalfluxerr1[i]=  irac1fluxerr[i]
      iracflux2_acs[i]=tflux2
      irac2fluxerr[i] = tfluxerr2*8.46      
      finalflux2[i] =  iracflux2_acs[i]
      finalfluxerr2[i]=  irac2fluxerr[i]
      iracflux3_acs[i]=tflux3
      irac3fluxerr[i] = tfluxerr3*8.46      
      finalflux3[i] =  iracflux3_acs[i]
      finalfluxerr3[i]=  irac3fluxerr[i]
      iracflux4_acs[i]=tflux4
      irac4fluxerr[i] = tfluxerr4*8.46      
      finalflux4[i] =  iracflux4_acs[i]
      finalfluxerr4[i]=  irac4fluxerr[i]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; do photometery on extended sources
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

      if ((semi_major_orig GT min_radius)) then begin
         print, 'extended acs '
         objectnew[i].extended = 1
; first mask the aperture itself
         for chanel=0, 3, do begin
            if chanel=0 then begin
               x = ch1x
               y = ch1y
               data = ch1data
               psf = psf1
               fiducial = fiducial1
            endif
            if chanel=1 then begin
               x = ch2x
               y = ch2y
               data = ch2data
               psf = psf2
               fiducial = fiducial2
            endif
            if chanel=2 then begin
               x = ch3x
               y = ch3y
               data = ch3data
               psf = psf3             
               fiducial = fiducial3
            endif
            if chanel=3 then begin
               x = ch4x
               y = ch4y
               data = ch4data
               psf = psf4
               fiducial = fiducial4
            endif

         dist_ellipse,mask_ellipse,[naxis1,naxis2],x[i],y[i],acsratio[i],90+object[i].acstheta
         mask[*]=0
         mask(where(mask_ellipse LE kron_factor*semi_major))=1
         mask_total=mask_total OR mask
        
         mask_object=where(mask EQ 1)
         mask_object_area=total(mask)
         mask_object_flux=total(data(mask_object))

; now process the sky
         mask[*]=0
         mask(where((mask_ellipse LE kron_factor*semi_major+15)and(mask_ellipse GE kron_factor*semi_major+5)))=1
         mask_sky=where(mask EQ 1)
         mmm,data(mask_sky),background,background_sigma,background_skew
         mask_sky_flux=background*mask_object_area

         mask_flux[i]=8.46*(mask_object_flux-mask_sky_flux)

; derive the appropriate aperture correction

; first, recreate the effective aperture
         mask_radius=sqrt(mask_object_area / 3.141)

         if (mask_radius LT 100) then begin
            aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,gain[chanel], mask_radius*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0

            apcor=pflux/fiducial

         endif else begin
            apcor=total(psf)/fiducial
         endelse

         mask_flux[i]=mask_flux[i]/apcor

;         print, i, object[i].ra, object[i].dec, object(i).acsisoarea, object[i].acsellip, semi_major_orig, semi_major, semi_minor , iracap_acs(i), iracflux_acs(i), iracfluxerr(i),  mask_flux[i], nearestdist, nearestflux, nearestisoarea
;         finalflux[i]=mask_flux[i]
            if finite(mask_flux[i]) gt 0 and chanel eq 0 then finalflux1[i]=mask_flux[i] else finalflux1[i] =iracflux1_acs(i)
            if finite(mask_flux[i]) gt 0 and chanel eq 1 then finalflux2[i]=mask_flux[i] else finalflux1[i] =iracflux2_acs(i)
            if finite(mask_flux[i]) gt 0 and chanel eq 2 then finalflux3[i]=mask_flux[i] else finalflux1[i] =iracflux3_acs(i)
            if finite(mask_flux[i]) gt 0 and chanel eq 3 then finalflux4[i]=mask_flux[i] else finalflux1[i] =iracflux4_acs(i)
         endfor   ;for all the chanells

      endif   ;if extended

   endif else begin
;-------------------------------------------------------------------------------------
;if the object is  NOT detected in acs but is on the irac frame
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;do all the same things as above.

      if object[i].irac1mag gt 0 then begin; and object[i].irac1mag lt 90 then begin
;         print, 'not detected in acs but in irac'
 ;match with SExtractor IRAC catalog.  
;if not in the catalog then must be a small object where a circular aperture would be a fine size.
         
         

         iracdist=sphdist( object[i].ra, object[i].dec,ch1x_world,ch1y_world,/degrees)        
         iracsep=min(iracdist,ind)
         sortediracdist = iracdist(sort(iracdist))
         iracb = sort(iracdist)

         if (iracsep ge 0.0002) then ap = 3.2
         if (iracsep LT 0.0002)  then begin
            ch1isoareanow = ch1isoarea_image(iracb(0))
            ch1rationow= ch1ratio(iracb(0))
            ch1thetanow = ch1theta(iracb(0))
; set the aperture size                                          ;

            ap=sqrt(ch1isoareanow/3.141)
            semi_minor=sqrt(ch1isoareanow/(3.141*ch1rationow))
            semi_major=ch1rationow*semi_minor


;instead try using SExtractor's a & b and not isoarea with ellipticity
           semi_major = ch1a(iracb(0))
           semi_minor = ch1b(iracb(0))
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
            if ch1flags(iracb(0))  ge 32  then ap = 3.2

         endif

;make sure the stars get the circular apertures
;???


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; for all objects, measure a circular aperture
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
         print, 'irac only ap', ap
         aper,ch1data,ch1x[i],ch1y[i],tflux,tfluxerr,tsky,tskyerr,phpadu1,ap,[ap+10,ap+20],[-100,1000],/nan,/exact,/flux,/silent

         tflux=tflux*8.46       ; convert to uJy


; derive the appropriate aperture correction
         apcor=fiducial_cor     ; default aperture
         if ((ap GT 3.2) and (ap LT 100))then begin
            aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,gain,ap*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
            apcor=pflux/fiducial
         endif

         if (ap GT 100) then apcor=total(psf)/fiducial

         tflux=tflux/apcor      ; apply aperture correction

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
            print, 'extended irac'
            objectnew[i].extended = 1

; first mask the aperture itself

            dist_ellipse,mask_ellipse,[naxis1,naxis2],ch1x[i],ch1y[i],ch1rationow,90+ch1thetanow
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
            mask_radius=sqrt(mask_object_area / 3.141)

            if (mask_radius LT 100) then begin
               aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,gain, mask_radius*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0
               apcor=pflux/fiducial

            endif else begin
               apcor=total(psf)/fiducial
            endelse

            mask_flux[i]=mask_flux[i]/apcor
            
;            print, i, iracb(0), object[i].ra, object[i].dec, ch1isoareanow, ch1rationow, semi_major_orig, semi_major, semi_minor , iracap_acs(i), iracflux_acs(i), iracfluxerr(i), mask_flux[i], nearestdist, nearestflux, nearestisoarea
            if finite(mask_flux[i]) gt 0 then finalflux[i]=mask_flux[i] else finalflux[i] =iracflux_acs(i)

;            testdeltamag[g] =  object[i].irac1mag - (8.926 - 2.5*alog10(finalflux[i]*1E-6))
;            testskyerr[g] =background_sigma
;            g = g + 1

         endif


      endif
      
   endelse


;update the catalog
objectnew[i].irac1flux = finalflux[i]
objectnew[i].irac1fluxerr=finalfluxerr[i]
objectnew[i].irac1mag = 8.926 - 2.5*alog10(finalflux[i]*1E-6)
objectnew[i].irac1magerr = objectnew[i].irac1mag - (8.926 - 2.5*alog10((finalflux[i] + finalfluxerr[i])*1E-6))

if object[i].irac1mag lt 90 and object[i].irac1mag gt 0 and object[i].irac1mag - objectnew[i].irac1mag gt 2 then print, 'new mag brighter by 2mags'

print, 'final', i, objectnew[i].irac1flux , objectnew[i].irac1fluxerr, objectnew[i].irac1mag, objectnew[i].irac1magerr, object[i].irac1flux, object[i].irac1mag, objectnew[i].extended

print, 'ratio', (object[i].irac1flux - objectnew[i].irac1flux) / object[i].irac1flux, objectnew[i].irac1mag
endif

endfor

;testdeltamag = testdeltamag[0:g-1]
;testskyerr=testskyerr[0:g-1]
;plot, testdeltamag, testskyerr, psym=8, symsize=0.2, xtitle='orig - new', ytitle='background_sigma', xrange=[0,6], yrange=[0.0009, 0.01], xstyle=1, ystyle=1



; write the aperture mask file
writefits,'mask_total.fits',mask_total,ch1head


; save some results
save,/variables,filename='extract_acs_shape_mask_v2.sav'
;save, objectnew, filename='/Users/jkrick/idlbin/objectnew.sav'

;plot, object[10000:17000].irac1mag, object[10000:17000].irac1mag - objectnew[10000:17000].irac1mag, psym = 3, xrange=[10,30], xtitle='original aperture magnitudes', ytitle='original - new',  charthick=3, xthick=3, ythick=3, title='non acs detected, extended objects'
;j = where(object.acsmag lt 0 or object.acsmag gt 90)
;i = where(object[j].irac1mag gt 0 and object[j].irac1mag lt 90)

;irac extended
;ei1 = where(objectnew[j[i]].extended gt 0)

;oplot, object[j[i[ei1]]].irac1mag, object[j[i[ei1]]].irac1mag - objectnew[j[i[ei1]]].irac1mag, psym = 2
;oplot, object[j[i[pi1]]].irac1mag, object[j[i[pi1]]].irac1mag - objectnew[j[i[pi1]]].irac1mag, psym = 2

print,"finishing at "+systime()


end
