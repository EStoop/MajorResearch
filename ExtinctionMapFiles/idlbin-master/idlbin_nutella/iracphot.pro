pro iracphot

readcol,'/Users/jkrick/hst/raw/wholeacs.cat', NUMBER,X_WORLD,Y_WORLD,X_IMAGE,Y_IMAGE,MAG_AUTO,FLUX_BEST,FLUXERR_BEST,MAG_BEST,MAGERR_BEST,BACKGROUND,FLUX_MAX,ISOAREA_IMAGE,ALPHA_J2000,DELTA_J2000,THETA_IMAGE,MU_THRESHOLD,MU_MAX,FLAGS,FWHM_IMAGE,CLASS_STAR,ELLIPTICITY, format="A"

acsratio=1./(1.-ellipticity) ; convert sextrator ellipticity to ratio of smaj/smin

readcol, '/Users/jkrick/spitzer/irac/SExtractor.1.cat', ch1NUMBER,ch1X_WORLD      ,ch1Y_WORLD      ,ch1X_IMAGE,ch1Y_IMAGE,ch1FLUX_AUTO,ch1MAG_AUTO,ch1MAGERR_AUTO,ch1FLUX_APER,ch1MAG_APER,ch1MAGERR_APER,ch1FLUX_BEST,ch1MAG_BEST,ch1MAGERR_BEST,ch1FWHM_IMAGE,ch1ISOAREA_IMAGE,ch1FLUX_MAX,ch1ELLIPTICITY,ch1CLASS_STAR,ch1FLAGS, format="A"

ch1ratio = 1./(1.-ch1ellipticity)

;read in irac image
fits_read, '/Users/jkrick/spitzer/irac/ch1/mosaic.fits', ch1data, ch1head  
fits_read, '/Users/jkrick/spitzer/irac/ch1/mosaic_cov.fits', ch1cov, ch1covhead  
naxis1=sxpar(ch1head,'NAXIS1')
naxis2=sxpar(ch1head,'NAXIS2')
psf=readfits('/Users/jason/psf/massimo/I1_psf_x4_best.fits',psf_hd)
restore, '/Users/jkrick/idlbin/object.sav'

;convert acs ra and dec to irac x, y
adxy, ch1head, x_world, y_world, ch1x, ch1y

; create arrays to hold output
npoints=n_elements(ch1x)
npoints=1000
iracflux_acs=fltarr(npoints)
iracap_acs=fltarr(npoints)
iracfluxerr = fltarr(npoints)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set the aperture size                                          ;
; ap size from gaia is sqrt(isoarea/3.141)), or, duh, the radius ;
; semi_major and minor axes derived from ellipse                  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ap=sqrt(isoarea_image/3.141)
semi_minor=sqrt(isoarea_image/(3.141*acsratio))
semi_major=acsratio*semi_minor
; convert that to IRAC mosaic units
ap=ap*0.05/0.6
semi_major=semi_major*0.05/0.6
semi_minor=semi_minor*0.05/0.6
; now set the minimum aperture
ap(where(semi_major LT 3.2))=3.2
;ap(where(finite(ap) NE 1))=3.2

;make sure the stars get the circular apertures
; define the bounding polygonal area for stars 
px=[12.5,25, 25, 10, 12.5] ; acsmag 
py=[6e5 ,60, 35, 3e5, 6e5] ; ascisoarea 
py=alog10(py) ; switch to log space 
roi = Obj_New('IDLanROI', px, py) 
star=roi->ContainsPoints(mag_best,alog10(isoarea_image)) 
temp=where(star GT 0) 
ap(temp)=3.2   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; compute the fiducial aperture correction ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,500,20,[40,50],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0

fiducial=pflux

aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,500,6.4,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0

fiducial_cor=pflux/fiducial
;print,total(psf)/fiducial
;print,fiducial_cor

;
;;;;;;;;;;;;;;;;;;;;;;;;;
; now for the craziness ;
;;;;;;;;;;;;;;;;;;;;;;;;;

kron_factor=2.5
min_radius=2.0 ; minimum radius for computing an ellipse (3.2 default)

mask=ch1data  ; initialize the mask to be the same as the image
mask[*]=0
mask_total=mask
mask_flux=fltarr(npoints)
mask_flux[*]=sqrt(-1)
resid_flux=mask_flux


; this fixes the ratio value to give a minimum semi_minor axis
; must be done _after_ semi-major axis is computed above
semi_major_orig=semi_major
semi_major(where(semi_major LT 3.2))=3.2
acsratio(where(semi_minor LT min_radius))=semi_major(where(semi_minor LT min_radius))/min_radius

gain = [3.3]; ch1, ch2-4=, 3.7, 3.8, 3.8]

for i=0, npoints-1 do begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; figure out distance to nearest neighbor
;and the brightness and isoarea of that neighbor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   dist=sphdist( x_world[i], y_world[i],object.ra,object.dec,/degrees)
   sorteddist = dist(sort(dist))
   b = sort(dist)
   nearestdist = sorteddist(1)
   nearestflux  = object[b(1)].acsflux
   nearestisoarea = object[b(1)].acsisoarea
   ;why isn't there an irac isoarea in the catalog?


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; calculate effective gain based on coverage map
;effective gain = N * gain
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   N= ch1cov[ch1x[i], ch1y[i]]
   phpadu = N*gain
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; for all objects, measure a circular aperture
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   aper,ch1data,ch1x[i],ch1y[i],tflux,tfluxerr,tsky,tskyerr,phpadu,ap[i],[ap[i]+10,ap[i]+20],[-100,1000],/nan,/exact,/flux,/silent

   tflux=tflux*8.46             ; convert to uJy


; derive the appropriate aperture correction
   apcor=fiducial_cor           ; default aperture
   if ((ap[i] GT 3.2) and (ap[i] LT 100))then begin
      aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,gain,ap[i]*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0

      apcor=pflux/fiducial

   endif

   if (ap[i] GT 100) then apcor=total(psf)/fiducial


   tflux=tflux/apcor            ; apply aperture correction

; register the fluxes and apertures
   iracflux_acs[i]=tflux
   iracap_acs[i]=ap[i]
   iracfluxerr[i] = tfluxerr

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; if extended, measure with mask
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   if ((semi_major_orig[i] GT min_radius)) then begin


; first mask the aperture itself

      dist_ellipse,mask_ellipse,[naxis1,naxis2],ch1x[i],ch1y[i],acsratio[i],90+theta_image[i]
      mask[*]=0
      mask(where(mask_ellipse LE kron_factor*semi_major[i]))=1
      mask_total=mask_total OR mask

      mask_object=where(mask EQ 1)
      mask_object_area=total(mask)
      mask_object_flux=total(ch1data(mask_object))
;      resid_object_flux=total(resid(mask_object))

; now process the sky
      mask[*]=0
      mask(where((mask_ellipse LE kron_factor*semi_major[i]+15)and(mask_ellipse GE kron_factor*semi_major[i]+5)))=1
      mask_sky=where(mask EQ 1)
      mmm,ch1data(mask_sky),background,background_sigma,background_skew
;      mmm,resid(mask_sky),resid_background,resid_background_sigma,resid_background_skew
;background=median(image(mask_sky))  ; median background
      mask_sky_flux=background*mask_object_area
;      resid_mask_sky_flux=resid_background*mask_object_area

      mask_flux[i]=8.46*(mask_object_flux-mask_sky_flux)
;      resid_flux[i]=8.46*(resid_object_flux-resid_mask_sky_flux)

; derive the appropriate aperture correction

; first, recreate the effective aperture
      mask_radius=sqrt(mask_object_area / 3.141)

      if (mask_radius LT 100) then begin
         aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,gain, mask_radius*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0

         apcor=pflux/fiducial

      endif else begin
         apcor=total(psf)/fiducial
      endelse

;print,apcor
      mask_flux[i]=mask_flux[i]/apcor
;      resid_flux[i]=resid_flux[i]/apcor

      print, i, x_world(i), y_world(i), isoarea_image(i), ellipticity(i), semi_major_orig(i), semi_major(i), semi_minor(i) , iracap_acs(i), iracflux_acs(i), iracfluxerr(i), semi_major[i], mask_flux[i], nearestdist, nearestflux, nearestisoarea

   endif


endfor

; write the aperture mask file
writefits,'mask_total.fits',mask_total,ch1head



 


; save some results
save,/variables,filename='extract_acs_shape_mask_v2.sav'
end
