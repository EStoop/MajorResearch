pro extract_acs_shape_mask_v2

; read in external data
; things to know:
;  daoflux and flux are on the same index system
; sex and acs have their own indices, matched via _match
;
restore,'catalog_nosystem.sav'
restore,'isostar.sav'
image=readfits('mosaic_middle.fits',hd)
psf=readfits('../massimo/I1_psf_x4_best.fits',psf_hd)
resid=readfits('substar9.fits',resid_hd)
resid=abs(resid)
naxis1=sxpar(hd,'NAXIS1')
naxis2=sxpar(hd,'NAXIS2')




; create arrays to hold output
npoints=n_elements(x)
iracflux_acs=fltarr(npoints)
iracap_acs=fltarr(npoints)

; read a bunch of variables in and match
acsfwhm=xcal(object.acsfwhm,acs_match)
isoarea=xcal(object.acsisoarea,acs_match)
acsclassstar=xcal(object.acsclassstar,acs_match)
acsellip=xcal(object.acsellip,acs_match)
acstheta=xcal(object.acstheta,acs_match)
ch1=xcal(object.irac1flux,acs_match)
ch2=xcal(object.irac2flux,acs_match)
ch1ch2=ch1/ch2
isostar_match=xcal(isostar,acs_match)
;star=where((acsclassstar EQ 1)and (ch1ch2 GT 1.4) and (ch1ch2 LT 1.8))
star=where(isostar_match GT 0)
acsratio=1/(1-acsellip) ; convert sextrator ellipticity to ratio of smaj/smin



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; set the aperture size                                          ;
; ap size from gaia is sqrt(isoarea/3.141)), or, duh, the radius ;
; semi_major and minor axes derived from ellipse                  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ap=sqrt(isoarea/3.141)
semi_minor=sqrt(isoarea/(3.141*acsratio))
semi_major=acsratio*semi_minor
; convert that to IRAC mosaic units
ap=ap*0.05/0.6
semi_major=semi_major*0.05/0.6
semi_minor=semi_minor*0.05/0.6
; now set the minimum aperture
ap(where(semi_major LT 3.2))=3.2
ap(where(finite(ap) NE 1))=3.2
ap(star)=3.2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; compute the fiducial aperture correction ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,500,20,[40,50],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0

fiducial=pflux

aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,500,6.4,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0

fiducial_cor=pflux/fiducial
print,total(psf)/fiducial
print,fiducial_cor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; loop through all points and measure ;
; for all using variable circular ap  ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for i=0,npoints-1 do begin

aper,image,x[i],y[i],tflux,tfluxerr,tsky,tskyerr,500,ap[i],[ap[i]+10,ap[i]+20],[-100,1000],/nan,/exact,/flux,/silent

tflux=tflux*8.46 ; convert to uJy


; derive the appropriate aperture correction
apcor=fiducial_cor ; default aperture
if ((ap[i] GT 3.2) and (ap[i] LT 100))then begin
     aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,500,ap[i]*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0

     apcor=pflux/fiducial

endif

if (ap[i] GT 100) then apcor=total(psf)/fiducial


tflux=tflux/apcor ; apply aperture correction

; register the fluxes and apertures
iracflux_acs[i]=tflux
iracap_acs[i]=ap[i]

endfor

plot,flux,iracflux_acs/flux,psym=3,/xlog,yrange=[0,2],xrange=[1,1e4]


;;;;;;;;;;;;;;;;;;;;;;;;;
; now for the craziness ;
;;;;;;;;;;;;;;;;;;;;;;;;;
print,'Processing '+'elliptical masks.'

kron_factor=2.5
min_radius=2.0 ; minimum radius for computing an ellipse (3.2 default)

mask=image  ; initialize the mask to be the same as the image
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


for i=0,npoints-1 do begin


; if extended, measure with mask
if ((semi_major_orig[i] GT min_radius)) then begin

; first mask the aperture itself

dist_ellipse,mask_ellipse,[naxis1,naxis2],x[i],y[i],acsratio[i],90+acstheta[i]
mask[*]=0
mask(where(mask_ellipse LE kron_factor*semi_major[i]))=1
mask_total=mask_total OR mask

mask_object=where(mask EQ 1)
mask_object_area=total(mask)
mask_object_flux=total(image(mask_object))
resid_object_flux=total(resid(mask_object))

; now process the sky
mask[*]=0
mask(where((mask_ellipse LE kron_factor*semi_major[i]+15)and(mask_ellipse GE kron_factor*semi_major[i]+5)))=1
mask_sky=where(mask EQ 1)
mmm,image(mask_sky),background,background_sigma,background_skew
mmm,resid(mask_sky),resid_background,resid_background_sigma,resid_background_skew
;background=median(image(mask_sky))  ; median background
mask_sky_flux=background*mask_object_area
resid_mask_sky_flux=resid_background*mask_object_area

mask_flux[i]=8.46*(mask_object_flux-mask_sky_flux)
resid_flux[i]=8.46*(resid_object_flux-resid_mask_sky_flux)

; derive the appropriate aperture correction

; first, recreate the effective aperture
mask_radius=sqrt(mask_object_area / 3.141)

if (mask_radius LT 100) then begin
      aper,psf,521,521,pflux,pfluxerr,psky,pskyerr,500,mask_radius*2,[10,20],[-100,1e5],/nan,/exact,/silent,/flux,setskyval=0

     apcor=pflux/fiducial

endif else begin
     apcor=total(psf)/fiducial
endelse

;print,apcor
mask_flux[i]=mask_flux[i]/apcor
resid_flux[i]=resid_flux[i]/apcor

endif


endfor

; write the aperture mask file
writefits,'mask_total.fits',mask_total,hd

window,xsize=320,ysize=240

sexbest=xcal(sex_best,sex_match)
ratio=mask_flux/sexbest
plot,flux,ratio,psym=3,/xlog,yrange=[0,2],xrange=[1,1e4]



; save some results
save,/variables,filename='extract_acs_shape_mask_v2.sav'
stop


end
