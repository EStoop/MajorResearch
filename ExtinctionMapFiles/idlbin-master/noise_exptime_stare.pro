pro noise_exptime_stare
;make a plot of the noise levels in an exoplanet staring mode observation
;some things are hard coded in so be careful to expand this to other observations

;maybe they should all be bands to show a range of eg. eclipse depths.

  bin_scale = findgen(1000) + 1
  root_n = sqrt(bin_scale)
 

;using WASP-14 data
  fits_read, '/Users/jkrick/irac_warm/KOI069/r41009920/ch2/bcd/SPITZER_I2_41009920_0023_0000_2_bcd.fits', im, h  ;2s sub
  fits_read, '/Users/jkrick/irac_warm/KOI069/r39438336/ch2/bcd/SPITZER_I2_39438336_0023_0000_2_bcd.fits', im, h  ;12s full

;assume aperture size of 2.25 pixels radius
;;code to plot SNR as a function of source magnitude and IRAC exptime

  ;;assume 2s subarray frametimes
  ;;assume ch 2
  ;;assum sp_type
  exptime = 1.92  ; 10.4
  ronoise =  9.4 ; 10.4 for 12s full
  ch = '2'
  sp_type = 'K0'
   
  ;;assume aperture size of 2.25 pixels radius
  apsize = 2.25
  aparea = !Pi*apsize^2
  ;; assuming background annulus [3,7]
  bkdgarea = !Pi*7^2 - !Pi*3^2  

  ;;assuming sky dark median from KOI69 data
  skydkmed = 1.3                ; 0.5 for 12s MJy/sr

  
  ;;basics of the detector
  pixel_scale = 1.22
  case ch of
     '1': begin
        gain =3.7
        flux_conv = .1253
        factor = 2.
     end
     '2': begin
        gain = 3.71
        flux_conv = .1469
        factor = 1.5
     end
  endcase
  
  SNR = fltarr(7)
  kmag = fltarr(7)
;;loop over a range of magnitudes
  for i = 0, 6, 1 do begin
     kmag(i) = i +6.
     
     ;;what is the SNR of this target?
     ;;----------------
     
     ;;calculate flux density of source in IRAC band
     source_mJy = convert_Kmag_IRAC(Kmag(i), ch, sp_type)
     
     ;;start with no binning
     bin_scale = exptime / 60.
     n_frames =  1              ; fix(bin_scale*60./exptime) 

     ;;calculate the number of electrons in the target
     source_electrons = mjy_to_electron( source_mjy, pixel_scale, gain, exptime, flux_conv)
     
     ;;poisson noise
     sigma_poisson = sqrt(source_electrons)
     print, 'sigma poisson', sigma_poisson
     
     ;;readnoise
     NR = sqrt(aparea*(1+aparea/bkdgarea)*ronoise^2)
     print, 'sigma readnoise', NR

     ;; background level
     ;;zodi level, CIB, dark level.
     skydkmed = skydkmed / flux_conv * exptime * gain
     sigma_dark = sqrt(skydkmed)
     ;;change this to be the worse case (high)  background expected for this exptime?
     ;;this does bin down
     B = 32                     ; electrons/s/pixel  from the data handbook
     B = B * exptime            ; now in electrons/pixel
     Nb = sqrt(B*aparea*(1+aparea/bkdgarea))
     print, 'sigma sky dark', sigma_dark, Nb


     ;;total noise
     noise = sqrt(sigma_poisson^2 + NR^2 + Nb^2)
     SNR(i) = source_electrons / noise
  endfor                        ;;loop over magnitudes


  p = plot(kmag, SNR, xtitle = 'Kmag', ytitle = 'SNR')
  
end
