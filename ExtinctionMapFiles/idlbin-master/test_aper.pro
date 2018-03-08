pro test_aper
  ;these two have latents in the center
  fits_read, '/Users/jkrick/irac_warm/darks/staring/diffim_44201728.fits', im, h
  fits_read, '/Users/jkrick/irac_warm/darks/staring/diffim_40152064.fits', im, h
  ;this one is blank
;  fits_read, '/Users/jkrick/irac_warm/darks/staring/diffim_42051584.fits', im, h

;this one has negative column band centered at ~x = 15.2
  fits_read, '/Users/jkrick/irac_warm/darks/staring/diffim_40151552.fits', im, h

;this one has positive column band centered at ~x = [28,32]
  fits_read, '/Users/jkrick/irac_warm/darks/staring/diffim_42506496.fits', im, h


  meanarr = fltarr(64)
  columnarr = fltarr(64)
for j = 0, 63 do begin  
   subframe = im[*,*,j]
  
  ronoise = sxpar(h, 'RONOISE')
  exptime = sxpar(h, 'EXPTIME')
  gain = sxpar(h, 'GAIN')
  fluxconv = sxpar(h, 'FLUXCONV')

  ronoise = 8.9 ; 0.4s
  exptime = 0.36
  gain = 3.7
  fluxconv = 0.1469

  subframe = (subframe / fluxconv) * exptime * gain  ; now in electrons
  badpix = [-9., 9.] * 1.D8

  aper, subframe, 29, 10.2, flux, err, back, backerr, 1.0, 3.0, [3,10], $
         /FLUX, /EXACT, /SILENT, /NAN,  /meanback,$
        READNOISE=ronoise

  meanarr(j) = flux


                                ;run 'photometry' on whole column set
                                ;this one has negative column band centered at ~x = 15.2
                                ;but this needs to be aperture photometry to see what the impact is on aperture photomtry of sources.
  on = total(subframe[13:16,*])
  off = total(subframe[23:26,*])
  
  columnarr(j) = on - off


endfor
meanclip, meanarr, meanmean, meansigma
print,'aperture mean', meanmean

;meanclip, columnarr, meanc, sigmacprint, 'column mean', meanc
end
