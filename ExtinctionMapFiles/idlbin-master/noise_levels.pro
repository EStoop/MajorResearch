pro noise_levels
;make a plot of the noise levels in an exoplanet staring mode observation
;some things are hard coded in so be careful to expand this to other observations

;maybe they should all be bands to show a range of eg. eclipse depths.

  bin_scale = findgen(1000) + 1
  root_n = sqrt(bin_scale)
 

;using WASP-14 data
  fits_read, '/Users/jkrick/irac_warm/pcrs_planets/wasp-14b/r48682752/ch2/bcd/SPITZER_I2_48682752_0013_0000_1_bcd.fits', im, h

;assume aperture size of 2.25 pixels radius
  apsize = 2.25
  aparea = !Pi*apsize^2

  bkdgarea = !Pi*7^2 - !Pi*3^2  ; assuming background annulus [3,7]

  ronoise = sxpar(h, 'RONOISE')  ; electrons
  exptime = sxpar(h, 'EXPTIME')  ;s
  gain = sxpar(h, 'GAIN')   ;e/DN
  fluxconv = sxpar(h, 'FLUXCONV')  ; MJY/sr / DN/s



;-----------------
;poisson noise of the source
;take the average over 64 subframes in a single image
  meane = fltarr(64)
  for j = 0, 63 do begin  
     subframe = im[*,*,j]
     subframe = (subframe / fluxconv) * exptime * gain ; now in electrons
     aper, subframe, 15.5, 15.5, flux, err, back, backerr, 1.0, 3.0, [3,7], $
           /FLUX, /EXACT, /SILENT, /NAN,  /meanback,$
           READNOISE=ronoise
     meane(j) = flux
  endfor
  meanclip, meane, avgelectrons, avgsigma
;  print, 'avgelectrons wasp 14', avgelectrons
  sigma_poisson = sqrt(avgelectrons)
;  print, 'poisson noise wasp 14', sigma_poisson

  
;instead of using a real source, I could just use half-well depth, which is what the SSC recommends
  full_well = 140000 ; electrons  according to the data handbook table 2.2
  avgelectrons = 0.5 * full_well ; in the peak pixel
;now looking at encircled energy plots I have made of HD7924, 54% of the total flux is in the central pixel whereas 87% is within 2.25 pixel aperture radius, so increase the total number of electrons to account for the aperture photometry.
  avgelectrons = avgelectrons*1.6
  sigma_poisson = sqrt(avgelectrons)
  print, 'poisson noise', sigma_poisson

  p = plot( bin_scale, sigma_poisson / root_n, thick = 3, xtitle = 'Binning Scale (N frames)', ytitle = 'Noise (Electrons)', axis_style = 1,/xlog,/ylog, xrange = [1, 1000], margin = 0.2, yrange = [1.0, 20000.])
  t = text(3., 200., 'Source Poisson', color = 'black', /current,/data)
  xaxis = axis('x', location = [1, 18100.], coord_transform = [0, 2/60.],target = p, textpos = 1, title = 'Binning Scale (Min.)')
  xaxis = axis('y', location = [925.,1.], coord_transform = [0, 1.],target = p, tickdir = 0, textpos = 0, showtext = 0)
;  xaxis = axis('x', location = [0, 4.25], coord_transform = [0, 2/60.],target = p, textpos = 1, title = 'Binning Scale (Min.)')
;  xaxis = axis('y', location = [2.9699,0], coord_transform = [0, 1.],target = p, tickdir = 0, textpos = 0, showtext = 0)

;-----------------
;readnoise
  ronoise = 9.4  ;2s
  NR = sqrt(aparea*(1+aparea/bkdgarea)*ronoise^2)
  print, 'sigma readnoise', NR
 
  p = plot(bin_scale, NR / root_n,  thick = 3, color = 'grey',/overplot, axis_style = 1)
  t = text(6., 15., 'Readnoise', color = 'grey', /current,/data)

;-----------------
; background level
;zodi level, CIB, dark level.
  skydkmed = sxpar(h, 'SKYDKMED') ; in MJy/sr
  skydkmed = skydkmed / fluxconv * exptime * gain
  sigma_dark = sqrt(skydkmed)
  print, 'sigma sky dark', sigma_dark
  
  
;change this to be the worse case (high)  background expected for this exptime?
;XXX eventually make it a stripe for low to high background
;this does bin down
  B = 32 ; electrons/s/pixel  from the data handbook
  B = B * exptime  ; now in electrons/pixel
  Nb = sqrt(B*aparea*(1+aparea/bkdgarea))
  p = plot( bin_scale, Nb / root_n, thick = 3, color = 'dark slate grey',/overplot, axis_style = 1)
  t = text(2, 30, 'Bkgd Poisson', color = 'dark slate grey', /current,/data)

;-----------------
;bias pattern = first frame effect 
  bias_pattern = 150.                 ; electrons
  print, 'bias pattern', bias_pattern
  
  bias_pattern = fltarr(n_elements(bin_scale)) + bias_pattern
  p = plot(bin_scale, bias_pattern, thick = 3, linestyle = 2, color = 'slate grey',/overplot, axis_style = 1)
  ;t = text(80., 177., 'Bias Pattern', color = 'slate grey', /current,/data)

;-----------------
;latent build-up, or left over from previous observation?
;from aperture photometry on staring mode darks diff imaging 
  latent = 175                   ; worse of the two examples we have (ch1) electrons
  print, 'latent signal', latent
  latent = fltarr(n_elements(bin_scale)) + latent
  p = plot(bin_scale, latent, thick = 3, linestyle = 2, color = 'slate grey',/overplot, axis_style = 1)
  t = text(50., 200., 'Latent & Bias Pattern', color = 'slate grey', /current,/data)

;-----------------
;pixel phase effect - sawtooth
; percent changes on sweet spot pixel for *warm* mission
  ch1_pixphase = 0.081
  ch2_pixphase = 0.021
  ch1_pixphase_signal = ch1_pixphase*avgelectrons
  ch2_pixphase_signal = ch2_pixphase*avgelectrons
  print, 'ch1 ch2 pixphasesignal', ch1_pixphase_signal, ch2_pixphase_signal

;XXXassuming not binning like root N, but probably does bin down somewhat
  ch2_pixphase_signal= fltarr(n_elements(bin_scale)) + ch2_pixphase_signal
  p = plot(bin_scale, ch2_pixphase_signal, thick = 3, color = 'dark grey',/overplot, axis_style = 1)
  t = text(100., 2800., 'ch2 Pixphase', color = 'dark grey', /current,/data)
  
  ch1_pixphase_signal= fltarr(n_elements(bin_scale)) + ch1_pixphase_signal
  p = plot(bin_scale, ch1_pixphase_signal, thick = 3, color = 'dark grey',/overplot, axis_style = 1)
  t = text(100., 10000., 'ch1 Pixphase', color = 'dark grey', /current,/data)

;-----------------
;residual non-linearity
; does this matter for the pmap dataset being used to correct a staring mode light curve?
; only if the detector effects are a function of flux (like maybe first frame effect) 
;keep track of it anyway.
  linearity = .01               ; we say we calibrate to 1%
;XXX do we have a better estimate of the level of residual nonlinarity?
  linearity_signal = linearity * avgelectrons
  print, 'linearity', linearity_signal


;;;---------------------------------------------------------
;astrophysical signals

;-----------------
;expected signal of the transit
  transit_depth = 0.01          ; assume 1%
  transit_signal = transit_depth*avgelectrons
  print, 'transit signa', transit_signal
  transit_signal= fltarr(n_elements(bin_scale)) + transit_signal
  p = plot(bin_scale, transit_signal, thick = 3, color = 'blue',/overplot, axis_style = 1)
  t = text(200., 600., '1% Transit', color = 'blue', /current,/data)

;-----------------
;expected signal of the eclipse
  eclipse_depth = 0.001         ; assume .1%
  eclipse_signal = eclipse_depth*avgelectrons
  print, 'eclipse signa', eclipse_signal
  eclipse_signal= fltarr(n_elements(bin_scale)) + eclipse_signal
  p = plot(bin_scale, eclipse_signal, thick = 3, color = 'red',/overplot, axis_style = 1)
  t = text(300., 70., 'Eclipse', color = 'red', /current,/data)


;-----------------
;expected signal of the phase variation
  phase_depth = 0.001           ; assume .1%
  phase_signal = phase_depth*avgelectrons
  print, 'phase signa', phase_signal
  phase_signal= fltarr(n_elements(bin_scale)) + phase_signal
  p = plot(bin_scale, phase_signal, thick = 3, color = 'red',/overplot, axis_style = 1)
  t = text(50., 70., '0.1% Phase,', color = 'red', /current,/data)


;-----------------
; stellar variability  and noise on the stellar variation
;there is a huge range here
;maybe I want to make this a big greyscale bar in the plot?
  s_vary = 25                   ; assume 0.1%
  svary_signal = s_vary*avgelectrons
  print, 'stellar variation', svary_signal

;-----------------
;do we have a prediction of an earth-like planet transitting in a habitable zone?
;Gillon et al 2014 for GJ1214 claim 70ppm
;I think that is 0.00007 or 
  superearth_depth = 0.00007          
  superearth_signal = superearth_depth*avgelectrons
  print, 'superearth signa', superearth_signal
  superearth_signal= fltarr(n_elements(bin_scale)) + superearth_signal
  p = plot(bin_scale, superearth_signal, thick = 3, color = 'maroon',/overplot, axis_style = 1)
  t = text(50., 9., '70 ppm Superearth', color = 'maroon', /current,/data)


;save
;p.save,'/Users/jkrick/conferences/IRACaas2014/noiseplot.png'
p.save, '/Users/jkrick/irac_warm/pcrs_planets/secondaries_paper/noiseplot.eps'


end
