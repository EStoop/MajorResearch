; PURPOSE:
;  This code calculates the predicted SNR of a transit given the 2MASS
;  K band magnitude of the source, the spectral type of the source, the
;  IRAC channel of interest, and the depth of the transit. The code
;  assumes a relation for how the noise of the source decreases with
;  increasing number of frames which is based on an empiracle relation
;  taken from the literature.  After data reduction, on average the
;  literature is finding that the poisson noise limit is routinely
;  achieved for ch1 on binning scales less than 30 frames and for ch2 on
;  binning scales less than 100 frames.  Ch1 data are able to achieve
;  2 times the poisson noise on binning scales greater than 30 frames
;  and a factor of 1.5 times poisson on binning scales greater than
;  100 frames for ch2.
;
; INPUTS:
;  - Kmag = Magnitude in the 2MASS K band
;  - sp_type = spectral type.  Options are ('A0', 'A5','F0','F5','G0','G5', 'K0','K5','M0','M2','M5')
;  - ch = the name of the channel for observations and in which
;         source_mJy is given.  Can be either the string '1' or '2' corresponding to 3.6
;         and 4.5 microns
;  - transit_depth in percentage of the stellar signal
;  - bin_scale = time on which to bin in minutes (can be int, float,
;    etc.)
;
; OPTIONAL INPUTS
;  - eclipse
;
; OUTPUTS:
;  SNR = Signal to Noise Ratio
; 
; EXAMPLE:
;  SNR = exoplanet_snr(12.0, 'M0','2', 1E-3)
;  SNR = exoplanet_snr(12.0, 'M0','2', 1E-3,/eclipse, AB = 1.0, Rp =
;  0.001, Rs = 1.0, semi_maj = .001)
;
; REQUIRED FUNCTIONS:
;  convert_Kmag_IRAC, find_exptime, mjy_to_electron, sigfig
;
; MODIFICATION HISTORY:
;  January 2015 Original Version  JK
;-
Function exoplanet_SNR, Kmag, sp_type, ch, transit_depth, bin_scale, eclipse = eclipse, AB = AB, Rp = Rp, Rs = Rs, semi_maj = semi_maj

  ;;error checking
  if ch ne '1' and ch ne '2' and ch ne 1 and ch ne 2 then begin
     print, 'ch should be either "1" or "2" where "1" is 3.6microns and "2" is 4.5 microns'
     return, 0
  endif
  if bin_scale gt 90 then print, 'Are you sure you want to bin on a', bin_scale, ' Minute timescale'

  ;;estimate eclipse depth or transit depth SNR?
  if keyword_set(eclipse) then begin

     ;;would need AB, Rp, Rs, semi_maj as inputs
     eclipse_depth = (AB/4)*((Rp/Rs)^2)*((Rs/semi_maj)^2)

     ;;if this keyword is set then I actually want to comput the SNR
     ;;of the predicted eclipse depth and not the transity depth                              
     transit_depth = eclipse_depth
  endif

  ;;calculate flux density of source in IRAC band
  source_mJy = convert_Kmag_IRAC(Kmag, ch, sp_type)


  ;;need to calculate frame time best suited for this target
  exptime = find_exptime(source_mJy, ch)
  print, 'exptime', exptime
  ;;given frame time; how many frames = input binning scale 
  n_frames = fix(bin_scale*60./exptime) 

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


  ;;calculate the number of electrons in the target
  source_electrons = mjy_to_electron( source_mjy, pixel_scale, gain, exptime, flux_conv)
  
  ;;poisson noise
  sigma_poisson = sqrt(source_electrons)

  ;;noise binns down like some factor times root n
  ;;this is an empiracle relation based on literature values
  root_n = sqrt(n_frames)
  if n_frames lt 100 then begin
     noise = sigma_poisson / root_n
  endif else begin
     noise = sigma_poisson / (factor * root_n)
  endelse

  ;;noise as percentage of the source flux
  noise = noise / source_electrons

  ;;signal to noise ratio
  SNR = transit_depth / noise
  SNR = sigfig(SNR,3) 

;-------------------------------------------


RETURN, SNR
end
