function convert_Kmag_IRAC, Kmag, ch, sp_type

;This is a function to convert 2MASS K magnitude into 3.6 or 4.5
;micron flux density given spectral type.

  ;;use online STAR-PET to determine colors as a function of spectral
  ;;type
  ;;http://ssc.spitzer.caltech.edu/warmmission/propkit/pet/starpet/
  if ch eq '1' then begin
     F0 = 280.9E5
     case sp_type of
        'A0': color = -0.01
        'A5': color = 0.01
        'F0': color = 0.02
        'F5': color = 0.03
        'G0': color = 0.04
        'G5': color = 0.04
        'K0': color = 0.05
        'K5': color = 0.08
        'M0': color = 0.12
        'M2': color = 0.15
        'M3': color = 0.15 ;;me extrapolating
        'M4': color = 0.16;;me extrapolating
        'M5': color = 0.16
     endcase
  endif
  if ch eq '2' then begin
     F0 = 179.7E5
     case sp_type of
        'A0': color = 0.03
        'A5': color = 0.04
        'F0': color = 0.06
        'F5': color = 0.05
        'G0': color = 0.03
        'G5': color = 0.01
        'K0': color = -0.04
        'K5': color = -0.07
        'M0': color = 0.01
        'M2': color = 0.04
        'M3': color = 0.053;;me extrapolating
        'M4': color = 0.066;;me extrapolating
        'M5': color = 0.08
     endcase
  endif
 
  ;change from Kmag into the IRAC bands
  source_mag = Kmag - color

  ;now go from mag to mJy
  source_mJy = F0*10.^(source_mag/(-2.5))

  return, source_mJy /100.
end

