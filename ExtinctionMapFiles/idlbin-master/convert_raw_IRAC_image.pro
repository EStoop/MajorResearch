function convert_raw_IRAC_image, raw_image, raw_header,$
     output_2D_for_subarray = output_2D_for_subarray,$
     transform_to_BCD_coordinates = transform_to_BCD_coordinates,$
     no_truncation_correction = no_truncation_correction,$
     minimum_possible_value = minimum_possible_value
; SUMMARY: converts ANY raw IRAC image (as read in by READFITS.PRO) to floating point or long integers, 
;           correcting for the wraparound or the input integers, with rectification of signals into the 
;	   positive domain for channels 1 and 2

; Inputs are the raw image (a DCE read in by READFITS.PRO) and the raw header, plus options
; Default output is a floating-point image (a 3-D stack for subarray DCEs) 

; CHANGE HISTORY:
;  BG 2010-09-17 
;
; CREATED by BG 2010-09-17 from dewrap3.pro
;
; Simplest Usage:
;   im_out=convert_raw_IRAC_image(raw_image,raw_header)   This will output a floating-point array: 2-D for Full-array DCEs
;                                          (including test patterns and signal or pedestal frames),
;                                          or 3-D if the DCE is in subarray mode. By default, the output
;                                          is in raw-array coordinates, and the truncation correction is made.
;
;
; *******************************************************************************************
; INPUTS:
;  raw_image     (required):	            A 256x256 array of integer values, read in as 16 bit signed integers by READFITS.PRO)
;
;  raw_header    (required):	            The DCE header (type STRARR) as read in by READFITS.PRO
;
;  output_2D_for_subarray (optional):       The DEFAULT is 0, which outputs an image stack for subarray DCEs. 
;					    Set it to 1 to output subarray images in the same 2-D arrangement of pixels as the raw image.
;
;  transform_to_BCD_coordinates (optional): If set, im_out will be in BCD coordinates (i.e., ch1, ch2 are row-reversed); 
;                                           except when the DCE is a test pattern, or when the DCE is in subarray mode
;                                           and output_2D_for_subarray is set.
;                                           The DEFAULT is to make the output in raw array coordinates for all but the two exceptions.
;
;  no_truncation_correction (optional):     If set, then the output will not have the theoretical average truncated value added in.
;                                           The DEFAULT is to have the average truncated value added to im_out, except when
;                                           output_long_integers is set.
;                                           This keyword has no effect on Fowler-0 frames, including test patterns; there 
;                                           is no truncation to correct for.
;                                            
;  minimum_possible_value (optional):       The most negative value possible for a normal Fowler frame. DEFAULT: -8000. 
;                                           For warm mission, use minimum_possible_value= -4000. or -5000. in the call.
;                                           For a frame from an IER where fowlnum NE 2^barrel, then this minimum_possible
;                                           values is multiplied by (2^barrel)/fowlnum.
;                                           
;**********************************************************
;
; USAGE:
;  nice_image = convert_raw_IRAC_image( raw_image ,raw_header, output_2D_for_subarray = output_2D_for_subarray, $
;     transform_to_BCD_coordinates = transform_to_BCD_coordinates, $
;     no_truncation_correction = no_truncation_correction, minimum_possible_value = minimum_possible_value
;
;  nice_image = call_function('convert_raw_IRAC_image', raw_image ,raw_header, output_2D_for_subarray = output_2D_for_subarray, $
;     transform_to_BCD_coordinates = transform_to_BCD_coordinates, $
;     no_truncation_correction = no_truncation_correction, minimum_possible_value = minimum_possible_value)
;
; EXAMPLE:
;  nice_image = call_function('convert_raw_IRAC_image', im, raw_header, output_2D_for_subarray = 0, /transform_to_BCD_coordinates, $
;     
; *************************************************************
; NOTE: Remember that when computing a Fowler-equivalent, which is
;   not done in this routine, for channels 1 and 2, it's P - S and
;   for channels 3 and 4 it's S - P.
;
; NOTE: The truncation correction is a correction for the average value of the bits that are lost when
;       the DSPs do a right-shift of the 24-bit values by FOWLNUM bits before the lower 16 bits of the data are 
;       written to the computer memory. The average error due to the truncation in a pixel raw DN is
;       
;          AVERAGE ERROR in DN CAUSED BY TRUNCATION = -0.5 * (1 - 2^(-ABARREL))
;
;       The truncation correction subtracts this average error from every pixel in Fowler or Pedestal or
;       Signal Frames when ABARREL > 0, after the raw image has been converted to floating point, and before the image is
;       rectified into the positive domain (i.e., channel 1 and 2 pixel values multiplied by -1), and
;       before any scaling for ABARREL and FOWLNUM. No correction is made for Fowler-0 frames or test patterns.
;
;       It is not necessary to correct for the truncation, but we must be consistent: In our pipeline processing,
;       we subtract a labdark from the image before flatfielding and scaling. The labdarks were made from DCEs
;       that were corrected for truncation, with a few exceptions (see below). And we have always used the same
;       relation between barrel shift and Fowler number. Therefore the pipeline must correct for truncation as well;
;       otherwise the skydarks (other than the excepted ones) will have mean offsets of between -1/4 (Fowler-1) and 
;       -31/64 (Fowler-32) DN. Since a skydark of the same Frame Type is subtracted in the process of making a BCD image,
;       BCDs are not affected by the truncation correction, as long as the skydark and the BCD are both processed
;       with or without the correction.
;
;       There are 2 AOR modes for which the labdarks have no truncation correction. They are are from the 
;       First-generation Labdarks (with version "v1.2.0" in their file names), which are still being
;       used in S18.14. There was insufficient data to make labdarks. The mean value for each pixel was
;       set to zero, and this is expected to cause a false pattern in the skydarks that is worse than the
;       lack of truncation correction.
;
;          1.  Subarray Mode, 0.02s frames: SUB_0.02s_0.02sf1d1r1_ch1_v1.2.0_dark.fits: The labdark is 0 for every
;              pixel in every subframe. (The pixel uncertainties are scaled (times sqrt(2))from the SUB_0.1s labdark,
;              and are duplicated in each subframe.)
;          2.  Subarray Mode, 0.4s frames: SUB_0.4s_0.4sf8d1r1_ch1_v1.2.0_dark.fits: The labdark is 0 for every
;              pixel in every subframe. (The pixel uncertainties are the SAME as the SUB_0.1s labdark, and are  
;              duplicated in each subframe.)
;
; The following are extracted from the header:
;  ichan	CHNLNUM IRAC channel number 1-4
;  barrel 	ABARREL barrel shift
;  fowlnum 	AFOWLNUM
;  readmode	AREADMOD  0=full-array, 1=subarray
;  pedsig 	APEDSIG 1=pedesal 2= signal 0 otherwise
;  pattern 	ADSPPATT  0=no pattern  1-4 = patterns 1-4
;
;***************************************************************

if arg_present(minimum_possible_value) $
 then min_possible_value=float(minimum_possible_value) $
 else min_possible_value= -8000. 	;lowest value a normal Fowler frame is allowed to have for Cold Mission --
;
; now extract the parameters from the header

ichan = sxpar(raw_header,'CHNLNUM')
barrel = sxpar(raw_header,'A0741D00')
fowlnum = sxpar(raw_header,'A0614D00')
readmode = sxpar(raw_header,'A0617D00')
pedsig = sxpar(raw_header,'A0742D00')
pattern = sxpar(raw_header,'A0744D00')
;
sub = readmode eq 1
;                                       ;For warm mission, use minimum_possible_value= -4000. or -5000. in the call.
if fowlnum ne 2^barrel then min_possible_value=min_possible_value*2.^barrel/fowlnum
;                 lower values will be assumed to be higher than
;                 32767*2.^barrel/fowlnum. The maximum possible value is 
;		  65535*2.^barrel/fowlnum. +  minimum_possible_value*2.^barrel/fowlnum
;		  This is done after correcting for truncation, barrel
;		  shift, and fowler number, and rectification of channels
;		  1 and 2.
;
if keyword_set(no_truncation_correction) then do_trunkcor=0 else do_trunkcor=1
im_out=float(raw_image)
value_resolution_step=1.0
;
if fowlnum eq 0 then begin ;FOWLER-0 and PATTERNs are in range 0-65535
;		patterns are ramp 0 to FFFFH and ramp FFFFH to 0 and
;               checkerboards 5555H (pos) and AAAAH (neg)
;		Fowler-0 or fixed patterns require only the correction for
;			signed integer to unsigned integer
;      
  a=where(im lt 0)
  if a[0] ne -1 then im_out[a]=im_out[a]+65536.
endif else begin
 if pedsig eq 1 then begin   ;PEDESTAL FRAME
   a=where(im gt 0)
   if a[0] ne -1 then im_out[a]=im_out[a]-65536.
   if do_trunkcor then begin	;  correct for bit truncation 
     corr1=0.5*(1.-2.^(-barrel))
     if barrel ne 0 then im_out=im_out+corr1
   endif
;  correct for barrel shift and fowler multiplication
  if fowlnum ne 2^barrel then begin
    im_out=im_out*(2.^barrel/fowlnum)
    value_resolution_step=2.^barrel/fowlnum
  endif
;   negate the pedestal -- it was loaded negative in DSP
  im_out= -im_out	;negate pedestal since it was loaded negative  
 endif else if pedsig eq 2 then begin    ;SIGNAL FRAME
   a=where(im lt 0)
   if a[0] ne -1 then im_out[a]=im_out[a]+65536.
   if do_trunkcor then begin	;  correct for bit truncation 
     corr1=0.5*(1.-2.^(-barrel))
     if barrel ne 0 then im_out=im_out+corr1
   endif
;  correct for barrel shift and fowler multiplication
   if fowlnum ne 2^barrel then begin
     im_out=im_out*(2.^barrel/fowlnum)
     value_resolution_step=2.^barrel/fowlnum
   endif  
 endif else begin	;FOWLER FRAME
   if do_trunkcor then begin	;  correct for bit truncation 
     corr1=0.5*(1.-2.^(-barrel))
     if barrel ne 0 then im_out=im_out+corr1
   endif
; negate the values for channels 1 and 2
   if ichan eq 1 or ichan eq 2 then im_out= -im_out
; now remove the assumed wraparound
   a=where(im_out lt min_possible_value)
   if a[0] ne -1 then im_out[a]=im_out[a]+65536.
;  correct for barrel shift and fowler multiplication
   if fowlnum ne 2^barrel then begin
     im_out=im_out*(2.^barrel/fowlnum)
     value_resolution_step=2.^barrel/fowlnum
   endif
 endelse
endelse
;
; Now convert a subarray image to a 32col x 32row x 64subframe stack unless output_2D_for_subarray is set.
;
if sub and not(keyword_set(output_2D_for_subarray)) then im_out=reform(im_out,32,32,64)
;
;
; Now convert the image or stack of subframes to BCD coordinates if desired.
;
if ( keyword_set(transform_to_BCD_coordinates) and (ichan eq 1 or ichan eq 2) ) and $
    pattern eq 0  and ( not(sub) or (sub and not(keyword_set(output_2D_for_subarray))) ) then begin
  temp=im_out
  if sub then begin
    for k=0,63 do begin
      for j=0,31 do im_out[*,j,k]=temp[*,31-j,k]  ;reverse rows in each subframe
    endfor
  endif else begin
    for j=0,255 do im_out[*,j]=temp[*,255-j]    ;reverse row order on imout (ch1, ch2 only)
  endelse
endif
;
return,im_out
end
