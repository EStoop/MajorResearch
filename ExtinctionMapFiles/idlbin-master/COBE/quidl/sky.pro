pro sky,vimage,skymode,skysig
;+
; NAME:
;   SKY
; PURPOSE:
;   Determine the sky level in an image by applying the procedure MMM
;   to approximately 4000 uniformly spaced pixels.  Adapted from the
;   DAOPHOT routine of the same name.
; CALLING SEQUENCE:
;   SKY,VIMAGE          ;Print the sky background and standard deviation
;   SKY,VIMAGE,SKYMODE,SKYSIG  ;Store output values in SKYMODE and SKYSIG
; INPUTS:
;    VIMAGE - One or two dimensional array
; OPTIONAL OUTPUT ARRAYS:
;    SKYMODE - Scalar, giving the mode of the sky pixel values of the 
;              array VIMAGE, as determined by the procedure MMM.
;    SKYSIG -  Scalar, giving standard deviation of sky brightness
; PROCEDURE:
;    A regular grid of points, not exceeding 4000 in number, is extracted
;    from the image array.  The mode of these pixel values is determined
;    by MMM.
; REVISION HISTORY:
;     Written, W. Landsman   STX Co.            September, 1987     
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;-
maxsky = 4000		;Maximum # of pixels to be used in sky calculation
err = string(7b) + 'SKY: ERROR - '
if n_params(0) eq 0 then begin
	print,string(7B),'CALLING SEQUENCE: sky,vimage,[skymode,skysig]
	return
endif
s = size(vimage)	                          
case s(0) of	;Number of dimensions in array? 
1: begin        ; 1-d vector
   istep = fix(s(1)/float(maxsky) + 0.5) > 1
   npts = s(1)/istep 
   ilow = (s(1) - istep*npts)/2
   index = ilow + indgen(npts)*istep
   end
2: begin
   fractn = float(maxsky)/s(4)   
   istep = fix(0.5 + 1/sqrt(fractn)) > 1
   ilow = (s(1) - istep*((s(1) - 2)/istep) -1) / 2
   jstep = s(2) *(s(1)/istep)/maxsky + 1
   jlow = (s(2) - jstep*((s(2) -2)/jstep) - 1) /2
   n_x = s(1)/istep &  n_y = s(2)/jstep
   npts = n_x*n_y
   nfirst = ilow + jlow*s(1)
   index = nfirst + s(1)*jstep*(indgen(npts)/n_x) + istep*(indgen(npts) mod n_x)
   end 
else: message,'Input array (first parameter) must be 1 or 2 dimensional'
endcase
mmm,vimage(index),skymode,skysig
if !debug then print,'Number of points used to find sky = ',npts
if not !quiet then begin
	print,'Approximate sky value for this frame = ',float(skymode)
	print,'Standard deviation of sky brightness = ',float(skysig)
endif
return
end
;DISCLAIMER:
;
;This software was written at the Cosmology Data Analysis Center in
;support of the Cosmic Background Explorer (COBE) Project under NASA
;contract number NAS5-30750.
;
;This software may be used, copied, modified or redistributed so long
;as it is not sold and this disclaimer is distributed along with the
;software.  If you modify the software please indicate your
;modifications in a prominent place in the source code.  
;
;All routines are provided "as is" without any express or implied
;warranties whatsoever.  All routines are distributed without guarantee
;of support.  If errors are found in this code it is requested that you
;contact us by sending email to the address below to report the errors
;but we make no claims regarding timely fixes.  This software has been 
;used for analysis of COBE data but has not been validated and has not 
;been used to create validated data sets of any type.
;
;Please send bug reports to CGIS@ZWICKY.GSFC.NASA.GOV.


