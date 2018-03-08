pro textclose,textout=textout
;+
; NAME:
;      TEXTCLOSE                   
; PURPOSE:
;      procedure to close file for text output as specifed
;      by !textout.
; CALLING SEQUENCE:
;      textclose
; KEYWORDS:
;      textout - Indicates output device that was used by
;                TEXTOPEN.
; SIDE EFFECTS:
;	if !textout is not equal to 5 and the textunit is
;	opened. Then unit !textunit is closed and released
; HISTORY:
;	D. Lindler  Dec. 1986  (Replaces PRTOPEN)
; spr 10448
;-
;-----------------------------------------------------------
; CLOSE PROPER UNIT
;
IF (n_elements(textout) eq 0) then textout = !textout ;use default

IF textout NE 5 THEN begin
	if !textunit ne 0 then begin
		free_lun,!textunit
                !textunit = 0
	end
end
RETURN
END
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


