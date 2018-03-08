PRO put_owner_id
;+
; NAME:
;     Put_owner_id		
;
; PURPOSE:
;     Put owner ID and system time on hard-copy plots
;
; CATEGORY:
;     Graphic display.
;
; CALLING SEQUENCE:
;     Put_owner_id
;
; INPUTS:
;     None
;
; KEYWORD PARAMETERS:
;     None
;
; OUTPUTS: 
;     No explicit outputs.
;
; COMMON BLOCKS:
;     None
;
; SIDE EFFECTS:
;     None
;
; RESTRICTIONS:
;     None
;
; PROCEDURE: 
;
;	IF the plot output device is Post Script THEN
;	   WRITE the user ID and system time at the 
;	         left hand bottom corner of the output page
;	ENDIF
;	RETURN
; 
; MODIFICATION HISTORY:
; 
; 	Jan 17, 1992 	Version 0; 
;		        Created by: V. Kumar/T. Hewagama (Hughes STX)
;
; SPR # 9416
;-
;
;IF !D.Name NE 'PS' then return
;
userid=getenv ('USER')
xyouts, 0,0,'Drawn on '+!stime+' by '+userid,/device,charsize=0.5
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


