pro box_erase
;+
; NAME:	
;	BOX_ERASE
; PURPOSE:
;	Erases boxes drawn by BOX_DRAW
; CATEGORY:
;	Display.
; CALLING SEQUENCE:
;	BOX_ERASE
; INPUTS:
;	None.
; OUTPUTS:
;	None.
; COMMON BLOCKS:
;	box
; SIDE EFFECTS:
; RESTRICTIONS:
;	For use with ZOOMW.
;	Intended for TV devices.
; PROCEDURE:
;	Corners of box and original data are obtained from common block: box
; MODIFICATION HISTORY:
;Frank Varosi NASA/GSFC 1989
;-

  common box, Lox,Loy,Hix,Hiy, sav_Horiz_B, sav_Horiz_T, sav_Vert_L, sav_Vert_R

	if N_elements( sav_Horiz_B ) LE 1 then return

	tv, sav_Horiz_B, Lox,Loy
	tv, sav_Vert_L, Lox,Loy
	tv, sav_Horiz_T, Lox,Hiy
	tv, sav_vert_R, Hix,Loy

	sav_Horiz_B = 0
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


