pro box_draw, x1,y1, x2,y2, color
;+
; NAME:	
;	BOX_DRAW
; PURPOSE:
;	Draws a box on the display.
; CATEGORY:
;	Display.
; CALLING SEQUENCE:
;	box_draw, x1, y1, x2, y2 [, color]
; INPUTS:
;	x1,y1 = lower left corner of the box
;	x2,y2 = upper right corner of the box
;	optional:
;	color = color used to draw box
; OUTPUTS:
;	none
; COMMON BLOCKS:
;	box
; SIDE EFFECTS:
; RESTRICTIONS:
;	For use with ZOOMW.
;	Intended for TV devices.
; PROCEDURE:
;	Straightforward.
; MODIFICATION HISTORY:	
;Frank Varosi NASA/GSFC 1989
;-

  common box, Lox,Loy,Hix,Hiy, sav_Horiz_B, sav_Horiz_T, sav_Vert_L, sav_Vert_R

	Lox = (x1 < x2)
	Loy = (y1 < y2)

	Hix = (x1 > x2)
	Hiy = (y1 > y2)

	xsiz = Hix - Lox +1
	ysiz = Hiy - Loy +1

	sav_Horiz_B = tvrd( Lox,Loy, xsiz,1 )
	sav_Horiz_T = tvrd( Lox,Hiy, xsiz,1 )

	sav_Vert_L = tvrd( Lox,Loy, 1,ysiz )
	sav_Vert_R = tvrd( Hix,Loy, 1,ysiz )

	if N_elements( color ) NE 1 then  color = !D.n_colors-1
	color = byte( color )

	horiz = replicate( color, xsiz, 1 )
	vert = replicate( color, 1, ysiz )

	tv, horiz, Lox,Loy
	tv, vert, Lox,Loy
	tv, horiz, Lox,Hiy
	tv, vert, Hix,Loy

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


