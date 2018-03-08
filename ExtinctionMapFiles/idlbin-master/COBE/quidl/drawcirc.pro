PRO drawcirc, fij_input, n_pix_arc, i_arc, j_arc, $
	      vec_arc, re_proj, w_pos, face_num, color, zoom=zoom
;+                                                                  
;  NAME:
;    drawcirc
;
;  PURPOSE: To draw a given great circle on a skycube or reprojection
;
;  CALLING SEQUENCE:
;    drawcirc, fij_input, n_pix_arc, i_arc, j_arc, $
;	       vec_arc, re_proj, w_pos, face_num, color
;
;  INPUT:
;    fij_input - (n,3) int array containing face, col (w/i face), and 
;                row (w/i face) of arc.  Only face number is used.
;                (Skycube Overlay Only)
;
;    n_pix_arc - number of points in arc
;
;    i_arc - int array containing column numbers within skycube 
;            (unfolded or sixpack) of arc.
;            (Skycube Overlay Only)
;
;    j_arc - int array containing row numbers within skycube 
;            (unfolded or sixpack) of arc.
;            (Skycube Overlay Only)
;
;    vec_arc - flt array containing unit vectors of arc
;              (Reprojection Overlay Only)
;
;    re_proj - reprojection flag (0 - skycube, 1 - reprojection input)
;
;    w_pos - offset of image within window (2 element int array)
;
;    face_num - face number if single face (-1 whole sky)
;
;    color - color number of overlay line (0-255)
;
;    zoom - image zoom factor
;
;  OUTPUT:
;    None
;
;  SUBROUTINES CALLED:
;    uv2proj
;
;  REVISION HISTORY
;
;  SPR 10476  Add Documentation  J.M. Gales  01/21/93
;-

COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
		    cube_side, proj, coord, sz_proj

IF(NOT KEYWORD_SET(zoom)) THEN zoom = 1

IF (re_proj EQ 0) THEN BEGIN	; skycube overlay

	n_seg = 0
	prev_face = fij_input(0,0)
	face_break = intarr(7)
	; get face number of 0th point in arc
	; generate face_break array

	FOR i=1,n_pix_arc-1 DO BEGIN
		; for all interior points of arc

		IF (fij_input(i,0) NE prev_face) THEN BEGIN
			n_seg = n_seg + 1
			face_break(n_seg) = i
			prev_face = fij_input(i,0)
			; if point in a new face then
			; increment number of segments
			; set point # of break
			; reset previous face #
		ENDIF

	ENDFOR
	face_break(n_seg+1)= n_pix_arc-1
			; set face break at last point

	FOR i=0,n_seg DO BEGIN
	; for all face segments

	IF ((face_num eq -1) OR (fij_input(face_break(i),0) EQ face_num)) $ 
	THEN BEGIN

      PLOTS,(i_arc(face_break(i):face_break(i+1)-1) + w_pos(0))*zoom + zoom/2, $
        (j_arc(face_break(i):face_break(i+1)-1) + w_pos(1))*zoom + zoom/2, $
        PSYM=-3,/DEVICE,LINESTYLE=0,THICK=1,COLOR=color
	; if whole sky or arc point within single face then
	; plot great arc by individual face segments
	; Note: This eliminates lines between disconnected faces

	ENDIF

	ENDFOR

ENDIF ELSE BEGIN	; reprojection overlay

	n_seg = 0
	quad_break = intarr(7)
	sign0 = fix(abs(vec_arc(0,1))/vec_arc(0,1))
	; generate quadrant break array
	; get sign of y-element of 0th point


	FOR i=1,n_pix_arc-1 DO BEGIN
	; for all interior points of arc

	sign1 = fix(abs(vec_arc(i,1))/vec_arc(i,1))
	; get sign of y-element of point

	IF (sign0 NE sign1) THEN BEGIN
		n_seg = n_seg + 1
		quad_break(n_seg) = i
		sign0 = sign1
			; if sign of quadrant changes then
			; increment number of segments
			; set point # of break
			; reset previous quadrant
	ENDIF

	ENDFOR
	quad_break(n_seg+1)= n_pix_arc-1
		; set quad break at last point

	gc = uv2proj(vec_arc,proj,sz_proj)
		; generate projection pixel #'s from arc unit vectors
		
	FOR i=0,n_seg DO BEGIN
	; for all quad segments

	IF (quad_break(i) LT quad_break(i+1)) THEN $

  PLOTS,(gc(quad_break(i):quad_break(i+1)-1,0) + w_pos(0))*zoom + zoom/2, $
    (gc(quad_break(i):quad_break(i+1)-1,1) + w_pos(1))*zoom + zoom/2, $
    PSYM=-3,/DEVICE,LINESTYLE=0,THICK=1,COLOR=color
	; plot great arc by individual quad segments
	; Note: This eliminates lines between disconnected quadrants

	ENDFOR

ENDELSE

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


