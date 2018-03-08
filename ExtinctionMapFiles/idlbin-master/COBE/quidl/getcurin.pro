FUNCTION getcurin, icur, jcur, mask, w_pos, color, zoom=zoom
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    getcurin gets cursor positions for skycut facility
;
;DESCRIPTION:
;    This IDL function returns the cursor position in image window
;    coordinates of the endpoints defining a great circle arc.
;
;CALLING SEQUENCE:
;    ret_stat = getcurin(icur, jcur, mask, w_pos, color, zoom=zoom)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    icur        O     int arr (2 elem)    cursor column positions
;    jcur        O     int arr (2 elem)    cursor row positions
;    mask        I     byte arr            reprojection mask
;                                          (if absent then skycube
;                                           input)
;    w_pos       I     int arr (2 elem)    window offset array
;    color       I     int                 cursor mark color
;
;    ret_stat    O     int                 return status
;                                          (0 - OK, 1 - error)
;    zoom       [I]    int                 zoom factor (1 by default)
;
;WARNINGS:
;
;    None.
;
;EXAMPLES: 
;
;    Not user routine.
;
;#
;COMMON BLOCKS:
;    imageparms
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Get first cursor position
;    If on active image then store
;    Get second cursor position
;    If on active image then store
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    None
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Jun 92
;
; SER 9955     Improved Documentation/Banner
; 9-SEP-92     J.M. Gales
;
; SPR 10285    Work with zoomed images (accept a ZOOM parameter).
; 30-NOV-92    J.A. Ewing
; SPR 10974  24 MAy 93 Added error message if point outside map. J Newmark.
;-
COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
		    cube_side, proj, coord, sz_proj
IF(NOT KEYWORD_SET(zoom)) THEN zoom = 1.
ret_stat = 0

sz_mask = SIZE(mask)

; Set cursor limits
; -----------------
IF (sz_mask(0) EQ 0) THEN BEGIN	; skycube

	i_max = input_l - 1
	j_max = input_h - 1

ENDIF ELSE BEGIN		; reprojection

	i_max = sz_mask(1) - 1
	j_max = sz_mask(2) - 1

ENDELSE

icur = intarr(2)
jcur = intarr(2)
	; allocate cursor position arrays

PRINT,'Mark first point with cursor'

first:
cursor,ic,jc,/down,/device
	; get cursor position

plots,ic,jc,psym=1,/device,COLOR=color
	; plot cross at cursor position

ic = (ic / zoom) - w_pos(0)
jc = (jc / zoom) - w_pos(1)
	; subtract window offsets

CASE !err OF

1:	BEGIN		; left button

	IF (ic LT 0 OR jc LT 0 OR ic GE i_max OR jc GE j_max) $
	THEN BEGIN
             PRINT,'Point selected is outside image limits'
             GOTO, first
	ENDIF
         	; if outside image limits then reject, get new pos

	IF (sz_mask(0) EQ 0) THEN BEGIN		; skycube

		IF (cur_2_face(ic/cube_side,jc/cube_side) EQ -1) $
		THEN BEGIN
                  PRINT,'Point selected is outside image limits'
                  GOTO, first
         	ENDIF
		; if not on image then reject, get new pos

	ENDIF ELSE IF (mask(ic,jc) EQ 255) THEN BEGIN
                    PRINT,'Point selected is outside image limits'
                    GOTO, first
	           ENDIF
		; if reprojection and outside image then get new pos

	icur(0) = ic
	jcur(0) = jc
		; store first point cursor positions

	END

2:	BEGIN		; middle button
	GOTO, first
	END

ELSE:	BEGIN		; right button
	ret_stat = 1
	RETURN, ret_stat
		; exit procedure
	END

ENDCASE


PRINT,'Mark second point with cursor'

second:
cursor,ic,jc,/down,/device

plots,ic,jc,psym=1,/device,COLOR=color

ic = (ic / zoom) - w_pos(0)
jc = (jc / zoom) - w_pos(1)

CASE !err OF

1:	BEGIN

	IF (ic LT 0 OR jc LT 0 OR ic GE i_max OR jc GE j_max) $
	THEN BEGIN
             PRINT,'Point selected is outside image limits'
             GOTO, second
	ENDIF

	IF (sz_mask(0) EQ 0) THEN BEGIN

		IF (cur_2_face(ic/cube_side,jc/cube_side) EQ -1) THEN $
		 BEGIN
                   PRINT,'Point selected is outside image limits'
                   GOTO, second
	         ENDIF

	ENDIF ELSE IF (mask(ic,jc) EQ 255) THEN BEGIN
                     PRINT,'Point selected is outside image limits'
                     GOTO, second
           	   ENDIF

	icur(1) = ic
	jcur(1) = jc

	END

2:	BEGIN

	icur = icur(0)
	jcur = jcur(0)

	END

ELSE:	BEGIN
	ret_stat = 1
	RETURN, ret_stat
	END

ENDCASE

exit:
PRINT,'Processing Cursor Entry'

RETURN, ret_stat
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


