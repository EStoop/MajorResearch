FUNCTION getinpar, input,nat_cube,res,re_proj,psize,mask,face_num
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    getinpar gets input data parameters
;
;DESCRIPTION:
;    This IDL function determines the input data parameters for use
;    in the SKYCUT facility.  It first checks the primary data input
;    (single face, unfolded skycube, sixpack cube, or reprojection).
;    In the last case, getinpar determines whether a native data cube 
;    auxiliary input exist or whether the cut will be taken directly 
;    from the reprojection.  The re_proj variable is zero in the first
;    three cases, one in the last if the auxiliary input exits and two
;    if it does not.  Also returned is the face_num (-1 if not single
;    face), the face offsets, the cursor to face translation matrix,
;    and the resolution.
;
;CALLING SEQUENCE:
;    ret_stat = getinpar(input,nat_cube,res,re_proj,psize,mask,face_num)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    input       I     flt arr             primary input array
;    [nat_cube]  I     flt arr             native cube input
;    res         O     int                 resolution
;    re_proj     O     int                 reprojection flag
;    psize       O     str                 '+' large, '-' small
;    mask        O     byte arr            projection mask
;    face_num    O     int                 face number
;    ret_stat    O     int                 return status
;                                          (0 - OK, 1 - error)
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
;    Determine input array type
;    If reprojection, determine whether aux native cube exists
;    Get cube resolution
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    None
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Jun 92
;
; SER 9955    Improved Documentation/Banner
; 9-SEP-92    J.M. Gales
;
; SPR 10317   Set resolution to 9 for reprojected input
; 07-DEC-92   J.M. Gales
;
; SPR 11162   Call to GET_PROJ_ARRAYS changed to GET_LUT
; 19-JUL-93   J.M. Gales
;
;-
COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
		    cube_side, proj, coord, sz_proj


ret_stat = 0

sz = SIZE(input)

input_l = sz(1)
input_h = sz(2)

IF (3*input_l EQ 4*input_h) THEN BEGIN	; unfolded skycube

	re_proj = 0
	cube_side = input_l/4
	face_num = -1
	offx = [0,0,1,2,3,0]
	offy = [2,1,1,1,1,0]

	cur_2_face = [[-1,-1,-1,5],[4,3,2,1],[-1,-1,-1,0]]

ENDIF


IF (2*input_l EQ 3*input_h) THEN BEGIN	; sixpacked skycube

	re_proj = 0
	cube_side = input_l/3
	face_num = -1
	offx = [0,0,1,2,2,1]
	offy = [1,0,0,0,1,1]

	cur_2_face = [[3,2,1],[4,5,0]]

ENDIF


IF (input_l EQ input_h) THEN BEGIN	; single face

	IF (N_ELEMENTS(face_num) EQ 0) THEN $
	   READ,'Please enter face number: ', face_num
	face_num = FIX(face_num)
	re_proj = 0
	cube_side = input_l

	offx = [0,0,0,0,0,0]
	offy = [0,0,0,0,0,0]
	cur_2_face = [face_num]

ENDIF


IF (input_l EQ 2*input_h) THEN BEGIN	; reprojected image

	sz_proj = SIZE(input)
	IF (sz_proj(1) EQ 512) THEN psize = '-' ELSE psize = '+'
	ret_stat = get_lut(proj,coord,psize,iarr,jarr,mask)
				; get mask


	IF (N_ELEMENTS(nat_cube) NE 0) THEN BEGIN  ; native cube present

		re_proj = 1	; native cube present
		sz_data = SIZE(nat_cube)
		input_l = sz_data(1)

		IF (sz_data(1) EQ sz_data(2)) THEN BEGIN
				; single face
			offx = [0,0,0,0,0,0]
			offy = [0,0,0,0,0,0]
			IF (N_ELEMENTS(face_num) EQ 0) THEN $
			READ,'Please enter face number: ', face_num
			cube_side = input_l

		ENDIF ELSE IF (3*sz_data(1) EQ 4*sz_data(2)) THEN BEGIN
				; unfolded sky cube
			face_num = -1
			offx = [0,0,1,2,3,0]
			offy = [2,1,1,1,1,0]
			cube_side = input_l / 4

		ENDIF ELSE IF (2*sz_data(1) EQ 3*sz_data(2)) THEN BEGIN
				; six pack
			face_num = -1
			offx = [0,0,1,2,2,1]
			offy = [1,0,0,0,1,1]
			cube_side = input_l / 3


		ENDIF

	ENDIF ELSE BEGIN	; skycut on reprojection
		re_proj = 2	; native cube not present
		face_num = -1
		res = 9
		GOTO, exit
	ENDELSE

ENDIF

; Determine resolution of quad cube
; ---------------------------------
FOR bit=0,15 DO BEGIN
	IF ((cube_side XOR 2^bit) EQ 0) THEN BEGIN
		res = bit+1
		GOTO, exit
	ENDIF
ENDFOR
MESSAGE,'Improper Image Size',/cont
ret_stat = 1

exit: 


RETURN,ret_stat
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


