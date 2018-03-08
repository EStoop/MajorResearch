FUNCTION projtran,input,i_arr,j_arr,mask,face_num,min
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    projtran - Builds projection image
;
;DESCRIPTION:
;    This IDL procedure actually creates reprojections of the native
;    skycubes.  It is called by JPRO.  There is a C module with the
;    same name with the same functionality.  If the linkimage is not
;    implemented, then this IDL module is compiled and run instead.
;
;CALLING SEQUENCE:
;    proj_img = projtran(input,i_arr,j_arr,mask,face_num,min)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    input       I/O   flt arr             Sky Cube/Face (Res 6-9)
;    proj_img    O     flt arr             Rasterized output array
;    i_arr       I     int arr             column # array
;    j_arr       I     int arr             row # array
;    mask        I     byte arr            mask array
;    face        I     int                 face number (for face inpt only)
;    min         I     flt                 Image Scale Minimum
;
;WARNINGS:
;    If displaying an image, an X-windows terminal must be used.
;
;EXAMPLE: 
;
; Not user routine
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Build projected image.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: None
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Oct 92
;    Delivered 13-NOV-1992  SPR 10212
;-

sz = SIZE(input)
input_l = FIX(sz(1))
input_h = FIX(sz(2))
		; get size OF input array


IF (3*sz(1) EQ 4*sz(2)) THEN input_type = 'CUBE'
IF (2*sz(1) EQ 3*sz(2)) THEN input_type = 'SIX'
IF (sz(1) EQ sz(2)) THEN input_type = 'FACE'
		; determine whether unfolded cube, sixpack or single face


IF (input_type EQ 'CUBE' OR input_type EQ 'SIX') THEN $

		proj_img = input(i_arr,j_arr) * (mask LT 255) + $
	        	   min * (mask EQ 255) $

ELSE $

		proj_img = input(i_arr,j_arr) * (mask EQ face_num) + $
		           min * (mask NE face_num)



RETURN,proj_img
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


