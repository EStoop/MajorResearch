FUNCTION get_lut, proj_type,coord_type,psize,i_arr,j_arr,mask,face
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    get_lut in projection lookup arrays
;
;DESCRIPTION:
;    This IDL function reads in the projection lookup arrays.  These
;    are the column array (i_arr), the row array (j_arr) and the mask
;    array.  It is called by JPRO and the skycut routine, getinpar and
;    is not a user routine.
;
;CALLING SEQUENCE:
;    get_lut, proj_type,coord_type,psize,i_arr,j_arr,mask
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    proj_type   I     str                 ('A','S','M')
;    coord_type  I     str                 ('E','G','Q')
;    psize       I     str                 ('L','S')  (from JPRO)
;                                          ('+','-')  (from getinpar)
;    i_arr       O     int arr             skycube column array
;    j_arr       O     int arr             skycube row array
;    mask        O     int arr             skycube mask array
;    face              int                 face number
;    retstat     O     int                 return status
;                                          (0-good, 1-error)
;
;WARNINGS:
;    None
;
;EXAMPLE: 
;
;    None.  Not user routine.
;
;#
;COMMON BLOCKS:
;    See code below
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Build lookup table array string.
;    Restore appropriate saveset if arrays not in memory.
;    Fill output arrays from saveset arrays if change in projection.
;    (Note: when called from getinpar then always read mask then exit)
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    None
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   May 92
;
; SER 9956   Improved Documentation/Banner
; 9-SEP-92   J.M. Gales
;
; SPR 11162  Change name to 'GET_LUT'
; 19-JUL-93  J.M. Gales
;
; SPR 11412  28 Oct 93  Add big_gs_ecl,gal reprojections. J. Newmark
;
; SPR 12144  16 Mar 95  Add face number to command line
;                       Add LAST_FACE variable to LAST common block
;                       Rebuild i_arr, j_arr, mask arrays if face changes
; J.M. Gales
;
;-
COMMON ait_ecl, ait_ecl_i,ait_ecl_j,ait_ecl_mask
COMMON ait_gal, ait_gal_i,ait_gal_j,ait_gal_mask
COMMON ait_equ, ait_equ_i,ait_equ_j,ait_equ_mask

COMMON gs_ecl, gs_ecl_i,gs_ecl_j,gs_ecl_mask
COMMON gs_gal, gs_gal_i,gs_gal_j,gs_gal_mask
COMMON gs_equ, gs_equ_i,gs_equ_j,gs_equ_mask

COMMON mwd_ecl, mwd_ecl_i,mwd_ecl_j,mwd_ecl_mask
COMMON mwd_gal, mwd_gal_i,mwd_gal_j,mwd_gal_mask
COMMON mwd_equ, mwd_equ_i,mwd_equ_j,mwd_equ_mask

COMMON big_ait_ecl, big_ait_ecl_i,big_ait_ecl_j,big_ait_ecl_msk
COMMON big_ait_gal, big_ait_gal_i,big_ait_gal_j,big_ait_gal_msk

COMMON big_mwd_ecl, big_mwd_ecl_i,big_mwd_ecl_j,big_mwd_ecl_msk
COMMON big_mwd_gal, big_mwd_gal_i,big_mwd_gal_j,big_mwd_gal_msk

COMMON big_gs_ecl, big_gs_ecl_i,big_gs_ecl_j,big_gs_ecl_msk
COMMON big_gs_gal, big_gs_gal_i,big_gs_gal_j,big_gs_gal_msk

COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
		    cube_side, proj, coord, sz_proj

COMMON last, last_l, last_h, last_proj, last_face, last_flag

ret_stat = 0
dir_str = getenv('CGIS_DATA')

CASE proj_type OF

	'A' :	p_str = 'ait'
	'S' :	p_str = 'gs'
	'M' :	p_str = 'mwd'

	ELSE :	BEGIN
		MESSAGE,'Unknown Projection Type: ' + proj_type,/CONT
		ret_stat = 1
		GOTO, exit
		END

ENDCASE

CASE coord_type OF

	'E' :	c_str = 'ecl'
	'G' :	c_str = 'gal'
	'Q' :	c_str = 'equ'

	ELSE :	BEGIN
		MESSAGE,'Unknown Coordinate System: ' + coord_type,/CONT
		ret_stat = 1
		GOTO, exit
		END

ENDCASE

IF ((psize EQ 'L' OR psize EQ '+') AND (coord_type EQ 'Q')) THEN BEGIN

	MESSAGE,'Lookup Table does not exist for this Projection/Coordinate',$
	/CONT
	ret_stat = 1
	GOTO, exit

ENDIF

tbl_str = p_str + "_" + c_str

IF (psize EQ 'L' OR psize EQ '+') THEN tbl_str = 'big_' + tbl_str

i = EXECUTE("j = N_ELEMENTS(" + tbl_str + "_i)")

IF (j EQ 0) THEN BEGIN

	PRINT, 'Restoring ' + STRUPCASE(tbl_str) + ' Lookup Table'
	i = EXECUTE("RESTORE,'" + dir_str + tbl_str + ".lut'")

ENDIF


IF (psize EQ '+' OR psize EQ '-') THEN BEGIN
	CASE psize OF	
		'-': i = EXECUTE("mask = " + tbl_str + "_mask")
		'+': i = EXECUTE("mask = " + tbl_str + "_msk")
	ENDCASE
	GOTO, exit
ENDIF
		; if call from skycut facility just get mask and exit


IF ((last_proj EQ tbl_str) AND $
    (last_face EQ face) AND $
    (last_h EQ input_h) AND $
    (last_l EQ input_l)) THEN BEGIN
	last_flag = 0
	GOTO, exit
ENDIF ELSE BEGIN
	last_proj = tbl_str
	last_face = face
	last_h = input_h
	last_l = input_l
	last_flag = -1
ENDELSE

i = EXECUTE("i_arr = " + tbl_str + "_i")
i = EXECUTE("j_arr = " + tbl_str + "_j")

CASE psize OF	
	'S': i = EXECUTE("mask = " + tbl_str + "_mask")
	'L': i = EXECUTE("mask = " + tbl_str + "_msk")
ENDCASE

exit:

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


