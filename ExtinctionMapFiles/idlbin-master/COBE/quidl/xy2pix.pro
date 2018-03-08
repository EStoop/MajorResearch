function xy2pix,x_in,y_in,res=res,face=face,sixpack=sixpack
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    XY2PIX calculates pixel numbers from input x-y raster coordinates.
;
;DESCRIPTION:                                   
;    This IDL function creates a pixel list given x and y raster 
;    coordinate arrays and the cube resolution.  The function assumes
;    the coordinate arrays are from right oriented skycubes.
;
;CALLING SEQUENCE:
;    pixel = xy2pix(x_in,y_in,res=res,[face=face],[/sixpack])
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    x_in        I     int arr         x-coordinate array
;    y_in        I     int arr         y-coordinate array
;    pixel       O     long arr        Output pixel list
;    res         I     int             Cube resolution
;    face       [I]    int             Face # for single face entry
;    sixpack    [I]    qualifier       Sixpack input format switch
;
;WARNINGS:
;    Only coordinate arrays of Right-T cubes are properly processed.
;    The pixel number is set to -1 if the coordinate is outside of
;    the T.
;
;EXAMPLE: 
;
; To create a FIRAS (res 6) pixel list from x & y raster (skycube/
; sixpack) coordinate arrays:
;
; firas_pixel = xy2pix(x_cube,y_cube,res=6) (unfolded skycube)
; firas_pixel = xy2pix(x_cube,y_cube,res=6,/sixpack) (sixpack)
;
;
; To create a DIRBE (res 9) pixel list from x & y raster (face)
; coordinate arrays:
;
; dirbe_pixel = xy2pix(x_face,y_face,res=9,face=3)
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Determines the face number of the raster coordinate.
;    Uses fij2pix to get the pixel number.
;    Sets pixel number to -1 if coordinate outside T.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: fij2pix
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Apr 92
;
; SPR 9741  Change PRINT statements to MESSAGE statements in
;           error situations.
; 2 JUN 1992
;
; SPR (UIMAGE) 067  Change KEYWORD_SET to N_ELEMENTS to determine
;                   whether FACE switch is present on command line.
;
; 24 JUL 1992
;
; SER 9857  Input x&y arrays can now refer to sixpack format.
; 2 AUG 1992
;
; SPR 9927  Fixed typos in help banner.
; 28 AUG 1992
;
;.TITLE
;Routine XY2PIX
;-
;

cube_side = 2^(res-1)

IF (KEYWORD_SET(sixpack)) THEN BEGIN

	offx = [0,0,1,2,2,1]
	offy = [1,0,0,0,1,1]
		; Note: These are packed Left_t offsets
	cur_2_face = [[3,2,1],[4,5,0]]
	t_len = 3 * cube_side

ENDIF ELSE BEGIN

	offx = [0,0,1,2,3,0]
	offy = [2,1,1,1,1,0]
		; Note: These are left_t offsets
	cur_2_face = [[-1,-1,-1,5],[4,3,2,1],[-1,-1,-1,0]]
		; Right-T map
	t_len = 4 * cube_side

ENDELSE


nx = N_ELEMENTS(x_in)
ny = N_ELEMENTS(y_in)

IF (nx NE ny) THEN BEGIN

	str = 'The number of elements in the input arrays must be the same.'
	MESSAGE,str,/CONT
	pixel = -1
	GOTO, exit
	pixel = -1

ENDIF

pixel = LONARR(nx)

IF (N_ELEMENTS(face) EQ 0) THEN BEGIN	; skycube

	fce = cur_2_face(x_in/cube_side,y_in/cube_side)

	i = x_in - offx(fce) * cube_side
	i = t_len - (i+1) 
		; flip to right-T

	j = y_in - offy(fce) * cube_side

ENDIF ELSE BEGIN			; face

	fce = REPLICATE(face,nx)
	i = cube_side - (x_in+1)
	j = y_in

ENDELSE


pixel = fij2pix([[fce],[i],[j]],res)

pixel = pixel * (fce NE -1) - (fce EQ -1)

exit:

RETURN,pixel
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


