PRO jpro,  input,proj=proj_type,coord=coord_type,orient=orient,win=win, $
	  face=face_num,min=min,max=max,para=para,merd=merd, $
	  gcoord=gcoord,noshow=noshow,psize=psize,proj_img
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    jpro - Produces reprojection arrays and coordinate overlays
;
;DESCRIPTION:
;    This IDL procedure actually creates reprojections of the native
;    skycubes.  It is called by REPROJ which acts as a shell around it.
;    For more information, consult the REPROJ UHELP banner.  The 
;    skycubes may be in either unfolded or sixpacked format.
;
;CALLING SEQUENCE:
;    Should not be called directly.  Use REPROJ.
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    input       I/O   flt arr             Sky Cube/Face (Res 6-9)
;    [proj_img]  O     flt arr             Rasterized output array
;    [proj]      I     string              Projection Type ('Q','A','S','M')
;    [coord]     I     string              Coordinate System ('E','G','Q')
;    [psize]     I     string              Projection Size ('S','L')
;    [face]      I     int                 face number (for face inpt only)
;    [win]       I     int                 Window # for display (default:0)
;    [merd]      I     flt arr             Meridian (longitude) array
;    [para]      I     flt arr             Parallel (latitude) array
;    [min]       I     flt                 Image Scale Minimum
;    [max]       I     flt                 Image Scale Maximum
;    [noshow]    I     qualifier           no display proj switch
;
;WARNINGS:
;    If displaying an image, an X-windows terminal must be used.
;
;EXAMPLE: 
;
; See examples in REPROJ.
;
;#
;COMMON BLOCKS:
;    IMAGE_PARMS, LAST, PROJ_ARRAY
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Get lookup tables using get_lut
;    Scale to proper resolution.
;    Flip to right T.
;    If sixpack then adjust lookup tables.
;    If face then adjust lookup tables.
;    Build projected image.
;    Draw parallels and meridians.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: get_lut, draw_pm, projtran
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Apr 92
;
; SPR 9771 (Enhancement)
; 
; If projection type, coordinate system, input size, and face
; number does not change for consecutive reprojections, mapping
; arrays, I_ARR & J_ARR, are not recalculated.
;
; 17 JUN 1992		JMG
;
; SER 9957    Allow sixpack input capabilities
; 9-SEP-92    J.M. Gales
;
; SPR 10043   '/no_c' command line switch to allow overide of
;             cproj call
; 24-SEP-92   J.M. Gales
;
; SPR 10212   Generate projections in call to PROJTRAN
;             PROJTRAN exists as both an IDL procedure and c module
; 13-NOV-92   J.M. Gales
;
; SPR 11162   Call to GET_PROJ_ARRAYS changed to GET_LUT
; 19-JUL-93   J.M. Gales
;
; SPR 12144   Add last_face variable to last COMMON
; 16-MAR-95   Fixes problem with lut rebuild for face change
;             J.M. Gales
;
;.TITLE
;Routine JPRO
;-
COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
		    cube_side, proj, coord, sz_proj

COMMON last, last_l, last_h, last_proj, last_face, last_flag

COMMON proj_array, i_arr, j_arr, mask


bim = 0


IF (n_elements(face_num) EQ 0) THEN face_num = -1
IF (keyword_set(orient) EQ 0) THEN orient = 'R'
IF (keyword_set(psize) EQ 0) THEN psize = 'S'
IF (keyword_set(noshow) EQ 0) THEN noshow = 0
IF (keyword_set(win) EQ 0) THEN win = 0

orient = STRUPCASE(strmid(orient,0,1))
psize = STRUPCASE(strmid(psize,0,1))
IF (N_ELEMENTS(coord_type) NE 0) THEN coord = STRUPCASE(strmid(coord_type,0,1))
IF (noshow EQ 0) THEN gcoord = STRUPCASE(strmid(gcoord,0,1))

color = 255

IF (N_ELEMENTS(last_l) EQ 0) THEN last_l = 0
IF (N_ELEMENTS(last_h) EQ 0) THEN last_h = 0
IF (N_ELEMENTS(last_face) EQ 0) THEN last_face = -1
IF (N_ELEMENTS(last_proj) EQ 0) THEN last_proj = ''

IF (N_ELEMENTS(min) EQ 0) THEN min = min(input)
IF (N_ELEMENTS(max) EQ 0) THEN max = max(input)
		; set min/max for projection

IF (N_ELEMENTS(merd) EQ 0) THEN merd=[0,60,120,180]
IF (N_ELEMENTS(para) EQ 0) THEN para=[0,30,60]


sz = SIZE(input)
input_l = FIX(sz(1))
input_h = FIX(sz(2))
		; get size OF input array


IF (3*sz(1) EQ 4*sz(2)) THEN BEGIN

	input_type = 'CUBE'
	reduce_factor = FIX(1024 / sz(1))
	cube_side = sz(1) / 4

	offx = [0,0,1,2,3,0]
	offy = [2,1,1,1,1,0]
		; left-T offsets

ENDIF ELSE IF (2*sz(1) EQ 3*sz(2)) THEN BEGIN

	input_type = 'SIX'
	reduce_factor = FIX(768 / sz(1))
	cube_side = sz(1) / 3

	offx = [0,0,1,2,2,1]
	offy = [1,0,0,0,1,1]
		; left-T offsets

ENDIF ELSE IF (sz(1) EQ sz(2)) THEN BEGIN

	IF (face_num EQ -1) THEN READ, 'Enter Face Number: ',face_num

	input_type = 'FACE'
	reduce_factor = FIX(256 / sz(1))
	cube_side = sz(1)

	offx = [0,0,0,0,0,0]
	offy = [0,0,0,0,0,0]

ENDIF
		; determine whether unfolded cube, sixpack or single face

FOR bit=0,15 DO BEGIN
	IF ((cube_side XOR 2^bit) EQ 0) THEN BEGIN
		res_input = bit+1
		GOTO, lbl_1
	ENDIF
ENDFOR
MESSAGE,'Improper Image Size',/cont
GOTO, exit

lbl_1:

IF (keyword_set(proj_type)) THEN BEGIN
	
	proj_type = STRUPCASE(strmid(proj_type,0,1))

	proj = proj_type
	re_proj = 1

	ret_stat = get_lut(proj_type,coord,psize,i_arr,j_arr,mask,face_num)
		; load projection lookup tables IF necessary

	IF (ret_stat NE 0) THEN RETURN
		; IF error in lookup, exit program

	sz_proj = SIZE(mask)

ENDIF ELSE BEGIN

	re_proj = 0

	IF (noshow NE 1) THEN BEGIN
		bim = BYTSCL(input,min = min,max = max)
		w_pos = [32,32]
		WINDOW,win,xsize=sz(1)+64,ysize=sz(2)+64,retain=2
		TV, bim, w_pos(0),w_pos(1)
		GOTO, draw
	ENDIF

ENDELSE


IF ((last_flag EQ -1) AND (reduce_factor GT 1)) THEN BEGIN

	i_arr = i_arr/reduce_factor
	j_arr = j_arr/reduce_factor

ENDIF
		; scale down lookup table values IF res < 9


IF ((last_flag EQ -1) AND (input_type EQ 'SIX')) THEN BEGIN

	PRINT,'Adjusting Lookup Tables to Sixpack Format'
	six_offx = [1,1,1,1,0,2]
	six_offy = [1,1,1,1,0,-1]

	FOR j=0,5 DO BEGIN
		i = where(mask EQ j)

		i_arr(i) = i_arr(i) - cube_side * six_offx(j)
		j_arr(i) = j_arr(i) - cube_side * six_offy(j)
	ENDFOR

ENDIF


IF ((last_flag EQ -1) AND (input_type EQ 'FACE')) THEN BEGIN

	p_offx = ((orient EQ 'L') * [0,0,1,2,3,0]) + $
		 ((orient EQ 'R') * [3,3,2,1,0,3])

	p_offy = [2,1,1,1,1,0]

	i_arr = i_arr - p_offx(face_num)*input_l
	j_arr = j_arr - p_offy(face_num)*input_l
			; subtract offsets

	i_arr = i_arr * (mask EQ face_num)
	j_arr = j_arr * (mask EQ face_num)
			; zero out all but desired face

ENDIF


PRINT, 'Building Projection'
proj_img = projtran(input,i_arr,j_arr,mask,face_num,min)

IF (noshow NE 1) THEN BEGIN
	bim = BYTSCL(proj_img,min = min,max = max)
	WINDOW,win,xsize=sz_proj(1)+64, ysize=sz_proj(2)+64,retain=2
	w_pos = [32,32]
	TV, bim, w_pos(0),w_pos(1)
ENDIF

draw:

IF (noshow NE 1) THEN draw_pm,para,merd,re_proj,face_num,w_pos, $
                              res_input,gcoord,color=color

exit:
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


