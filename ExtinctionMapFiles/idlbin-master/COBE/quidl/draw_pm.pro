pro draw_pm,para,merd,re_proj,face_num,w_pos,res_input,gcoord, $
	    color=color,zoom=zoom,subset=subset
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    draw_pm calculates and draws the coordinates of the grid lines
;
;DESCRIPTION:
;    This IDL procedure calculates the coordinate grid lines 
;    (meridians and parallels) on either unfolded skycubes or 
;    reprojected images.  It is called by JPRO and is not a user
;    routine.
;
;CALLING SEQUENCE:
;    draw_pm,para,merd,re_proj,face_num,w_pos,res_input, $
;	     zoom=zoom,subset=subset
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    para        I     int arr             parallels
;    merd        I     int arr             meridians
;    re_proj     I     int                 reprojection flag
;    face_num    I     int                 face number (-1 if cube)
;    w_pos       I     int arr             pos of proj in window
;    res_input   I     int                 input resolution
;    color       I     int                 grid line color
;    zoom        I     int                 zoom factor
;    subset      I     int arr             subsetting indices [x0,x1,y0,y1]
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
;    image_parms
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Calculate sky coordinates of parallels.
;    If sky cube projection, determine device pixel coordinates.
;    Call draw_grid_line to draw parallels.
;    Repeat procedure for meridians.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    draw_grid_line
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   May 92
;  SPR 10325  Mar 05 93  Coordinate grids on zoomed images.  J Ewing
;
;  SPR 12144  Mar 16 95  Change call to DRAW_GRID_LINE  J.M. Gales
;-

COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
		    cube_side, proj, coord, sz_proj

d2r = !dpi / 180
minus_1 = -1
out_res = STRCOMPRESS('R' + STRING(res_input),REMOVE_ALL=1)

IF (KEYWORD_SET(color) EQ 0) THEN color = 255

; Draw Parallels
; --------------
vec_arc = fltarr(361,3)

sz = SIZE(para)
IF (sz(0) EQ 0) THEN no_para = 1 ELSE BEGIN 
	no_para = 0
	para = FLOAT(para)
	lon_rad = d2r*(indgen(361)-180+0.5)
	lon_rad(0)   = d2r * (-180 + 1.e-3)
	lon_rad(360) = d2r * (+180 - 1.e-3)
ENDELSE

FOR j=1,2*n_elements(para)*(no_para EQ 0) DO BEGIN

	k = (j-1) mod 2

	IF (para((j-1)/2) EQ 0) THEN lat_rad = (-1)^k * 1.e-6 ELSE $
	lat_rad = d2r * para((j-1) / 2)

	IF ((k EQ 0) OR (para((j-1)/2) EQ 0)) THEN BEGIN

		vec_arc(*,0) = cos(lon_rad)*cos(lat_rad)
		vec_arc(*,1) = sin(lon_rad)*cos(lat_rad)
		vec_arc(*,2) = replicate(sin(lat_rad),361)

	ENDIF ELSE vec_arc(*,2) = -vec_arc(*,2)

	n_pix_arc = 361

	IF (re_proj eq 0) THEN BEGIN

		pixel_arc = coorconv(vec_arc,infmt='U',inco=gcoord, $
				     outfmt='P',outco=out_res)

		fij = pix2fij(pixel_arc,res_input)

		i_pix = offx(fij(*,0)) * cube_side + fij(*,1)
		i_pix = input_l - (i_pix+1)
		j_pix = offy(fij(*,0)) * cube_side + fij(*,2)

		IF (face_num NE -1) THEN BEGIN
			i_pix = i_pix * (fij(*,0) EQ face_num)
			j_pix = j_pix * (fij(*,0) EQ face_num)
			i_pix = i_pix - w_pos(0) * (fij(*,0) NE face_num)
			j_pix = j_pix - w_pos(1) * (fij(*,0) NE face_num)
		ENDIF

          IF(KEYWORD_SET(subset)) THEN BEGIN
            w = WHERE((i_pix LT subset(0)) or (i_pix GT subset(1)), nw)
            IF(nw GT 0) THEN i_pix(w) = -1000
            w = WHERE((j_pix LT subset(2)) or (j_pix GT subset(3)), nw)
            IF(nw GT 0) THEN j_pix(w) = -1000
          ENDIF

		cube_input = [[fij(*,0)],[i_pix],[j_pix]]

	ENDIF ELSE $

		vec_arc2 = coorconv(vec_arc,infmt='U',inco=gcoord,outco=coord)

	draw_grid_line,cube_input,vec_arc2,proj,sz_proj,w_pos,color, $
		       zoom=zoom

ENDFOR



; Draw meridians
; --------------
vec_arc = fltarr(181,3)

sz = SIZE(merd)
IF (sz(0) EQ 0) THEN no_merd = 1 ELSE BEGIN 
	no_merd = 0
	merd = FLOAT(merd)
	lat_rad = d2r*(indgen(181)-90+0.5)
	lat_rad(0)   = d2r * (-90 + 1.e-6)
	lat_rad(180) = d2r * (+90 - 1.e-6)
ENDELSE

FOR j=1,2*n_elements(merd)*(no_merd EQ 0) DO BEGIN

	k = (j-1) mod 2

	IF (merd((j-1)/2) EQ 0) THEN lon_rad = 1.e-6 ELSE $
	lon_rad = d2r * merd((j-1) / 2)

	IF ((k EQ 0) OR (merd((j-1)/2) EQ 0)) THEN BEGIN

		vec_arc(*,0) = cos(lon_rad)*cos(lat_rad)
		vec_arc(*,1) = sin(lon_rad)*cos(lat_rad)
		vec_arc(*,2) = sin(lat_rad)

	ENDIF ELSE vec_arc(*,1) = -vec_arc(*,1)

	n_pix_arc = 181

	n = INDGEN(n_pix_arc/2 + 1)
	s = n + (n_pix_arc/2)


	IF (re_proj eq 0) THEN BEGIN
	
		pixel_arc = coorconv(vec_arc,infmt='U',inco=gcoord, $
				     outfmt='P',outco=out_res)

		fij = pix2fij(pixel_arc,res_input)

		i_pix = offx(fij(*,0)) * cube_side + fij(*,1)
		i_pix = input_l - (i_pix+1)
		j_pix = offy(fij(*,0)) * cube_side + fij(*,2)

		IF (face_num NE -1) THEN BEGIN
			i_pix = i_pix * (fij(*,0) EQ face_num)
			j_pix = j_pix * (fij(*,0) EQ face_num)
			i_pix = i_pix - w_pos(0) * (fij(*,0) NE face_num)
			j_pix = j_pix - w_pos(1) * (fij(*,0) NE face_num)
		ENDIF

          IF(KEYWORD_SET(subset)) THEN BEGIN
            w = WHERE((i_pix LT subset(0)) or (i_pix GT subset(1)), nw)
            IF(nw GT 0) THEN i_pix(w) = -1000
            w = WHERE((j_pix LT subset(2)) or (j_pix GT subset(3)), nw)
            IF(nw GT 0) THEN j_pix(w) = -1000
          ENDIF

		cube_input = [[fij(*,0)],[i_pix],[j_pix]]

	ENDIF ELSE $

		vec_arc2 = coorconv(vec_arc,infmt='U',inco=gcoord,outco=coord)


	IF (re_proj eq 0) THEN BEGIN

		draw_grid_line,cube_input(n,*),vec_arc2,proj,sz_proj, $
			       w_pos,color,zoom=zoom
		draw_grid_line,cube_input(s,*),vec_arc2,proj,sz_proj, $
			       w_pos,color,zoom=zoom

	ENDIF ELSE BEGIN

		draw_grid_line,cube_input,vec_arc2(n,*),proj,sz_proj, $
			       w_pos,color,zoom=zoom
		draw_grid_line,cube_input,vec_arc2(s,*),proj,sz_proj, $
			       w_pos,color,zoom=zoom

	ENDELSE

ENDFOR

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


