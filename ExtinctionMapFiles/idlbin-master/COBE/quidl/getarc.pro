FUNCTION getarc, entry_mode,input,nat_cube,mask,res_input,face_num,w_pos, $
                 re_proj,arc_type,lon_flag,lat_flag,noinpt,color,$
                 zoom=zoom
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    getarc retrieves the coordinates along a skycut arc
;
;DESCRIPTION:
;    This IDL function queries the user for the endpoints of the sky
;    cut, either via the cursor or the long/lat, then determines
;    the coordinate of this arc in both pixel and unit vector format.
;
;CALLING SEQUENCE:
;    ret_stat = getarc(entry_mode,input,nat_cube,mask,res_input, $
;                      face_num,w_pos,re_proj,arc_type,lon_flag, $
;                      lat_flag,noinpt,color,zoom=1)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    entry_mode  I     int                 Entry mode (1-cur, 2-lon/lat)
;    input       I     flt arr             skycube/face/projection input
;    nat_cube    I     flt arr             Native skycube/face input
;                                          (unfolded or sixpack)
;    mask        I     byte arr            reprojection mask
;    res_input   I     int                 skymap resolution
;    face_num    I     int                 face number (face inpt only)
;    w_pos       I     int arr (2 elem)    window offset array
;    re_proj     I     int                 Reprojection value
;                                          (0-cube input, 1-reproj
;                                           input w/cube aux,
;                                           2-direct reproj cut)
;    arc_type    I     str                 'S' - short arc < 180
;                                          'L' - long arc > 180
;    lon_flag    O     int                 -1 if longitude cut
;    lat_flag    O     int                 -1 if latitude cut
;    noinpt      I     int                 no input flag
;                                          -1 if no user query
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
;    imageparms, arc
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Get arc endpoints
;    Get arc coordinates using GRTCIRC
;    Build window title
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    coorconv,proj2uv,fij2pix,grtcirc,validnum
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Jun 92
;
; SER 9955     Improved Documentation/Banner
; 9-SEP-92     J.M. Gales
;
; SPR 10254    Fix misleading prompt for long/lat
; 18-NOV-1992  JM Gales
;
; SPR 10262    Call VALIDNUM to catch invalid numerical input-strings.
; 20-NOV-1992  J.A. Ewing
;
; SPR 10285    Pass a zoom parameter to GETCURIN.
; 30-NOV-1992  J.A. Ewing
;
;-
COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
		    cube_side, proj, coord, sz_proj

COMMON arc, x_plot,y_plot,win_title,fij_input,n_pix_arc,i_arc,j_arc, $
            vec_arc,end_pixs
IF(NOT KEYWORD_SET(zoom)) THEN zoom = 1
arc_deg = 0.5
d2r = !dpi / 180
ret_stat = 0

out_res = STRCOMPRESS('R' + STRING(res_input),REMOVE_ALL=1)

IF (noinpt EQ 1) THEN BEGIN
	noinpt = 0
	GOTO, lbl
ENDIF			; skip entry section


; Cursor Entry Section
; --------------------
IF (entry_mode EQ 1) THEN BEGIN

	ret_stat = getcurin(icur,jcur,mask,w_pos,color,zoom=zoom)
	IF (ret_stat NE 0) THEN GOTO, exit
				; get cursor input

	lat_flag = 0
	lon_flag = 0

	IF (re_proj EQ 0) THEN BEGIN	; skycube/face

		cur_face = cur_2_face(icur/cube_side,jcur/cube_side)
		i_input = icur - offx(cur_face) * cube_side
		i_input = input_l - (i_input+1) 	; flip to RT
		j_input = jcur - offy(cur_face) * cube_side
		end_pixs = fij2pix([[cur_face],[i_input],[j_input]],res_input)
				; determine face number
				; determine x&y raster coordinate
				; determine pixel # of endpoints

	ENDIF ELSE BEGIN		; reprojection

		first_pnt  = proj2uv([[icur(0)],[jcur(0)]],proj,sz_proj)
		second_pnt  = proj2uv([[icur(1)],[jcur(1)]],proj,sz_proj)
				; determine unit vectors

		first_pix = coorconv(first_pnt,infmt='U',outfmt='P', $
				     inco=coord,outco=out_res)

		second_pix = coorconv(second_pnt,infmt='U',outfmt='P', $
				      inco=coord,outco=out_res)
				; convert to pixel number

		end_pixs = [first_pix,second_pix]

	ENDELSE

ENDIF ELSE BEGIN


; Lon/Lat Entry Section
; ---------------------
	first_pnt = ' '
	second_pnt = ' '
	lat_flag = 0
	lon_flag = 0

readfirst:
	read,'Enter longitude & latitude of first point (Ex: 30,25):  ', $
	first_pnt

	first_pnt = STRTRIM(first_pnt,2)
	len = STRLEN(first_pnt)
			; trim leading and trailing blanks
			; get string length

	FOR i=0,len-1 DO IF (STRMID(first_pnt,i,1) EQ ',') THEN $
	first_pnt = STRMID(first_pnt,0,i) + ' ' + $
		    STRMID(first_pnt,i+1,len-i-1)
			; replace ',' with ' '

	spc = STRPOS(first_pnt,' ')
			; get space position

	IF (spc NE -1) THEN BEGIN
                str1 = STRMID(first_pnt,0,spc)
                str2 = STRTRIM(STRMID(first_pnt,spc,len-spc),2)
                IF(NOT (validnum(str1) AND validnum(str2))) THEN BEGIN
                  MESSAGE, 'Invalid numerical input was supplied.', /CONT
                  GOTO, readfirst
                ENDIF
		first_pnt = [FLOAT(str1), FLOAT(str2)]

		IF (first_pnt(0) EQ 0.0) THEN first_pnt(0) = 1.e-6
		IF (first_pnt(1) EQ 0.0) THEN first_pnt(1) = 1.e-6

readsecond:
		read,'Enter longitude & latitude of second point:  ',second_pnt
	ENDIF ELSE BEGIN	; longitude cut
                IF(NOT validnum(first_pnt)) THEN BEGIN
                  MESSAGE, 'Invalid numerical input was supplied.', /CONT
                  GOTO, readfirst
                ENDIF
		vec_arc = fltarr(361,3)
		lon_rad = d2r * FLOAT(first_pnt)
		IF (lon_rad EQ 0.0) THEN lon_rad = d2r*1.e-6

		lat_rad = d2r*(0.5*indgen(361)-90)
		lat_rad(0)   = d2r * (-90 + 1.e-6)
		lat_rad(180) = d2r * (1.e-6)
		lat_rad(360) = d2r * (+90 - 1.e-6)

		vec_arc(*,0) = cos(lon_rad)*cos(lat_rad)
		vec_arc(*,1) = sin(lon_rad)*cos(lat_rad)
		vec_arc(*,2) = sin(lat_rad)

		arc_actl = 0.5	
		pixel_arc = coorconv(vec_arc,infmt='U',inco=coord, $
				     outfmt='P',outco=out_res)

		lon_flag = -1
		GOTO, lbl
	ENDELSE


	second_pnt = STRTRIM(second_pnt,2)
	len = STRLEN(second_pnt)

	FOR i=0,len-1 DO IF (STRMID(second_pnt,i,1) EQ ',') THEN $
	second_pnt = STRMID(second_pnt,0,i) + ' ' + $
		    STRMID(second_pnt,i+1,len-i-1)
				; replace ',' with ' '

	spc = STRPOS(second_pnt,' ')

	IF (spc NE -1) THEN BEGIN
                str1 = STRMID(second_pnt,0,spc)
                str2 = STRTRIM(STRMID(second_pnt,spc,len-spc),2)
                IF(NOT (validnum(str1) AND validnum(str2))) THEN BEGIN
                  MESSAGE, 'Invalid numerical input was supplied.', /CONT
                  GOTO, readsecond
                ENDIF
		second_pnt = [FLOAT(str1), FLOAT(str2)]

		IF (second_pnt(0) EQ 0.0) THEN second_pnt(0) = 1.e-6
		IF (second_pnt(1) EQ 0.0) THEN second_pnt(1) = 1.e-6

	ENDIF ELSE BEGIN	; latitude cut
                IF(NOT validnum(second_pnt)) THEN BEGIN
                  MESSAGE, 'Invalid numerical input was supplied.', /CONT
                  GOTO, readsecond
                ENDIF
		second_pnt = [FLOAT(second_pnt), $
			      first_pnt(1)]
		long_1 = first_pnt(0)
		long_2 = second_pnt(0)
		IF (long_1 EQ 0.0) THEN long_1 = 1.e-6
		IF (long_2 EQ 0.0) THEN long_2 = 1.e-6
		IF (long_2 LT long_1) THEN long_2 = long_2 + 360
		lat_rad = d2r*first_pnt(1)
		IF (lat_rad EQ 0.0) THEN lat_rad = 1.e-6
		lat_len = (long_2-long_1) * cos(lat_rad)
		n_arc = fix(abs(lat_len) / arc_deg) + 1
		del_lat = (long_2-long_1) / n_arc
		vec_arc = dblarr(n_arc+1,3)
		arc_actl = del_lat

		FOR i=0,n_arc DO BEGIN
		  vec_arc(i,0) = cos(d2r*(long_1+(i*del_lat)))*cos(lat_rad)
		  vec_arc(i,1) = sin(d2r*(long_1+(i*del_lat)))*cos(lat_rad)
		  vec_arc(i,2) = sin(lat_rad)
		ENDFOR

		pixel_arc = coorconv(vec_arc,infmt='U',inco=coord, $
				     outfmt='P',outco=out_res)

		lat_flag = -1
		GOTO, lbl

	ENDELSE

	first_pnt_tran = TRANSPOSE(first_pnt)
	second_pnt_tran = TRANSPOSE(second_pnt)

	first_pix = coorconv(first_pnt_tran,infmt='L',outfmt='P', $
			     inco=coord,outco=out_res)
				; convert from lon/lat to pixel #

	second_pix = coorconv(second_pnt_tran,infmt='L',outfmt='P', $
			      inco=coord,outco=out_res)

	end_pixs = [first_pix,second_pix]

ENDELSE


lbl:

IF ((lat_flag EQ 0) AND (lon_flag EQ 0)) THEN BEGIN
				; if not long or lat arc

	pixel_arc = grtcirc(end_pixs,vec_arc,input_type='P', $
		    res=res_input,arc_inc=arc_deg,arc_type=arc_type, $
		    arc_actl)
				; get pixel #'s along great circle

	IF (arc_actl EQ -1) THEN BEGIN
		PRINT,'Great Circle not uniquely defined'
		ret_stat = 1
		GOTO, exit
	ENDIF

	vec_arc = coorconv(vec_arc,infmt='U',inco='E',outco=coord)
				; convert path in unit vectors 
				; from ecliptic (skycube system)
				; to user specified system

ENDIF


n_pix_arc = n_elements(pixel_arc)

IF (re_proj EQ 2) THEN BEGIN	; reprojection cut

	proj_xy = uv2proj(vec_arc,proj,sz_proj)
	i_arc = proj_xy(*,0)
	j_arc = proj_xy(*,1)
	win_title = 'PR  '
			; determine projection coordinates along arc

ENDIF ELSE BEGIN		; cube cut

	fij_input = pix2fij(pixel_arc,res_input)

	i_arc = offx(fij_input(*,0)) * cube_side + fij_input(*,1)
	i_arc = input_l - (i_arc+1) 
	j_arc = offy(fij_input(*,0)) * cube_side + fij_input(*,2)

	win_title = 'SC  '
			; determine skycube coordinates along arc

ENDELSE


; Build window title
; ------------------
IF (lat_flag NE 0) THEN win_title = win_title + ' (lat) '

ll_pos = coorconv([vec_arc(0,*),vec_arc(n_pix_arc-1,*)], $
		  infmt='U',outfmt='L')

IF(lon_flag NE 0) THEN 	ll_pos(*,0) = first_pnt

w_str = string(ll_pos,format='(f6.1)')
win_title = win_title + $
  STRCOMPRESS('(' +w_str(0)+ ',' +w_str(2)+ ')',REMOVE_ALL=1) + '  to  ' + $
  STRCOMPRESS('(' +w_str(1)+ ',' +w_str(3)+ ')',REMOVE_ALL=1) + ' ' + coord



d_mask = intarr(n_pix_arc)
IF (face_num EQ -1) THEN d_mask(*) = 1 ELSE d_mask = (fij_input(*,0) $
					    EQ face_num)
				; mask out points not in active face


x_plot = FINDGEN(n_pix_arc)*arc_actl
IF (re_proj EQ 1) THEN y_plot = nat_cube(i_arc,j_arc)*d_mask $
		  ELSE y_plot = input(i_arc,j_arc)*d_mask
				; extract data along path

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


