FUNCTION define_area, name_sel, poly_pix
;+
;  USKYCUT - a UIMAGE-specific routine.  This is UIMAGE's driver for the
;  SKYCUT procedure.
;#
;  Written by J Ewing & J Gales
;  SPR 10385  Dec 22 92  Set XSTART, XSTOP, etc for graph-object.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;  SPR 11174  Jun 14 94  Upgrade for use of PXINPOLY, J. Newmark
;--------------------------------------------------------------------------
  COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
    cube_side, proj, coord, sz_proj
  COMMON arc, x_plot,y_plot,win_title,fij_input,n_pix_arc,i_arc,j_arc, $
    vec_arc,end_pixs
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON nonxwindow,term_type
  COMMON journal,journal_on,luj
  COMMON color_values,c_badpix,c_draw,c_scalemin
  COMMON history,uimage_version
  COMMON zback,zbgr,zdsrcmap,zbsub
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  desc = ' '
  mask = 0
  area_avg = 0
  first = 1
  noinpt = 0
  n_plot = 0
  nblocks_x = 5
  nblocks_y = 5
  color = 1
;

have_name:
;
;  At this point, NAME_SEL identifies the image on which the user will
;  mark the two arc-endpoints.
;  -------------------------------------------------------------------
  j = EXECUTE('title_sel = ' + name_sel + '.title')
  IF(journal_on) THEN BEGIN
    IF(first EQ 1) THEN PRINTF, luj, 'Define Area'
    PRINTF, luj, '  operand:  ' + title_sel
  ENDIF
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN zoom = 2.^zoom_index $
    ELSE zoom = 1
;
;  Define NAME_BASE to be equal to NAME_SEL if NAME_SEL does not identify
;  a zoomed image, else set NAME_BASE so it identifies the image which
;  had been zoomed in on.
;  ----------------------------------------------------------------------
  type_sel = STRMID(name_sel, 0, STRPOS(name_sel, '('))
  IF(type_sel EQ 'ZOOMED') THEN BEGIN
    j = EXECUTE('win_orig = ' + name_sel + '.win_orig')
    bflag=0
    IF (win_orig LT 0) THEN BEGIN
      win_orig=ABS(win_orig)
      bflag=1
    ENDIF
    j=EXECUTE('zoomflag = ' + name_sel + '.zoomflag')
    IF (zoomflag NE 0) THEN bflag=1
    name_base = get_name(win_orig)
  ENDIF ELSE name_base = name_sel

  type_base = STRMID(name_base, 0, STRPOS(name_base, '('))
  j = EXECUTE('coordinate = ' + name_base + '.coordinate_system')
  coord = STRMID(coordinate, 0, 1)
  IF(coordinate EQ 'EQUATORIAL') THEN coord = 'Q'
;
;  If NAME_BASE identifies a projected image, then set NAME_ORIG so it
;  identifies the image which had been reprojected.
;  -------------------------------------------------------------------
  IF((type_base EQ 'PROJ_MAP') OR (type_base EQ 'PROJ2_MAP'))THEN BEGIN
    j = EXECUTE('projection = ' + name_base + '.projection')
    CASE projection OF
      'AITOFF'            : proj = 'A'
      'GLOBAL SINUSOIDAL' : proj = 'S'
      'MOLLWEIDE'         : proj = 'M'
       ELSE               : BEGIN
                              PRINT, 'Unsupported projection.'
                              RETURN,-1
                            END
    ENDCASE
    j = EXECUTE('win_orig =' + name_base + '.win_orig')
    name_orig = get_name(win_orig)
    type_orig = STRMID(name_orig, 0, STRPOS(name_orig, '('))
    IF(type_orig EQ 'FACE') THEN $
      j = EXECUTE('face_num = ' + name_orig + '.faceno')
    j = EXECUTE('nat_cube = ' + name_orig + '.data')
    re_proj = 1
  ENDIF ELSE re_proj = 0
  j = EXECUTE('img_win = ' + name_sel + '.window')
  j = EXECUTE('orient = ' + name_base + '.orient')
  j = EXECUTE('pos_x = ' + name_base + '.pos_x')
  j = EXECUTE('pos_y = ' + name_base + '.pos_y')
  j = EXECUTE('badpixval = ' + name_base + '.badpixval')
  j = EXECUTE('min = ' + name_base + '.scale_min')
  j = EXECUTE('max = ' + name_base + '.scale_max')
  j = EXECUTE('y_title = ' + name_base + '.units')
  j = EXECUTE('input = ' + name_base + '.data')
  IF(type_sel EQ 'ZOOMED') THEN BEGIN
    j = EXECUTE('start_x = ' + name_sel + '.start_x')
    j = EXECUTE('start_y = ' + name_sel + '.start_y')
    j = EXECUTE('stop_x = ' + name_sel + '.stop_x')
    j = EXECUTE('stop_y = ' + name_sel + '.stop_y')
    j = EXECUTE('specific_zoom = ' + name_sel + '.specific_zoom')
    pos_x = -start_x * zoom
    pos_y = -start_y * zoom
    IF (bflag EQ 1) THEN BEGIN
      IF (zoomflag EQ 1) THEN input(start_x:stop_x,start_y:stop_y)=zbgr
      IF (zoomflag EQ 2) THEN input(start_x:stop_x,start_y:stop_y)=zdsrcmap
      IF (zoomflag EQ 3) THEN input(start_x:stop_x,start_y:stop_y)=zbsub
    ENDIF
  ENDIF ELSE specific_zoom = 1.
  w_pos = [pos_x, pos_y]
  sz = SIZE(input)
  dim1 = sz(1)
  dim2 = sz(2)
  input_l = zoom * dim1
  input_h = zoom * dim2
  IF((type_base EQ 'PROJ_MAP') OR (type_base EQ 'PROJ2_MAP')) THEN BEGIN
    re_proj = 1
    sz_proj = SIZE(input)
    sz_proj(1) = sz_proj(1) * zoom
    sz_proj(2) = sz_proj(2) * zoom
    res_input = 0
    init_res=res_input
  ENDIF ELSE BEGIN
    re_proj = 0
    ic = STRPOS(name_base, '(')
    res_input = FIX(STRMID(name_base, ic-1, 1))
    init_res=res_input
    IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) $
      THEN res_input = res_input + zoom_index
    cube_side = zoom * dim2 / 3.
  ENDELSE
  face_num = -1
  IF(STRMID(type_base, 0, 4) EQ 'FACE') THEN BEGIN
    j = EXECUTE('face_num = ' + name_base + '.faceno')
    offx = [0,0,0,0,0,0]
    offy = [0,0,0,0,0,0]
    cube_side = zoom * dim2
  ENDIF ELSE BEGIN
    offx = [0,0,1,2,3,0]
    offy = [2,1,1,1,1,0]
  ENDELSE
  input = CONGRID(input, sz(1)*zoom, sz(2)*zoom)
  IF(type_base EQ 'FACE') THEN j = EXECUTE('face_num = '+name_base+'.faceno')
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN noshow = 1 ELSE noshow = 0
  nostore = 1
;
;  Preparatory set-up is completed, now call Joel's primitives.
;  ------------------------------------------------------------
  coord0 = coord 
  ret_stat = getinpar(input, nat_cube, res, re_proj, psize, mask, face_num)
e_mode:
  arc_type = 'S'
  getemode, entry_mode, noshow, coord, infmt, re_proj
  out_res = STRCOMPRESS('R' + STRING(res),REMOVE_ALL=1)
  sz_mask = SIZE(mask)
loop:
  IF(noshow EQ 0) THEN BEGIN
    PRINT, 'Please mark points on the window titled:  ' + bold(title_sel)
    WSET, img_win
    WSHOW, img_win
;    TVSCL, max<input>min,w_pos(0),w_pos(1) 
  ENDIF

n_vert = 0
vert = [-1]

; Set cursor limits
; -----------------
IF (sz_mask(0) EQ 0) THEN BEGIN	; skycube

	i_max = input_l - 1
	j_max = input_h - 1

ENDIF ELSE BEGIN		; reprojection

	i_max = sz_mask(1) - 1
	j_max = sz_mask(2) - 1

ENDELSE


; Cursor Entry Section
; --------------------
IF (entry_mode EQ 1) THEN BEGIN

	icur = intarr(100)
	jcur = intarr(100)
		; allocate cursor position arrays

PRINT,'Press left button to mark vertices'
PRINT,'To exit press middle or right button'
PRINT, ' '

mark:

	cursor,ic,jc,/down,/device
		; get cursor position
        ic=FLOAT(ic)
        jc=FLOAT(jc)
	ic = (ic / specific_zoom) - w_pos(0)
	jc = (jc / specific_zoom) - w_pos(1)
		; subtract window offsets

	CASE !err OF

	1:	BEGIN		; left button

		IF (ic LT 0 OR jc LT 0 OR ic GE i_max OR jc GE j_max) $
		THEN GOTO, mark
			; if outside image limits then reject, get new pos

		IF (sz_mask(0) EQ 0) THEN BEGIN		; skycube

			IF (cur_2_face(ic/cube_side,jc/cube_side) EQ -1) $
			THEN GOTO, mark
			; if not on image then reject, get new pos

		ENDIF ELSE IF (mask(ic,jc) EQ 255) THEN GOTO, mark
			; if reprojection and outside image then get new pos

		plots,specific_zoom*(ic+w_pos(0)),specific_zoom*(jc+w_pos(1)),$
                   psym=1, /DEVICE,COLOR=color
			; plot cross at cursor position

		icur(n_vert) = FIX(ic)
		jcur(n_vert) = FIX(jc)
		n_vert = n_vert + 1
		GOTO, mark
			; store cursor positions

		END

	ELSE:	BEGIN		; middle or right button
                     IF (n_vert LE 2) THEN BEGIN
                 	MESSAGE,'Polygon not completely defined',/CONT
           	        GOTO, mark
                     ENDIF
                 	; Polygon must be at least triangle
		     icur = icur(0:n_vert-1)
		     jcur = jcur(0:n_vert-1)
		     GOTO, extcur
		END

	ENDCASE

extcur:

	IF (re_proj EQ 0) THEN BEGIN	; skycube/face

		cur_face = cur_2_face(icur/cube_side,jcur/cube_side)
		i_input = icur - offx(cur_face) * cube_side
		i_input = input_l - (i_input+1) 	; flip to RT
		j_input = jcur - offy(cur_face) * cube_side
		vert = fij2pix([[cur_face],[i_input],[j_input]],res)
				; determine face number
				; determine x&y raster coordinate
				; determine pixel # of endpoints

	ENDIF ELSE BEGIN		; reprojection

		FOR i=0,n_vert-1 DO BEGIN
			pnt  = proj2uv([[icur(i)],[jcur(i)]],proj,sz_proj)
				; determine unit vectors
			pix = coorconv(pnt,infmt='U',outfmt='P', $
				       inco=coord0,outco=out_res)
			vert = [vert,pix]
		ENDFOR
		vert = vert(1:n_vert)

	ENDELSE

ENDIF ELSE BEGIN


; Lon/Lat Entry Section
; ---------------------

lonlat:
    PRINT,' '
   PRINT,'You will be prompted for points that define the vertices of an area.'
   PRINT,'You need to define at least 3 points. Enter the Longitude & Latitude'
   PRINT,'in degrees followed by a <CR>. After entering all the points, hit a <CR>.'  
   PRINT,' '
;   PRINT,'Please identify the coordinate system you will use.'
;   menu = ['Coordinate system:','Ecliptic','Equatorial','Galactic']
;   sel = one_line_menu(menu)
;   icoord = STRMID('EQG',sel-1,1)
   PRINT,'Next enter the lon/lat pair in degrees (example: '+$
        bold('20. -30.')+').  [<CR> to exit]'
keyin:
	pnt = ' '
	read,'Enter longitude & latitude of vertex :  ', pnt

	IF (STRLEN(pnt)EQ 0) THEN BEGIN
		vert = vert(1:n_vert)
		GOTO, outvert
	ENDIF
		; Exit section if <CR>


	pnt = STRTRIM(pnt,2)
	len = STRLEN(pnt)
			; trim leading and trailing blanks
			; get string length

	FOR i=0,len-1 DO IF (STRMID(pnt,i,1) EQ ',') THEN $
	pnt = STRMID(pnt,0,i) + ' ' + STRMID(pnt,i+1,len-i-1)
			; replace ',' with ' '

	spc = STRPOS(pnt,' ')
			; get space position

	IF (spc NE -1) THEN BEGIN
		pnt = [FLOAT(STRMID(pnt,0,spc)),FLOAT(STRMID(pnt,spc,len-spc))]

		IF (pnt(0) EQ 0.0) THEN pnt(0) = 1.e-6
		IF (pnt(1) EQ 0.0) THEN pnt(1) = 1.e-6
		pnt_tran = TRANSPOSE(pnt)
		pix = coorconv(pnt_tran,infmt='L',outfmt='P', $
			       inco=coord,outco=out_res)
				; convert from lon/lat to pixel #

		vert = [vert,pix]
		n_vert = n_vert + 1

		GOTO, keyin

	ENDIF ELSE BEGIN

		MESSAGE,'Latitude must be specified',/CONT
		GOTO, keyin

	ENDELSE

outvert:

ENDELSE


IF (n_vert LE 2) THEN BEGIN
	MESSAGE,'Polygon not completely defined',/CONT
	GOTO, mark
ENDIF
	; Polygon must be at least triangle


; Draw polygon if map displayed
; -----------------------------
IF (noshow EQ 0) THEN BEGIN

	FOR i=0,n_vert-1 DO BEGIN

		IF (i NE n_vert-1) THEN	end_pixs = vert(i:i+1) $
				   ELSE end_pixs = [vert(n_vert-1),vert(0)]


		pixel_arc = grtcirc(end_pixs,vec_arc,input_type='P', $
			    res=res,arc_inc=1.0,arc_type='S',arc_actl)
				; get pixel #'s along great circle


		IF (re_proj EQ 1) THEN $
		vec_arc = coorconv(vec_arc,infmt='U',inco='E',outco=coord0)
				; convert path in unit vectors 
				; from ecliptic (skycube system)
				; to user specified system


		fij_input = pix2fij(pixel_arc,res)

		i_arc = offx(fij_input(*,0)) * cube_side + fij_input(*,1)
		i_arc = input_l - (i_arc+1) 
		j_arc = offy(fij_input(*,0)) * cube_side + fij_input(*,2)

		n_pix_arc = n_elements(pixel_arc)

		drawcirc, pix2fij(pixel_arc,res), n_pix_arc, i_arc, j_arc, $
		    vec_arc, re_proj, w_pos, face_num, color,$
 		    zoom=specific_zoom

	ENDFOR

ENDIF
                                                                  
; Get pixels within polygon
; -------------------------
   diffres=init_res-res
   IF((type_base NE 'PROJ_MAP') AND (type_base NE 'PROJ2_MAP') AND $
      (diffres NE 0)) THEN BEGIN
      IF (NOT DEFINED(zoom_index)) THEN zoom_index=0
      input = DOWNRES(input, zoom_index)
      inco='R'+STRTRIM(STRING(res),1)
      outco='R'+STRTRIM(STRING(init_res),1)
      vert = coorconv(vert,infmt='P', outfmt='P',inco=inco,outco=outco)
   ENDIF ELSE init_res=res
   poly_pix = pxinpoly(vert=vert,res=init_res)
   IF (poly_pix(0) EQ -1) THEN RETURN,-1

   IF (N_ELEMENTS(nat_cube) GT 0) THEN $
     area=pix2dat(pixel=poly_pix,raster=nat_cube) ELSE $
     area=pix2dat(pixel=poly_pix,raster=input)

   RETURN,area
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


