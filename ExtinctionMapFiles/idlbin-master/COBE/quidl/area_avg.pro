FUNCTION area_avg, input, nat_cube=nat_cube, weight=weight, $
         poly_pix=poly_pix, proj=proj, coord=coord, win=win, $
         w_pos=w_pos, min=min, max=max, face=face_num, color=color, $
         sent_val=sent_val, desc_reg=desc_reg, zero=zero, uimage=uimage
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    AREA_AVG returns the pixels and area ave of a polygonal region.
;
;DESCRIPTION:
;    This procedure determines the pixels contained within a user-
;    defined polygon and calculates the area average for this region.
;    If the pixel list is supplied on the command line, then an average
;    is computed.  The polygon can be specified using the cursor on a
;    skycube or reprojection or by specifying the longitude and latitude
;    of the vertices.
;
;
;CALLING SEQUENCE:
;    avg = area_avg(input,[nat_cube=nat_cube],[weight=weight], $
;                   [poly_pix=poly_pix],[proj=proj],[coord=coord], $
;                   [win=win],[w_pos=w_pos],[min=min],[max=max], $
;                   [face=face_num],[color=color],[sent_val=sent_val], $
;                   [desc_reg=desc_reg],[/zero])
;
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    input       I     flt arr             skycube/face/projection input
;    [nat_cube]  I     flt arr             Native skycube/face input
;                                          (unfolded or sixpack)
;    [weight]    I     int/flt arr         Weighting skycube
;    [poly_pix]  I/O   int vector          Pixel list    
;    [proj]      I     string              Projection Type:
;                                          ('A' - Aitoff
;                                           'S' - Global Sinusoidal
;                                           'M' - Mollweide)
;    [coord]     I     string              Coordinate System:
;                                          ('E' - Ecliptic
;                                           'G' - Galactic
;                                           'Q' - Equatorial)
;    [win]       I     int                 Window # of image
;    [w_pos]     I     int arr             Image offset within window
;    [min]       I     flt                 Image Scale Minimum
;    [max]       I     flt                 Image Scale Maximum
;    [face]      I     int                 face number (face inpt only)
;    [color]     I     int                 Region boundary color
;    [sent_val]  I     int/flt             Sentinal value
;    [desc_reg]  O     string              Region Description Vector
;    [/zero]                               Suppress zero values in
;                                          average qualifier
;    avg         O     flt                 Area average vector
;
;WARNINGS:
;    If displaying an image, an X-windows terminal must be used.
;    It is recommended when accessing DIRBE data to use either a
;    single face or a small (512 by 256) projection.  Use of an
;    entire DIRBE sky cube or large projection will slow the
;    routine.
;
;EXAMPLES: 
;
;    If the user already has a pixel list, 'pix_list', from which they
;    want to extract an average for a skycube, 'incube', they should
;    enter:
;
;    	avg = area_avg(incube,poly_pix=pix_list).
;
;
;    In the calculation of the average, each pixel is given equal 
;    weight.  Pixels with value 0 are included in the average.  To
;    exclude these pixels the user should enter:
;
;
;    	avg = area_avg(incube,poly_pix=pix_list,/zero).
;
;    In order to include a weighting function, the user should supply
;    a weight skymap of the same dimensions as the input data.
;
;    	avg = area_avg(incube,poly_pix=pix_list,weight=weight).
;
;
;    If the skycube contains a sentinal value this can be supplied
;    using then 'sent_val' keyword.  Any pixel value less than this
;    is ignored when calculating the average.
;
;    	avg = area_avg(incube,poly_pix=pix_list,sent_val=-32000, $
;                      weight=weight).
;
;
;
;    If a pixel list is not supplied on the command line, then the
;    user can define a region.  When using an X-windows terminal, this
;    region can be defined by marking the vertices with the cursor.
;    Alternatively, the user can define the vertices by their longitudes
;    and latitudes.  This is the only option of non-X-windows terminals.
;    If the input array is a reprojection then the original native
;    skycube (unpacked or packed) must be supplied using the 'nat_cube'
;    keyword.  If the user desires the pixel lists for these region
;    then they should supply the 'poly_pix' keyword on the command
;    line which will they contain the pixels numbers of each of the
;    regions defined.  This will be an integer vector with the different
;    pixel lists for each region separated by a '-1' entry.
;
;    For example, if the user has a galactic Aitoff projection, 'repro',
;    of the original 'incube' skycube displayed in X-window 1, and 
;    offset from the window origin by 32,32 then the call will be:
;
;    	avg = area_avg(repro,nat_cube=incube,proj='a',coord='g',win=1, $
;                      w_pos=[32,32],min=0,max=30,poly_pix=pix_list, $
;                      desc_reg=desc_reg,weight=weight).
;
;
;     The 'min' and 'max' keywords are supplied so that the reprojection
;     is properly displayed in the window each time the user is prompted
;     to define a region.  The 'desc_reg' keyword is user to define a 
;     string vector with as many elements as regions defined and 
;     containing a short description of the region supplied by the user.
;     
;     If the skycube itself is displayed, say in window 0, then the call
;     would be:
;
;    	avg = area_avg(incube,win=0,w_pos=[32,32],min=0,max=30, $
;                      poly_pix=pix_list,desc_reg=desc_reg,weight=weight).
;
;
;     If a single cube face is supplied then the call is:
;
;    	avg = area_avg(cube_face,win=0,w_pos=[32,32],min=0,max=30, $
;                      face=1,poly_pix=pix_list,desc_reg=desc_reg, $
;                      weight=face_weight).
;
;#
;COMMON BLOCKS:
;    image_parms
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Get input parameters.
;    Get entry mode (cursor of lon/lat).
;    Specify polygon.
;    Extract pixels contained in polygon.
;    Compute area average.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: getinpar, getemode, grtcirc, drawcirc,
;                        pxinpoly
;
;MODIFICATION HISTORY
;
;   SPR 11174 Written by J.M. Gales, Applied Research Corp.   July 93
;
;
;.TITLE
;Routine AREA_AVG
;-
;

COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
		    cube_side, proj0, coord0, sz_proj


desc = ' '
mask = 0
first = 1B
area_avg = 0


If (N_ELEMENTS(min) EQ 0) THEN min = min(input)
If (N_ELEMENTS(max) EQ 0) THEN max = max(input)
	; get image min and max if not specified

IF (KEYWORD_SET(color) EQ 0) THEN color = 255
	; set boundary color to max if not specified


IF ((!d.name NE 'X') OR (N_ELEMENTS(win) EQ 0)) THEN noshow = 1 $
                                                ELSE noshow = 0
	; If non-X-windows terminal or 'win' keyword not specified
	; then noshow is true



; Make sure skymap and weights map are same size
; ----------------------------------------------
str_wgt = 'Skymap and Weight Arrays not compatible'
IF (N_ELEMENTS(weight) GT 0) THEN BEGIN
	sz_wgt = SIZE(weight)
	IF (N_ELEMENTS(nat_cube) GT 0) THEN BEGIN
		sz0 = SIZE(nat_cube)
		IF (sz0(1) NE sz_wgt(1) OR sz0(2) NE sz_wgt(2)) THEN BEGIN
			MESSAGE,str_wgt,/CONT
			GOTO, exit
		ENDIF
	ENDIF ELSE BEGIN
		sz0 = SIZE(input)
		IF (sz0(1) NE sz_wgt(1) OR sz0(2) NE sz_wgt(2)) THEN BEGIN
			MESSAGE,str_wgt,/CONT
			GOTO, exit
		ENDIF
	ENDELSE
ENDIF
		
	
IF ((noshow EQ 0) AND (N_ELEMENTS(w_pos) EQ 0)) THEN BEGIN
	w_pos = intarr(2)
	READ,'Please enter image offsets within window. (Ex: 32 32): ',w_pos
ENDIF
	; If image offset in window not specified then query



; Set sentinal value to min(skymap)-1 if not specified
; ----------------------------------------------------
IF (N_ELEMENTS(sent_val) EQ 0 AND N_ELEMENTS(nat_cube) EQ 0) $
	THEN sent_val = MIN(input) - 1

IF (N_ELEMENTS(sent_val) EQ 0 AND N_ELEMENTS(nat_cube) GT 0) $
	THEN sent_val = MIN(nat_cube) - 1




; Get projection type and coordinate system if not specified
; ----------------------------------------------------------
sz = SIZE(input)
IF (sz(1) EQ 2*sz(2)) THEN BEGIN

	IF (KEYWORD_SET(proj) EQ 0) THEN BEGIN
		entry = umenu(['Projection Type','Aitoff', $
        	               'Global Sinusoidal', 'Mollweide'], $
                	        title=0,init=1)

		proj = STRMID('ASM',entry-1,1)
	ENDIF

	IF (KEYWORD_SET(coord) EQ 0) THEN BEGIN
		entry = umenu(['Projection coordinate System', $
			       'Ecliptic', 'Galactic', 'Equatorial'], $
                	        title=0,init=1)

		coord = STRMID('EGQ',entry-1,1)
	ENDIF

ENDIF


IF (KEYWORD_SET(coord)) THEN coord0 = STRUPCASE(STRMID(coord,0,1))
IF (KEYWORD_SET(proj)) THEN proj0 = STRUPCASE(STRMID(proj,0,1))



; Get input parameters
; --------------------
IF (getinpar(input,nat_cube,res,re_proj,psize,mask,face_num) NE 0) $
THEN GOTO, exit
out_res = STRCOMPRESS('R' + STRING(res),REMOVE_ALL=1)
sz_mask = SIZE(mask)




; Check that native skympap is present if input is projection
; -----------------------------------------------------------
IF (re_proj EQ 2) THEN BEGIN
	MESSAGE,'Native Skycube must be present',/CONT
	GOTO, exit
ENDIF
	; Skycube must be present when input is reprojection




; If pixel number not specified on command line then get entry mode
; -----------------------------------------------------------------
IF (N_ELEMENTS(poly_pix) EQ 0) THEN BEGIN
	getemode, entry_mode,noshow,coord0,infmt,re_proj
	pp_pres = 0B
ENDIF ELSE BEGIN
	pix = poly_pix
	pp_pres = 1B
	GOTO, get_data
ENDELSE





; Main Loop
; ---------

loop:

IF (noshow EQ 0) THEN BEGIN
	WSET, win
	TVSCL, max<input>min,w_pos(0),w_pos(1) 
ENDIF
	; Redisplay input if X-windows terminal



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

	zoom = 1

PRINT,'Press left button to mark vertices'
PRINT,'To exit press middle or right button'
PRINT, ' '

mark:

	cursor,ic,jc,/down,/device
		; get cursor position

	ic = (ic / zoom) - w_pos(0)
	jc = (jc / zoom) - w_pos(1)
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


		plots,zoom*ic+w_pos(0),zoom*jc+w_pos(1),psym=1, $
		      /DEVICE,COLOR=color
			; plot cross at cursor position

		icur(n_vert) = ic
		jcur(n_vert) = jc
		n_vert = n_vert + 1
		GOTO, mark
			; store cursor positions

		END

	ELSE:	BEGIN		; middle or right button
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

	pnt = ' '

	read,'Enter longitude & latitude of vertex (Ex: 30,25):  ', pnt

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
			       inco=coord0,outco=out_res)
				; convert from lon/lat to pixel #

		vert = [vert,pix]
		n_vert = n_vert + 1

		GOTO, lonlat

	ENDIF ELSE BEGIN

		MESSAGE,'Latitude must be specified',/CONT
		GOTO, lonlat

	ENDELSE

outvert:

ENDELSE


IF (n_vert LE 2) THEN BEGIN
	MESSAGE,'Polygon not completely defined',/CONT
	GOTO, wish
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
			  vec_arc, re_proj, w_pos, face_num, color

	ENDFOR

ENDIF





; Get pixels within polygon
; -------------------------
print,'vert=',vert
IF (N_ELEMENTS(poly_pix) EQ 0 OR first EQ 0B) THEN BEGIN
	PRINT,'Extracting pixels numbers within polygon'
	pix = pxinpoly(vert=vert,res=res)
ENDIF

IF (KEYWORD_SET(uimage)) THEN BEGIN
	poly_pix = pix
	GOTO, exit
ENDIF

IF (pix(0) EQ -1) THEN GOTO, wish


; REMOVE THIS SECTION AFTER DEBUGGING
; -----------------------------------
;goto, get_data
; ----------------------------------------------------------------


get_data:



; Extract data from skymap and calculate average
; ----------------------------------------------
IF (pix(0) NE -1) THEN BEGIN

	IF (N_ELEMENTS(nat_cube) GT 0) THEN $
		data = pix2dat(pixel=pix,raster=nat_cube) ELSE $
		data = pix2dat(pixel=pix,raster=input)

	i = WHERE(data GT sent_val)
	data = data(i)
	pix = pix(i)
		; exclude sentinal values

	IF (KEYWORD_SET(zero) EQ 1) THEN BEGIN
		i = WHERE(data NE 0.0)
		data = data(i)
		pix = pix(i)
	ENDIF
		; exclude zero values if desired

	IF (N_ELEMENTS(weight) NE 0) THEN $
		wgt = pix2dat(pixel=pix,raster=weight) ELSE $
		wgt = 0*pix + 1
		; get weights (default = 1 for all pixels)

	avg = TOTAL(data*wgt) / TOTAL(wgt)
		; calculate weighted average


	PRINT,'Average Value:',avg
	PRINT,'Number of Pixels:', N_ELEMENTS(pix)
	PRINT,' '


	IF (pp_pres EQ 0B) THEN READ,'Description: ',desc


	IF (first EQ 1B) THEN BEGIN
		first = 0B
		poly_pix = [pix]
		area_avg = [avg]
		desc_reg = [desc]
	ENDIF ELSE BEGIN
		poly_pix = [poly_pix,-1,pix]
		area_avg = [area_avg,avg]
		desc_reg = [desc_reg,desc]
	ENDELSE

ENDIF


IF (pp_pres EQ 1B) THEN GOTO, exit





; Get next move
; -------------
wish:	wish = umenu(['Do you wish to:', $
		      'Specify new area?', $
		      'Change Entry Mode?', $
		      'Exit procedure?', $
		      'HELP'],title=0,init=1)

	CASE wish OF

	1 :	GOTO, loop

	2 :	BEGIN

		IF (noshow EQ 1) THEN BEGIN 
			PRINT,'Cursor entry not allowed'
			GOTO, wish
		ENDIF ELSE BEGIN
			getemode, entry_mode,noshow,coord0,infmt,re_proj
			GOTO, loop
		ENDELSE

		END


	3 :	GOTO, exit

	4 :	BEGIN
		PRINT,'No help available at this time'
		GOTO, wish
		END

	ENDCASE


exit:



RETURN,area_avg

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


