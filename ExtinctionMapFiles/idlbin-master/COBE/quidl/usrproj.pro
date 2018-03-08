pro usrproj,input,proj_img,ilut=ilut,jlut=jlut,mask=mask, $
            proj=proj_type,coord=coord_type,gcoord=gcoord, $
            euler=euler,face=face_num,win=win,merd=merd,para=para, $
            min=min,max=max,noshow=noshow
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    USRPROJ generates projections from user-generated lookup tables.
;
;DESCRIPTION:
;    This IDL procedure creates reprojections of skycubes using user-
;    defined lookup tables created by the PROJBLD facility.  The user
;    supplies the input skycube (in either unfolded or sixpack format,
;    and the the three projection lookup tables built previously by 
;    PROJBLD.  (In contrast, the REPROJ facility restores the needed
;    lookup files, transparently, from the IDL data directory.)  If
;    default coordinate grid overlays are desired, the user must provide
;    the projection type.  The default grid is simply the standard
;    (symmetric) grid drawn by REPROJ when the projection and grid
;    coordinate systems are the same.  If the user wants a non-standard
;    grid, for example, a galactic grid over a projection centered on
;    Orion, they must also provide the reference coordinate system and 
;    Euler angle vector used to define the projection in PROJBLD and the
;    grid coordinate system ('E','G','Q').  If the Euler angle vector
;    is not supplied on the command line, the default grid is used. 
;
;CALLING SEQUENCE:
;    pro usrproj,input,[proj_img],ilut=ilut,jlut=jlut,mask=mask, $
;                [proj=proj_type],[coord=coord_type],[gcoord=gcoord], $
;	         [euler=euler],[face=face_num],[win=win],[merd=merd], $
;                [para=para],[min=min],[max=max],[/noshow]
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    input       I     flt arr             Sky Cube/Face (Res 6-9)
;                                          (Unfolded or Sixpack)
;    [proj_img]  O     flt arr             Rasterized output array
;    ilut        I     int arr             I*2 array containing skycube
;                                          column numbers for reproj.
;                                          Dimensions: (4*n,2*n) 
;    jlut        I     int arr             I*2 array containing skycube
;                                          row numbers for reproj.
;                                          Dimensions: (4*n,2*n) 
;    mask        I     byte arr            Face number / Boundary mask
;    [proj]      I     string              Projection Type:
;                                          ('A','S','M')
;    [coord]     I     string              Reference Coordinate System:
;                                          ('E','G','Q')
;    [gcoord]    I     string              Grid coordinate system
;                                          (
;    [euler]     I     flt arr             Euler angle vector (3 elem)
;    [face]      I     int                 face number (face inpt only)
;    [win]       I     int                 Window # (default:0)
;    [merd]      I     flt arr             Meridian (longitude) array
;                                          Default: [0,60,120,180]
;    [para]      I     flt arr             Parallel (latitude) array
;                                          Default: [0,30,60]
;    [min]       I     flt                 Image Scale Minimum
;    [max]       I     flt                 Image Scale Maximum
;    [noshow]    I     qualifier           no display proj switch
;
;WARNINGS:
;    If displaying an image, an X-windows terminal must be used.
;
;EXAMPLE: 
;
; To generate and display an "anti-ecliptic" Mollweide projection
; of a skycube called 'cube' using the default coordinate overlay and
; default meridian and parallel values and return the image in the
; 'outimg' array:
;
; usrproj,cube,outimg,proj='m',ilut=ae_i,jlut=ae_j,mask=ae_mask
;
; ____
;
; To generate and display in window 2, a "CMB dipole" Global Sinusoidal
; projection with a ecliptic coordinate overlay with meridians at
; 0,45,90,135,180 degrees and no parallels:
;
; usrproj,cube,outimg,proj='s',coord='g',gcoord='e', $
; euler=[267.7,43,0],merd=[0,45,90,135,180],para=-1,ilut=dp_i, $
; jlut=dp_j,mask=dp_mask,win=2
;
; ____
;
; To generate a "Leo-centered" Aitoff projection from a skycube called
; 'incube' and store the rasterized image in 'leo_proj' without 
; displaying it:
;
; usrproj, incube,leo_proj,proj='a',/noshow
;
; ____
;
;#
;COMMON BLOCKS:
;    image_parms, last_usr, proj_array
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Prompt for necessary information.
;    Create and display reprojection.
;    Draw coordinate grids.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: None
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Jun 92
;
; SER 9957   Sixpack input capabilities
; 9-SEP-92   J.M. Gales
;
; SPR 11033  Fix sixpack lookup table adjustments
; 9-JUN-93   J.M. Gales
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;
; SPR 12144  Add last face to LAST_USR common
;            Add polar projection
; 16-MAR-95  J.M. Gales
;
;.TITLE
;Routine USRPROJ
;-
;
COMMON image_parms_usr, input_l,input_h, cur_2_face, offx, offy, $
		        cube_side, proj, coord, sz_proj

COMMON last_usr, last_l, last_h, last_face, last_proj

COMMON proj_array_usr, i_arr, j_arr

IF (N_ELEMENTS(input) EQ 0) THEN BEGIN
	MESSAGE,'Input Array Does Not Exist',/CONT
	GOTO, exit
ENDIF

if (N_ELEMENTS(face_num) eq 0) then face_num = -1
if (keyword_set(win) eq 0) then win=0
if (keyword_set(noshow) eq 0) then noshow=0
IF ((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN noshow = 1
if (keyword_set(merd) eq 0) then merd = [0,60,120,180]
if (keyword_set(para) eq 0) then para = [0,30,60]

min_in = ' '
max_in = ' '

color = 255
re_proj = 1
d2r = !dpi / 180

sz = SIZE(input)
		; get size OF input array

IF (sz(0) NE 2) THEN BEGIN
	MESSAGE,'Input Array Must be Two-Dimensional',/CONT
	GOTO, exit
ENDIF

orient = 'R'
sz_proj = SIZE(ilut)


sz_m = SIZE(merd)
sz_p = SIZE(para)


default = 0
IF (KEYWORD_SET(euler) EQ 0) THEN BEGIN
	euler = [0.,0.,0.]
	default = 1
ENDIF

cos_z1 = cos(d2r*euler(0))
sin_z1 = sin(d2r*euler(0))

cos_y1 = cos(d2r*euler(1))
sin_y1 = sin(d2r*euler(1))

cos_z2 = cos(d2r*euler(2))
sin_z2 = sin(d2r*euler(2))

m_z1 = [[cos_z1,sin_z1,0],[-sin_z1,cos_z1,0],[0,0,1]]
m_y1 = [[cos_y1,0,-sin_y1],[0,1,0],[sin_y1,0,cos_y1]]
m_z2 = [[cos_z2,sin_z2,0],[-sin_z2,cos_z2,0],[0,0,1]]

rot_mat = m_z1 # m_y1 # m_z2


proj_table = ['A','S','M','P']

IF ((keyword_set(proj_type)) EQ 0 AND (noshow EQ 0) AND $
    (sz_m(0)+sz_p(0) GT 0)) THEN BEGIN
	i_proj = umenu(['Enter Projection Type', $
                	'Aitoff', $
	                'Global Sinusoidal', $
        	        'Mollweide','Polar'],title=0,init=1) 

	proj = proj_table(i_proj-1)

ENDIF ELSE IF (KEYWORD_SET(proj_type) NE 0) THEN $
		proj = STRUPCASE(STRMID(proj_type,0,1))


coor_table = ['E','G','Q']

IF ((keyword_set(coord_type) EQ 0)  AND (noshow EQ 0) AND $
    (sz_m(0)+sz_p(0) GT 0) AND (default EQ 0)) THEN BEGIN

	i_coor = umenu(['Enter Reference Coordinate System', $
		        'Ecliptic', $
                	'Galactic', $
	                'Equatorial'], title=0,init=1)

	coord = coor_table(i_coor-1)

ENDIF ELSE IF (KEYWORD_SET(coord_type) NE 0) THEN $
		coord = STRUPCASE(STRMID(coord_type,0,1))


IF ((keyword_set(gcoord) EQ 0) AND (noshow EQ 0) AND $
    (sz_m(0)+sz_p(0) GT 0) AND (default EQ 0)) THEN BEGIN

	i_coor = umenu(['Enter Grid Coordinate System', $
			'Default', $
		        'Ecliptic', $
                	'Galactic', $
	                'Equatorial'], title=0,init=1)

	IF (i_coor GT 1) THEN $
		gcoord = coor_table(i_coor-2) ELSE BEGIN
		gcoord = coord
		rot_mat = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
	ENDELSE

ENDIF


IF (noshow NE 1) THEN BEGIN

	IF (N_ELEMENTS(min) EQ 0) THEN BEGIN

	print, 'Minimum value for Image Scaling <Sky Min>'
	read,min_in
	if (min_in ne '') then min=float(min_in) else min = min(input)

	ENDIF

	IF (N_ELEMENTS(max) EQ 0) THEN BEGIN

	print, 'Maximum value for Image Scaling <Sky Max>'
	read,max_in
	if (max_in ne '') then max=float(max_in) else max = max(input)

	ENDIF

ENDIF	; noshow

IF ((noshow EQ 1) AND (N_ELEMENTS(min) EQ 0)) THEN min = min(input)
IF ((noshow EQ 1) AND (N_ELEMENTS(max) EQ 0)) THEN max = max(input)

input_l = FIX(sz(1))
input_h = FIX(sz(2))


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
		; determine whether unfolded cube, six pack or single face


FOR bit=0,15 DO BEGIN
	IF ((cube_side XOR 2^bit) EQ 0) THEN BEGIN
		res_input = bit+1
		GOTO, lbl_1
	ENDIF
ENDFOR
MESSAGE,'Improper Image Size',/cont
GOTO, exit

lbl_1:

sz = size(ilut)
new_proj = ilut(*,sz(2)/2)
IF (N_ELEMENTS(i_arr) NE 0) THEN BEGIN

	check = ABS(MIN(last_proj-new_proj)) + ABS(MAX(last_proj-new_proj))
	IF ((check NE 0.0) OR (last_face NE face_num) OR $
	    (last_l NE input_l) OR (last_h NE input_h)) THEN BEGIN

		last_flag = -1
		last_l = input_l
		last_h = input_h
		last_face = face_num
		last_proj = new_proj

	ENDIF ELSE last_flag = 0

ENDIF ELSE BEGIN

	last_flag = -1
	last_l = input_l
	last_h = input_h
	last_face = face_num
	last_proj = new_proj

ENDELSE


IF (last_flag EQ -1) THEN BEGIN

	i_arr = ilut/reduce_factor
	j_arr = jlut/reduce_factor

ENDIF
		; scale down lookup table values


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

IF (noshow EQ 1) THEN GOTO, exit


; Draw Parallels
; --------------
vec_arc = fltarr(361,3)

IF (sz_p(0) EQ 0) THEN no_para = 1 ELSE BEGIN 
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

	IF (default EQ 0) THEN $
	vec_arc2 = coorconv(vec_arc,infmt='U',inco=gcoord,outco=coord) $
	ELSE vec_arc2 = vec_arc

	vec_arc2 = vec_arc2 # rot_mat

	draw_grid_line,x,vec_arc2,proj,sz_proj,w_pos,color

ENDFOR



; Draw meridians
; --------------
vec_arc = fltarr(181,3)

IF (sz_m(0) EQ 0) THEN no_merd = 1 ELSE BEGIN 
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
	s = n + (n_pix_arc/2 + 1)

	IF (default EQ 0) THEN $
	vec_arc2 = coorconv(vec_arc,infmt='U',inco=gcoord,outco=coord) $
	ELSE vec_arc2 = vec_arc

	vec_arc2 = vec_arc2 # rot_mat

	draw_grid_line,x,vec_arc2(n,*),proj,sz_proj,w_pos,color
	draw_grid_line,x,vec_arc2(s,*),proj,sz_proj,w_pos,color


ENDFOR


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


