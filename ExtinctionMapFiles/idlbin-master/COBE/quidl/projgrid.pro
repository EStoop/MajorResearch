pro projgrid,vimage,proj=proj_type,coord=coord_type,gcoord=gcoord, $
             face=face_num,w_pos=w_pos,merd=merd,para=para,color=color
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    PROJGRID overlays coordinate grids on projections.
;
;DESCRIPTION:
;    This IDL procedure generates a coordinate grid that is overlaid
;    on a currently displayed image.  The user must provide the
;    projection, which may be an unfolded skycube or an Aitoff, Global
;    Sinusoidal, or Mollweide projection, in an existing window.  In
;    the latter case, the image array can be of any size but it must
;    have a 2:1 aspect ratio.  In addition, the grid coordinate system
;    does not have to be the same as the projection coordinate system.
;
;CALLING SEQUENCE:
;    pro projgrid,vimage,[proj=proj_type],[coord=coord_type], $
;                 [gcoord=gcoord],[w_pos=w_pos],$
;                 [merd=merd],[para=para],[color=color]
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    vimage      I     flt arr             Image
;    [proj]      I     string              Projection Type:
;                                          ('A','S','M')
;    [coord]     I     string              Projection Coordinate System:
;                                          ('E','G','Q')
;    [gcoord]    I     string              Grid Coordinate System:
;                                          ('E','G','Q')
;    [w_pos]     I     int arr             Offset of image in window
;                                          (Default: [0,0])
;    [merd]      I     flt arr             Meridian (longitude) array
;                                          (Default: [0,60,120,180])
;                                          (merd = -1 if none)
;    [para]      I     flt arr             Parallel (latitude) array
;                                          (Default: [0,30,60])
;                                          (para = -1 if none)
;    [color]     I     int                 Grid color (0-255)
;
; Note: The meridian and parallel entries are "symmetrized", ie,
;       the default parallels are drawn at -60, -30, 0, 30, 60.
;
;WARNINGS:
;    The user must set the window (with the WSET command) before
;    running PROJGRID.
;
;EXAMPLE: 
;
; To display an equatorial coordinate overlay on a skycube called, 
; 'cube', with no offset in the display window:
;
; projgrid,cube,gcoord='q'
;
; In this case it is not necessary to give the projection type as
; the program can recognize a skycube by its 4:3 aspect ratio.
;
; ----
;
; To display an ecliptic coordinate overlay with the default meridians
; (0,60,120,180) and parallels (0,30,60) on an Galactic Aitoff 
; reprojection, 'ait_img', offset within the window by [32,32].
;
; projgrid, ait_img,proj='a',coord='g',gcoord='e',w_pos=[32,32]
; ____
;
; To display a galactic coordinate overlay with the meridians at
; (0,45,90,135,180) and no parallels on a Equatorial Mollweide 
; reprojection, 'mol_img', with zero window offset and grid color 128.
;
; projgrid, mol_img,proj='m',coord='q',gcoord='g', $
;           merd=[0,45,90,135,180),para=-1,color=128
;
;
; If the 'proj', 'coord', or 'gcoord' keywords are not provided, the user
; is prompted for them.
; ____
;
;#
;COMMON BLOCKS:
;    image_parms
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Essentially a routing routine.  Querys for info and uses 
;    if-then-else and case statements to make the proper call
;    to DRAW_PM which actually generates and displays the grid
;    overlay.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: draw_pm
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   May 92
;
;    SER 9756
;
;    SPR 9834	When the 'proj' and 'coord' keywords are used, the
;		variables, 'proj_type' and 'coord_type' are defined.
;		These values must then be assigned to 'proj' and
;		'coord'.
; 24JUL 1992
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;
;.TITLE
;Routine PROJGRID
;-
;

COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
		    cube_side, proj, coord, sz_proj

sz_proj = SIZE(VIMAGE)

IF (sz_proj(0) NE 2) THEN BEGIN
	MESSAGE,'Image Array Must be Two-Dimensional',/CONT
	GOTO, exit
ENDIF

IF (n_elements(face_num) EQ 0) THEN face_num = -1
if (keyword_set(merd) eq 0) then merd = [0,60,120,180]
if (keyword_set(para) eq 0) then para = [0,30,60]
IF (N_ELEMENTS(w_pos) EQ 0) THEN w_pos = [0,0]
IF (keyword_set(color) EQ 0) THEN color = 255

proj_table = ['a','s','m']
coor_table = ['e','g','q']
res_input = 0

IF (3*sz_proj(1) EQ 4*sz_proj(2)) THEN BEGIN

	cube_side = sz_proj(1) / 4
	offx = [0,0,1,2,3,0]
	offy = [2,1,1,1,1,0]
	input_l = sz_proj(1)
	input_h = sz_proj(2)

	FOR bit=0,15 DO BEGIN
		IF ((cube_side XOR 2^bit) EQ 0) THEN BEGIN
		res_input = bit+1
		GOTO, lbl_1
		ENDIF
	ENDFOR
	MESSAGE,'Improper Image Size',/cont
	GOTO, exit

lbl_1:

	proj = 'Q'
	coord = ' '

	GOTO, gcd

ENDIF

IF (sz_proj(1) EQ sz_proj(2)) THEN BEGIN

	cube_side = sz_proj(1)
	offx = [0,0,0,0,0,0]
	offy = [0,0,0,0,0,0]
	input_l = sz_proj(1)
	input_h = sz_proj(2)

	FOR bit=0,15 DO BEGIN
		IF ((cube_side XOR 2^bit) EQ 0) THEN BEGIN
		res_input = bit+1
		GOTO, lbl_2
		ENDIF
	ENDFOR
	MESSAGE,'Improper Image Size',/cont
	GOTO, exit

lbl_2:

	proj = 'Q'
	coord = ' '

	IF (face_num EQ -1) THEN READ, 'Enter Face Number: ',face_num

	GOTO, gcd

ENDIF


IF (sz_proj(1) NE sz_proj(2)*2) THEN BEGIN
	MESSAGE,'Reprojected Image Array Must Have 2:1 Aspect Ratio',/CONT
	GOTO, exit
ENDIF


IF (keyword_set(proj_type) EQ 0) THEN BEGIN
	i_proj = umenu(['Enter Projection Type', $
                	'Aitoff', $
	                'Global Sinusoidal', $
        	        'Mollweide'],title=0,init=1) 

	proj = proj_table(i_proj-1)

ENDIF ELSE proj = proj_type

proj = STRUPCASE(STRMID(proj,0,1))

IF (keyword_set(coord_type) EQ 0) THEN BEGIN

	i_coor = umenu(['Enter Projection Coordinate System', $
		        'Ecliptic', $
                	'Galactic', $
	                'Equatorial'], title=0,init=1)

	coord = coor_table(i_coor-1)

ENDIF ELSE coord = coord_type

coord = STRUPCASE(STRMID(coord,0,1))


gcd:

IF (keyword_set(gcoord) EQ 0) THEN BEGIN

	i_coor = umenu(['Enter Grid Coordinate System', $
		        'Ecliptic', $
                	'Galactic', $
	                'Equatorial'], title=0,init=1)

	gcoord = coor_table(i_coor-1)

ENDIF

gcoord = STRUPCASE(STRMID(gcoord,0,1))

IF (proj EQ 'Q') THEN re_proj = 0 ELSE re_proj = 1

draw_pm,para,merd,re_proj,face_num,w_pos,res_input,gcoord,color=color

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


