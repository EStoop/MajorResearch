FUNCTION pxinpoly,vertices=vertices,coord=coord,res=res,core,boundary
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    PXINPOLY generates pixels within a given spherical polygon
;
;DESCRIPTION:
;    This IDL function generates the pixel numbers of those pixels
;    whose centers fall within a given spherical polygon, that is,
;    points on a sphere that are connected with great circle segments.
;
;CALLING SEQUENCE:
;    px = pxinpoly(vertices=vertices,res=res,[coord=coord])
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    vertices    I     int/flt vector      Vertices of polygon
;                                          (n) for pixel
;                                          (n,2) for lon/lat
;                                          (n,3) for unit vector
;    coord       I     string              Center coordinate system
;                                          ['E' (default), 'G', 'Q']
;                                          (Ignored for pixel input)
;    res         I     int                 Resolution
;    px          O     int/long arr        Pixels within polygon
;
;WARNINGS:
;
; The great circular arcs connecting the vertices are, by default,
; less than 180 degrees.
;
; Polygons must be simply connected, that is, the edges must not
; cross.  PXINPOLY checks for intersecting segments and returns
; [-1] if this is the case.
;
;EXAMPLE: 
;
; To extract the pixels within a polygon defined by the vertices in
; the array 'vert' for a DIRBE (res=9) skymap:
; 
; vert = [45,1010,2032,27]
;
; px = pxinpoly(vertices=vert,res=9)
;
;
; Longitude/Latitude
; ------------------
; vert = [[10,12,45,42],[-5,15,25,0]]  ; "Galactic coordinates"
;
; Note: Longitude are specified in the first row, latitude in the second
;
; px = pxinpoly(vertices=vert,res=9,coord='g')
;
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Find points within circle about polygon
;    Find points within polygon at low res
;    Get points along polygon edges
;    Remove boundary pixels from core
;    Expand core pixels to output resolution
;    Find boundary pixels (at output resolution) within polygon
;    Combine core and boundary pixels together
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: coorconv, get_tri, inpoly, pxincirc, grtcirc
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Apr 93
;    Initial Delivery   SPR 11174
;
; SPR 12146  Fix bug when FLAG variable from inpoly are all zeros, ie
;            there are no res 6 "core" pixels.
;
;
;.TITLE
;Routine PXINPOLY
;-
;

rs = STRCOMPRESS('r' + STRING(res),/REMOVE_ALL)
	; build resolution string


IF (KEYWORD_SET(coord) EQ 0) THEN coord = 'E'
coord = STRUPCASE(STRMID(coord,0,1))
	; set default coordinate system to ecliptic


IF (res LE 6) THEN BEGIN
	res_low = STRCOMPRESS('r' + STRING(res),/REMOVE_ALL)
	num_pix = 6 * 4^(res-1)
ENDIF ELSE BEGIN
	res_low = 'r6'
	num_pix = 6144
ENDELSE

diff_res = res - 6

uv = coorconv(INDGEN(num_pix),infmt='p',inco=res_low,outfmt='u',outco='e')
	; get resolution of initial search
	; if res gt 6 then make initial search res 6
	; get complete skycube unit vectors for initial search



; Get polygon vertices in unit vector format
; ------------------------------------------
sz = SIZE(vertices)

IF (sz(0) EQ 1) THEN BEGIN
	uv_poly = coorconv(vertices,infmt='p',inco=rs,outfmt='u',outco='E')
ENDIF
	; pixels to uv

IF (sz(0) EQ 2 AND sz(2) EQ 2) THEN $
	uv_poly = coorconv(vertices,infmt='l',inco=coord,outfmt='u',outco='E')
	; lon/lat to uv

IF (sz(0) EQ 2 AND sz(2) EQ 3) THEN BEGIN
	uv_poly = coorconv(vertices,infmt='u',inco=coord,outco='E')
	uv_poly = REFORM(uv_poly/SQRT(TOTAL(uv_poly*uv_poly)))
ENDIF
	; uv to uv ecliptic



; Find points within circle about polygon
; ---------------------------------------
avg = TOTAL(uv_poly,1)
avg = REFORM(avg/SQRT(TOTAL(avg*avg)))
	; get average of vertices

cs = MIN(uv_poly # avg)
	; get largest angular separation from average

IF (cs LT 0) THEN BEGIN
	PRINT,STRING([7B]),'Angular Separation of vertices too large'
	px0 = [-1]
	GOTO, exit
ENDIF

ang0 = (ACOS(cs)*180/!dpi) + 10
	; get largest angular separation of vertex from "center"
	; add 10 degrees to be safe

PRINT,'Extracting pixels within polygon'

px0 = pxincirc(cen=avg,ang=ang0,res=MIN([res,6]))
	; get pixels within circle around polygon

vec_in_circ = coorconv(px0,infmt='p',inco=res_low,outfmt='u',outco='e')
	; convert to unit vectors

rot = [[cos(1.e-4),sin(1.e-4),0],[-sin(1.e-4),cos(1.e-4),0],[0,0,1]]
vec_in_circ = vec_in_circ # rot
	; rotate a bit to get vertices off pixel centers


; Break polygon into triangles and find points within at low res
; --------------------------------------------------------------
get_tri,uv_poly,edge_vec,plane,index,type
IF (type(0) EQ 0) THEN BEGIN
	px0 = [-1]
	GOTO, exit
ENDIF
	; Polygon edges cross

n_tri = N_ELEMENTS(index)/3

flag = inpoly(uv_poly,edge_vec,plane,index,vec_in_circ)

k = WHERE(flag EQ 1)

IF (k(0) NE -1) THEN px0 = px0(k)

IF (diff_res LE 0) THEN GOTO, exit
	; if res le 6 then exit




; Get points along polygon edges (at output res)
; ----------------------------------------------
n_vert = N_ELEMENTS(uv_poly)/3
vert_px = coorconv(uv_poly,infmt='u',inco='e',outfmt='p',outco=rs)
		; get vertices in pixel format	

bnd = grtcirc(vert_px(0:1),input_type='P',res=res,arc_inc=0.01,arc_type='S')
FOR i=1,n_vert-1 DO BEGIN
	p = vert_px([i,(i+1) MOD n_vert])
	bnd0 = grtcirc(p,input_type='P',res=res,arc_inc=0.01,arc_type='S')
	bnd = [bnd,bnd0]
ENDFOR


bnd = coorconv(bnd,infmt='p',inco=rs,outfmt='p',outco=res_low)
	; convert to lower resolution




; Remove boundary pixels from existing interior pixels
; ----------------------------------------------------
h = HISTOGRAM(MIN=0,MAX=num_pix-1,bnd)
boundary = WHERE(h GT 0)

core = [-1]
IF (k(0) NE -1) THEN BEGIN
	h = HISTOGRAM(MIN=0,MAX=num_pix-1,px0)
	h(boundary) = 0
	core = WHERE(h GT 0)
ENDIF



; Expand core pixels to output resolution
; --------------------------------------
IF (core(0) NE -1) THEN BEGIN
	no_core = 0
	FOR i=1,2*diff_res DO core = [[core],[core]]
	core = core * (4^diff_res)
	FOR i=1,(4^diff_res)-1 DO core(*,i) = core(*,i) + i
	core = REFORM(core,N_ELEMENTS(core))
		; expand all core pixels to higher resolution
ENDIF ELSE no_core = 1




; Find boundary pixels (at output resolution) within polygon
; ---------------------------------------------------------
FOR i=1,2*diff_res DO boundary = [[boundary],[boundary]]
boundary = boundary * (4^diff_res)
FOR i=1,(4^diff_res)-1 DO boundary(*,i) = boundary(*,i) + i
boundary = REFORM(boundary,N_ELEMENTS(boundary))
	; expand to higher resolution


vec = coorconv(boundary,infmt='p',inco=rs,outfmt='u',outco='e')
		; convert to unit vectors

flag = inpoly(uv_poly,edge_vec,plane,index,vec)
boundary = boundary(WHERE(flag EQ 1))




; Combine core and boundary pixels together
; -----------------------------------------
IF (no_core EQ 0) THEN px0 = [core,boundary] ELSE px0 = [boundary]
		; add these boundary pixels to interior

px0 = px0(SORT(px0))
		; sort by pixel number



exit:

RETURN,px0
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


