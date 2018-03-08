FUNCTION LATCIRC, pole=pole,angle=angle,vec_lc,res=res, $
		  arc_inc=arc_inc,arc_inc_actl
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    latcirc generates a 'latitude' circle
;
;DESCRIPTION:
;    This IDL function generates a circle a about a given 'polar'
;    point.  The circle center is defined by the pole and the boundary
;    by the number of degrees, along a great circular arc, away.
;
;CALLING SEQUENCE:
;    circ = latcirc(pole=pole,angle=angle,vec_lc,[res=res], $
;                   [arc_inc=arc_inc],arc_inc_actl)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    pole        I     int/flt vec         n element vector where n=1 for 
;                                          pixel entry, n=2 for lon/lat
;                                          entry, and n=3 for unit 
;                                          vector entry.
;    angle       I     flt                 
;    [vec_lc]    O     flt arr             circular boundary in
;                                          unit vector format
;    [res]       I     int                 Pixel resolution
;    [arc_inc]   I     flt                 Arc length increment in
;                                          degrees (Default: 1 deg)
;    arc_inc_actl
;                O     flt                 Actual arc length incr
;
;    circ        O     int/flt arr         circular boundary in same
;                                          format as input
;
;WARNINGS:
;
;    None.
;
;EXAMPLES: 
;
;    To draw a circle five degress about the point longitude 30,
;    latitude 25, in two degree increments:
;
;    circ = latcirc(pole=[30,25],angle=5,arc_inc=2.0)
;
;
;    To draw a circle 15 degress about the point centered on DIRBE
;    pixel 1002, in one degree increments, returning the unit vector
;    values in the array uv_circ:
;
;    circ = latcirc(pole=1002,res=9,angle=15,uv_circ)
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Convert polar input to unit vector format (if necessary)
;    Compute rotation matrix to system with z-axis along pole
;    Build boundary vector and rotate about pole
;    Transform boundary vectors back to original coord. system
;    Convert to final output format
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    coorconv
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Dec 92
;
;-
;
IF (KEYWORD_SET(pole) EQ 0) THEN BEGIN
	MESSAGE,'Pole vector not defined',/CONT
	output = 0
	GOTO, exit
ENDIF

IF (KEYWORD_SET(angle) EQ 0) THEN BEGIN
	MESSAGE,'Angle not defined',/CONT
	output = 0
	GOTO, exit
ENDIF

IF ((KEYWORD_SET(res) EQ 0) AND (N_ELEMENTS(pole) EQ 1))THEN BEGIN
	MESSAGE,'Resolution not defined',/CONT
	output = 0
	GOTO, exit
ENDIF


CASE N_ELEMENTS(pole) OF

	3  :	BEGIN
		nrml = pole
		END

	2  :	BEGIN
		nrml = coorconv(pole*1.,infmt='L', outfmt='U')
		END

	1  :	BEGIN

		in_res = STRCOMPRESS('R' + STRING(res),REMOVE_ALL=1)
		nrml = coorconv(pole,infmt='P', inco=in_res, $
			             outfmt='U',outco='E')
		END

	else :	BEGIN
		MESSAGE,'Improper input type'
		output = pole
		GOTO,exit
		END

ENDCASE



nrml = nrml / SQRT(TOTAL(nrml * nrml))
			; normalize normal vector

th_y = ACOS(nrml(2))
cs = SQRT(1 - nrml(2)^2)

IF (cs ne 0.0) THEN BEGIN

	m_z = [[nrml(0)/cs,nrml(1)/cs,0],[-nrml(1)/cs,nrml(0)/cs,0],[0,0,1]]
	m_y = [[cos(th_y),0,-sin(th_y)],[0,1,0],[sin(th_y),0,cos(th_y)]]
	s = (m_z # m_y)

ENDIF ELSE BEGIN

	s = [[1.,0,0],[0,1.,0],[0,0,1.]]

ENDELSE

s_inv = transpose(s)	; orthogonal transformation (inverse = transpose)

IF (N_ELEMENTS(arc_inc) EQ 0) THEN arc_inc = 1.
IF ((360./arc_inc) EQ FIX(360./arc_inc)) THEN n = FIX(360./arc_inc) ELSE $
n = FIX(360./arc_inc) + 1
		; number of positions along circle

arc_inc_actl = 2 * !dpi / n
		; angle increment

vec_lc = fltarr(n+1,3)
		; define output array

cs = cos(arc_inc_actl)
sn = sin(arc_inc_actl)
rot = [[cs,-sn,0],[sn,cs,0],[0,0,1]]
		; build 'equatorial' rotation matrix

arc_inc_actl = arc_inc_actl * 180 / !dpi

vec_lc(0,0) = sin(angle*!dpi/180)
vec_lc(0,1) = 0
vec_lc(0,2) = cos(angle*!dpi/180)
		; build first boundary point

FOR i=1,n DO vec_lc(i,*) = vec_lc(i-1,*) # rot
		; rotate boundary vector about normal vector

vec_lc = vec_lc # s_inv
		; rotate back to original system


CASE N_ELEMENTS(pole) OF

	3  :	BEGIN
		output = vec_lc
		END

	2  :	BEGIN
		output = coorconv(vec_lc,infmt='U', outfmt='L')
		END

	1  :	BEGIN
		output = coorconv(vec_lc,infmt='U', inco='E', $
			          outfmt='P',outco=in_res)
		END

ENDCASE


exit:

RETURN, output
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


