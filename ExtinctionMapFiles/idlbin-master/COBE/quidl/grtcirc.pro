FUNCTION GRTCIRC, input,vec_gc,input_type=input_type,res=res, $
		  arc_inc=arc_inc,arc_type=arc_type,arc_inc_actl
;+                                                                  
;  NAME:
;    grtcirc
;
;  PURPOSE:                                   
;    Generate great circle arc in pixel, radec, lonlat, or unit vector format
;
;  CALLING SEQUENCE:
;    output = grtcirc(input,arc_inc=arc_inc,res=res,arc_type=arc_type)
;
;  INPUT:
;    input - coordinate(s) defining great circle.  if one coordinate is 
;	     given, it is interpreted as normal to the great circle, 
;	     two coordinates are interpreted as the ENDpoints of the 
;	     great circle.
;
;    arc_inc - desired arc length increment in degrees.  Because this 
;	       may not be exactly comiserate with the total arc length,
;	       the actual arc length increment is derived by dividing 
;	       the total length by arc_inc, truncating to the nearest 
;	       integer, and then dividing the total length by this 
;	       number.  This value is then returned in arc_inc.
;
;     res - resolution of pixel number (if pixel input)
;
;     arc type - 'S' short arc (<180 deg), 'L' long arc (>180 deg)
;
;  OUTPUT:
;    Coordinates of great circle
;
;  SUBROUTINES CALLED:
;
;  REVISION HISTORY
;    J.M. Gales
;    Jan 92
;-
;
pi = 4 * ATAN(1.0d0)

input_type = STRUPCASE(STRMID(input_type,0,1))

CASE input_type OF

	'U'  :	BEGIN
		inp = input
		END

	'H'  :	BEGIN
		inp = [[input(*,0)*15],[input(*,1)]]
		inp = coorconv(inp,infmt='L', outfmt='U')
		END

	'L'  :	BEGIN
		inp = coorconv(input,infmt='L', outfmt='U')
		END

	'P'  :	BEGIN
		in_res = STRCOMPRESS('R' + STRING(res),REMOVE_ALL=1)
		inp = coorconv(input,infmt='P', inco=in_res, $
			             outfmt='U',outco='E')
		END

	else :	BEGIN
		MESSAGE,'Improper input type'
		output = input
		GOTO,exit
		END

ENDCASE


sz = SIZE(input)

IF ((sz(1) EQ 2) AND (N_ELEMENTS(input) GE 2)) THEN BEGIN

	desig = 'A'
	arc_type = STRUPCASE(STRMID(arc_type,0,1))
	vec = inp
	;print,desig


ENDIF ELSE IF (sz(0) eq 0 or sz(1) eq 1) THEN BEGIN

	desig = 'N'
	;print,desig
	vec_start = [[1.],[0],[0]]
	theta = 2 * pi
	nrml = inp

ENDIF ELSE BEGIN

	MESSAGE, 'Improper input',/CONT
	output = input
	GOTO, exit

ENDELSE


IF (desig eq 'A') THEN BEGIN

nrml = [vec(0,1)*vec(1,2)-vec(1,1)*vec(0,2), $
	vec(1,0)*vec(0,2)-vec(0,0)*vec(1,2), $
	vec(0,0)*vec(1,1)-vec(1,0)*vec(0,1)]
			; calculate cross product
			; (vector normal to plane containing vectors)

ENDIF

IF (TOTAL(nrml * nrml) EQ 0) THEN BEGIN
	output = -1
	arc_inc_actl = -1
	GOTO, exit
ENDIF

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

IF (desig eq 'A') THEN BEGIN

	vec = vec # s
		; rotate vectors to 'equatorial' plane

	theta = acos(total(vec(0,*)*vec(1,*)))
		; get angle between vectors

	vec_start = vec(0,*)

	IF (arc_type eq 'L') THEN BEGIN

		theta = 2*pi - theta
			; If long arc then get "complement"

		vec_start = vec(1,*)

	ENDIF

ENDIF

n = fix(theta/(arc_inc*pi/180)) + 1
		; number of positions along g.c.

arc_inc_actl = theta/n
		; angle increment

vec_gc = fltarr(n+1,3)
		; define output array

cs = cos(arc_inc_actl)
sn = sin(arc_inc_actl)
rot = [[cs,-sn,0],[sn,cs,0],[0,0,1]]
		; build 'equatorial' rotation matrix

arc_inc_actl = arc_inc_actl * 180 / !dpi

vec_gc(0,*) = vec_start

FOR i=1,n DO vec_gc(i,*) = vec_gc(i-1,*) # rot


vec_gc = vec_gc # s_inv
		; rotate back to original system
;print,vec_gc(0,*)
;print,vec_gc(n,*)


CASE input_type OF

	'U'  :	BEGIN
		output = vec_gc
		END

	'H'  :	BEGIN
		output = coorconv(vec_gc,infmt='U', outfmt='L')
		output(*,0) = output(*,0)/15
		END

	'L'  :	BEGIN
		output = coorconv(vec_gc,infmt='U', outfmt='L')
		END

	'P'  :	BEGIN
		output = coorconv(vec_gc,infmt='U', inco='E', $
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


