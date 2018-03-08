FUNCTION proj2uv, ij, proj, sz_proj
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    proj2uv - Converts projection coodinates to unit vectors.
;
;DESCRIPTION:
;    This IDL function converts an array of projection screen coordinates 
;    to an array of unit vectors.  The input array has dimensions (n,2) 
;    where n is the number of points to be converted.  For a 512 by 256 
;    the numbers in this array range from 0 to 511 for the elements in 
;    the first row and 0 to 255 for elements in the second.  The output
;    unit vector array has dimensions (n,3).  The coordinate system for 
;    the output is the same as that of the projection.  This routine is 
;    used by SKY_CUT/CROSS_SECTION and the pixel identification routines
;    to determine the sky coordinates of a projection pixel.
;
;CALLING SEQUENCE:
;    uvec = proj2uv, ij, proj, sz_proj
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    ij          I     int/long arr        Projection screen coord. array
;    proj        I     string              projection string ('A','S','M')
;    sz_proj     I     long arr            Array = SIZE(projection array)
;    uvec        O     flt arr             Output unit vector array
;
;WARNINGS:
;    Only works for Aitoff, Global Sinusoidal, and Molweide projections.
;    Also input coordinates outside the projection image will bomb.  
;    This is NOT a fully supported user program.
;
;EXAMPLE: 
;
; To convert from the screen coordinates (100,50) for an Aitoff 
; projection (stored in the IDL array, 'proj_image') to the unit
; vector 'uv':
;
; uv = proj2uv([[100],[50]],'a',SIZE(proj_image))
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Uses CASE statement to direct to proper code for each projection.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Uses coordinate transformation formuli that can be found in any 
;    complete book on mapping/projections.
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Apr 92
;
; SPR 10000   Error in y-axis input for Mollweide projection
; 17-SEP-92   J.M. Gales
;
;.TITLE
;Routine PROJ2UV
;-
;
n = sz_proj(2)/2
invc = 1 / (4.0 * n*n)
pi = 4 * atan(1.d0)
half_pi = pi / 2

half_l = sz_proj(1)/2
half_h = sz_proj(2)/2

CASE STRUPCASE(STRMID(proj,0,1)) OF

'A' :	BEGIN


	xa = half_l - ij(*,0) - 0.5
	ya = ij(*,1) - half_h + 0.5


	r = ya / xa
	s_prime = 1 - invc * (xa*xa + 4*ya*ya)

	c1_sqr = (1 + (2*r*s_prime)^2) / (1 + (2*r)^2)
	c1 = SQRT(c1_sqr)
	c2 = s_prime / c1
	c3 = 2 * c2*c2 - 1

	x = c1 * c3
	y = 2 * s_prime * SQRT(ABS(1-c2*c2)) * (xa/ABS(xa))
	z = SQRT(ABS(1-c1_sqr)) * (ya/ABS(ya))

	END


'S' :	BEGIN


	xs = pi * (half_l - ij(*,0) - 0.5) / (2*n)
	ys = half_pi * (ij(*,1) - half_h + 0.5) / n

	cs = COS(ys)
	sn = SQRT(1 - cs*cs)

	x = cs * COS(xs/cs)
	y = cs * SIN(xs/cs)
	z = sn * (ys/abs(ys))


	END


'M' :	BEGIN


	xm = (ij(*,0) - half_l + 0.5) / n
	ym = (ij(*,1) - half_h + 0.5) / n

	sn_beta = ym
	cs_beta = SQRT(1 - sn_beta^2)
	beta = ASIN(sn_beta)
	sn_2beta = 2 * sn_beta * cs_beta
	sn_del = (2*beta + sn_2beta) / pi
	cs_del = SQRT(1 - sn_del^2)

	s_prime = 1 - 0.25 * (xm*xm + 4*ym*ym)


	alph = 0.5 * xm * pi / cs_beta
	cs_alph = COS(alph)
	sn_alph = SQRT(1 - cs_alph^2)

	x = cs_del * cs_alph
	y = cs_del * sn_alph * (-xm/abs(xm))
	z = sn_del


	END


ENDCASE

RETURN,[[x],[y],[z]]
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


