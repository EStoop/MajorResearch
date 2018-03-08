FUNCTION inpoly, vert_vec,edge_vec,plane,index,uv
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    INPOLY determines whether a point is inside a spherical polygon.
;
;DESCRIPTION:
;    This IDL divides a polygon into separate triagonal faces.  It is
;    used by the PXINPOLY routine.
;
;CALLING SEQUENCE:
;    flag = inpoly(vert_vec,edge_vec,plane,index,uv)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    vert_vec    I     flt array           Vertices of polygon
;                                          (n_vert,3) in unit vector fmt
;    edge_vec    I     flt array           Edge vectors of triangular
;                                          facets (3,3*n_tri) array
;    plane       I     flt array           Coefficients of plane
;                                          defined by facets (4,n_tri)
;                                          ['E' (default), 'G', 'Q']
;    index       I     int array           Polygon vertex indicies of
;                                          each facet  (3,n_tri)
;    uv          I     flt array           Vertices of polygon
;                                          (n_pnts,3) in unit vector fmt
;    flag        O     byte vec            Flag value (1 if inside,
;                                          0 if outside)
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    See code documentation
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: crsprod
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   May 93
;    Initial Delivery   SPR 11174
;
;
;.TITLE
;Routine INPOLY
;-
;
n_pnts = N_ELEMENTS(uv)/3
	; # of points to test
n_tri = (N_ELEMENTS(vert_vec)/3) - 2
	; # of triangular facets in polygon
flag = BYTARR(n_pnts)
	; allocate flag array


FOR tst=0L,n_pnts-1 DO BEGIN

	x = uv(tst,0)
	y = uv(tst,1)
	z = uv(tst,2)
	IF (z NE 0) THEN BEGIN
		a0 = x/z
		a1 = y/z
	ENDIF ELSE IF (y NE 0) THEN a0 = x/y


	FOR tri=0,n_tri-1 DO BEGIN

		tri_vec = vert_vec(index(*,tri),*)
		pl0 = plane(0,tri)
		pl1 = plane(1,tri)
		pl2 = plane(2,tri)
		pl3 = plane(3,tri)

		IF (z NE 0) THEN BEGIN
			v3 = pl3 / (pl0 * a0 + pl1 * a1 + pl2)
			uv_pl = v3 * [a0,a1,1]
		ENDIF ELSE IF (y NE 0) THEN BEGIN
			v2 = pl3 / (pl0 * a0 + pl1)
			uv_pl = v2 * [a0,1,0]
		ENDIF ELSE uv_pl = [pl3/pl0,0,0]
		; get intersection of plane defined by triangular
		; facet and line defined by the sphere center and
		; the test point.

		diff = tri_vec - TRANSPOSE([[uv_pl],[uv_pl],[uv_pl]])
		; get difference between this point and vertices
		; of triangular facet

		FOR j=0,2 DO BEGIN
			crs = crsprod(edge_vec(*,3*tri+j),diff(j,*))
			test = TOTAL(crs*plane(0:2,tri))
			IF (test LT 0) THEN GOTO, next_tri
		ENDFOR
		; for each of these difference vectors
		; take cross product of edge vector to next vertex
		; and this difference
		; if in opposite dirction of normal outward from
		; triangular face then outside
		; check next facet

		flag(tst) = 1B
		GOTO, next_pix
		; Inside facet therefore set flag and go on to next point

		next_tri:

	ENDFOR

next_pix:

ENDFOR

RETURN, flag
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


