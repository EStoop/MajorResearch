PRO get_tri, vert_vec, edge_vec, plane, index, type
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    GET_TRI divides a polygon into separate triagonal faces.
;
;DESCRIPTION:
;    This IDL divides a polygon into separate triagonal faces.  It is
;    used by the PXINPOLY routine.
;
;CALLING SEQUENCE:
;    get_tri, vert_vec, edge_vec, plane, index, type
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    vert_vec    I     flt array           Vertices of polygon
;                                          (n,3) in unit vector fmt
;    edge_vec    O     flt array           Edge vectors of triangular
;                                          facets (3,3*n_tri) array
;    plane       O     flt array           Coefficients of plane
;                                          defined by facets (4,n_tri)
;                                          ['E' (default), 'G', 'Q']
;    index       O     int array           Polygon vertex indicies of
;                                          each facet  (3,n_tri)
;
;    type        O     int array           Type of vertex
;                                          (1 = convex, -1 = concave)
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Determine whether each vertex is convex or concave
;    Check whether test triangle is outside polygon
;    Check whether any vertices are inside a test facet
;    If enough facets found then quit
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: inpoly,crsprod
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   July 93
;    Initial Delivery   SPR 11174
;
;
;.TITLE
;Routine GET_TRI
;-
;

n_vert = N_ELEMENTS(vert_vec)/3
	; # of vertices
n_tri = n_vert - 2
	; # of triangular facets


; Initialize Arrays
; -----------------
plane = FLTARR(4,n_tri)
index = INTARR(3,n_tri)
type = INTARR(n_vert) + 1
angle = DBLARR(n_vert)

edge_vec = FLTARR(3,3*n_tri)


poly_edge = TRANSPOSE(vert_vec - vert_vec(SHIFT(INDGEN(n_vert),1),*))
	; get vectors between consecutive vertices
	; jth poly_edge() points to jth vertex



; Determine whether vertex is concave (1 = convex, -1 = concave)
; --------------------------------------------------------------
FOR cur_vert=0,n_vert-1 DO BEGIN

	IF (cur_vert EQ 0) THEN prv_vert = n_vert-1 $
			   ELSE prv_vert = cur_vert-1
		; get index of previous vertex

	nxt_vert = (cur_vert+1) MOD n_vert
		; get index of next vertex



	out = vert_vec(nxt_vert,*)-vert_vec(cur_vert,*)
		; get vector from current vertex to next vertex

	in =  vert_vec(cur_vert,*)-vert_vec(prv_vert,*)
		; get vector from previous vertex to current vertex

	crs = crsprod(in,out)
		; get cross product


	IF (TOTAL(crs*vert_vec(cur_vert,*)) LT 0) THEN type(cur_vert) = -1
		; if cross product is in opposite direction 
		; of vector from center of sphere to the vertex 
		; then vertex is concave


	in =  poly_edge(*,cur_vert)
	out = poly_edge(*,nxt_vert)
	dot = TOTAL(in*out)
	norm = SQRT(TOTAL(in*in) * TOTAL(out*out))
	angle(cur_vert) = ACOS(dot/norm) * 180 / !dpi
		; get angle (in degrees) between incoming and 
		; outgoing edge vectors

ENDFOR


tot_ang = TOTAL(type*angle + (1-type) * 180)
IF (2*tot_ang GT 360*n_vert) THEN type = -type
		; if two times sum of angles > proper sum then switch types
		; tot_ang = Sum(ang)
		; tot_ang' = Sum(360-ang) = 360*n_vert-Sum(ang)
		; If Sum(ang) > 360*n_vert-Sum(ang) =>
		;  2*Sum(ang) > 360*n_vert then switch types





; Shift 0th vertex to concave point if it exists
; ----------------------------------------------
i = WHERE(type EQ -1)
IF (i(0) NE -1) THEN BEGIN
	i = min(i)
	type = SHIFT(type,-i)
	vert_vec = SHIFT(vert_vec,-i,0)
	poly_edge = TRANSPOSE(vert_vec - vert_vec(SHIFT(INDGEN(n_vert),1),*))
ENDIF





; Get rotational sense of polygon
; -------------------------------
FOR i=0,n_vert-1 DO BEGIN
	nxt = (i+1) MOD n_vert
	lst = (i+2) MOD n_vert
	IF (type(nxt) EQ 1) THEN BEGIN
		mat = [[REFORM(vert_vec(i,*))], $
		       [REFORM(vert_vec(nxt,*))], $
		       [REFORM(vert_vec(lst,*))]]
		norm = [1,1,1] # INVERT(mat)
		crs = crsprod(poly_edge(*,nxt),poly_edge(*,lst))
		rot = TOTAL(crs) / TOTAL(norm)
		rot = rot/ABS(rot)
		GOTO, find_facets

	ENDIF
ENDFOR
	; Pick 3 consecutive vertices
	; If central vertex is convex then compute normal to
	; triangle
	; Compute cross product of edges
	; Determine whether parallel or anti-parallel
	; rot = +1 for clockwise, -1 for counterclockwise




; Divide polygon into triangular facets
; -------------------------------------
find_facets:

n_facet = 0
edge_used = BYTARR(n_vert)
pnts = INDGEN(n_vert)
next_pnts = INTARR(n_vert)


loop:	i = 0
	l = 0
	n_pnts = N_ELEMENTS(pnts)

WHILE (i LT n_pnts) DO BEGIN

	beg_vert = pnts(i)

	FOR k=1,n_pnts-2 DO BEGIN

		cen_vert = pnts((i+k) MOD n_pnts)
		end_vert = pnts((i+k+1) MOD n_pnts)
			; get central and last vertex

		IF ((cen_vert-beg_vert EQ 1 OR $
		     cen_vert+n_pnts-beg_vert EQ 1) AND $
		     edge_used(beg_vert) EQ 1) THEN BEGIN
			next_pnts(l) = beg_vert
			l = l + 1
			i = i + k
			GOTO, next_vert
		ENDIF
			; If 1st and 2nd vertices are neighbors
			; and this edge belongs to a previously
			; found facet then go to next vertex



		IF ((end_vert-cen_vert EQ 1 OR $
		     end_vert+n_pnts-cen_vert EQ 1) AND $
		     edge_used(cen_vert) EQ 1) THEN BEGIN
			next_pnts(l) = beg_vert
			l = l + 1
			i = i + k
			GOTO, next_vert
		ENDIF
			; If 2nd and 3rd vertices are neighbors
			; and this edge belongs to a previously
			; found facet then go to next vertex



		; Compute triangle info
		; ---------------------
		idx = [beg_vert,cen_vert,end_vert]
		tri_vec = vert_vec(idx,*)
		prev = [2,0,1]
		tri_edge = TRANSPOSE(tri_vec-tri_vec(prev,*))
		nrml = crsprod(tri_edge(*,1),tri_edge(*,0))
		pl = [nrml,TOTAL(nrml*tri_vec(0,*))]



		; Check if exterior triangle 
		; --------------------------
		mat = [[REFORM(vert_vec(idx(0),*))], $
		       [REFORM(vert_vec(idx(1),*))], $
		       [REFORM(vert_vec(idx(2),*))]]
		norm = [1,1,1] # INVERT(mat)
		crs = crsprod(tri_edge(*,0),tri_edge(*,1))
		tri_rot = TOTAL(crs) / TOTAL(norm)
		tri_rot = tri_rot/ABS(tri_rot)

		IF (tri_rot*rot LT 0) THEN BEGIN
			IF (n_pnts EQ 3) THEN BEGIN
				type(0) = 0
				PRINT,STRING([7B])
				PRINT,'Polygon Edges Cross'
				PRINT,' '
				GOTO, exit
			ENDIF
			next_pnts(l) = beg_vert
			l = l + 1
			i = i + k
			GOTO, next_vert
		ENDIF
		; If sense of triangular facet is opposite from
		; overall sense of polygon then this is an
		; exterior triangle




		; Check if any vertices are inside triangle
		; -----------------------------------------
		IF (n_tri GT 1) THEN BEGIN

			check = INTARR(n_pnts) + 1
			check(idx) = 0
			check = WHERE(check EQ 1)

			FOR j=0,N_ELEMENTS(check)-1 DO BEGIN
				v = vert_vec(check(j),*)
				flag = inpoly(tri_vec,tri_edge,pl,[0,1,2],v)
				IF (flag(0) EQ 1) THEN BEGIN
					next_pnts(l) = beg_vert
					l = l + 1
					i = i + k
					GOTO, next_vert
				ENDIF
			ENDFOR
		ENDIF




		; New Triangular Facet Found
		; --------------------------
		IF (cen_vert-beg_vert EQ 1 OR $
		    cen_vert+n_pnts-beg_vert EQ 1) THEN $
			edge_used(beg_vert) = 1
			; If 1st and 2nd vertices are neighbors
			; then mark this edge as used

		IF (end_vert-cen_vert EQ 1 OR $
		    end_vert+n_pnts-cen_vert EQ 1) THEN $
			edge_used(cen_vert) = 1
			; If 2nd and 3rd vertices are neighbors
			; then mark this edge as used

		plane(*,n_facet) = pl
		edge_vec(*,3*n_facet:3*n_facet+2) = tri_edge
		index(*,n_facet) = idx
			; Store facet info for later use in INPOLY

		n_facet = n_facet + 1
		IF (n_facet EQ n_tri) THEN GOTO, exit
			; Exit if enough facets found	

	ENDFOR	; k-loop

i = i + 1

next_vert:

ENDWHILE


pnts = next_pnts(0:l-1)
GOTO, loop

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


