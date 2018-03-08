FUNCTION pxincirc,cen=cen,coord=coord,ang=ang,res=res,ii=ii
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    PXINCIRC generates pixels within a given circle
;
;DESCRIPTION:
;    This IDL function generates the pixel numbers of those pixels
;    whose centers fall within a given angular radius of a given
;    point.  The angle is expressed in degrees and the center can
;    be given in pixel number, unit vector, or long/lat.
;
;CALLING SEQUENCE:
;    px = pxincirc(cen=cen,ang=ang,res=res,[coord=coord])
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    cen         I     int/flt vector      Circle center
;    coord       I     string              Center coordinate system
;                                          ['E' (default), 'G', 'Q']
;    ang         I     flt                 Angular radius (in degrees)
;    res         I     int                 Resolution
;    px          O     int/long arr        Pixels within circle
;
;WARNINGS:
;
;EXAMPLE: 
;
; To extract the pixels within 10 degrees of a circle centered on 
; pixel 3457 for a DIRBE (res=9) skymap:
; 
; px = pxincirc(cen=3457,ang=10,res=9)
;
; To extract the pixels within 30 degrees of a circle centered on 
; [30,-23] galactic coordinates for a FIRAS (res=6) skymap:
; 
; px = pxincirc(cen=[30,-23],coord='g',ang=30,res=9)
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Extracts pixels for res 6 or smaller skymap directly
;    If larger res then check corners of pixels
;    If all corners within circle than expand pixels to higher res
;    If any corner outside circle than expand pixels to higher res
;    and check within center directly.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: coorconv
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Apr 93
;    Initial Delivery   SPR 10776
;
;
;.TITLE
;Routine PXINCIRC
;-
;

sq2 = 1/SQRT(2.d0)
n_vec_index = [2,0,1,0,1,2]
	; define index of face normals

corn_a = [[sq2,sq2,0], $
          [0,sq2,sq2], $
          [sq2,0,sq2], $
          [0,sq2,sq2], $
          [sq2,0,sq2], $
          [sq2,sq2,0]]

corn_b = [[sq2,-sq2,0], $
          [0,sq2,-sq2], $
          [sq2,0,-sq2], $
          [0,sq2,-sq2], $
          [sq2,0,-sq2], $
          [sq2,-sq2,0]]
	; define vectors from pixel centers to corners


rs = STRCOMPRESS('r' + STRING(res),/REMOVE_ALL)
	; build resolution string

IF (res LE 6) THEN BEGIN
	rs0 = STRCOMPRESS('r' + STRING(res),/REMOVE_ALL)
	num_pix = 6 * 4^(res-1)
	ang0 = ang
ENDIF ELSE BEGIN
	rs0 = 'r6'
	num_pix = 6144
	ang0 = ang + 10
ENDELSE
	; get resolution of initial search
	; if res gt 6 then make initial res 6 search
	; increase angular radius of search circle just to be safe

cs0 = COS(ang0*!dpi/180)
	; compute minimun cosine for inclusion for initial search

IF (KEYWORD_SET(coord) EQ 0) THEN coord = 'E'
coord = STRUPCASE(STRMID(coord,0,1))

n = N_ELEMENTS(cen)
CASE n OF
	1: BEGIN
		uv0 = coorconv([cen],infmt='p',inco=rs,outfmt='u',outco='e')
		uv0 = REFORM(uv0/SQRT(TOTAL(uv0*uv0)))
	   END

	2: BEGIN
		uv0 = coorconv(cen,infmt='l',inco=coord,outfmt='u',outco='E')
		uv0 = REFORM(uv0/SQRT(TOTAL(uv0*uv0)))
	   END

	3 : BEGIN
		uv0 = coorconv(cen,infmt='u',inco=coord,outco='E')
		uv0 = REFORM(uv0/SQRT(TOTAL(uv0*uv0)))
	    END
ENDCASE
	; get circle center in unit vector format


uv = coorconv(INDGEN(num_pix),infmt='p',inco=rs0,outfmt='u',outco='e')
	; get complete skycube unit vectors for initial search


cs = uv # uv0
	; calculate dot product
interior = WHERE(cs GE cs0)
	; pixel in interior if dot product >= limit



IF (res GT 6) THEN BEGIN

	diff_res = res - 6
	cs1 = COS(ang*!dpi/180)
		; compute minimun cosine for inclusion for final search
	face = interior / 1024
		; get pixel face
	n_vec = n_vec_index(face)
		; get index of normal vector to face

	vec0 = coorconv(interior,infmt='p',inco='r6',outfmt='u',outco='e')
		; get complete skycube unit vectors of interior

	FOR i=0,N_ELEMENTS(interior)-1 DO BEGIN
		v_to_face = vec0(i,n_vec(i))
		vec0(i,*) = vec0(i,*) * 16 / v_to_face
		vec0(i,*) = vec0(i,*) * (v_to_face / abs(v_to_face))
	ENDFOR
		; scale so that distance along face normal is +/-16

	vec = FLTARR(N_ELEMENTS(interior),3)
	flag = BYTARR(N_ELEMENTS(interior),4)
		; allocate corner and flag arrays


	; Corner 1
	; --------
	vec(*,0) = vec0(*,0) + corn_a(0,face)
	vec(*,1) = vec0(*,1) + corn_a(1,face)
	vec(*,2) = vec0(*,2) + corn_a(2,face)
		; get vector to first corner

	cs = vec # uv0
	norm = SQRT(TOTAL(vec*vec,2))
	flag(*,0) = (cs gt cs1*norm)
		; calculate dot product
		; calculate norm of vector
		; check if within circle


	; Corner 2
	; --------
	vec(*,0) = vec0(*,0) - corn_a(0,face)
	vec(*,1) = vec0(*,1) - corn_a(1,face)
	vec(*,2) = vec0(*,2) - corn_a(2,face)

	cs = vec # uv0
	norm = SQRT(TOTAL(vec*vec,2))
	flag(*,1) = (cs gt cs1*norm)


	; Corner 3
	; --------
	vec(*,0) = vec0(*,0) + corn_b(0,face)
	vec(*,1) = vec0(*,1) + corn_b(1,face)
	vec(*,2) = vec0(*,2) + corn_b(2,face)

	cs = vec # uv0
	norm = SQRT(TOTAL(vec*vec,2))
	flag(*,2) = (cs gt cs1*norm)


	; Corner 4
	; --------
	vec(*,0) = vec0(*,0) - corn_b(0,face)
	vec(*,1) = vec0(*,1) - corn_b(1,face)
	vec(*,2) = vec0(*,2) - corn_b(2,face)

	cs = vec # uv0
	norm = SQRT(TOTAL(vec*vec,2))
	flag(*,3) = (cs gt cs1*norm)


	flag0 = flag(*,0) * flag(*,1) * flag(*,2) * flag(*,3)

	i = WHERE(flag0 EQ 1)
		; get pixels where all corners within circle

	no_core = 0
	IF (i(0) NE -1) THEN core = interior(i) ELSE BEGIN
		no_core = -1
		GOTO, bdry
	ENDELSE
		; these are the "core" pixels


	FOR i=1,2*diff_res DO core = [[core],[core]]
	core = core * (4^diff_res)
	FOR i=1,(4^diff_res)-1 DO core(*,i) = core(*,i) + i
	core = REFORM(core,N_ELEMENTS(core))
		; expand all these pixels to higher resolution



bdry:	boundary = interior(WHERE(flag0 EQ 0))
		; get pixels where at least one corner is outside

	FOR i=1,2*diff_res DO boundary = [[boundary],[boundary]]
	boundary = boundary * (4^diff_res)
	FOR i=1,(4^diff_res)-1 DO boundary(*,i) = boundary(*,i) + i
	boundary = REFORM(boundary,N_ELEMENTS(boundary))
		; expand to higher resolution

	vec = coorconv(boundary,infmt='p',inco=rs,outfmt='u',outco='e')
		; convert to unit vectors


	cs = vec # uv0
	norm = SQRT(TOTAL(vec*vec,2))
	flag = (cs gt cs1*norm)
	boundary = boundary(WHERE(flag EQ 1))
		; calculate dot product
		; calculate norm of vector
		; check if within circle

	IF (no_core EQ 0) THEN interior = [core,boundary] ELSE $
	interior = [boundary]
		; add these boundary pixels to interior



ENDIF	; res > 6

interior = interior(SORT(interior))
		; sort by pixel number


RETURN,interior
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


