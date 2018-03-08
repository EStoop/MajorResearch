pro projbld, n,proj=proj,coord=coord,euler=euler,ilut,jlut,mask
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    PROJBLD generates projection lookup tables.
;
;DESCRIPTION:
;    This IDL procedure generates the projection lookup tables that are
;    used by the USRPROJ facility to create and display reprojections
;    of the native skycube.  The user enters half the dimension of the
;    vertical component of the projection array (all such arrays have a
;    2:1 aspect ratio), the projection type (Aitoff, Global Sinusoidal,
;    Mollweide), the reference coordinate system (Ecliptic, Galactic,
;    Equatorial), and the set of Euler angles defining the user 
;    coordinate system.  The procedure returns the three lookup tables
;    (in unfolded skycube coordinates) to be used by the USRPROJ facility.
;
;CALLING SEQUENCE:
;    pro projbld, n,proj=proj,coord=coord,[euler=euler],ilut,jlut,mask
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    n           I     int                 half dimension of vertical
;                                          components of proj arrays
;    proj        I     string              Projection Type:
;                                          ('A','S','M','P')
;    coord       I     string              Reference Coordinate System:
;                                          ('E','G','Q')
;    [euler]     I     flt arr             Euler angle vector (3 elem)
;                                          Default: [0,0,0]
;    [ilut]      O     int arr             I*2 array containing skycube
;                                          column numbers for reproj.
;                                          Dimensions: (4*n,2*n) 
;    [jlut]      O     int arr             I*2 array containing skycube
;                                          row numbers for reproj.
;                                          Dimensions: (4*n,2*n) 
;    [mask]      O     byte arr            Face number / Boundary mask
;                                          (face # within boundary
;                                           of projection, 255 outside
;                                           boundary.)
;
;WARNINGS:
;
;    To generate a 512 by 256 set of lookup tables takes about
;    90 sec of CPU time and 1283 blocks of storage (IDL save set).
;
;    A 1024 by 512 set takes about 380 CPU sec and 5123 blocks.
;
;EXAMPLES: 
;
; The Euler angles are defined in terms of a left-handed coordinate
; system with positive angles specifying clockwise rotations of the
; coordinate axes.  The first element gives the rotation about the
; z-axis of the reference coordinate system, the second, the rotation
; about the new y-axis and the third, the rotation about the z-axis
; of the current system.
;
;
; To generate a 256 by 128 Aitoff galactic-centered projection:
;
; projbld,64,proj='a',coord='g',ia,ja,ma
;
; ____
;
; To generate a 512 by 256 "anti-ecliptic" Mollweide projection:
;
; projbld,128,proj='m',coord='e',euler=[180,0,0],ae_i,ae_j,ae_mask
;
; ____
;
; To generate a 512 by 256 "CMB dipole (l = 267.7, b = 47)" Global
; Sinusoidal projection:
;
; projbld,128,proj='s',coord='g',euler=[267.7,43,0],dp_i,dp_j, $
; dp_mask
;
; ____
;
; To generate a 1024 by 512 "Leo-centered (RA 10h30m dec 20)" Aitoff
; projection:
;
; projbld,128,proj='a',coord='q',euler=[157.5,-20,0],leo_i,leo_j, $
; leo_mask
;
; ____
;
;#
;COMMON BLOCKS:
;    cmn
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Generate rotation matrix for Euler angle vector.
;    For given projection type generate unit vectors for each
;    projection pixel (row by row) and calculate skycube face
;    number,column and row.
;    Reflect column number to give Right-T tables.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: xyz2ijm
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Jun 92
;
; SER 9956    Improved Documentation
; 9-SEP-92    J.M. Gales
;
;.TITLE
;Routine PROJBLD
;-
;
COMMON cmn, col,row,cube_side,offx,offy,left,right

ilut = INTARR(4*n,2*n)
jlut = INTARR(4*n,2*n)
mask = BYTARR(4*n,2*n)
bound = BYTARR(4*n,2*n)

proj = STRUPCASE(STRMID(proj,0,1))

invc = 1 / (4.0 * n*n)
pi = 4 * atan(1.d0)
half_pi = pi / 2
d2r = !dpi / 180

offx = [0,0,1,2,3,0]
offy = [2,1,1,1,1,0]

cube_side = 256
	; default resolution 9

col = indgen(2*n)
	; generate column numbers
left = 2*n-1 - col
	; reflect about y-axis
right = col + (2*n)
	; shift to right of y-axis


IF (KEYWORD_SET(euler) EQ 0) THEN euler = [0.,0.,0.]

; Generate rotation matrix
; ------------------------
cos_z1 = cos(d2r*euler(0))
sin_z1 = sin(d2r*euler(0))

cos_y1 = cos(d2r*euler(1))
sin_y1 = sin(d2r*euler(1))

cos_z2 = cos(d2r*euler(2))
sin_z2 = sin(d2r*euler(2))

m_z1 = [[cos_z1,sin_z1,0],[-sin_z1,cos_z1,0],[0,0,1]]
m_y1 = [[cos_y1,0,-sin_y1],[0,1,0],[sin_y1,0,cos_y1]]
m_z2 = [[cos_z2,sin_z2,0],[-sin_z2,cos_z2,0],[0,0,1]]

rot_mat = TRANSPOSE(m_z1 # m_y1 # m_z2)


CASE proj OF

'A' :	BEGIN

	FOR row=0l,2*n-1 DO BEGIN
		; for each row of the projection

	IF (row/10 EQ float(row)/10) THEN print,row

	ya = row - n + 0.5
		; shift y_axis origin to middle of projection
		; offset to middle of pixel

	xa = col + 0.5
		; offset to middle of pixel

	r = ya / xa
	s_prime = 1 - invc * (xa*xa + 4*ya*ya)

	c1_sqr = (1 + (2*r*s_prime)^2) / (1 + (2*r)^2)
	c1 = SQRT(c1_sqr)
	c2 = s_prime / c1
	c3 = 2 * c2*c2 - 1

	x = c1 * c3
	y = 2 * s_prime * SQRT(ABS(1-c2*c2))
	z = SQRT(ABS(1-c1_sqr)) * (ya/abs(ya))
		; generate unit vectors

	bound(left,row) = (s_prime GE 0.0)
	bound(right,row) = (s_prime GE 0.0)
		; fill in boundary mask

	xyz2ijm, n,coord,x,y,z,rot_mat,ilut,jlut,mask

	ENDFOR


	END


'S' :	BEGIN


	FOR row=0l,2*n-1 DO BEGIN

	IF (row/10 EQ float(row)/10) THEN print,row

	ys = half_pi * (row - n + 0.5) / n

	cs = COS(ys)
	sn = SQRT(1 - cs*cs)
	max_x = pi * cs

	xs = pi * (col + 0.5) / (2*n)

	x = cs * COS(xs/cs)
	y = cs * SIN(xs/cs)
	z = replicate(sn,2*n)* (ys/abs(ys))

	bound(left,row) = (xs LE max_x)
	bound(right,row) = (xs LE max_x)

	xyz2ijm, n,coord,x,y,z,rot_mat,ilut,jlut,mask

	ENDFOR

	end


'M' :	BEGIN

	FOR row=0l,2*n-1 DO BEGIN

	IF (row/10 EQ float(row)/10) THEN print,row

	ym = (row - n + 0.5) / n

	sn_beta = ym
	cs_beta = SQRT(1 - sn_beta^2)
	beta = ASIN(sn_beta)
	sn_2beta = 2 * sn_beta * cs_beta
	sn_del = (2*beta + sn_2beta) / pi
	cs_del = SQRT(1 - sn_del^2)

	xm = (col + 0.5) / n

	s_prime = 1 - 0.25 * (xm*xm + 4*ym*ym)

	alph = 0.5 * xm * pi / cs_beta
	cs_alph = COS(alph)
	sn_alph = SQRT(1 - cs_alph^2)

	x = cs_del * cs_alph
	y = cs_del * sn_alph
	z = replicate(sn_del,2*n)

	bound(left,row) = (s_prime GE 0.0)
	bound(right,row) = (s_prime GE 0.0)

	xyz2ijm, n,coord,x,y,z,rot_mat,ilut,jlut,mask


	ENDFOR


	END


'P' :	BEGIN

	col = indgen(n)
		; generate column numbers

	flip = n - 1 - col
		; flip about lon=90 & 270

	FOR row=0l,n-1 DO BEGIN
		; for each row of the projection

		IF (row/10 EQ float(row)/10) THEN print,row

		xp = n - row - 0.5
		yp = n - col - 0.5
			; offset to middle of pixel

		theta = atan(yp/xp)
		phi = asin(1 - ((xp/n) / cos(theta))^2)

		x = cos(theta) * cos(phi)
		y = sin(theta) * cos(phi)
		z = sin(phi)
			; generate unit vectors

		one = WHERE(z GT 0)

		x = [x, x(flip), x, x(flip)]
		y = [y,-y(flip),-y, y(flip)]
		z = [z, z(flip),-z,-z(flip)]
			; extend along row

		bnd = BYTARR(n)
		bnd(one) = 1B
		bound(*,row) = [bnd,bnd(flip),bnd,bnd(flip)]
		bound(*,2*n-1-row) = [bnd,bnd(flip),bnd,bnd(flip)]
			; fill in boundary mask


		; Bottom Part of Map
		; ------------------
		vec = [[x],[y],[z]]
		vec = vec # rot_mat
		vec = coorconv(vec,infmt='u',inco=coord,outfmt='u',outco='e')

		axisxy, vec,face,x0,y0
		icube = fix(x0 * cube_side) < 255
		jcube = fix(y0 * cube_side) < 255

		ilut(*,row) = offx(face) * cube_side + icube
		jlut(*,row) = offy(face) * cube_side + jcube
		mask(*,row) = face	



		; Top Part of Map
		; ---------------
		vec = [[-x],[y],[z]]
		vec = vec # rot_mat
		vec = coorconv(vec,infmt='u',inco=coord,outfmt='u',outco='e')

		axisxy, vec,face,x0,y0
		icube = fix(x0 * cube_side) < 255
		jcube = fix(y0 * cube_side) < 255

		ilut(*,2*n-1-row) = offx(face) * cube_side + icube
		jlut(*,2*n-1-row) = offy(face) * cube_side + jcube
		mask(*,2*n-1-row) = face	

	ENDFOR

	END

ENDCASE

	ilut = 1023 - ilut	; switch to right-T

	ilut = ilut * bound
	jlut = jlut * bound
	mask = mask * bound + 255B * (1B - bound)

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


