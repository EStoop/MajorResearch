pro xyz2ijm, n,coord,x,y,z,rot_mat,i_arr,j_arr,mask
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    xyz2ijm calculates projection lut values from unit vectors.
;
;DESCRIPTION:
;    This IDL procedure is called by PROJBLD to generate the projection
;    lookup tables values.  The set of unit vectors are first rotated
;    back to the reference coordinate system (ecliptic, galactic, or
;    equatorial), then rotated to the native ecliptic coordinate
;    system of the quadcube.  The AXISXY procedure returns the cor-
;    responding face, face column, and face row number which are
;    then transformed into unfolded skycube coordinates.  The above
;    procedure is repeated for the input vectors reflected about
;    the y-axis.
;
;CALLING SEQUENCE:
;    xyz2ijm, n,coord,x,y,z,rot_mat,i_arr,j_arr,mask
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    n           I     int                 half dimension of vertical
;                                          components of proj arrays
;    coord       I     string              Reference Coordinate System:
;                                          ('E','G','Q')
;
;    x           I     flt arr             x coordinate array
;    y           I     flt arr             y coordinate array
;    z           I     flt arr             z coordinate array
;
;    rot_mat     I     flt arr             3 by 3 euler angle rotation
;                                          matrix
;
;    i_arr       O     int arr             I*2 array containing skycube
;                                          column numbers for reproj.
;                                          Dimensions: (4*n,2*n) 
;    j_arr       O     int arr             I*2 array containing skycube
;                                          row numbers for reproj.
;                                          Dimensions: (4*n,2*n) 
;    mask        O     byte arr            Face number / Boundary mask
;                                          (face # within boundary
;                                           of projection, 255 outside
;                                           boundary.)
;
;WARNINGS:
;
;    None.
;
;EXAMPLES: 
;
;    Not user routine.
;
;#
;COMMON BLOCKS:
;    cmn
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Rotate unit vector via Euler angle rotation matrix
;    Convert to ecliptic coordinates
;    Get face, face column, & face row # from axisxy
;    Generate unfolded cube column and row #
;    Load into lookup table arrays
;    Repeat above steps for unit vector reflected about y-axis
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Calls coorconv, axisxy
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Jun 92
;
; SER 9956    Improved Documentation/Banner
; 9-SEP-92    J.M. Gales
;
;-
COMMON cmn, col,row,cube_side,offx,offy,left,right

	vec = [[x],[y],[z]]
	vec = vec # rot_mat
	vec = coorconv(vec,infmt='u',inco=coord,outfmt='u',outco='e')

	axisxy, vec,face,x0,y0
	icube = fix(x0 * cube_side) < 255
	jcube = fix(y0 * cube_side) < 255

	i_arr(left,row) = offx(face) * cube_side + icube
	j_arr(left,row) = offy(face) * cube_side + jcube
	mask(left,row) = face	

	vec = [[x],[-y],[z]]
	vec = vec # rot_mat
	vec = coorconv(vec,infmt='u',inco=coord,outfmt='u',outco='e')

	axisxy, vec,face,x0,y0
	icube = fix(x0 * cube_side) < 255
	jcube = fix(y0 * cube_side) < 255

	i_arr(right,row) = offx(face) * cube_side + icube
	j_arr(right,row) = offy(face) * cube_side + jcube
	mask(right,row) = face	

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


