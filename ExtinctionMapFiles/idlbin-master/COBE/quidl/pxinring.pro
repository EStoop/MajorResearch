FUNCTION pxinring,cen=cen,coord=coord,out_ang=out_ang,in_ang=in_ang,res=res
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    PXINCIRC generates pixels within a given ring
;
;DESCRIPTION:
;    This IDL function generates the pixel numbers of those pixels
;    whose centers fall between two angular radii of a given
;    point.  The angles are expressed in degrees and the center can
;    be given in pixel number, unit vector, or long/lat.
;
;CALLING SEQUENCE:
;    px = pxincirc(cen=cen,in_ang=in_ang,out_ang=out_ang,res=res, $
;                  [coord=coord])
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    cen         I     int/flt vector      Circle center
;    coord       I     string              Center coordinate system
;                                          ['E' (default), 'G', 'Q']
;    in_ang      I     flt                 Inner angular radius (in degrees)
;    out_ang     I     flt                 Outer angular radius (in degrees)
;    res         I     int                 Resolution
;    px          O     int/long arr        Pixels within ring
;
;WARNINGS:
;
;EXAMPLE: 
;
; To extract the pixels within a ring defined by two circles, one 
; 10 degrees, the other 20 degrees about a point centered on 
; pixel 3457 for a DIRBE (res=9) skymap:
; 
; px = pxinring(cen=3457,in_ang=10,out_ang=20,res=9)
;
; To extract the pixels within a ring defined by two circles, one 
; 20 degrees, the other 25 degrees about a point at 30 degrees galactic
; longitude, -23 degrees galactic latitude for a FIRAS (res=6) skymap:
; 
; px = pxincirc(cen=[30,-23],coord='g',in_ang=20,out_ang=25,res=6)
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Call pxincirc for outer circle
;    Call pxincirc for inner circle
;    Generate long int array
;    Set outer pixels to 1
;    Reset inner pixels to 0
;    Ring pixels where long array = 1
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: pxincirc
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Apr 93
;    Initial Delivery   SPR 10825
;
;
;.TITLE
;Routine PXINRING
;-
;

IF (KEYWORD_SET(res) EQ 0) THEN BEGIN
	MESSAGE,'Resolution must be supplied.',/CONT
	ring = -1
	GOTO, exit
ENDIF

IF (KEYWORD_SET(coord)) THEN BEGIN
	out = pxincirc(cen=cen,ang=out_ang,res=res,coord=coord)
	in = pxincirc(cen=cen,ang=in_ang,res=res,coord=coord)
ENDIF ELSE BEGIN
	out = pxincirc(cen=cen,ang=out_ang,res=res)
	in = pxincirc(cen=cen,ang=in_ang,res=res)
ENDELSE


min_val = MIN(out)
ring = LONARR(MAX(out)-min_val+1)

ring(out-min_val) = 1
ring(in-min_val) = 0

ring = WHERE(ring EQ 1) + min_val

exit:

RETURN, ring
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


