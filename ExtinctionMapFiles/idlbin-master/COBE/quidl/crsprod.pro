FUNCTION crsprod,a,b
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    CRDPROD computes the cross product of two 3-element vectors
;
;DESCRIPTION:
;    This function returns the cross product of two 3-element vectors.
;    These vectors can be either row or column vectors.
;
;CALLING SEQUENCE:
;    c = crsprod(a,b)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    a           I     3 elem vector
;    b           I     3 elem vector
;    c           O     3 elem vector       c = a X b
;
;#
;COMMON BLOCKS:
;    None
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    None
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   May 93
;    Initial Delivery   SPR 11174
;
;
;.TITLE
;Routine CRSPROD
;-
;
n = N_ELEMENTS(a)
IF (n NE 3) THEN BEGIN
	MESSAGE, 'First Input Vector must have 3 elements',/CONT
	crs = -1
	GOTO, exit
ENDIF

n = N_ELEMENTS(b)
IF (n NE 3) THEN BEGIN
	MESSAGE, 'Second Input Vector must have 3 elements',/CONT
	crs = -1
	GOTO, exit
ENDIF

sz = SIZE(a)
IF (sz(0) EQ 2) THEN aa = TRANSPOSE(a) ELSE aa = a

sz = SIZE(b)
IF (sz(0) EQ 2) THEN bb = TRANSPOSE(b) ELSE bb = b

crs = [aa(1)*bb(2)-aa(2)*bb(1), $
       aa(2)*bb(0)-aa(0)*bb(2), $
       aa(0)*bb(1)-aa(1)*bb(0)]
		; calculate cross product

exit:

RETURN,crs
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


