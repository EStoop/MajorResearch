FUNCTION merge,vec1,vec2
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    MERGE combines two vectors into one, eliminating duplicates.
;
;DESCRIPTION:
;    This function combines two integer vectors together.  Each value
;    occurs only once.  The second vector can be empty in which case
;    duplicates in the first vector are eliminated.  The output vector
;    is sorted in ascending order.
;
;CALLING SEQUENCE:
;    vec = merge(vec1,vec2)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    vec1        I     integer vector
;    vec2        I     integer vector
;    vec         O     integer vector
;
;
;EXAMPLES:
;
;    a = [2,1,5,4],  b = [3,2,3,5,7]
;
;    c = merge(a,b)
;
;    c = [1,2,3,4,5,7]
;
;    d = merge(b)
;
;    d = [2,3,5,7]
;
; Note: The merged vector will always be sorted in ascending order.
;       Therefore, DO NOT USE MERGE if order must be preserved.
;
;#
;COMMON BLOCKS:
;    None
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    None
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   July 93
;    Initial Delivery   SPR 11199
;
;
;.TITLE
;Routine MERGE
;-
;
n1 = N_ELEMENTS(vec1)
n2 = N_ELEMENTS(vec2)

IF (n1 EQ 0 AND n2 EQ 0) THEN BEGIN
	MESSAGE,'Vectors undefined',/CONT
	vec = [-1]
	GOTO, exit
ENDIF

IF (n1 EQ 0) THEN BEGIN
	MESSAGE,'Vector 1 undefined',/CONT
	vec = [-1]
	GOTO, exit
ENDIF

IF (n2 EQ 0) THEN BEGIN
	sz = SIZE(vec1)
	type = sz(sz(0)+1)
	IF (type GE 2 AND type LE 3) THEN BEGIN
	 	vec = WHERE(HISTOGRAM(vec1,MIN=0) GT 0)
	ENDIF ELSE BEGIN
		MESSAGE,'Vector 1 must be of integer type',/CONT
		vec = [-1]
		GOTO, exit
	ENDELSE
ENDIF ELSE BEGIN
	sz1 = SIZE(vec1)
	type1 = sz1(sz1(0)+1)
	sz2 = SIZE(vec2)
	type2 = sz2(sz2(0)+1)
	IF (type1 GE 2 AND type1 LE 3 AND type2 GE 2 AND type2 LE 3) THEN BEGIN
		vec = WHERE(HISTOGRAM([vec1,vec2],MIN=0) GT 0)
	ENDIF ELSE BEGIN
		MESSAGE,'Vectors must be of integer type',/CONT
		vec = [-1]
		GOTO, exit
	ENDELSE
ENDELSE

exit:

RETURN,vec
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


