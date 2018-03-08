	FUNCTION sum, array, index
;
;+NAME/ONE-LINE DESCRIPTION:
;    SUM totals an array over one of its dimensions
;
; DESCRIPTION:
;    This function performs a summation over either any one dimension of 
;    an N-dimensional array, or of the whole array (i.e. summing all 
;    elements).  When called with an array index 0,1,2... it will sum 
;    over the first, second, third, etc., index of the array; when 
;    called with no index number, it totals the whole array. 
;
; CALLING SEQUENCE:
;    result = SUM (array, [index])
;
; ARGUMENTS (I=input, O=output, []=optional):
;    RESULT    O    any data type    Output array or scalar
;    ARRAY     I    any data type    Input array, any dimensions
;   [INDEX]   [I]   integer          Dimension over which to sum
;
; WARNINGS:
;    Note that the numbering of the dimensions differs from that 
;    in the TOTAL routine!  In SUM, the first dimension of the array 
;    has index 0 (not 1!) for consistency with the rest of IDL index 
;    notation.
;
; EXAMPLE:
;    Suppose ARR is an array of dimensions (3,4,5).  The call
;
;                 UIDL> OUT = SUM(ARR,0)
;
;    will yield an array OUT of dimensions (4,5).  Similarly, 
;
;                 UIDL> OUT = SUM(ARR,2)
;
;    will cause OUT to have dimensions (3,4).  Finally, note that
;
;                 UIDL> OUT = SUM(ARR)
;
;    will sum all elements of ARR, so that OUT will be a scalar.
;#
; COMMON BLOCKS:
;    None.
;
; PROCEDURE:
;    Uses the enhanced TOTAL routine (IDL 2.3 and later).  For 
;    some awful reason TOTAL departs from IDL convention by numbering 
;    the dimension indices from 1 instead of zero, but SUM "fixes" 
;    this.  This is also consistent with the linkimage SUM.FOR, 
;    which it supercedes.  The calling sequence is identical.
;
; ALGORITHMS, LIBRARY CALLS:
;    None.
;
; MODIFICATION HISTORY:
;    Written 16 July 92 by Rich Isaacman, General Sciences Corp.
;
;.title
;.Function SUM
;-
	ON_ERROR, 2
	ndims = SIZE(array)
	IF (index+1 GT ndims(0)) THEN MESSAGE, $
		'Specified array index is greater than array dimension.'
;
;   If index is negative or nonexistent, sum over all elements
;
	IF ((index LT 0) or (N_ELEMENTS(index) EQ 0)) THEN $
		RETURN, TOTAL(array)
	RETURN, TOTAL(array,index+1)
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


