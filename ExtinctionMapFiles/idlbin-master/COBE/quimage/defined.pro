FUNCTION defined, var
;+NAME/ONE LINE DESCRIPTION:
;    DEFINED tells whether or not the supplied variable is defined.
;
;DESCRIPTION:
;    DEFINED returns a value of 1 if the supplied variable is defined,
;    else it returns a 0.
;
;CALLING SEQUENCE:
;    IF(defined(variable)) THEN ...
;
;ARGUMENTS (I=input, o=output, []=optional):
;    var        I   any type            Arbitrary input variable.
;  <returned value> scalar integer      Boolean value (0 or 1).
;
;WARNINGS: none.
;
;EXAMPLE:
;    IF(defined(variable)) THEN ...
;#
;COMMON BLOCKS: none
;
;PROCEDURE:
;    Call the SIZE command, giving it the supplied argument.  Look at
;    the first two values in the array returned by SIZE.  If both those
;    values are zero, then that means the argument is undefined, else
;    it is defined.
;
;LIBRARY CALLS:  None.
;
;REVISION HISTORY:
;    Written by John Ewing.
; SPR 10383  Jan 06 93  Add more comments.  J Ewing
;.TITLE
;Routine DEFINED
;-
;
  sz = SIZE(var)
  IF(sz(0)+sz(1) EQ 0) THEN RETURN, 0 ELSE RETURN, 1
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


