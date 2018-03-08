FUNCTION hilight, string
;+NAME/ONE LINE DESCRIPTION:
;    HILIGHT allows a string to be printed in reverse video.
;
;DESCRIPTION:
;    This function returns a string which consists of a certain escape
;    sequence appended by the input string, and then appended by another
;    escape sequence.  If this output string is then PRINTed on a
;    Tektronix-compatible terminal, then the input string will appear
;    in reverse video.
;
;CALLING SEQUENCE:
;    PRINT, hilight('hello')
;
;ARGUMENTS (I=input, o=output, []=optional): 
;    string     I   scalar   str        Arbitrary character string.
; <returned quantity> scalar str        A string.
;
;WARNINGS: none.
;
;EXAMPLE:
;    PRINT, hilight('hello')
;#
;COMMON BLOCKS: none
;
;PROCEDURE:
;    RETURN a string which is the bold-text escape sequence concatenated
;    with the input string concatenated with the normal-text escape
;    sequence.
;
;LIBRARY CALLS:  None.
;
;REVISION HISTORY:
;    Written by John Ewing.
;    SPR 10383  Jan 07 93  Add more comments.  J Ewing
;.TITLE
;Routine HILIGHT
;-
  IF(N_PARAMS(0) NE 1) THEN BEGIN
    MESSAGE, 'One argument (a string) must be supplied.'
    RETURN, ''
  ENDIF
  esc = STRING(27b)
  wr = esc + '[7m'
  wz = esc + '[m'
  RETURN, wr + string + wz
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


