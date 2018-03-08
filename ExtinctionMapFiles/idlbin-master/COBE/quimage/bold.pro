FUNCTION bold, string
;
;+NAME/ONE LINE DESCRIPTION:
;    BOLD modifies a string so that, if printed, it will appear bolded.
;
;DESCRIPTION:
;    For IDL environments exclusive of PCs and the SunOs X-environment,
;    the supplied string argument will have certain escape sequences
;    appended to its head and tail so that the string will appear
;    bolded if it is printed out.
;
;CALLING SEQUENCE:
;    newstr = bold(oldstr)
;
;ARGUMENTS (I=input, o=output, []=optional): 
;    oldstr     I   string
;    newstr     O   string
;
;WARNINGS:
;    Bolding doesn't work on PC IDL and on text windows on a SPARCstation.
;
;EXAMPLE:
;    PRINT, bold('hi')
;#
;COMMON BLOCKS: none
;
;PROCEDURE:
;    For the appropriate environments, add a bolding escape sequence to the
;    head of the string, and add a reset (normal text) escape sequence to
;    the tail of the string.
;
;LIBRARY CALLS:  None.
;
;REVISION HISTORY:
;    Written by John Ewing (ARC)                       January 1992
;    SPR 11003 2 Jun 93   Change !version.os to !cgis_os. J Newmark
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;
;.TITLE
;Routine BOLD
;-
  IF(N_PARAMS(0) NE 1) THEN BEGIN
    MESSAGE, 'One argument (a string) must be supplied.', /CONT
    RETURN, ''
  ENDIF

  IF !cgis_os EQ 'windows' THEN RETURN, string

  esc = STRING(27b)
  wb = esc + '[1m'
  wz = esc + '[m'
  RETURN, wb + string  + wz

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


