PRO text
;+
; NAME:	
;	TEXT
; PURPOSE:
;       To switch back to the text window for Graph-On terminals.
; CATEGORY:
;	User interface.
; CALLING SEQUENCE:
;       TEXT
; INPUTS:
;       None.
; OUTPUTS:
;       None.
; COMMON BLOCKS:
;       None.
; RESTRICTIONS:
;       Verified for Graph-On terminals only.
; PROCEDURE:
;       Print out appropriate escape sequence.
; SUBROUTINES CALLED:
;       None.
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, February 1992.
;-
  PRINT, '', FORMAT = '(a,$)'
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


