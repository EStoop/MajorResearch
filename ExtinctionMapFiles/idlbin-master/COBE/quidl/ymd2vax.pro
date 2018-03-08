FUNCTION YMD2VAX, Y, M, D
;+
; NAME:
;     YMD2VAX
;
; PURPOSE:
;      Convert from year, month, day numbers to date string.
;
; CALLING SEQUENCE:
;          DATE = YMD2VAX(Y,M,D)
;
; INPUTS:
;          Y = year number (like 1986 or 86)
;          M = month number (1 - 12)
;          D = day of month number (like 1 to 31).
;
; OUTPUTS:
;          Date
;
; EXAMPLES:
;          Y=YMD2VAX(1985,9,24)
;          Where Y will contain 24-SEP-1985
;
; COMMON BLOCKS:
;          None
;
; RESTRICTIONS:
;          None
;
; MODIFICATION HISTORY:
;          R. Sterner.
;          Urmila Prasad (ARC) August 1991
;-


	IF Y LT 0 THEN BEGIN
	  PRINT,'Error in YMD2VAX: invalid year.'
	  RETURN, -1
	ENDIF
	IF Y LT 100 THEN Y = Y + 1900
	IF (M LT 1) OR (M GT 12) THEN BEGIN
	  PRINT,'Error in YMD2VAX: invalid month.'
	  RETURN, -1
	ENDIF
	IF (D LT 1) OR (D GT MONTHDAYS(Y,M)) THEN BEGIN
	  PRINT,'Error in YMD2VAX: invalid month day.'
	  RETURN, -1
	ENDIF

	DATE = STRTRIM( D, 2) + '-' + $
	  STRMID( 'JanFebMarAprMayJunJulAugSepOctNovDec', 3*(M-1), 3) + $
	  '-' + STRTRIM( Y, 2)

	RETURN, DATE

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


