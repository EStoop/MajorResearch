FUNCTION YDN2VAX,YEAR,DAY
;+
; NAME:
;	YDN2VAX
;
; PURPOSE:
;	YDN2VAX converts day-of-year to ordinary Gregorian date;
;	works on cumulative day from 1-JAN-1979 (for instance),
;	as well as more well-mannered days-of-year.
;
; CALLING SEQUENCE:
;       D_STRING = YDN2VAX(Year,day)
;
; INPUTS:
;	Year - Years after 1900 may either be written out in full
;	(e.g. 1986) or with just the last two digits (e.g. 86)
;       Day - number of days after Jan 1 of the specified year
;
; OUTPUTS:
;	D_STRING - String giving date in format '13-MAR-1986'
;
; EXANPLES:
;       Y=YDN2VAX(1990,23)
;       Will return Y with the value equal to 23-Jan-1990.
;
; COMMON BLOCKS:
;       None
;
; RESTRICTIONS:
;	Will only work for years after 1900.
;
; MODIFICATION HISTORY:
;       D.M. fecit  24 October,1983
;       Urmila Prasad (ARC) July 1991
;-
IF DAY LE 0 THEN BEGIN
	D_STRING = '%DATE-F-DAY.LE.ZERO'
ENDIF ELSE BEGIN
	LAST_DAY = [31,59,90,120,151,181,212,243,273,304,334,365]
	LD = [0,INTARR(11)+1]
	DAY_OF_YEAR = DAY
	MONTHS = 'JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC'
	IF YEAR GT 1899 THEN YR = YEAR - 1900 ELSE YR = YEAR
	N_DAYS = 365 + ((YR MOD 4) EQ 0)
	WHILE DAY_OF_YEAR GT N_DAYS DO BEGIN
		DAY_OF_YEAR = DAY_OF_YEAR - N_DAYS
		YR = YR + 1
		N_DAYS = 365 + ((YR MOD 4) EQ 0)
	END
	END_DATE = '-19' + STRTRIM(YR,2)
	IF (YR MOD 4) EQ 0 THEN LAST_DAY = LAST_DAY + LD
	LAST_MONTH = DAY_OF_YEAR LE LAST_DAY
	WHERE_LD = WHERE(LAST_MONTH)
	N_MONTH = !ERR
	IF N_MONTH EQ 12 THEN BEGIN
		D_STRING = STRTRIM(DAY_OF_YEAR,2) + '-JAN' + END_DATE
	ENDIF ELSE BEGIN
		LAST_MONTH = WHERE_LD(0)
		MONTH = STRMID(MONTHS,3*LAST_MONTH,3)
		DAY_OF_MONTH = DAY_OF_YEAR - LAST_DAY(LAST_MONTH-1)
		D_STRING = STRTRIM(DAY_OF_MONTH,2) + '-' + MONTH + END_DATE
	END
END
RETURN,D_STRING
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


