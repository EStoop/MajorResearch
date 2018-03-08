FUNCTION MONTHDAYS,YR,M
;+
; NAME:
;    MONTHDAYS
;
; PURPOSE:
;    Gives number of days in a month given year and month.
;
; CALLING SEQUENCE:
;      DAYS = MONTHDAYS(YR,MON)
;
; INPUTS:
;      YR = year,  MON = month number (jan = 1).
;      If MON is 0 then return month days for all of year as array.
;
; OUTPUTS:
;      DAYS = number of days in month (like 31).
;      Or array for MON = 0.
;
; EXAMPLES:
;         DAYS = MONTHDAYS(84,2)  
;         Returns the value 29 in the variable DAYS.
;
; COMMON BLOCKS:
;      None
;
; RESTRICTIONS:
;      Only for 1901 to 2099.
;
; MODIFICATION HISTORY:
;      RES  14 Aug, 1985.
;      Urmila Prasad (ARC) August 1991
;
;-


	DYS = [0,31,28,31,30,31,30,31,31,30,31,30,31]
	IF (YR MOD 4) EQ 0 THEN DYS(2) = 29

	IF M EQ 0 THEN RETURN, DYS
	RETURN, DYS(M)

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


