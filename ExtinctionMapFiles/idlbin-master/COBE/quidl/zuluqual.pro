pro zuluqual, zulu_time, qual_flag
;+                                                                  
;  NAME:
;    zuluqual
;
;  PURPOSE:                                   
;    Tests "data quality" of zulu date
;    (Proper day of year, hour, minute, sec values)
;
;  CALLING SEQUENCE:
;    zuluqual,zulu_time,qual_flag
;
;  INPUT:
;    Date - The format being 'yydddhhmmsssss'.
;
;  OUTPUT:
;    quality flag ( -1 if OK, 0 if bad)
;    (Error message is printed to screen for each "failure")
;
;  SUBROUTINES CALLED:
;    leapchk
;
;  REVISION HISTORY
;    J.M Gales
;    Dec 91
;
; SPR 11331	Fix problem with maximum date check for leap years
; J. M. Gales	09/24/93
;
;-
;
COMMON zulu, yr,dy,hr,mn,sc

qual_flag = -1				; Assume good time


IF (STRPOS(zulu_time,'.') NE -1) THEN BEGIN
	str = 'Zulu Times must not include decimal points: ' + zulu_time
	MESSAGE,str,/CONT
	qual_flag = 0
	GOTO, exit
ENDIF

yr = fix(strmid(zulu_time,0,4))			; get year
dy = fix(strmid(zulu_time,4,3))			; get day of year
hr = fix(strmid(zulu_time,7,2))			; get hour
mn = fix(strmid(zulu_time,9,2))			; get minute
sc = long(strmid(zulu_time,11,5))/1000.0	; get seconds


if ((leapchk(yr) eq -1) and (dy gt 366)) then begin
	str = 'Improper day of year in ZULU date ' + STRING(dy)
	MESSAGE,str,/INFO
	qual_flag = 0			; bad leap year day
endif

if ((leapchk(yr) eq 0) and (dy gt 365)) then begin
	str = 'Improper day of year in ZULU date ' + STRING(dy)
	MESSAGE,str,/INFO
	qual_flag = 0			; bad non-leap year day
endif

if (dy lt 1) then begin
	str = 'Improper day of year in ZULU date ' + STRING(dy)
	MESSAGE,str,/INFO
	qual_flag = 0			; bad day
endif

if (hr gt 23) then begin
	str = 'Improper hour in ZULU date ' + STRING(hr)
	MESSAGE,str,/INFO
	qual_flag = 0			; bad hour
endif

if (mn gt 59) then begin
	str = 'Improper minute in ZULU date ' + STRING(mn)
	MESSAGE,str,/INFO
	qual_flag = 0			; bad minute
endif

if (sc gt 59.999) then begin
	str = 'Improper second in ZULU date ' + STRING(sc)
	MESSAGE,str,/INFO
	qual_flag = 0			; bad second
endif

exit:

return
end
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


