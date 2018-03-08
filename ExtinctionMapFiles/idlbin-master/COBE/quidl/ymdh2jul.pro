PRO YMDH2JUL,YR,MN,DAY,HR,VAX2ZULU
;+
; NAME:
;       YMDH2JUL
;
; PURPOSE:
;       Converts Gregorian dates to VAX2ZULU days   
;
; CALLING SEQUENCE:
;       YMDH2JUL,YR,MN,DAY,HR,VAX2ZULU
;
; INPUTS:
;       YR = Year (integer)  
;       MN = Month (integer 1-12)
;       DAY = Day  (integer 1-31) 
;       HR  = Hours and fractions of hours of universal time (U.T.)
;
; OUTPUTS:
;       VAX2ZULU = VAX2ZULU date (double precision) 
;
; EXAMPLE:
;       To find the VAX2ZULU Date at 1978 January 1, 0h (U.T.)
;       YMDH2JUL,1978,1,1,0.,VAX2ZULU
;       will give VAX2ZULU = 2443509.5
;      Notes:
;       (1) YMDH2JUL will accept vector arguments 
;       (2) JULDATE is an alternate procedure to perform the same function
;
; COMMON BLOCKS:
;        None
;
; RESTRICTIONS:
;        None
;
; MODIFICATION HISTORY:
;       Converted to IDL from Don Yeomans Comet Ephemeris Generator,
;       B. Pfarr, STX, 6/15/88
;       Urmila Prasad (ARC) Aug 1991
;-
if n_params(0) lt 4 then begin
	print,string(7B),'CALLING SEQUENCE: YMDH2JUL,YR,MN,DAY,HR,VAX2ZULU
        return
endif
yr = long(yr) & mn = long(mn) &  day = long(day)	;Make sure integral
L = (mn-14)/12		;In leap years, -1 for Jan, Feb, else 0
VAX2ZULU = day - 32075l + 1461l*(yr+4800l+L)/4 + $
         367l*(mn - 2-L*12)/12 - 3*((yr+4900l+L)/100)/4
VAX2ZULU = double(VAX2ZULU) + (HR/24.0D) - 0.5D
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


