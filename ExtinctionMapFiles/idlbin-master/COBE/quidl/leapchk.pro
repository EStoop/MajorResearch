function leapchk, yr
;+                                                                  
;  NAME:
;    leapchk
;
;  PURPOSE:                                   
;    Determines where year is a leap year
;
;  CALLING SEQUENCE:
;    leap = leapchk(yr)
;
;  INPUT:
;    yr - year to check
;
;  OUTPUT:
;    leap - leap year flag (-1 if leap, 0 if non-leap)
;
;  SUBROUTINES CALLED:
;    None
;
;  REVISION HISTORY
;    J.M Gales
;    Jan 92
;-
;
leap = 0
;
if (yr mod 4 eq 0) then leap = -1
;
if (yr mod 100 eq 0) then leap = 0
;
if (yr mod 400 eq 0) then leap = -1
;
return,leap
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


