pro strform, num, str, l_str
;+                                                                  
;  NAME:
;    strform
;
;  PURPOSE:                                   
;    Form right-justified string of desired length from numeric input
;
;  CALLING SEQUENCE:
;    strform, num, str, l_str
;
;  INPUT:
;    num - numeric value
;    l_str - length of output string
;
;  OUTPUT:
;    str - output string
;
;  SUBROUTINES CALLED:
;    None
;
;  REVISION HISTORY
;    J.M Gales
;    Jan 92
;-
;
z_str = '00000'

	str = strtrim(string(num),2)	; convert number to string format
					; trim leading and trailing blanks

	str_len = strlen(str)		; get new string length

	if (str_len lt l_str) then $
	   str = strmid(z_str,1,l_str-str_len) + str
					; left pad with blanks to desired
					; string length
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


