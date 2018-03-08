function addx, a, b, n
;+                                                                  
;  NAME:
;    addx
;
;  PURPOSE:                                   
;    Performs extended word integer addition
;
;  CALLING SEQUENCE:
;    sum = addx(addend1, addend2, length)
;
;  INPUT:
;    addend1 - extended integer (32 * n bits)
;    addend2 - extended integer (32 * n bits)
;    length - length of extended integer in longwords
;
;  OUTPUT:
;    sum - extended integer (32 * n bits)
;
;  SUBROUTINES CALLED:
;    None
;
;  REVISION HISTORY
;    J.M. Gales
;    Delivered 13-NOV-1992  SPR 10212
;-
;
c = lonarr(n)
carry = bytarr(n+1)
a0 = a
sign_mask = '80000000'x	; sign bit set, all others zero

for i=0,n-1 do begin

	c(i) = a(i) + b(i) + carry(i)

	if ( (a(i) lt 0 and b(i) lt 0) or $
	     ((a(i) lt 0 and b(i) gt 0) and c(i) ge 0) or $
	     ((a(i) gt 0 and b(i) lt 0) and c(i) ge 0) ) then $
	carry(i+1) = 1b

endfor

return, c
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


