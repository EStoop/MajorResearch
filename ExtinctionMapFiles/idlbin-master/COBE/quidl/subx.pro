function subx, a, b, n, borrow_bit
;+                                                                  
;  NAME:
;    subx
;
;  PURPOSE:                                   
;    Performs extended word integer subtraction
;
;  CALLING SEQUENCE:
;    difference = subx(minuend, subtrahend, length, borrow_bit)
;
;  INPUT:
;    minuend - extended integer (32 * n bits)
;    subtrahend - extended integer (32 * n bits)
;    length - length of extended integer in longwords
;    borrow bit - one if minuend lt subtrahend, zero otherwise
;
;  OUTPUT:
;    difference - extended integer (32 * n bits)
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
borrow = bytarr(n+1)
a0 = a
sign_mask = '80000000'x	; sign bit set, all others zero

sign_bit_a = lonarr(n)
sign_bit_b = b and sign_mask

for i=0,n-1 do begin

	if (a0(i) eq 0) then borrow(i+1) = borrow(i)
			; to propagate borrow bit through string of 0's
			; to handle a=[0,0], b=[-1,1] properly
			; otherwise, borrow(n) not set

	a0(i) = a0(i) - borrow(i)

	c(i) = a0(i) - b(i)

	sign_bit_a(i) = a0(i) and sign_mask

	if ( (sign_bit_a(i) eq sign_bit_b(i)) and $
	     (a0(i) lt b(i)) ) then borrow(i+1) = 1
			; if sign bits are equal and minuend < subtrahend

	if (sign_bit_a(i) eq 0) and (sign_bit_b(i) eq sign_mask) then $
		borrow(i+1) = 1
			; if msb(a)=0 and msb(b)=1

endfor

borrow_bit = borrow(n)

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


