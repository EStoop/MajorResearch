pro emul, a,b,c
;+                                                                  
;  NAME:
;    emul
;
;  PURPOSE:                                   
;    Performs quadword integer multiplication
;
;  CALLING SEQUENCE:
;    emul, multiplicand1, multiplicand2, product
;
;  INPUT:
;    multiplicand1 - quadword (2 dim longword array)
;    multiplicand2 - longword
;
;  OUTPUT:
;    product - quadword (2 dim longword array)
;
;  SUBROUTINES CALLED:
;    qs1_l
;    addx
;
;  REVISION HISTORY
;    J.M. Gales
;    Delivered 13-NOV-1992  SPR 10212
;-
;
c = lonarr(2)
a0 = a
mask = long(1)

for i=0,31 do begin

	if ((b and mask) eq mask) then c = addx(a0,c,2)

	mask = ishft(mask,1)	; build new mask

	a0 = qs1_l(a0)		; quadword shift left 1 bit

endfor

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


