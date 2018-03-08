pro ediv,a,b,quot,rem
;+                                                                  
;  NAME:
;    idl$ediv
;
;  PURPOSE:                                   
;    Performs quadword integer division
;
;  CALLING SEQUENCE:
;    idl$ediv, dividend, divisor, quotient, remainder
;
;  INPUT:
;    dividend - quadword (2 dim longword array)
;    divisor - quadword (2 dim longword array)
;
;  OUTPUT:
;    quotient - quadword (2 dim longword array)
;    remainder - quadword (2 dim longword array)
;
;  SUBROUTINES CALLED:
;    qs1_l
;    qs1_r
;    subx
;
;  REVISION HISTORY
;    J.M. Gales
;    Delivered 13-NOV-1992  SPR 10212
;-
;
a0 = a				; temporary dividend
b0 = b				; temporary divisor
quot = lonarr(2)		; define quotient array
bit0 = 0			; set bit counter to zero


	if (b(0) eq 0) and (b(1) eq 0) then begin	; divisor = 0

		print,'Division by zero'
		rem = a
		goto, exit

	endif


	if (a(1) eq 0) and (b(1) eq 0) 	 then begin	; standard integer
							; division

	first_1 = 31			; start at msb
	mask = '80000000'x		; generate mask

	while (a0(0) and mask) ne mask do begin		; while not found
		first_1 = first_1 - 1			; decrement bit
		mask = ishft(mask,-1)			; build new mask
	endwhile


;	*Shift divisor to left until first nonzero bits align*

	while (b0(0) and mask) ne mask do begin
		b0(0) = ishft(b0(0),1)	; quadword shift left 1 bit
		bit0 = bit0 + 1		; increment bit counter
	endwhile


	for i=bit0,0,-1 do begin

		j = i / 32		; get longword index
		bit = i mod 32		; get bit number within longword

		rem = subx(a0,b0,2,borrow)
					; subtract shifted divisor from
					; remaining dividend

		if (borrow eq 0) then begin	; if a0 > b0
			a0 = rem		; new remaining dividend
			quot(j) = quot(j) or (2l)^bit
						; set bit in quotient
		endif

		b0(0) = ishft(b0(0),-1)		; shift 'divisor' right
						; one bit

	endfor

	rem = a0				; get final remainder

		goto, exit

	endif


;	if (a(1) lt b(1)) then begin			; a < b
;
;		rem = a
;		goto, exit
;
;	endif


;	-----------------------
;	start quadword division
;	-----------------------

;	*Find first nonzero bit in dividend starting from left*

	first_1 = 31			; start at msb
	mask = '80000000'x		; generate mask

	while (a0(1) and mask) ne mask do begin		; while not found
		first_1 = first_1 - 1			; decrement bit
		mask = ishft(mask,-1)			; build new mask
	endwhile


;	*Shift divisor to left until first nonzero bits align*

	while (b0(1) and mask) ne mask do begin
		b0 = qs1_l(b0)		; quadword shift left 1 bit
		bit0 = bit0 + 1		; increment bit counter
	endwhile


	for i=bit0,0,-1 do begin

		j = i / 32		; get longword index
		bit = i mod 32		; get bit number within longword

		rem = subx(a0,b0,2,borrow)
					; subtract shifted divisor from
					; remaining dividend

		if (borrow eq 0) then begin	; if a0 > b0
			a0 = rem		; new remaining dividend
			quot(j) = quot(j) or (2l)^bit
						; set bit in quotient
		endif

		b0 = qs1_r(b0)			; shift 'divisor' right
						; one bit

	endfor

	rem = a0				; get final remainder

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


