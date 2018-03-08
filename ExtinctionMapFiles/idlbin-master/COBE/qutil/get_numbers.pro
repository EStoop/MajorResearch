;+
; NAME: GET-NUMBERS
;
; PURPOSE: To extract numbers from a string with various delimiters
; between the numbers.  The delimiters can be anything and in any
; quantity except plus and minus sign and decimal points (+ - .)'s.
;
; CATEGORY: Utility
;
; CALLING SEQUENCE: get_numbers,input,numbers,ct
;
; INPUTS:
;          input        a string, scalar or array with numbers in it.
;
; OUTPUTS:
;	   numbers	a floating point array with just the numbers.
;
;          ct           a scalar value of how many numbers are in the array.
;#
; WARNINGS: None 
;
; COMMON BLOCKS: None
;
; INCLUDE FILES: None
;
; RESTRICTIONS: None
;
; SUBROUTINES CALLED:
;
; MODIFICATION HISTORY: Created J. Newmark May 27, 1993
;---------------------------------------------------------------------
PRO get_numbers,input,numbers,ct

a = input(0)
n = strlen(a)
b = byte(a)

pos = where(b gt 47 and b lt 58,c)

if (c eq 0) then begin
	numbers = [0]
	ct = 0
	return
endif

if (c eq n) then begin
	numbers = [float(a)]
	ct = 1
	return
endif

pos = where((b ne 43) and (b ne 45) and (b ne 46) and $
	not ((b gt 47) and (b lt 58)),c)

if (c eq 0) then begin
	numbers = [float(a)]
	ct = 1
	return
endif else b(pos) = 32

b = byte(strtrim(strcompress(string(b)),2))
pos = [0,where(b eq 32,c),strlen(string(b))-1]
ct = c + 1
st = strarr(ct)
numbers = intarr(ct)
for i = 0,c do numbers(i) = fix(string(b(pos(i):pos(i+1))))

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


