function qs1_l, a
;+                                                                  
;  NAME:
;    qs1_l
;
;  PURPOSE:                                   
;    Shift quadword integer 1 bit to left
;
;  CALLING SEQUENCE:
;    a = qs1_l(a)
;
;  INPUT:
;    a - quadword integer (2 dim longword array)
;
;  OUTPUT:
;    a - shifted quadword integer
;
;  SUBROUTINES CALLED:
;    None
;
;  REVISION HISTORY
;    J.M Gales
;    Delivered 13-NOV-1992  SPR 10212
;-
;
a(1) = ishft(a(1),1)
a(1) = a(1) + ishft((a(0) and '80000000'x),-31)
a(0) = ishft(a(0),1)

return,a
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

