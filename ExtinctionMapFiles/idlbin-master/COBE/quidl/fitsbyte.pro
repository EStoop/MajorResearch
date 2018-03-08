PRO fitsbyte, im, HOST=host
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     FITSBYTE byteswaps between FITS and internal machine formats.
;
;DESCRIPTION:
;     To convert an image array from FITS byte ordering to the ordering
;     used by the computer, and vice-versa
;
;CALLING SEQUENCE:
;     fitsbyte, im [, /host]
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     im    I,O  arr   The image array (both input and output).
;     HOST [I]   key   If present and non-zero, the image is converted 
;                      to "host" (machine internal) form.  
;                      Otherwise, it is converted to "network", or FITS,
;                      form.
;
;WARNINGS:
;     1.  This routine alters the input array values on output.
;
;EXAMPLE:
;     To convert data array DATA from FITS byte ordering to 
;     machine-internal ordering:
;         fitsbyte,data,/host
;
;#
;COMMON BLOCKS:
;     None.
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     The byteorder procedure is used to perform the conversion.  The
;     form of this conversion depends upon the datatype of the array.
;     If it is of I*2 format, a short integer swap is performed.  If it 
;     is of I*4 format, a long integer swap is performed.  Otherwise, 
;     nothing is done.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS,ETC.:
;     Uses IDL routine BYTEORDER.     
;
;MODIFICATION HISTORY:
;     Written by Michael R. Greason, STX, 21 August 1990.
;
; SPR 9616
;.TITLE
;Routine FITSBYTE
;-
;			Initialize.
;
s = size(im)
dtype = s(s(0) + 1)			; The datatype.
h = keyword_set(host)			; Is the host keyword set?
i2 = 2					; Size type code for I*2.
i4 = 3					; Size type code for I*4.
;
;			If I*2, do a short integer swap.
;
IF (dtype EQ i2) THEN BEGIN
	IF (h) THEN byteorder, im, /ntohs $	; Convert to Host format.
	       ELSE byteorder, im, /htons	; Convert to Network format.
ENDIF ELSE BEGIN
;
;			If I*4, do a long integer swap.
;
	IF (dtype EQ i4) THEN BEGIN
		IF (h) THEN byteorder, im, /ntohl $ ; Convert to Host format.
		       ELSE byteorder, im, /htonl ; Convert to Network format.
	ENDIF
ENDELSE
;
;			Done.
;
RETURN
END

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


