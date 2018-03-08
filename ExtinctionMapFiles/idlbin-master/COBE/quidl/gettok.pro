FUNCTION GETTOK,ST,CHAR
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     GETTOK retrieves a segment of a character string.
;
;DESCRIPTION:
;     IDL function to retrieve the first part of a character string
;     until the character CHAR is encountered.
;
;CALLING SEQUENCE:
;     Y = GETTOK(ST,CHAR)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     Y             O   str        output character string which
;                                  terminates at CHAR.
;
;     ST            I   str        input character string
;
;     CHAR          I   str        character at which output character
;                                  string is to be terminated.
;
;WARNINGS:
;     NONE
;
;EXAMPLE:
;     ST='ABC=999'
;     Y=GETTOK(ST,'=')
;     PRINT,Y
;       ABC
;
;#
;COMMON BLOCKS:
;     NONE
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     GETTOK uses the IDL function STRPOS to determine the position of
;     CHAR in the character string ST, and the IDL funtion STRMID to
;     extract the character string segment which occurs before CHAR.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     NONE
;
;MODIFICATION HISTORY
;     Written by D. Lindler, April 1986
;     Revised by Vidya Sagar, Applied Research Corp.   Aug 1991
;     Improved documentation
;
; SPR 9616
;.TITLE
;ROUTINE GETTOK
;
;-
; if char is a blank treat tabs as blanks
;
	tab='	'
	while strpos(st,tab) ge 0 do begin
		pos=strpos(st,tab)
		strput,st,' ',pos
	end

	;
	; find character in string
	;
	pos=strpos(st,char)
	if pos eq -1 then begin	;char not found?
		token=st
		st=''
		return,token
	endif

	;
	; extract token
	;
	token=strmid(st,0,pos)
	len=strlen(st)
	if pos eq (len-1) then st='' else st=strmid(st,pos+1,len-pos-1)

	;
	;  Return the result.
	;
	return,token
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


