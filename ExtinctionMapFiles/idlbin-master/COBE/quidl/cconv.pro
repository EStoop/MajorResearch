FUNCTION cconv,input,code
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    CCONV is a shell around pixel to unit vector to lon/lat
;
;DESCRIPTION:                                   
;    IDL function that acts as shell around the pixel <--> unit
;    vector, pixel <--> lon/lat, and unit vector <--> lon/lat
;    routines.  It is called by COORCONV.
;
;CALLING SEQUENCE:
;    output = convert(input,code)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    input       I     arr         Input coordinate array
;    code        I     int         If > 0 then res else
;                                  uv <--> ll conversion
;    output      O     [arr]       Output coordinate array
;
;
;WARNINGS:
;    There is a C module linkimage by the same name that is a shell 
;    around C functions that also perform conversions between pixels
;    and unit vectors and lon/lat.  If this linkimage is not implemented
;    then this IDL procedure is compiled and run instead.
;
;EXAMPLE: 
;
; Not user routine
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Uses CASE statements to route the input through the proper 
;    procedures.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: pix2uv, uv2pix, ll2uv, uv2ll
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Oct 92
;    Delivered 13-NOV-1992  SPR 10212
;
;    SPR 10316   Correct call to UV2LL
;    07-dec-1992 JM Gales
;-
;
sz = SIZE(input)

IF (sz(0) EQ 1) THEN in_fmt = 'P'

IF (sz(0) EQ 2) THEN BEGIN
	CASE sz(2) OF
	  2: in_fmt = 'L'
	  3: in_fmt = 'U'
	ENDCASE
ENDIF

IF (code GT 0) THEN out_fmt = 'P' ELSE BEGIN
	CASE in_fmt OF
	  'L': out_fmt = 'U'
	  'U': out_fmt = 'L'
	ENDCASE
ENDELSE

CASE in_fmt OF

'P' :	begin

	output = pix2uv(input,code)
	IF (out_fmt EQ 'L') THEN output = uv2ll(output)

	end


'L' :	begin

	output = ll2uv(input)
	IF (out_fmt EQ 'P') THEN output = uv2pix(output,code)

	end


'U' :	begin

	IF (out_fmt EQ 'P') THEN output = uv2pix(input,code)
	IF (out_fmt EQ 'L') THEN output = uv2ll(input)

	end

ENDCASE


RETURN, output
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


