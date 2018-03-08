FUNCTION coorconv, input, infmt=in_fmt0, inco=in_coord0, $
	           outfmt=out_fmt0, outco=out_coord0
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    COORCONV is a shell routine for coordinate conversion procedures
;
;DESCRIPTION:                                   
;    IDL function that acts as shell around the various coordinate
;    conversion routines.  It is the program that the user actually
;    calls.  Input and output formats are pixel number, decimal 
;    hour/latitude, longitude/latitude, and unit vector.
;
;CALLING SEQUENCE:
;    outcoor = coorconv(in_coor,infmt=in_fmt,outfmt=out_fmt, $
;                       inco=in_co, outco=out_co)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    in_coor     I     [arr]       Input coordinate array
;    infmt       I     str         Input coordinate format
;    outfmt      I     str         Output coordinate format
;    inco        I     str         Input coord. system or resolution
;    outo        I     str         Output coord. system or resolution
;    out_coor    O     [arr]       Output coordinate array
;
; where in_coor & out_coor arrays have the following dimensions
; for n points:
;
;          Pixel          - input_arr(n)
;          Decimal Hr/deg - input_arr(n,2)
;          lon/lat        - input_arr(n,2)
;          Unit vector    - input_arr(n,3)
;
;
; infmt & outfmt are specified as follows:
;
;          Pixel          - 'P'
;          Decimal Hr/dec - 'H'
;          Lon/Lat        - 'L'
;          Unit vector    - 'U'
;
;			 
; If pixel number is the input format then inco (outco) is specified
; as follows:
;
;          Firas (Res 6) - 'F' or 'R6'
;          DMR   (Res 6) - 'D' or 'R6'
;          Dirbe (Res 9) - 'B' or 'R9'
;          Resolution m  - 'Rm'  where m is <= 15
;
; (Note: The default coordinate system for pixels is ecliptic.)
;
; otherwise use:
;
;          Ecliptic   - 'E'
;          Galactic   - 'G'
;          Equatorial - 'Q'  (Epoch 2000)
;
; 
; Note: Lower case may also be used.
;
;WARNINGS:
;********** Input of (0,0) in ECLIPTIC coordinates will result in ******
;********** a division by zero error. Additionally, input of      ******
;********** coordinates that lie exactly on a pixel boundary will ******
;********** yield UNPREDICTABLE results.                          ******
;
;    In order to optimize throughput, this routine does minimal data 
;    checking.  Improper input such as "unit" vectors with norm not 
;    equal to one or pixel numbers with values outside the range
;    defined by the resolution will be processed without crashing,
;    however the output might lead to unexpected results. Therefore it
;    is up to the user to take proper care when generating and/or
;    passing the input.
;
;    The forward and inverse operations between pixel numbers and
;    pixel unit vector coordinates have been verified up to and
;    including resolution 12.  (The forward procedure is based on
;    the forward_cube.for program in the upx directory.)
;
;EXAMPLE: 
;
; To convert an array of sky coordinates in 'L' format from ecliptic
; coordinates to an array in 'U' format in galactic coordinates:
;
; output = coorconv(input,infmt='L',outfmt='U',inco='E',outco='G')
;
;
; To convert an array of sky coordinates in 'L' format to an array in
; 'U' format without changing coordinate systems:
;
; output = coorconv(input,infmt='L',outfmt='U')
;
;
; To convert an array of sky coordinates in 'H' format in equatorial
; coordinates to an array of Firas pixel numbers:
;
; output = coorconv(input,infmt='H',outfmt='P',inco='Q',outco='F') or
;
; output = coorconv(input,infmt='H',outfmt='P',inco='Q',outco='R7')
;
;
; To convert an array of resolution 8 pixel numbers to sky coordinates
; in 'L' format in galactic coordinates:
;
; output = coorconv(input,infmt='P',outfmt='L',inco='R8',outco='G')
;
;
; To convert an array of Dirbe pixel numbers to an array of DMR pixels:
;
; output = coorconv(input,infmt='P',outfmt='p',inco='P',outco='D')
;
; Note: You cannot "up res" pixel numbers from a lower resolution to a
; higher one.
;
;#
;COMMON BLOCKS:
;    sky_conv_com
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Uses CASE statements to route the input through the proper 
;    procedures.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: loadsky, skyconv, cconv
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Mar 92
;
; SPR 9741  Change PRINT statements to MESSAGE statements in
;           error situations.
; 2 JUN 1992
;
; SER 9984    Implement cconv (c module for p <-> uv <-> ll conversions)
; 9-SEP-1992  JMG
;
; SPR 10212   Generate conversions in call to CCONV module
;             Note: New CCONV exists as both IDL procedure and C module.
; 13-NOV-1992 JMG
;
; SPR 10318   Fix problem with scalar pixel entry
; 07-DEC-1992 JMG
;
; SPR 10348   Fix problem with integer format lon/lat and RA/dec input
;             Fix incorrect arithmetic in lon/lat <=> RA/dec conversion
; 14-DEC-1992 JMG
;
; SPR 10842   Fix problem with conversion of long integers to short int
; 28-SEP-1993 JMG
;
;.TITLE
;Routine COORCONV
;-
;
COMMON sky_conv_com, e2g,e2q,g2e,g2q,q2e,q2g,load_flag


in_fmt = strupcase(strmid(in_fmt0,0,1))

IF (keyword_set(out_fmt0) eq 0) THEN out_fmt0 = in_fmt0

out_fmt = strupcase(strmid(out_fmt0,0,1))

IF (keyword_set(in_coord0)) THEN $
	in_coord = strupcase(in_coord0) ELSE $
	in_coord = ' '

IF (keyword_set(out_coord0)) THEN $
	out_coord = strupcase(out_coord0) ELSE $
	out_coord = ' '

row = 0
long_flag = 0
scalar_flag = 0
		; initialize row, long_flag, and scalar_flag flags

IF (N_ELEMENTS(input) EQ 0) THEN BEGIN
	MESSAGE,'Input is undefined',/cont
	output = 0
	GOTO, exit
ENDIF

sz = SIZE(input)


CASE in_fmt OF

'P' :	begin
	n = n_elements(input)

	IF (N EQ 1) THEN BEGIN
		scalar_flag = -1
		input = [input]
		sz = SIZE(input)
	ENDIF

	IF (sz(2) NE 3) THEN BEGIN
		long_flag = -1
		input = LONG(input)
	ENDIF

	end

'L' :	begin
	n = n_elements(input) / 2
	end

'H' :	begin
	n = n_elements(input) / 2
	end

'U' :	begin
	n = n_elements(input) / 3
	end

else :	begin
	str = 'Unknown or improper input format: '
	str = str + '"' + in_fmt + '"'
	message,str,/cont
	output = input
	goto, exit
	end

endcase

IF ((out_fmt NE 'P') AND (out_fmt NE 'L') AND (out_fmt NE 'H') AND $
    (out_fmt NE 'U')) THEN BEGIN
	str = 'Unknown or improper output format: '
	str = str + '"' + out_fmt + '"'
	message,str,/cont
	output = input
	goto, exit
ENDIF


sz = SIZE(input)
row = 0
IF (((in_fmt EQ 'L') OR (in_fmt EQ 'H') OR (in_fmt EQ 'U')) AND $
    (sz(0) EQ 1)) THEN BEGIN
	row = 1
	temp_input = input
	input = TRANSPOSE(input)
ENDIF			; flip array if row vector input



IF (in_fmt eq 'P') THEN BEGIN

	case strmid(in_coord,0,1) of

	'B' :	begin
		in_res = 9
		end	

	'D' :	begin
		in_res = 6
		end	

	'F' :	begin
		in_res = 6
		end	

	'R' :	begin
		in_res = fix(strmid(in_coord,1,2))
		end

	else :	begin
		str = 'Improper input resolution designation: '
		str = str + '"' + strmid(in_coord,0,1) + '"'
		message,str,/cont
		output = input
		goto, exit
		end

	endcase

ENDIF


IF (out_fmt eq 'P') THEN BEGIN

	case strmid(out_coord,0,1) of

	'B' :	begin
		out_res = 9
		end	

	'D' :	begin
		out_res = 6
		end	

	'F' :	begin
		out_res = 6
		end	

	'R' :	begin
		out_res = fix(strmid(out_coord,1,2))
		end

	else :	begin
		str = 'Improper output resolution designation: '
		str = str + '"' + strmid(out_coord,0,1) + '"'
		message,str,/cont
		output = input
		goto, exit
		end

	endcase

ENDIF


IF (in_fmt eq 'P' AND out_fmt eq 'P') THEN BEGIN

	diff_res = out_res - in_res

	IF (diff_res ge 0) THEN BEGIN

		str = 'Input resolution must be larger than output resolution'
		message,str,/cont
		output = input
		goto, exit

	ENDIF ELSE BEGIN

		output = ishft(input,2*diff_res)
		goto, exit

	ENDELSE

ENDIF

	
if (n_elements(load_flag) eq 0) then loadsky
		; load coordinate transformation matrices (E, G, Q)
		; if not done so already


case in_fmt of

'P' :	begin

	case out_fmt of

	'L' :	begin

		uvec = cconv(input,in_res)
		uvec = skyconv(uvec, inco = 'E', outco=out_coord)
		output = cconv(uvec,-1)

		end

	'H' :	begin

		uvec = cconv(input,in_res)
		uvec = skyconv(uvec, inco = 'E', outco=out_coord)
		output = cconv(uvec,-1)

		output(*,0) = output(*,0) / 15

		end

	'U'  :	begin

		uvec = cconv(input,in_res)

		output = skyconv(uvec, inco = 'E', outco=out_coord)

		end

	endcase

	end



'L' :	begin

	case out_fmt of

	'H' :	begin

		IF (in_coord NE out_coord) THEN BEGIN

		uvec = cconv(input*1.,-1)
		uvec = skyconv(uvec, inco=in_coord, outco=out_coord)
		output = cconv(uvec,-1)
		output(*,0) = output(*,0) / 15.

		ENDIF ELSE BEGIN

		output = input * 1.
		output(*,0) = input(*,0) / 15.

		ENDELSE

		end

	'U' :	begin

		uvec = cconv(input*1.,-1)
		output = skyconv(uvec, inco=in_coord, outco=out_coord)

		end

	'P' :	begin

		uvec = cconv(input*1.,-1)
		uvec = skyconv(uvec, inco=in_coord, outco='E')
		output = cconv(uvec,out_res)

		end

	'L' :	begin

		IF (in_coord NE out_coord) THEN BEGIN

		uvec = cconv(input*1.,-1)
		uvec = skyconv(uvec, inco=in_coord, outco=out_coord)
		output = cconv(uvec,-1)

		ENDIF ELSE output = input * 1.

		end

	endcase

	end


'H' :	begin

	case out_fmt of

	'U' :	begin

		uvec = cconv([[input(*,0)*15.],[input(*,1)]],-1)
		output = skyconv(uvec, inco=in_coord, outco=out_coord)

		end

	'P' :	begin

		uvec = cconv([[input(*,0)*15.],[input(*,1)]],-1)
		uvec = skyconv(uvec, inco=in_coord, outco='E')
		output = cconv(uvec,out_res)

		end

	'L' :	begin

		IF (in_coord NE out_coord) THEN BEGIN

		uvec = cconv([[input(*,0)*15.],[input(*,1)]],-1)
		uvec = skyconv(uvec, inco=in_coord, outco=out_coord)
		output = cconv(uvec,-1)

		ENDIF ELSE BEGIN

		output = input * 1.
		output(*,0) = input(*,0) * 15.

		ENDELSE

		end

	'H' :	begin

		IF (in_coord NE out_coord) THEN BEGIN

		uvec = cconv([[input(*,0)*15.],[input(*,1)]],-1)
		uvec = skyconv(uvec, inco=in_coord, outco=out_coord)
		output = cconv(uvec,-1)
		output(*,0) = output(*,0) / 15.

		ENDIF ELSE output = input * 1.

		end


	endcase

	end


'U' :	begin

	case out_fmt of

	'P' :	begin

		uvec = skyconv(input, inco=in_coord, outco='E')
		output = cconv(uvec,out_res)

		end

	'L' :	begin

		uvec = skyconv(input, inco=in_coord, outco=out_coord)

		output = cconv(uvec,-1)

		end

	'H' :	begin

		uvec = skyconv(input, inco=in_coord, outco=out_coord)

		output = cconv(uvec,-1)

		output(*,0) = output(*,0) / 15

		end

	'U' :	begin
		output = skyconv(input, inco=in_coord, outco=out_coord)
		end

	endcase

	end

endcase

exit:

IF (row EQ 1) THEN input = TRANSPOSE(input)
		; flip input back to row vector if input is row vector


IF (in_fmt EQ 'P' AND long_flag EQ -1) THEN input = FIX(input)
		; convert pixel input back to short int if necessary

IF (in_fmt EQ 'P' AND scalar_flag EQ -1) THEN input = input(0)
		; convert pixel input back to scalar if necessary

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


