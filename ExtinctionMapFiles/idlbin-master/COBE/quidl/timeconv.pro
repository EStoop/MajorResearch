function timeconv, in_time, ret_stat, infmt=in_fmt, outfmt=out_fmt
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    TIMECONV is the shell routine for time conversion procedures.
;
;DESCRIPTION:                                   
;    IDL function that acts as shell around the various time conversion
;    routines.  It is the program that the user actually calls.
;
;CALLING SEQUENCE:
;    out_time = timeconv(in_time,ret_stat, infmt=in_fmt, outfmt=out_fmt)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    in_time     I     [arr/scl]   Input time array
;    infmt       I     str         Input time format (see table below)
;    outfmt      I     str         Output time format
;    out_time    O     [arr/scl]   Output time array
;    ret_stat    O     int [arr]   return status array [0 = good, 
;                                                      -1 = bad]
;
;
;
;                            Time Formats
;                            ============
;
;
;  'Zulu' (Format Designation: 'z')
;   -------------------------------
;   A 16 character string consisting of the year (1858-9999), the day
;   of the year (001-365/366), the hour (00-23), the minute (00-59), 
;   the second (00-59), and the number of milliseconds (000-999).  If
;   the string is less than 16 characters long, it is assumed that the
;   first two characters give the year within the 20th century (00-99)
;   unless an 'e' (for extended) is appended to the end in which case 
;   it is assumed that the first FOUR characters give the year. 
;
;
;   Example:	'2090123091530053'	(123rd day of 2090)
;		'86264134550416'	(264th day of 1986)
;		'1895002e'		(2nd day of 1895)
;
;
;
;  'Julian' (Format Designation: 'j')
;   ---------------------------------
;   A (non-negative) real*8 floating point variable giving the number
;   of days (with fractional part) from 00hr 17-Nov-1858.  (This is
;   actually the "reduced" Julian date; "true" Julian is computed from
;   12hr 01-Jan-4713 BC.) 
;
;
;   Example:	40235.25
;
;
;
;  'Vax' (Format Designation: 'v')
;   ------------------------------
;   A character string consisting of the day of the month, the month 
;   (using the standard 3 character abrieviation), the (full) year,
;   and the hour-minute-second information as in 'Zulu'.  This last
;   portion is optional in which case '00:00:00.000' is assumed. 
;   The month may also be specified in lower case. 
;
;
;   Examples:	'03-MAR-1989:14:07:23.952'
;		'10-oct-1993'
;
;
;  'Sec' (Format Designation: 's')
;   ------------------------------
;   A (non-negative) real*8 floating point variable giving the number
;   of seconds (with fractional part) since 0hr 01-JAN-1968.  It must
;   be lessthan 2.14748365e9, corresponding to:
;   19-JAN-2036:03:14:00.000. 
;
;
;   Example:	2.035456675e+8
;
;
;  'ADT' (Format Designation: 'a')
;   ------------------------------
;   VMS system time, giving the number of 100 nanosecond intervals
;   since 00:00 hours, 17-NOV-1858.  It is stored in a two-dimensional
;   longword (integer*4) array.  Only times after 01-JAN-1985
;   00:00hr should be used.
;
;   Example:	[46556L,9642143L]  (09-FEB-1990:08:28:04.84)
;
;
;WARNINGS:
;    The user should check the return status array to see whether 
;    element in output array was processed without an error.
;
;EXAMPLE: 
;    To convert from Julian to zulu time:
;
;    out_time = timeconv([34500.35,36783.9],r_stat,infmt='j',outfmt='z')
;
;    out_time will be a string array with 2 elements
;
;    To convert from adt to 'sec' time:
;
;    in_time = [[345345l,9350000l],[27l,9372260l],[-74634l,9367600l]]
;
;    out_time = timeconv(in_time,ret_stat,infmt='a',outfmt='s')
;
;    out_time will be a double array with 3 elements.
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Uses CASE statements to route the input through the proper 
;    procedures.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: zulu2jul, zulu2vax, zulu2sec, jul2zulu, 
;                        vax2zulu, sec2zulu, adt2zulu
;                        zulu2adt
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Dec 91
;
;    SPR 10414  Change loop index (i) from short to long integer
;    01/05/93   J.M. Gales
;
; SPR 9734  - TIMECONV expects array input
;             Solution: Force all single element input to arrays
;             Fixed: 28 May 1992 JMG
;
; SPR 9750  - Check that input arrays are of the correct type.
;             Solution: Use SIZE statement to determine input data type.
;             Fixed: 5 Jun 1992 JMG
;
; SPR 9759  - Error messages are now displayed using MESSAGE statement
;             rather than PRINT statement.
;             Fixed: 11 Jun 1992 JMG
;
; SER 10039   Add ADT output capabilities (calls zulu2adt c module)
; 23-SEP-1992 J.M. Gales
;
; SER 10212   Call ADT2ZULU and ZULU2ADT routines to perform conversion
;             between ADT and ZULU formats
;             ADT2ZULU and ZULU2ADT exist as both IDL procedures and
;             c modules.
; 13-NOV-1992 J.M. Gales
; SPR 10845  2-May-1994  if input is a scalar then do not alter it to array.
;             J. Newmark
;
;.TITLE
;Routine TIMECONV
;-
;
n = n_elements(in_time)
scflag=0
if (n eq 1) then begin
  sz=size(in_time)
  if sz(0) eq 0 then scflag=1
  in_time = [in_time]	; force scalar to array
endif

in_fmt = strupcase(strmid(in_fmt,0,1))
out_fmt = strupcase(strmid(out_fmt,0,1))

if (n eq 1) and (in_fmt eq 'A') then begin
	MESSAGE,'Insufficient number of elements',/CONT
	out_time = ''
	goto, exit
endif

if (in_fmt eq 'A') then n = n / 2

ret_stat = intarr(n)
sz = SIZE(in_time)
data_type = sz(sz(0)+1)

case in_fmt of

'Z' :	begin

	IF (data_type NE 7) THEN BEGIN
		MESSAGE,'Zulu Times must be strings.',/CONT
		out_time = in_time
		GOTO, exit
	ENDIF
	
	case out_fmt of

	'J' :	begin
		zulu2jul,in_time,out_time,ret_stat
		end

	'V' :	begin
		zulu2vax,in_time,out_time,ret_stat
		end

	'S' :	begin
		zulu2sec,in_time,out_time,ret_stat
		end

	'A' :	begin
		zulu2adt,in_time,out_time,ret_stat
		end

	else :	begin
		str = 'Unknown or improper output format ' + out_fmt
		MESSAGE,str,/CONT
		out_time = in_time
		end

	endcase

	end

'J' :	begin

	IF (data_type EQ 7) THEN BEGIN
		MESSAGE,'Julian Times must be numeric.',/CONT
		out_time = in_time
		GOTO, exit
	ENDIF

	case out_fmt of

	'Z' :	begin
	
		jul2zulu,in_time,out_time,ret_stat

		end

	'V' :	begin

		jul2zulu,in_time,temp_time,ret_stat
		zulu2vax,temp_time,out_time,ret_stat

		end

	'S' :	begin

		jul2zulu,in_time,temp_time,ret_stat
		zulu2sec,temp_time,out_time,ret_stat

		end

	'A' :	begin

		jul2zulu,in_time,temp_time,ret_stat
		zulu2adt,temp_time,out_time,ret_stat

		end

	else :	begin
		str = 'Unknown or improper output format ' + out_fmt
		MESSAGE,str,/CONT
		out_time = in_time
		end

	endcase

	end

'V' :	begin

	IF (data_type NE 7) THEN BEGIN
		MESSAGE,'Vax Times must be strings.',/CONT
		out_time = in_time
		GOTO, exit
	ENDIF

	case out_fmt of

	'Z' :	begin

		vax2zulu,in_time,out_time,ret_stat

		end

	'J' :	begin

		vax2zulu,in_time,temp_time,ret_stat
		zulu2jul,temp_time,out_time,ret_stat

		end

	'S' :	begin

		vax2zulu,in_time,temp_time,ret_stat
		zulu2sec,temp_time,out_time,ret_stat

		end

	'A' :	begin

		vax2zulu,in_time,temp_time,ret_stat
		zulu2adt,temp_time,out_time,ret_stat

		end

	else :	begin
		str = 'Unknown or improper output format ' + out_fmt
		MESSAGE,str,/CONT
		out_time = in_time
		end

	endcase

	end

'S' :	begin

	IF (data_type EQ 7) THEN BEGIN
		MESSAGE,'Second Times must be numeric.',/CONT
		out_time = in_time
		GOTO, exit
	ENDIF

	case out_fmt of

	'Z' :	begin

		sec2zulu,in_time,out_time,ret_stat

		end

	'J' :	begin

		sec2zulu,in_time,temp_time,ret_stat
		zulu2jul,temp_time,out_time,ret_stat

		end

	'V' :	begin

		sec2zulu,in_time,temp_time,ret_stat
		zulu2vax,temp_time,out_time,ret_stat

		end

	'A' :	begin

		sec2zulu,in_time,temp_time,ret_stat
		zulu2adt,temp_time,out_time,ret_stat

		end

	else :	begin
		str = 'Unknown or improper output format ' + out_fmt
		MESSAGE,str,/CONT
		out_time = in_time
		end

	endcase

	end


'A' :	begin

	case out_fmt of

	'Z' :	begin

		adt2zulu,in_time,out_time,ret_stat

		end

	'J' :	begin

		adt2zulu,in_time,temp_time,ret_stat
		zulu2jul,temp_time,out_time,ret_stat

		end

	'V' :	begin

		adt2zulu,in_time,temp_time,ret_stat
		zulu2vax,temp_time,out_time,ret_stat

		end

	'S' :	begin

		adt2zulu,in_time,temp_time,ret_stat
		zulu2sec,temp_time,out_time,ret_stat

		end

	else :	begin
		str = 'Unknown or improper output format ' + out_fmt
		MESSAGE,str,/CONT
		out_time = in_time
		end

	endcase

	end


else :	begin 
	str = 'Unknown or improper input format ' + in_fmt
	MESSAGE,str,/CONT
	out_time = in_time
	end

endcase

if (out_fmt eq 'Z') then begin
for i=0l,n-1 do begin
;
	if (strmid(out_time(i),0,2) eq '19') then $
	    out_time(i) = strmid(out_time(i),2,14)
			; trim to 14 character if 20th century
endfor
endif

if n_elements(out_time) eq 1 then begin
        if scflag eq 1 then in_time = in_time(0)
	out_scalar = out_time(0)
	out_time = out_scalar
			; convert to scalar if 1-dim array
endif

exit:
;
return, out_time
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


