pro sec2zulu, sec_time, zulu_time, istat
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    sec2zulu - Convert from # of seconds from fiduicial date to zulu time
;
;DESCRIPTION:                                   
;    IDL routine called by TIMECONV to convert from number of seconds 
;    from 01-JAN-1968:00:00:00.000 to zulu time.
;
;CALLING SEQUENCE:
;    sec2zulu, sec_time, zulu_time, istat
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    sec_time    I     dbl [arr]   seconds from 01/01/68 - 00hr
;    zulu_time   O     str [arr]   zulu date - "yydddhhmmsssss"
;    istat       I/O   int [arr]   status array
;
;WARNINGS:
;    This routine should not be called directly.
;
;EXAMPLE: 
;    Not a user-callable routine
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES)
;    Subtract years (in seconds) from 1968 until remainder less than sec/yr
;    Get day of year, hour, minute, and secs in sequence by dividing 
;    by appropriate factor and passing remainder to next operation.
;    Convert to strings and concatanate.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: leapchk, strform
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Dec 91
;-
;
n = n_elements(sec_time)			; get input array size
;
zulu_time = strarr(n)				; define output array
;
num_sec_n = 3600l * 24 * 365
num_sec_l = 3600l * 24 * 366
;
;
	for j = 0l,n-1 do begin			; loop over all elements

	if (sec_time(j) lt 0 or sec_time(j) gt 2.14748365e9) then begin

	zulu_time(j) = '0000000000000000'
	istat(j) = -1
	str$ = STRCOMPRESS('Number of seconds ' + STRING(sec_time(j)) + $
			   ' out of range')
	MESSAGE,str,/INFO

	endif else begin

	sec_rem = sec_time(j)
	for yr=1968,9999 do begin

	if (leapchk(yr)) then num_sec_ln = num_sec_l else $
				 num_sec_ln = num_sec_n
			; decide whether leap year or not

	if (sec_rem ge num_sec_ln) then begin

	sec_rem = sec_rem - num_sec_ln

	endif else begin

	strform, yr, yr_str, 4			; form 4 char string
	goto, days

	endelse

	endfor

days:	d_of_y = long(sec_rem/86400l) + 1
	strform, d_of_y, doy_str, 3		; form 3 char string

	sec_rem = sec_rem - (d_of_y-1) * 86400l

	hr = long(sec_rem/3600)			; calculate nearest hour
	strform, hr, hr_str, 2			; form 2 char string

	sec_rem = sec_rem - 3600*hr		; calculate remaining

	mn = long(sec_rem/60)			; calculate nearest minute
	strform, mn, mn_str, 2			; form 2 char string

	sec_rem = sec_rem - 60*mn		; calculate rem. sec

	sec_rem = sec_rem * 1000		; convert to microsec
	strform, long(sec_rem+0.5), sc_str, 5	; form 5 char string

	zulu_time(j) = yr_str + doy_str + hr_str + mn_str + sc_str
						; concatenate strings

	zulu_time(j) = strcompress(zulu_time(j))
						; delete white space

	endelse

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


