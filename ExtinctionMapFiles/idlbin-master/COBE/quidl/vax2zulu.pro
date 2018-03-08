pro vax2zulu, vax_time, zulu_time, istat
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    vax2zulu - Convert VAX time to zulu time
;
;DESCRIPTION:                                   
;    IDL routine called by TIMECONV to convert from VAX time
;    'dd-mmm-yy:hh:mm:ss.sss' to zulu time.
;
;CALLING SEQUENCE:
;    vax2zulu, vax_time, zulu_time, istat
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    vax_time    I     str [arr]   VAX time string
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
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Extract year, month, and day from VAX string.
;    Convert month and day to day of year.
;    Extract hour, minute, and second from VAX string.
;    Convert to strings and concatanate.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: leapchk, strform
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Dec 91
;-
;
n = n_elements(vax_time)			; get input array size
;
zulu_time = strarr(n)				; define output array
;
day_n = [1,32,60,91,121,152,182,213,244,274,305,335,366]
day_l = [1,32,61,92,122,153,183,214,245,275,306,336,367]
			; day_n - day of year of 1st of each month (non-leap)
			; day_l - day of year of 1st of each month (leap)

month = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT', $
	 'NOV','DEC']

;
	for j = 0l,n-1 do begin			; loop over all elements

	vax_time(j) = strcompress(strupcase(vax_time(j)))
	first_dash = strpos(vax_time(j),'-')
	sec_dash = strpos(vax_time(j),'-',first_dash+1)

	if ((first_dash eq -1) or (sec_dash eq -1) or $
	    (sec_dash-first_dash ne 4)) then begin
		str = 'Improper Vax time: ' + STRING(vax_time(j))
		MESSAGE,STRCOMPRESS(str),/INFO
		zulu_time(j) = '0000000000000000'
		istat(j) = -1
		goto, nxt_el
	endif
				; check for two dashes 3 spaces apart


;		 ----------------
;		| Determine year |
;		 ----------------

	yr = fix(strmid(vax_time(j),sec_dash+1,4))	; get year

	if (leapchk(yr)) then day_ln = day_l else day_ln = day_n
			; decide whether leap year or not

	strform, yr, yr_str, 4			; form 4 char string



;		 -----------------
;		| Determine month |
;		 -----------------

	mn_str = strmid(vax_time(j),first_dash+1,3)
	for i=0,11 do begin

	if (mn_str eq month(i)) then begin
		mon = i
		max_day = day_ln(i+1) - day_ln(i)
		goto, fndmth
	endif

	endfor

	str = 'Improper Vax time: ' + STRING(vax_time(j))
	MESSAGE,STRCOMPRESS(str),/INFO
	zulu_time(j) = '0000000000000000'
	istat(j) = -1
	goto, nxt_el


;		 ---------------
;		| Determine day |
;		 ---------------

fndmth:	dy = fix(strmid(vax_time(j),0,first_dash))	; get day

	if ((dy lt 0) or (dy gt max_day)) then begin
		str = 'Improper Vax time: ' + STRING(vax_time(j))
		MESSAGE,STRCOMPRESS(str),/INFO
		zulu_time(j) = '0000000000000000'
		istat(j) = -1
		goto, nxt_el
	endif

	d_of_y = day_ln(mon) + dy - 1		; get day of year from
						; month
	strform, d_of_y, doy_str, 3		; form 3 char string



;		 ---------------------------
;		| Determine hour:minute:sec |
;		 ---------------------------

	hr = fix(strmid(vax_time(j),12,2))
	mn = fix(strmid(vax_time(j),15,2))
	sc = float(strmid(vax_time(j),18,6)) * 1000

	if (hr gt 23) then begin
		str = 'Hour greater than 23: ' + STRING(vax_time(j))
		MESSAGE,STRCOMPREESS(str),/INFO
		zulu_time(j) = '0000000000000000'
		istat(j) = -1
		goto, nxt_el
	endif

	if (mn gt 59) then begin
		str = 'Minute greater than 59: ' + STRING(vax_time(j))
		MESSAGE,STRCOMPRESS(str),/INFO
		zulu_time(j) = '0000000000000000'
		istat(j) = -1
		goto, nxt_el
	endif

	strform, hr, hr_str, 2
	strform, mn, mn_str, 2

	strform, long(sc+0.5), sc_str, 5	; form 5 char string

	if (fix(strmid(sc_str,0,2) gt 59)) then begin
		str = 'Second greater than 59: ' + STRING(vax_time(j))
		MESSAGE,STRCOMPRESS(str),/INFO
		zulu_time(j) = '0000000000000000'
		istat(j) = -1
		goto, nxt_el
	endif


	zulu_time(j) = yr_str + doy_str + hr_str + mn_str + sc_str
						; concatenate strings

nxt_el:	zulu_time(j) = strcompress(zulu_time(j))
						; delete white space


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


