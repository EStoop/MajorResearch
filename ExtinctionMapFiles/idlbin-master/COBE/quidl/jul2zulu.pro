pro jul2zulu, jul_time, zulu_time, istat
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    jul2zulu - Converts from reduced Julian Date to zulu time
;
;DESCRIPTION:                                   
;    IDL routine called by TIMECONV to convert from reduced Julian Date
;    to zulu time.  The reduced Julian Date can be converted to the full
;    Julian Date by adding 2400000.5 to it.
;
;CALLING SEQUENCE:
;    jul2zulu, jul_time, zulu_time, istat
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    jul_time    I     dbl [arr]   Reduced julian date
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
;    Convert to full Julian Date
;    Extract year, month, and day of month
;    Calculate day of year from month and day of month
;    Get fractional part of day and calculate hour, minute, and
;    second (up to millisec accuracy).
;    Convert to strings and concatanate.
;
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;
;    Uses a conversion algorithm from Fliegel, H.F. & van Flandern, 
;    T.C., Comm. ACM, Vol 11, pg. 657 (1968)
;
;    Subroutines called: leapchk, strform
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Dec 91
;-
;
n = n_elements(jul_time)			; get input array size
jd = long(jul_time + 2400001.d0)
;
zulu_time = strarr(n)				; define output array
;
day_n = [1,32,60,91,121,152,182,213,244,274,305,335]
day_l = [1,32,61,92,122,153,183,214,245,275,306,336]
			; day_n - day of year of 1st of each month (non-leap)
			; day_l - day of year of 1st of each month (leap)

;
	for j = 0l,n-1 do begin			; loop over all elements

	if (jul_time(j) lt 0.0d0) then begin

	istat(j) = -1
	zulu_time = '0000000000000000'
	MESSAGE,'Julian date must be non-negative',/INFO

	endif else begin

	l = jd(j) + 68569
	m = 4*l/146097
	l = l - (146097*m + 3)/4
	yr = 4000*(l+1)/1461001
	l = l - 1461*yr/4 + 31			; magic formuli !
	mon = 80*l/2447
	dy = l - 2447*mon/80
	l = mon/11
	mon = mon + 2 - 12*l
	yr = 100*(m-49) + yr + l

	if (leapchk(yr)) then day_ln = day_l else day_ln = day_n
			; decide whether leap year or not

	strform, yr, yr_str, 4			; form 4 char string

	d_of_y = day_ln(mon-1) + dy - 1		; get day of year from
						; month
	strform, d_of_y, doy_str, 3		; form 3 char string

	frac = double(jul_time(j) - long(jul_time(j)))
						; get fractional part of
						; julian time (in days)

	sc = long(frac * 86400000l + 0.5)	; convert from days to
						; milliseconds

	hr = sc/3600000l			; calculate nearest hour
	strform, hr, hr_str, 2			; form 2 char string

	sc = sc - 3600000l*hr			; calculate remaining

	mn = long(sc/60000l)			; calculate nearest minute
	strform, mn, mn_str, 2			; form 2 char string

	sc = sc - 60000l*mn			; calculate rem. millisec

	strform, sc, sc_str, 5			; form 5 char string

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


