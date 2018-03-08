pro adt2zulu, bin_time, zulu_time, istat
;+                                                                  
;  NAME:
;    adt2zulu
;
;  PURPOSE:                                   
;    Convert from VAX system time (100 nanosecond units since 
;    17-NOV-1858) to calendar date
;
;  CALLING SEQUENCE:
;    adt2zulu, bin_time, zulu_time
;
;  INPUT:
;    VAX system time - (2n) dimensional i*4 array
;
;  OUTPUT:
;    Date - The format being 'yydddhhmmsssss'.
;
;  SUBROUTINES CALLED:
;    subx
;    ediv
;    leapchk
;    strform
;
;  REVISION HISTORY
;    J.M. Gales
;    Jan 92
;
;    SPR 10414  Change loop index (j) from short to long integer
;    01/05/93   J.M. Gales
;-
;
n = n_elements(bin_time)/2			; get input array size
;
zulu_time = strarr(n)				; define output array
;
base_time = [154632192l,9266898l]	; # of intervals to 01-JAN-1985
num_sec_n = [2026291200l,73425l]	; # of intervals in non-leap year
num_sec_l = [-1557102592l,73626l]	; # of intervals in leap year
sec_in_dy = [711573504l,201l]		; # of intervals in day
sec_in_hr = [1640261632l,8l]		; # of intervals in hour
sec_in_mn = [600000000l,0l]		; # of intervals in minute
;
;
;
	for j = 0l,n-1 do begin			; loop over all elements

	sec_rem = subx([bin_time(2*j),bin_time(2*j+1)],base_time,2,borrow)
						; subtract base time

	if (borrow eq 1) then begin

		zulu_time(j) = '0000000000000000'
		istat(j) = -1
		MESSAGE,'Binary time before 01-JAN-1985',/INFO
						; before base time

	endif else begin

	for yr=1985,9999 do begin		; loop over years

	if (leapchk(yr)) then num_sec_ln = num_sec_l else $
				 num_sec_ln = num_sec_n
			; decide whether leap year or not

	sec_rem_0 = subx(sec_rem,num_sec_ln,2,borrow)
						; subtract year intervals

	if (borrow eq 0) then begin		; difference > 0

		sec_rem = sec_rem_0		; new remainder

	endif else begin

		strform, yr, yr_str, 4		; form 4 char string
		goto, days			; exit to day determination

	endelse

	endfor


days:	ediv,sec_rem,sec_in_dy,quot,sec_rem	; divide remaining intervals
						; by intervals per day

	d_of_y = quot(0) + 1			; day = quotient plus 1
	strform, d_of_y, doy_str, 3		; form 3 char string

	ediv,sec_rem,sec_in_hr,quot,sec_rem	; divide remaining intervals
						; by intervals per hour

	hr = quot(0)				; hour = quotient
	strform, hr, hr_str, 2			; form 2 char string

	ediv,sec_rem,sec_in_mn,quot,sec_rem	; divide remaining intervals
						; by intervals per minute

	mn = quot(0)				; minute = quotient
	strform, mn, mn_str, 2			; form 2 char string

	sec = sec_rem(0) / 10000l		; convert to microsec
	strform, long(sec+0.5), sc_str, 5	; form 5 char string

	zulu_time(j) = yr_str + doy_str + hr_str + mn_str + sc_str
						; concatenate strings

	zulu_time(j) = strcompress(zulu_time(j))
						; delete white space

	endelse

	endfor

exit:

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


