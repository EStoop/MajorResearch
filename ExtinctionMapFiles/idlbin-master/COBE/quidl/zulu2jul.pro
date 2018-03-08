pro zulu2jul, zulu_in, jul_time, istat
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    zulu2jul - Converts from zulu time to reduced Julian Date
;
;DESCRIPTION:                                   
;    IDL routine called by TIMECONV to convert from zulu time
;    to the reduced Julian Date.  The reduced Julian Date can be
;    converted to the full Julian Date by adding 2400000.5 to it.
;
;CALLING SEQUENCE:
;    zulu2jul, zulu_in, jul_time, istat
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    zulu_in     I     str [arr]   zulu date - "yydddhhmmsssss"
;    jul_time    O     dbl [arr]   Reduced julian date
;    istat       I/O   int [arr]   status array
;
;WARNINGS:
;    This routine should not be called directly.
;
;EXAMPLE: 
;    Not a user-callable routine
;#
;COMMON BLOCKS:
;    zulu
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Convert to "extended" Zulu format if necessary.
;    Extract year, day of year, hour, minute, and second.
;    Calculate month and day of month from day of year
;    Get full Julian date (integer part).
;    Convert hour, minute, second to portion on day.
;    Add this to Julian date and subtract offset to get reduced Julian
;
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;
;    Uses a conversion algorithm from Fliegel, H.F. & van Flandern, 
;    T.C., Comm. ACM, Vol 11, pg. 657 (1968)
;
;    Subroutines called: leapchk, zuluqual
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Dec 91
;
;    SPR 10414  Change loop index (j) from short to long integer
;    01/05/93   J.M. Gales
;
; SPR 9750  - Blanks in zulu times can lead to conversion errors.
;             Solution: Compress zulu input.
;             Fixed: 5 Jun 1992 JMG
;
;-
;
COMMON zulu, yr,dy,hr,mn,sc

zulu_time = zulu_in				; 'insulate' input

n = n_elements(zulu_time)			; get input array size

zulu_time = STRCOMPRESS(zulu_time,/REMOVE_ALL)
;
z_str = '0000000000000000'

for i=0l,n-1 do begin
;
l_str = strlen(zulu_time(i))			; get string length
if (strupcase(strmid(zulu_time(i),l_str-1,1)) eq 'E') then begin
;
zulu_time(i) = strmid(zulu_time(i),0,l_str-1)	; trim off 'E'
zulu_time(i) = zulu_time(i) + z_str		; pad with zeros
zulu_time(i) = strmid(zulu_time(i),0,16)	; trim to 16 char
l_str = strlen(zulu_time(i))			; get string length
;
endif
;
if (l_str lt 16) then begin
;
zulu_time(i) = '19' + zulu_time(i)		; assume 20th century
zulu_time(i) = zulu_time(i) + z_str		; pad with zeros
zulu_time(i) = strmid(zulu_time(i),0,16)	; trim to 16 char
;
endif
;
endfor
;
;
jul_time = dblarr(n)				; define output array
;						  (double precison)
;
;
day_n = [1,32,60,91,121,152,182,213,244,274,305,335]
day_l = [1,32,61,92,122,153,183,214,245,275,306,336]
			; day_n - day of year of 1st of each month (non-leap)
			; day_l - day of year of 1st of each month (leap)

rdc_off = double(2400001l)
;
for j = 0l,n-1 do begin

zuluqual,zulu_time(j), qual_flag

if (qual_flag eq 0) or (istat(j) ne 0) then begin

	jul_time(j) = -1
	istat(j) = -1

endif else begin

if leapchk(yr) then day_ln = day_l else day_ln = day_n
			; decide whether leap year or not

	done = 0			; done flag
	i = 11				; initialize month counter
	while (done eq 0) do begin

	if (day_ln(i) le dy) then begin	; find month starting 
						; before (or at)
						; day of year
						; (backwards loop)

		month = i + 1
		dy = dy - day_ln(i) + 1	; subtract offset
		done = -1			; set done flag

	endif

	i = i-1		; decrement month counter

	endwhile


mm = long((month-14)/12)
yr_long = long(yr)		; long integer yr
dy_long = long(dy)		; long integer day


jt = dy_long - 32075 + 1461*(yr_long + 4800 + mm)/4  $
                     +  367*(month - 2 - mm*12)/12  $
                     -    3*((yr_long + 4900 + mm)/100)/4

				; magic formula !


hr_day = double(hr)/24		; portion of day for hour
mn_day = double(mn)/1440		; portion of day for minute
sc_day = double(sc)/86400		; portion of day for second

jul_time(j) = double(jt) + hr_day + mn_day + sc_day - rdc_off

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


