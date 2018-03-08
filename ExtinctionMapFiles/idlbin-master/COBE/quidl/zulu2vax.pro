pro zulu2vax, zulu_in, vax_time, istat
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    zulu2vax - Convert zulu time to VAX time
;
;DESCRIPTION:                                   
;    IDL routine called by TIMECONV to convert from zulu time to
;    VAX time 'dd-mmm-yy:hh:mm:ss.sss'.
;
;CALLING SEQUENCE:
;    zulu2vax, zulu_in, vax_time, istat
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    zulu_in     I     str [arr]   zulu date - "yydddhhmmsssss"
;    vax_time    O     str [arr]   VAX time string

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
;    Extract year, day of year, hour-minute-second string.
;    Calculate month and day of month form day of year
;    Form strings from year and month and day of month
;    Concatanate with hour-minute-second string to get VAX format
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: zuluqual, leapchk, strform
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

zulu_time = zulu_in				; 'insulate input'

n = n_elements(zulu_time)			; get input array size
;
zulu_time = STRCOMPRESS(zulu_time,/REMOVE_ALL)

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
hr_mn_sc = strmid(zulu_time,7,9)		; get hour-min-sec string
;
vax_time = strarr(n)				; define output array
;
;
day_n = [1,32,60,91,121,152,182,213,244,274,305,335,365]
day_l = [1,32,61,92,122,153,183,214,245,275,306,336,366]
			; day_n - day of year of 1st of each month (non-leap)
			; day_l - day of year of 1st of each month (leap)
;
month = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT', $
	 'NOV','DEC']
;
for j = 0l,n-1 do begin

zuluqual, zulu_time(j), qual_flag

if (qual_flag eq 0) or (istat(j) ne 0) then begin

	vax_time(j) = '00-XXX-0000:00:00:00.000'
	istat(j) = -1

endif else begin

if (leapchk(yr)) then day_ln = day_l else day_ln = day_n
			; decide whether leap year or not


	done = 0			; done flag
	i = 11				; initialize month counter
	while (done eq 0) do begin

	if (day_ln(i) le dy) then begin	; find month starting 
						; before (or at)
						; day of year
						; (backwards loop)

		mnth_str = month(i)
		dy = dy - day_ln(i) + 1	; subtract offset
		done = -1			; set done flag

	endif

	i = i-1		; decrement month counter

	endwhile


strform, yr, yr_str, 4

strform, dy, dy_str, 2

hr_mn_sc(j) = strmid(hr_mn_sc(j),0,2) + ':' + strmid(hr_mn_sc(j),2,2) + ':' + $
	      strmid(hr_mn_sc(j),4,2) + '.' + strmid(hr_mn_sc(j),6,3)

vax_time(j) = dy_str + '-' + mnth_str + '-' + yr_str + ':' + hr_mn_sc(j)

vax_time(j) = strcompress(vax_time(j))

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


