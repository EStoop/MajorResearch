pro zul2adt, zulu_in, adt_time, istat
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    zulu2adt - Convert from zulu time to ADT time
;
;DESCRIPTION:                                   
;    IDL routine called by TIMECONV to convert from zulu time to
;    ADT time
;
;CALLING SEQUENCE:
;    zulu2adt, zulu_in, adt_time, istat
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    zulu_in     I     str [arr]   zulu date - "yydddhhmmsssss"
;    adt_time    O     dbl [arr]   seconds from 01/01/68 - 00hr
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
;    Add number of secs from year, day of year, hour, minute, and sec
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: leapchk, zuluqual, emul, addx
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

base_time = [154632192l,9266898l]	; # of intervals to 01-JAN-1985
num_sec_n = [2026291200l,73425l]	; # of intervals in non-leap year
num_sec_l = [-1557102592l,73626l]	; # of intervals in leap year
sec_in_dy = [711573504l,201l]		; # of intervals in day
sec_in_hr = [1640261632l,8l]		; # of intervals in hour
sec_in_mn = [600000000l,0l]		; # of intervals in minute

zulu_time = zulu_in				; 'insulate' input

n = n_elements(zulu_time)			; get input array size

zulu_time0 = STRCOMPRESS(zulu_time,/REMOVE_ALL)
;
z_str = '0000000000000000'

for i=0l,n-1 do begin
;
l_str = strlen(zulu_time0(i))			; get string length
if (strupcase(strmid(zulu_time0(i),l_str-1,1)) eq 'E') then begin
;
zulu_time0(i) = strmid(zulu_time0(i),0,l_str-1)	; trim off 'E'
zulu_time0(i) = zulu_time0(i) + z_str		; pad with zeros
zulu_time0(i) = strmid(zulu_time0(i),0,16)	; trim to 16 char
l_str = strlen(zulu_time0(i))			; get string length
;
endif
;
if (l_str lt 16) then begin
;
zulu_time0(i) = '19' + zulu_time0(i)		; assume 20th century
zulu_time0(i) = zulu_time0(i) + z_str		; pad with zeros
zulu_time0(i) = strmid(zulu_time0(i),0,16)	; trim to 16 char
;
endif
;
endfor
;
adt_time = lonarr(2*n)				; define output array
;		
;
for j = 0l,n-1 do begin

zuluqual, zulu_time0(j), qual_flag		; check qual of zulu time

if (qual_flag eq 0) or (istat(j) ne 0) then begin

	adt_time(2*j) = -1		; if bad qual or bad input
	adt_time(2*j+1) = -1
	istat(j) = -1			; then bail

endif else begin


	bin_time = base_time
			; initialize adt to base time 

	if (yr lt 1985) then begin

		istat(j) = -1
		message,'Binary time before 01-JAN-1985',/cont
			; ADT time before base time  
			; set return time to all zeros 
			; set return status to bad & exit 

	endif else begin

		for i=1985,yr-1 do begin

		leap = 0
		if (i / 4 eq i / 4.) then leap = -1
		if (i / 100 eq i / 100.) then leap = 0
		if (i / 400 eq i / 400.) then leap = -1

		if (leap eq -1) then $
			sec_in_yr = num_sec_l $
		else $
			sec_in_yr = num_sec_n
			; decide whether leap year or not
			; get corresponding secs in year 

		bin_time = addx(bin_time,sec_in_yr,2)
			; add year intervals 

		endfor  ; end year loop 

	endelse		; end year if 


	emul,sec_in_dy,dy-1,prod
		; multiply # of days times units per day 
	bin_time = addx(bin_time,prod,2)
		; add days 


	emul,sec_in_hr,hr,prod
		; multiply # of hours times units per hour 
	bin_time = addx(bin_time,prod,2)
		; add hours 


	emul,sec_in_mn,mn,prod
		; multiply # of minutes times units per minute 
	bin_time = addx(bin_time,prod,2)
		;add minutes 

	prod(0) = sc * 10000000l 
		; convert to 100 nanosecond units 
	prod(1) = 0
	bin_time = addx(bin_time,prod,2)
		; add seconds 


	adt_time(2*j) = bin_time(0)
	adt_time(2*j+1) = bin_time(1)

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


