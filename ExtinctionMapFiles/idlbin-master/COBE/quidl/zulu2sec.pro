pro zulu2sec, zulu_in, sec_time, istat
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    zulu2sec - Convert from zulu time to # of seconds from fiduicial date
;
;DESCRIPTION:                                   
;    IDL routine called by TIMECONV to convert from zulu time to
;    number of seconds from 01-JAN-1968:00:00:00.000.
;
;CALLING SEQUENCE:
;    zulu2sec, zulu_in, sec_time, istat
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    zulu_in     I     str [arr]   zulu date - "yydddhhmmsssss"
;    sec_time    O     dbl [arr]   seconds from 01/01/68 - 00hr
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
sec_time = dblarr(n)				; define output array
;						  (double precison)
;
;
num_sec_n = 3600l * 24 * 365			; num sec in non-leap year
num_sec_l = 3600l * 24 * 366			; num sec in leap year
;
;
for j = 0l,n-1 do begin

zuluqual, zulu_time(j), qual_flag		; check qual of zulu time

if (qual_flag eq 0) or (istat(j) ne 0) then begin

	sec_time(j) = -1		; if bad qual or bad input
	istat(j) = -1			; then bail

endif else begin

if (yr lt 1968 or yr gt 2035) then begin

str = 'Improper year - must be between 1968 and 2035 (inclusive): ' + $
      STRING(yr)
MESSAGE,STRCOMPRESS(str),/INFO
istat(j) = -1

endif

	tot_sec = 0l			; define tot_sec as long int
	for i=1968,yr-1 do begin	; loop through years from 1968

	if (leapchk(i)) then $
	tot_sec = tot_sec + num_sec_l  $
	else $
	tot_sec = tot_sec + num_sec_n

	endfor				; end year loop

	tot_sec = tot_sec + (dy-1) * 86400l		; add day secs
	tot_sec = tot_sec + hr *      3600l		; add hour secs
	tot_sec = tot_sec + mn *        60l		; add min secs

	tot_sec = tot_sec + fix(sc)
	frac = sc - fix(sc)

	sec_time(j) = double(tot_sec) + frac

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


