	pro verify_range, signal_range, start_range, end_range, $
	    delim=delim

	if n_elements(delim) eq 0 then delim = ':'

	!err = 0
;	handle no input given
	if signal_range eq ' ' then begin
	    !err = 10
	    return
	endif

;	if this is true then we have a range
	pos = strpos(signal_range,delim)

;	first handle if no beginning range value was given
	if pos eq 0 then begin
	   !err = 10
	   return
	endif

;	now handle if no ending range value was given
	if pos eq strlen(signal_range)-1 then begin
	   !err = 10
	   return
	endif

	if pos gt 0 then begin
	    num1 = fix(strmid(signal_range,0,pos))
	    if num1 ge start_range and num1 le end_range then begin
		len = strlen(signal_range) - (pos+1)
		num2 = fix(strmid(signal_range,pos+1,len))
		if num2 ge start_range and num2 le end_range then $
		   !err = 0 else !err = 10
	    endif else !err = 10

	    if !err eq 1 then if num1 ge num2 then !err = 10

;	else we have a scalar value
	endif else begin
	    num = fix(signal_range)
	    if num ge start_range and num le end_range then !err = 0 $
		else !err = 10
	endelse

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


