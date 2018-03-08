    pro wave_dirbe, header, field_name, signal_range, fsize

    if strupcase(field_name) eq 'STOKES  ' or $
        strupcase(field_name) eq 'STOKESSD' or $
        strupcase(field_name) eq 'STOKQUAL' or $
        strupcase(field_name) eq 'FRACUSED' or $
        strupcase(field_name) eq 'DISPLACE' then begin
menu1:
        signal_range = ''
 	print, ' '
	print, 'This field contains '+strtrim(fsize,2)+' elements.'
	print, ' '
	print, 'Enter an element or a range of elements'
	print, 'To specify a range use the format start:end, example 1:3'
	print, 'Use ALL to specify all.  Use BACK to go back a menu.
	read, '----> ', signal_range
	if strupcase(signal_range) eq 'ALL' then $
               signal_range='1:'+strtrim(fsize,2)
	if strupcase(signal_range) eq 'BACK' then begin
	   signal_range = 'back'
	   return
	endif
        start_range = 1
        end_range = fsize
    
;       check to see that the numbers input are within the valid range
        verify_range, signal_range, start_range, end_range
        if !err eq  10 then goto, menu1
    endif else begin
	waves = sxpar(header, 'WAVE*')
	waves = waves(where(waves ne ''))
	n_waves = n_elements(waves)

	; handle there only being one wave value in the file
	if n_waves eq 1 then begin
	    signal_range = ''
	    return
	endif

menu:
	signal_range = ''
 	print, ' '
	print, 'The following DIRBE wave values are in the file'
	print, ' '
	for i = 0, n_waves-1 do print, format='("Wave: ",i2,"  ",a)', $
            i+1,waves(i)
	print, ' '
	print, 'Enter a DIRBE wave value or a range of values'
	print, 'To specify a range use the format start:end, example 1:3'
	print, 'Use ALL to specify all.  Use BACK to go back a menu.
	read, '----> ', signal_range

	if strupcase(signal_range) eq 'ALL' then $
            signal_range = '1:'+strtrim(n_waves,2)
	if strupcase(signal_range) eq 'BACK' then begin
	   signal_range = 'back'
	   return
	endif
        start_range = 1
        end_range = n_waves
    
;       check to see that the numbers input are within the valid range
        verify_range, signal_range, start_range, end_range
        if !err eq  10 then goto, menu
    ENDELSE	

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


