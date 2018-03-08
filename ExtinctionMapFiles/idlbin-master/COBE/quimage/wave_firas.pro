	pro wave_firas, header, field_name, signal_range, fsize

	signal_range = ''

	waves = sxpar(header, 'LINFRQ*')
	waves = waves(where(waves ne ''))
	n_waves = n_elements(waves)

	; handle there only being one wave value in the file
	if n_waves eq 1 then begin
	    signal_range = ''
	    return
	endif
        
        ; handle the case when choosing a non-standard field
        if n_waves ne fsize(1) then begin
          signal_range='1:'+strtrim(fsize(1),2) 
          goto, get_all
        endif
menu:
	signal_range = ''
	print, ' '
	print, 'The following FIRAS Line Map frequency values are in the file'
	print, ' '
	for i = 0, n_waves-1 do print, format='("Wave: ",i2,"  ",a)', $
            i+1,waves(i)
	print, ' '
	print, 'Enter a FIRAS line frequency value or a range of values'
	print, 'To specify a range use the format start:end, example 1:3'
	print, 'Use ALL to specify all.  Use BACK to go back a menu.
	read, '----> ', signal_range

	if strupcase(signal_range) eq 'ALL' then signal_range = '1:'+ $
               strtrim(n_waves,2)
	if strupcase(signal_range) eq 'BACK' then begin
	   signal_range = 'back'
	   return
	endif
	
	start_range = 1
	end_range = n_waves

;	check to see that the numbers input are within the valid range
	verify_range, signal_range, start_range, end_range
	if !err eq  10 then goto, menu

get_all:
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


