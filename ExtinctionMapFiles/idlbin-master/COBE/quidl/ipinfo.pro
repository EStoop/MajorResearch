;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    IPINFO reads the PDS/IP FITS header and displays field and frequency info
;
;DESCRIPTION:
;    IPINFO reads in the FITS headers and displays some useful information
;    on what's in the file.  This can be useful when working at the UIDL
;    level with the PDSs/IPs
;
;CALLING SEQUENCE:
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    param		I/O   	type		description
;    dsname              I      string		filename to read
;    freq_only		 I 	switch		if /freq_only set then only
;                                               the frequency data will be
;                                               displayed
;
;WARNINGS:
;
;EXAMPLE: 
;
;#
;COMMON BLOCKS:
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;
;MODIFICATION HISTORY
;     19-Aug-1993   Dalroy Ward, GSC,   initial routine
; SPR 11905 07-Sep-1994   Ingest DIRBE PDS's. J. Newmark
;
;
;.TITLE
;IPINFO
;-

	pro ipinfo, dsname, freq_only=freq_only

	files = findfile(dsname,count=fcount)
	if fcount eq 0 then begin 
	   print, 'IPINFO: file not found '+dsname
	   return
	endif

	result = is_fits(dsname)
	if result eq 0 then begin
	   print, 'IPINFO: input file must be a FITS format file'
	   return
	endif

	; get the fits header
	header = headfits(dsname)

	telescop = sxpar(header, 'TELESCOP')
	if strtrim(telescop) ne 'COBE' then begin
	   print, 'IPINFO is designed to read the COBE PDSs/IPs only'
	   return
	endif

	instrume = sxpar(header, 'INSTRUME')

	if strtrim(instrume) eq 'FIRAS' then begin
	   num_freq  = sxpar(header, 'NUM_FREQ')
	   nu_zero	 = sxpar(header, 'NU_ZERO')
	   delta_nu  = sxpar(header, 'DELTA_NU')
	   freq = fltarr(num_freq)
	   for i = 0, num_freq-1 do $
	       freq(i) = nu_zero + (delta_nu * i)
	endif

	if strtrim(instrume) eq 'DIRBE' then begin
	   waves = sxpar(header,'WAVE*')
	endif

	if keyword_set(freq_only) then goto, printfreq

	; read in the header for the fits BINARY TABLE extension
	fxbopen,unit,dsname,1,bin_header
	fxbclose, unit

	fieldlist,bin_header,instrume,fields,dims,fieldcnt

	print, ' '
	print, 'Information on file: '+dsname
	print, 'The file contains data from the '+instrume+' instrument'
	print, ' '
	print, 'And contains datafields:'
	print, fields
	print, ' '

printfreq:
	if strtrim(instrume) eq 'FIRAS' then begin
	   print, 'The file contains the frequencies (in GHz):'

	   askcr = 15
	   if keyword_set(freq_only) then askcr = 22
	   str = ''
	   indices = indgen(num_freq)
	   last = fix((num_freq-1) / 5)
	   rest = num_freq - (last*5)
	   for i = 0, last-1 do begin
	       j = (i) * 5
	       if i eq askcr then $
		  read, 'Press carriage return to continue: ', str
	       print, indices(j), freq(j), indices(j+1), freq(j+1),$
	              indices(j+2), freq(j+2), indices(j+3), freq(j+3),$
	  	      indices(j+4), freq(j+4), FORMAT='(5(i3,2x,f8.3,3x))'
	   endfor

	   for i = 0, rest-1 do begin
	       j = 5 * last + i
	       print, indices(j),freq(j), FORMAT='(i3,2x,f8.3,3x,$)'
	   endfor
	endif

	if strtrim(instrume) eq 'DIRBE' then begin
	   print, 'The file contains the wavelengths (in microns):'
	   print, waves
	endif
	     
	End
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


