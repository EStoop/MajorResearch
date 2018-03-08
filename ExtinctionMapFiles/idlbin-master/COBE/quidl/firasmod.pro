;+NAME/ONE-LINE DESCRIPTION:
;      FIRASMOD: used to read in the FIRAS Calibration Model IPs, 
;                FIRAS Line Profile PDSs, and other FIRAS PDS.
;
;DESCRIPTION: UIDL routine which reads in fields of the FIRAS IP/PDS
;             files. This is mainly used for non-skymap (i.e. not
;             pixilzed) fields.
;
;CALLING SEQUENCE:
;        firasmod, dsname, dsfield, data=data, units=units, message=message
;
;        dsname: the name of the data file to be read
;        dsfield: the name of the field to be read in (TTYPE field in
;                 the fits extension header.  Pass this parameter in as
;                 a question mark ('?') in order to get a list of
;                 fieldnames.
;        data     array for the data that was read in to be returned to
;        units    var to hold the units of the returned data
;        message  Y (or not set) for messages to be printed as operations
;                 are performed, N to be silent
;
;#
; SUBROUTINES CALLED: most of the FXB routines for handling the FITS
;                     I/O.
;
; COMMON BLOCKS:  None
;
; LIBRARY CALLS:  None
;
; WARNINGS:
;
; PROGRAMMING NOTES:
;
; MODIFICATION HISTORY: D. Ward   GSC   June 1993  Orginal program
; SPR 11759 17-May-1994   Ingest FIRAS PDS's. J. Newmark
;
;.TITLE
;Routine FIRASMOD
;-
	Pro Firasmod, dsname, dsfield, data=data, units=units, message=message

	if n_params() lt 1 then begin
	   print, 'Syntax: firasmod, dsname, dsfield, data=data, units=units'
	   goto, abort
	endif

	if n_elements(dsname) lt 1 then begin
	   print, 'FIRMOD: You must supply the name of the data file'
	   goto, abort
	endif

	if n_elements(message) lt 1 then message = 'Y'

	if message ne 'N' then Print, 'Reading FITS format file, '+dsname

        if n_elements(dsfield) eq 0 then dsfield = '?'
        dsfield = strupcase(dsfield)

	if dsfield ne '?' then $
	   if strlen(dsfield) lt 8 then $
		for i = 0,7-strlen(dsfield) do dsfield = dsfield+' '

	; get the primary FITS header
	header = headfits(dsname)

        instrume  = sxpar(header, 'INSTRUME') 
	if strtrim(instrume) ne 'FIRAS' then begin
	  print, 'FIRMOD: This is not a FIRAS IP/PDS File'
	  goto, abort
	endif

        ;------------------------------------------------------------
	; get the keywords which define the frequencies in the file
        ;------------------------------------------------------------
	nu_zero	 = sxpar(header, 'NU_ZERO')
	delta_nu  = sxpar(header, 'DELTA_NU')
	num_freq  = sxpar(header, 'NUM_FREQ')
        if num_freq eq 0 then num_freq = 1
	frequency = fltarr(num_freq)
	for i = 0, num_freq-1 do $
	    frequency(i) = nu_zero + (delta_nu * i)

	; read in the header for the fits BINARY TABLE extension
;	fxbopen,unit,dsname,1,bin_header
        all_data = readfits(dsname,bin_header,/ext)

	fields = fxpar(bin_header,'TTYPE*',count=fieldcnt)
	unitlist  = fxpar(bin_header,'TUNIT*',count=unitcnt)

	if strtrim(dsfield,2) eq '?' then begin
           print,'This file contains the following fields:'
	   for i = 0,fieldcnt-1 do print, fields(i)
           dsfield = ''
           read,dsfield,prompt='Enter Field Name as shown: '
 	   while strlen(dsfield) lt 8 do dsfield = dsfield + ' '
	endif

	choice = where(strupcase(dsfield) Eq strupcase(fields))
	units =  unitlist(choice(0))

        ;------------------------------------------------------------
	;IDL arrays start at 0, fits columns start at 1 so we have
	;to add 1 to the variable choice to point to the correct one
	;------------------------------------------------------------
	readfield = choice(0) + 1
	fieldpointer = choice(0)

	if message ne 'N' then print,'reading in data ',systime(0)
;	fxbread,unit,data,readfield
        data = tbget(bin_header,all_data,readfield)
;	fxbclose,unit
        all_data = 0
        dum = temporary(all_data)

        ;------------------------------------------------------------
	; here we re-scale the data arrays based on the actual number
	; of array positions which have data
        ;------------------------------------------------------------
	if fields(fieldpointer) eq 'SIGNAL  ' or $
        fields(fieldpointer) eq 'SERROR  ' or $
        fields(fieldpointer) eq 'SPECTRUM' or $
        fields(fieldpointer) eq 'SIGMAS  ' or $
        fields(fieldpointer) eq 'REAL_SPE' or $
        fields(fieldpointer) eq 'IMAG_SPE' or $
        fields(fieldpointer) eq 'REAL_VAR' or $
        fields(fieldpointer) eq 'IMAG_VAR' or $
        fields(fieldpointer) eq 'REAL_IMA' or $
	fields(fieldpointer) eq 'IM_SPEC ' then $
	   data = data(*,*,0:num_freq-1) $
        else begin
            sz=size(data)
            if (sz(1) eq 210) or (sz(1) eq 182) then data=data(0:num_freq-1)
        endelse

	Return

Abort:

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


