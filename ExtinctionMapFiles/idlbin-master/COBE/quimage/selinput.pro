;+
; NAME:
;  SELINPUT 
; PURPOSE:
;       This function gets all the input data set information from
;	the user. 
; CATEGORY:
;  User interface, Menu, Data I/O
; CALLING SEQUENCE:
;  SelInput ...
; INPUTS:
;   none
;
; OUTPUTS:
;   input_format	format, CSM, FITS, CISS
;   data_set_type	type, FIRAS, DIRBE, DMR, External
;   data_set_name	full file specification
;   fits_exten_type	the FITS extension type (if any)
;   field_name		field name
;   subscript		subscripts for field name
;   weights 		name of weighting array file
;   face_number  	face number, 7 for all
;
;#
; COMMON BLOCKS:
;  None.
; RESTRICTIONS:
;  None.
; MODIFICATION HISTORY:
;       Creation:  Pete Kryszak, GSC, March 1992.
;
;	16-Dec-1992	10330	Kryszak - added BACK and Quit functionality
;	23-Dec-1992	10387	Kryszak - added help (this module had to
;					  change)
;       23-Feb-1993	10819	Ward    - added extension variable to
;					  dsname, to support FITS extensions
;       15-Apr-1993     10829   Ward,Newmark - properly handle FITS fields.
;  SPR 11759 17-May-1994   Ingest FIRAS PDS's. J. Newmark
; SPR 11905 07-Sep-1994   Ingest DIRBE PDS's. J. Newmark
;  SPR 11981 08-Nov-1994  Read in skymap_info files.  J. Newmark 
; SPR 12140 16-Mar-1995   Ingest FIRAS line maps. J. Newmark
;
;.TITLE
;Routine SELINPUT
;-
;
	Pro SelInput,	$ 
		input_format,		$ 
		data_set_type,		$
		data_set_name,		$ 
                fits_exten_type,        $
		field_name,		$
		subscript,		$
		weights,		$
		face_number 

	;
	; Initialize outputs
	;
	!err = 1
	!error = 1
        old_quiet=!quiet

	input_format = ''
	data_set_type = ''
	data_set_name = ''
	fits_exten_type = ''
	field_name = ''
	subscript = ''
	weights = ''
	face_number = 7

inpform:
	; Get data set format
	inpform, input_format 
	IF !err LT 0 THEN Begin
            !quiet=0
	    message, 'Input format invalid.', /informational 
            !quiet=old_quiet
	    GOTO, abort ; return if errors
	EndIf
	IF !err EQ 99 THEN Goto, abort
	IF !err EQ 88 THEN Goto, abort


; get the dataset type (DMR, DIRBE, FIRAS)
intype:
	if (input_format eq 'CSM') then begin
	   intype, data_set_type
	   IF !err LT 0 THEN Begin
              !quiet=0
	      message, 'Data set type invalid.', /informational 
              !quiet=old_quiet
	      GOTO, abort ; return if errors
	   EndIf
	   IF !err EQ 99 THEN Goto, abort
	   IF !err EQ 88 THEN Goto, inpform
	EndIf

;
; user format data
;
If data_set_type Eq 'U' Then Begin
	; use Restore, etc, then UPutData in DConvert
	If input_format Eq 'CISS' Then goto, done
	If input_format Eq 'CSM' Then Goto, done
EndIf
	
dsname:
	; Get data set name
	dsname, input_format, data_set_type, data_set_name, fits_exten_type
	IF !err LT 0 THEN Begin
            !quiet=0
	    message, 'Data set name invalid.', /informational 
            !quiet=old_quiet
	    GOTO, abort ; return if errors
	EndIf

	IF !err EQ 99 THEN Goto, abort

	; handle back...
	IF !err EQ 88 THEN Begin
	    If (input_format Eq 'CISS') or (input_format Eq 'FITS') Then Begin
		GOTO, inpform
		EndIf $
	    Else Begin
		GOTO, intype
	    EndElse
	EndIf

	; if save set then skip all this...
	If input_format Eq 'CISS' Then Goto, done

If input_format Eq 'FITS' Then Begin
	subscript = ''
	; get the fits header
	header = headfits(data_set_name)

	data_set_type = sxpar(header, 'INSTRUME')
	product_type = sxpar(header, 'PRODUCT')
        if string(product_type) eq 'PIXINFO6' or $
           string(product_type) eq 'PIXINFO9' then data_set_type='PXINFO'
        
        frequency_type = sxpar(header, 'CHANNEL')
	weights = sxpar(header, 'WEIGHTS')
	if !err lt 0 then weights = ''

	if weights ne '' then begin
	   print, ' '
	   print, 'The file that you have selected has a weighting array'
	   question = 'Would you like to read in the weighting array? '
	   wans = one_line_menu([question,'yes','no'],init=1)
	   if wans eq 1 then weights = 'yes'
	endif

	; if this is an image written out by UImage earlier, we don't allow
	; subsetting, we just bring in the entire thing.
	index = where(strpos(header,'Written by UImage DATA IO') gt 0)
	if index(0) ne -1 then goto, leaveok
        index = sxpar(header,'EXTEND')
        if (index ne 1) then goto, leaveok

        if strtrim(data_set_type) eq 'FIRAS' then begin
;
; check to see if one of the FIRAS PDS which contain no fields that
; are in quad-cube format. If one of these file types then abort.
          if strtrim(product_type) eq 'HLINE-PR' or $
               strtrim(product_type) eq 'LLINE-PR' or $
               strtrim(product_type) eq 'CAL-ERR' or $
               strtrim(product_type) eq 'CCAL-ERR' or $
               strtrim(product_type) eq 'COV-MAT' or $
               strtrim(product_type) eq 'C-VEC' or $
               strtrim(product_type) eq 'A-VEC' or $
               strtrim(product_type) eq 'PHYS-ST' or $
               strtrim(product_type) eq 'ORTHO-ST' or $
               strtrim(product_type) eq 'MOD-ERRS' or $
               strtrim(product_type) eq 'DS-QMAT' or $
               strtrim(product_type) eq 'DS-CSPEC' or $
               strtrim(product_type) eq 'DS-SCMAT' or $
               strtrim(product_type) eq 'CAL-MOD' then begin
                  !quiet=0
         	  message, 'Selected data set contains no ' + $
                    'Quad-cube fields.', /informational 
                  !quiet=old_quiet
                  goto, abort
           endif
;
           if string(product_type) eq 'HLIN-MAP' then $
             data_set_type='HLIN-MAP'
           if string(product_type) eq 'LLIN-MAP' then $
             data_set_type='LLIN-MAP'
;
	   num_freq  = sxpar(header, 'NUM_FREQ')
	   nu_zero	 = sxpar(header, 'NU_ZERO')
	   delta_nu  = sxpar(header, 'DELTA_NU')
           if num_freq eq 0 then num_freq=1
	   frequency = fltarr(num_freq)
	   for i = 0, num_freq-1 do $
	       frequency(i) = nu_zero + (delta_nu * i)
	endif

	; read in the header for the fits BINARY TABLE extension
	fxbopen,unit,data_set_name,1,bin_header
	fxbclose, unit

    fitsfldname:
	fieldlist,bin_header,data_set_type,fields,dims,fieldcnt

	choice = fitsflds(bin_header, fields, fieldcnt)
	IF !err LT 0 THEN Begin
            !quiet=0
	    message, 'Field name invalid.', /informational 
            !quiet=old_quiet
	    GOTO, abort ; return if errors
	EndIf
	IF !err EQ 99 THEN Goto, abort
	IF !err EQ 88 THEN Goto, dsname

	field_name = fields(choice-1)

   subscripts:
	if strtrim(data_set_type) eq 'DIRBE' then begin
	   fsize=fix(strmid(dims(choice-1),0,strlen(dims(choice-1))-1))
	   if fsize gt 1 then begin
               wave_dirbe,header,field_name,subscript,fsize
	       if strupcase(subscript) eq 'BACK' then goto, fitsfldname
           endif
	endif

	if strtrim(data_set_type) eq 'FIRAS' then begin
	   fsize=intarr(2)
	   fsize(0) = 1
	   fsize(1)=fix(strmid(dims(choice-1),0,strlen(dims(choice-1))-1))
	   if fsize(1) eq 210 or fsize(1) eq 182 then fsize(1) = num_freq
	   if fsize(1) gt 1 and fsize(1) eq num_freq then begin
	      firas_sub, field_name,frequency,fsize,subscript
	      if strupcase(subscript) eq 'BACK' then goto, fitsfldname
	   endif else begin
; handle the case for non-standard (non-frequency) fields: get all.
              if fsize(1) gt 1 then subscript='1:'+strtrim(fsize(1),2) 
           endelse
	endif
	if data_set_type eq 'HLIN-MAP' or data_set_type eq 'LLIN-MAP' then begin
	   fsize=intarr(2)
	   fsize(0) = 1
	   fsize(1) = $
	      fix(strmid(dims(choice-1),1,strpos(dims(choice-1),')')-1))
	   if fsize(1) gt 1 then begin
             wave_firas,header,field_name,subscript,fsize
	     if strupcase(subscript) eq 'BACK' then goto, fitsfldname
           endif
	endif

	face, face_number
	IF !err LT 0 THEN Begin
            !quiet=0
	    message, 'Face number invalid.', /informational 
            !quiet=old_quiet
	    GOTO, abort ; return if errors
	EndIf
	IF !err EQ 99 THEN Goto, abort
	IF !err EQ 88 THEN Goto, fitsfldname
Endif

If input_format Eq 'CSM' Then Begin
    csmfldname:
	fldname, data_set_type, data_set_name, field_name, subscript
	IF !err LT 0 THEN Begin
            !quiet=0
	    message, 'Field name or subscript invalid.', /informational 
            !quiet=old_quiet
	    GOTO, abort ; return if errors
	EndIf
	IF !err EQ 99 THEN Goto, abort
	IF !err EQ 88 THEN Goto, dsname

	face, face_number
	IF !err LT 0 THEN Begin
            !quiet=0
	    message, 'Face number invalid.', /informational 
            !quiet=old_quiet
	    GOTO, abort ; return if errors
	EndIf
	IF !err EQ 99 THEN Goto, abort
	IF !err EQ 88 THEN Goto, csmfldname
Endif

leaveok:
	!err = 1
	!error = 1

done:
abort:  
	return

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


