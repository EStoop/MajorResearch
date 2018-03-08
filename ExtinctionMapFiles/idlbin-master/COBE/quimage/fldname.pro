;+
; NAME:
;	FLDNAME 
; PURPOSE:
;       This routine allows the user to specify which field is desired
;	from the specified data set name.  The array subscript is also
;	selected in this routine.
;	
; CATEGORY:
;	User interface, Menu, DATA I/O
;
; CALLING SEQUENCE:
;	fldname, dstype, dsname, fldname, subscripts
;	
; INPUTS:
;	data set type
;	data set name
;       
; OUTPUTS:
;	The field name in text string format and any subscripts.
;
;       !Err = -1 if this function is called incorrectly
;	     = 88 is returned if "Return to previous menu" is chosen.
;
; COMMON BLOCKS:
;	None.
; RESTRICTIONS:
;	None.
; MODIFICATION HISTORY:
;       Creation:  Pete Kryszak, GSC, March 1992.
;
;	16-Dec-1992	10330	Kryszak - added/edit menu options
;	18-Dec-1992	10311	Kryszak - 'ALL' to specify all channels
;	23-Dec-1992	10387	Kryszak - added help
;	05-Jan-1993	10177	Kryszak - subscript check
;       11-Mar-1993     10672   Ewing - change backward control options.
;
;-
;
Pro fldname, dstype, dsname, name, subscripts
;
;
old_quiet=!quiet
; get the full list for 'Other' choice
; and subscript checking 
tfield = getflds( dsname )

menu:
; Get available field names - This gets only STANDARD FIELD NAMES
; this should be removed if ADBs and STANDARD FLDs go away (soon I hope)
status =  lstflds( dstype, dsname, flist, fcount )

IF ( status LE 0 ) THEN BEGIN
        !quiet=0
	Message,/informational, 'Valid field list not retrieved.'
	!err = -1
	!error = -1
        !quiet=old_quiet
	goto, abort
	ENDIF

field_list = strarr( fcount + 5)
field_list( 0 ) = 'Field List'
field_list( fcount+1 ) = 'Other'
field_list( fcount+2 ) = '    '
field_list( fcount+3 ) = 'HELP'
field_list( fcount+4 ) = 'Return to previous menu'

for i = 0, fcount-1 do begin
	field_list(i+1) = flist(i)
	endfor 

if (fcount Eq 0 ) then begin
	message, /continue, 'No fields found.'
	!ERR = -1
	!ERROR = -1
	goto, abort
	EndIf



jmenu:		;just the menu, see also menu (above)
status = umenu( field_list, title = 0, init=1, xpos=0)
if ( status LE 0 ) then Begin
	!ERR = -1
        !ERROR = -1
	goto, abort ; return if errors
	EndIf

if ( field_list(status) Eq 'HELP' ) Then Begin
	action, 'fldname.hlp'
	goto, jmenu
	EndIf

if ( field_list(status) Eq 'Return to previous menu' ) Then Begin
	!ERR = 88
	!ERROR = 88
	goto, abort 
	EndIf

other = 0

other:
if (other eq 1) or (STRUPCASE(field_list(status)) eq 'OTHER') Then Begin
   
	fcount = n_elements( tfield )
	field_list = strarr( fcount + 4)
	field_list( 0 ) = 'Select Field'   ;first
	field_list( 1:fcount) = tfield
	field_list( fcount + 1 ) = '    '
	field_list( fcount + 2 ) = 'HELP'
	field_list( fcount + 3 ) = 'Return to previous menu'

	other_menu:
	status = umenu( field_list, title = 0, init=1, xpos=0)

	if ( status LE 0 ) then Begin                   
                !ERR = -1                               
                !ERROR = -1                             
                goto, abort ; return if errors          
                EndIf                                   
                                                        
        if ( field_list(status) Eq 'HELP' ) Then Begin  
                action, 'fldothr.hlp' 
                goto, other_menu                              
                EndIf                                   
                                                        
	if ( field_list(status) Eq 'Return to previous menu' ) Then Begin
		goto, menu
		EndIf

	EndIf

; set data field name according to selection by user
spc = STRPOS( field_list(status), ' ' )
fldname = STRMID(field_list( status ),0, spc )
;print, fldname

;
; get field array subscripts

; this will emphasize the difference between ADB formats 
; with STD FIELDS and other formats which don't have these fields
validfld = where( strpos(tfield,fldname) ne -1 )
if (validfld(0) eq -1) then begin
	message, /continue, fldname + ' is not in this data set.'
	message, /continue, 'Choose "Other" to see the valid list.'
	other = 1
	goto, other
	endif

fldname = tfield( validfld(0) )
getfinf, fldname, name, fsize


If fsize(0) gt 0 then begin
   subscr:        
   print, 'You have selected '+ name
   print, 'Now select a subscript within the array.'
   print, 'The array is a '+strtrim(string(fsize(0)),2)+'-D array,' + $
   		' with FORTRAN subscripts:'
   for i = 1, fsize(0) do begin
   	print, '(1:'+strtrim(string(fsize(i)),2)+')'
	endfor
   print, '3D images require that you ' + $
		'use a colon as separator for ranges, e.g. "5:7"'
   print, 'For multidimensional arrays, separate ranges with commas.'
   print, 'Type "ALL" if you want the entire array.'
   read,  'Enter subscript(s) > ', subscripts
   if subscripts eq '' then begin
	print, 'Type "ALL" to indicate all channels,'
	print, ' "BACK" to backup, "EXIT" to quit Import,'
	print, 'or type a subscript or range.'
	goto, subscr
	EndIf

   ; if ALL specified then make it blank to get all channels
   ; from data server
   if STRUPCASE(STRTRIM(subscripts,2)) eq 'ALL' then begin
	subscripts = ''
	goto, nocheck
	endif

   if STRUPCASE(STRTRIM(subscripts,2)) eq 'EXIT' then begin
		!ERR = 99
 		!ERROR = 99
		goto, abort 
		EndIf
   if STRUPCASE(STRTRIM(subscripts,2)) eq 'BACK' then begin
		!ERR = 88
		!ERROR = 88
		goto, abort 
		EndIf

   ; check subscripts against rdf info
   groups = strparse( subscripts, delim=',' )
   if n_elements(groups) ne fsize(0) then begin
	message, /continue, 'Too many or too few subscripts given.'
	goto, subscr
	endif

   for i=1,fsize(0) do begin
	dims = strparse(groups(i-1),delim=':')
	if n_elements(dims) gt 2 then begin
		message,/continue, 'Format of ranges not correct.'
		goto, subscr
		endif
	
	if n_elements(dims) eq 1 then begin
		if fix(dims) gt fsize(i) then begin
			message, /continue, 'Subscript out of range.'
			goto, subscr
			endif 
		goto, nocheck
		endif

	if fix(dims(0)) le 0 then begin
		message,/continue, 'Subscript cannot be <= 0'
		goto, subscr
		endif

	if fix(dims(1)) gt fsize(i) then begin
		message, /continue, 'Subscript out of range.'
		goto, subscr
		endif

	endfor

   EndIf

nocheck:
 


; everything a.o.k.
status = 1
!err = 1
!error = 1

return

abort:
; return to calling program
RETURN
;
END
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


