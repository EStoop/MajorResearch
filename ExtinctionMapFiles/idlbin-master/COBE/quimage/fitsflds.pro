;+
; NAME:
;	FITSFLDS
; PURPOSE:
;       This routine presents the user with the fields which are available
;	in the selected file to be read in.  If a FIRAS file (with multiple
;	freqs / pixel) then the freqs to read in will also be selected here.
;	
; CATEGORY:
;	User interface, Menu, DATA I/O
;
; CALLING SEQUENCE:
;	fitsflds, bin_header, flist, fcount
;	
; INPUTS:
;	bin_header: binary extension header, to be displayed if help selected
;	flist:  field list of TTYPE field in binary table file
;	fcount: number of fields passed
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
;       Creation:  Dalroy Ward, GSC, March 1993. SPR 10819
;-
;
	Function fitsflds, bin_header, flist, fcount

menu:
	; setup the menu
	field_list = strarr( fcount + 4)
	field_list( 0 ) = 'Field List'
	field_list( fcount+1 ) = '    '
	field_list( fcount+2 ) = 'HELP'
	field_list( fcount+3 ) = 'Return to previous menu'

	for i = 0, fcount-1 do begin
	    field_list(i+1) = flist(i)
	endfor 

jmenu:
	; present the menu
	status = umenu( field_list, title = 0, init=1, xpos=0)
	if ( status LE 0 ) then Begin
	    !ERR = -1
            !ERROR = -1
	    goto, abort		; return if errors
	EndIf

	field = status

	if ( field_list(status) Eq 'HELP' ) Then Begin
	    action, 'fldname.hlp'
	    uscroll, bin_header, title='Extension Header',quit=tquit
	    goto, jmenu
	EndIf

	if ( field_list(status) Eq 'Return to previous menu' ) Then Begin
	    !ERR = 88
	    !ERROR = 88
	    goto, abort 
	EndIf

nocheck:

	; set return status to ok
	status = 1
	!err = 1
	!error = 1
	return, field

abort:
	return, -1

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


