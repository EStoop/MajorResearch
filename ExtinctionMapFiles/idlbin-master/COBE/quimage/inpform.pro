;+
; NAME:
;	INPFORM
; PURPOSE:
;       This function gets and displays the input format menu 
; CATEGORY:
;	User interface, Menu, Data I/O
; CALLING SEQUENCE:
;	choice=INPFORM, input_format 
; INPUTS:
;       none
; OUTPUTS:
;	The input format in text string format.
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
;	16-Dec-1992	10330	Kryszak - added/edited help options
;	23-Dec-1992	10387	Kryszak - added help
;       11-Mar-1993     10672   Ewing - changed backward control options.
;
;-
;
Pro inpform, input_format
;
status = 0
title = 'Data Set Format'

mnuoptions = [	$
	'Select Input Format',	$
	'FITS Files',		$
	'COBE IDL Save Sets',	$
	'COBE Sky Maps',	$
;	'Pixel List FITS',	$
	'    ',			$
	'HELP',			$
	'Return to previous menu' $
	]

mnufunctions = [	$
	'title', $
	'FITS',	$
	'CISS',	$
	'CSM',	$
;	'CFITS',$
	'    ', $
	'HELP', $
	'Return to previous menu' $
	]

!err = 1
!error = 1
;
; display and get the user's input
menu:
status = umenu( mnuoptions, title=0, init=1, xpos=0 )
if status LE 0 then goto, abort ; return if errors
;

if mnufunctions(status) eq 'Return to previous menu' then begin
	!ERR = 88
	!ERROR = 88
	goto, abort
	endif

if ( strupcase(mnufunctions(status)) Eq 'HELP' ) Then Begin
	action, 'inpform.hlp'
	!err = 1
	!error = 1
	goto, menu
	EndIf


;
; set the return value 
input_format = mnufunctions(status)
;
;
abort:
IF status LE 0 THEN !ERR = -1 ;not ok
;
; return to calling program
RETURN 
;
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


