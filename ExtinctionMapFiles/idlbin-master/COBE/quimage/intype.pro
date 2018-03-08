;+
; NAME:
;	INTYPE
; PURPOSE:
;       This function gets and displays the input format menu 
; CATEGORY:
;	User interface, Menu, Data I/O
; CALLING SEQUENCE:
;	choice=INTYPE, ds_type 
; INPUTS:
;       None. 
; OUTPUTS:
;	The input format in text string format.
;       !Err = -1 if this function is called incorrectly
;	     = 88 is returned if "Return to previous menu" is chosen.
; COMMON BLOCKS:
;	None.
; RESTRICTIONS:
;	None.
; MODIFICATION HISTORY:
;       Creation:  Pete Kryszak, GSC, March 1992.
;
;	16-Dec-1992	10330	Kryszak - added/edited menu options
;	23-Dec-1992	10387	Kryszak - added help
;       11-Mar-1993     10672   Ewing - changed backward control options
;       13-apr-1993     10819   Ward  - made changes to handle fits files
;                                       with binary extensions, cleaned
;                                       up backward control options
;
;-
;
Pro intype, ds_type

;print, 'in intype.pro'

!err = 1
!err = 1

status = 0

mnuoptions = [	$
	'Select Data Type', $
	'FIRAS',	$
	'DIRBE',	$
	'DMR',	$
	'ADB',	$
	'User Defined', $
	'    ', $
	'HELP', $
	'Return to previous menu' $
	]

mnufunctions = [	$
	'title', $
	'F',	$
	'B',	$
	'D',	$
	'X',	$
	'U',    $
	' ',	$
	'HELP', $
	'Return to previous menu' $
	]

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
	action, 'intype.hlp'
	!err = 0
	!error = 0
	goto, menu
	EndIf


;
; set the return value 
ds_type = mnufunctions(status)
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


