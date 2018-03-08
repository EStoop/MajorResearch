;+
; NAME:
;	FACE 
; PURPOSE:
;       This function allows the user to specify the face number. 
; CATEGORY:
;	User interface, Menu, Data I/O
; CALLING SEQUENCE:
;	FACE, face_number 
; INPUTS:
;       None. 
; OUTPUTS:
;	The face_number, 0-5, or 7 meaning all faces.
;       !Err = -1 if this function is called incorrectly
;	     = 88 is returned "Return to previous menu" is chosen.
;
; COMMON BLOCKS:
;	None.
; RESTRICTIONS:
;	None.
; MODIFICATION HISTORY:
;       Creation:  Pete Kryszak, GSC, March 1992.
;	
;	16-Dec-1992 	10330	Kryszak - added/edited menus
;	21-Dec-1992	10387	Kryszak - added help screen
;       11-Mar-1993     10672   Ewing - changed backward control options
;
;-
;
Pro face, face_number 
;
status = 0

mnuoptions = [	$
	'Select Map Face', 		$
	'All -   Right-handed T', 	$
	'Zero -  North Ecliptic Pole', 	$
	'One -   South Galactic Pole',  $
	'Two -   Galactic Anti-Center',	$
	'Three - North Galactic Pole',	$
	'Four -  Galactic Center', 	$
	'Five -  South Ecliptic Pole', 	$
	'    ',				$
	'HELP',				$
	'Return to previous menu' 	$
	]

mnufunctions = [	$
	'title',$
	'7',	$
	'0',	$
	'1',	$
	'2',	$
	'3',	$
	'4',	$
	'5',	$
	' ',	$
	'HELP', $
	'Return to previous menu' $
	]

!ERR = 1
!ERROR = 1
	
face_number = 7 ; means all faces (default)

menu:
status = umenu( mnuoptions, title=0, init=1, xpos=0 )
if (status le 0)  then Begin
	!ERR = -1
	!ERROR = -1
	goto, abort
	EndIf

If ( strupcase(mnufunctions(status)) Eq 'HELP' ) Then Begin
    action, 'face.hlp'
    Goto, menu
    EndIf

if (mnufunctions(status) eq 'Return to previous menu') then begin
	!ERR = 88
	!ERROR = 88
	goto, abort
	endif


;
fnumber = mnufunctions(status)
face_number = FIX(fnumber)
;
;
;
abort:
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


