;+
; NAME:
;  SELOUTPT 
; PURPOSE:
;       This function gets all the output data set information from
;	the user. 
; CATEGORY:
;  User interface, Menu, Data I/O
; CALLING SEQUENCE:
;  seloutpt ...
; INPUTS:
;   none
;
; OUTPUTS:
;   outform 		format, FITS, CISS
;   outfname		user entered filename
;
; COMMON BLOCKS:
;  None.
; RESTRICTIONS:
;  None.
; MODIFICATION HISTORY:
;       Creation:  Pete Kryszak, GSC, March 1992.
;
;	16-Dec-1992	10330	Kryszak - Quit functionality
;			10310	Kryszak - user is now required to enter
;				an extension
;	23-Dec-1992	10387	Kryszak - added help, needed changes here too
;	30-Dec-1992	10403   Kryszak - shows default directory
;       22-feb-1993     10819   Ward - modifications for fits io (binary
;                                      tables)
;
;-
;
Pro SelOutpt, outform, outfname

outform = ''
outfname = ''
old_quiet=!quiet
outform, outform 

IF !err LT 0 THEN Begin
        !quiet=0
	message, 'Output format invalid.', /informational 
        !quiet=old_quiet
	GOTO, abort ; return if errors
	EndIf
IF (!err EQ 99) or (!err EQ 88) THEN Begin
        !quiet=0
	message, 'The user selected Quit.' , /informational
        !quiet=old_quiet
	GOTO, abort; return if exit
	EndIf

; get the current directory
cd, current=default

print, ' '
print, 'The file will be written to ' + default 
print, 'or a directory that you specify in the file specification.'
print, 'Include an appropriate file name extension, e.g. FITS, or CISS.'

read, 'Data set name? ', outfname

return

abort:
	outform = ''
	outfname = ''
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


