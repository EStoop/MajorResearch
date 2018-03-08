FUNCTION umenu,strings,title=title,init=init,valid=valid, $
               xpos=xpos, ypos=ypos, demowait=demowait
;+
; NAME:	
;       UMENU
; PURPOSE:
;       This function displays a menu whose choices are given by the
;       elements of a string array and then returns the index of that
;       choice which the user selects.  If the program is being run in
;       an X-window environment, then WMENU is called to accomplish this
;       objective, otherwise TMENU is called.
; CATEGORY:
;	User interface, Menu.
; CALLING SEQUENCE:
;	i=UMENU(['Select one','Item_1','Item_2'],title=0,init=1,xpos=10)
; INPUTS:
;       STRINGS  = A string array with each element containing the text
;       of one menu choice.
;       TITLE    = The index of the STRINGS element that is the title,
;       normally 0.  The title element is not selectable and is displayed
;       bolded and underlined.  If this parameter is omitted, all items
;       are selectable.
;       INIT     = The index of the initial selection.  If this parameter
;       is specified and within the range of STRINGS indices (excluding
;       the TITLE index), then the initial menu display is made with the
;       designated item selected.
;       VALID    = List of 0's and 1's, which designate valid choices.
;       (optional)  If supplied, then valid items will be shown in
;       bold type.  Invalid items will not be selectable.
; OUTPUTS:
;       The index number of the selected menu item.  A value of -1 is
;       returned if this function is called incorrectly.
; COMMON BLOCKS:
;	None.
; RESTRICTIONS:
;       (the following restrictions are for a TEKTRONIX environment only)
;       Menu choices must be less than 52 characters.
;       There must be less than 21 menu choices.
;       The maximum allowable value of XPOS is limited for different cases.
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, November 1991.
;       Eliminated a call to EXECUTE:  John A. Ewing, ARC, January 1992.
;       Date       SPR	  Name	   Description
;	29-Dec-92  10336  Kryszak  Uses QMENU instead of BIGWMENU
;       30-Dec-92  10394  Turpie   Added ypos keyword.
;-
;
;  Check to see is an appropriate array of strings was passed in for "STRINGS"
;  --------------------------------------------------------------------------
  sz = SIZE(strings)
  IF(sz(0) NE 1) OR (sz(2) NE 7) THEN BEGIN
    PRINT, '%UMENU:  An inappropriate argument was supplied.'
    RETURN, -1
  ENDIF
  old_quiet = !QUIET
  !QUIET = 1
;
;  Call the appropriate function
;  -----------------------------
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
    IF(!D.WINDOW NE -1) THEN DEVICE,CURSOR_STANDARD=22
    select = qmenu(strings,title=title,init=init,demowait=demowait)
    IF(!D.WINDOW NE -1) THEN DEVICE,CURSOR_STANDARD=33
  ENDIF ELSE BEGIN
    select = tmenu(strings,title=title,init=init,valid=valid, $
                   xpos=xpos,ypos=ypos,demowait=demowait)
  ENDELSE
  !QUIET = old_quiet
  RETURN, select
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


