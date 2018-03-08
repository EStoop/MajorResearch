FUNCTION one_line_menu, strings, init = init, demowait = demowait
;+
; NAME:	
;	one_line_menu
; PURPOSE:
;       This function displays a menu whose choices are given by the
;       elements of a string array, and then returns the index of that
;       choice which the user selects.  The menu is laid out horizontally.
;       The first element in the STRINGS array is interpreted as a
;       prompt string.
; CATEGORY:
;	Menu, User interface.
; CALLING SEQUENCE:
;	i = one_line_menu(strings, init = init)
; INPUTS:
;       STRINGS  = A string array.  The first element is interpreted
;       as a prompt string, and other elemnents are options.  There
;       must be at least 3 elements altogether.
;       INIT = index of the element in STRINGS which will be the initial
;       choice.
;       DEMOWAIT = delay (in seconds) after which the program will exit,
;       returning the value of INIT (for demo purposes only).
; OUTPUTS:
;       The index number of the selected menu item.  A value of -1 is
;       returned if this function is called incorrectly.
;
; EXAMPLE:
;       sel = one_line_menu(['{B}Question 1:{N}  Which option do you' + $
;         'prefer?','Portrait','Landscape'])
; COMMON BLOCKS:
;	None.
; RESTRICTIONS:
;       There is a limit on the number of characters allowable.
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, March 1992.
;       Modified so cursor can jump from last item to first and vice versa
;          John Ewing  Oct. 23, 1992.
;       SPR 10260 - changed GET_KBRD calls to use 1 instead of 0 (JAE)
;       SPR 11003  Jun 02 93  Changed !version.os to !cgis_os. J Newmark
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;  SPR 11296 09 Sep 1993  Call UMENU for actual menu display. J. Newmark
;-
  IF (!cgis_os EQ 'vms') THEN cr = STRING(13b) ELSE cr = STRING(10b)
  IF (!cgis_os EQ 'windows') THEN cr=''
;
;  Check the validity of what was supplied for the "STRINGS" argument.
;  -------------------------------------------------------------------
  IF(N_PARAMS(0) NE 1) THEN BEGIN
    MESSAGE, 'One argument (an array of strings) must be supplied.', /CONT
    RETURN, -1
  ENDIF
  sz = SIZE(strings)
  IF(sz(0) NE 1) OR (sz(2) NE 7) THEN BEGIN
    MESSAGE, 'An inappropriate argument was supplied.', /CONT
    RETURN, -1
  ENDIF
  nstr = sz(1)
  IF(nstr LE 2) THEN BEGIN
    MESSAGE, 'Not enough strings were supplied (should be at least 3).', /CONT
    RETURN, -1
  ENDIF
  total_len = STRLEN(strings(0)) + 5
  IF(total_len GT 80) THEN BEGIN
    MESSAGE, 'Too many characters.', /CONT
    RETURN, -1
  ENDIF
  IF(NOT KEYWORD_SET(init)) THEN init = 1
  init = (init>1)<(nstr-1)
  strlist = strings
  IF(KEYWORD_SET(demowait)) THEN BEGIN
      WAIT, demowait
      sel=init
      GOTO, exit
  ENDIF
;
; Call UMENU to put up a Widget for a X window or else a vertical menu.
;
  first=init
  sel=UMENU(strlist,init=first,title=0)
exit:
  RETURN, sel
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


