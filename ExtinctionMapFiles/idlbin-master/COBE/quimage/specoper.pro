PRO specoper
;+
;  SPECOPER - a UIMAGE-specific routine.  This routine puts up the
;  menu labeled 'Spectrum operations', and takes appropriate
;  action depending on selections made by the user.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 06 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;--------------------------------------------------------------------------
;
;  Define routines which can be called from here, as well as the
;  phrases, to appear in a menu, which identify those routines.
;  [The VT field below tells whether that routine should be accessible
;  from a VT200 (non-X-window) terminal.]
;  -------------------------------------------------------------------
  n_routines = 9
  routine = REPLICATE({rs, vt: 0, name: '', title: ''}, n_routines)
  routine(00) = {rs, 1, 'EXTRACT_SPEC',  'Extract spectrum from a pixel'}
  routine(01) = {rs, 1, 'AREA',          'Average spectra in an area'}
  routine(02) = {rs, 0, '',              ''}
  routine(03) = {rs, 1, 'SLICE_EXTRACT', 'Extract a frequency/wavelength slice'}
  routine(04) = {rs, 1, 'INTEG_POWER',   'Integrate over frequency/wavelength'}
  routine(05) = {rs, 1, 'DISPLAY_FREQ',  'Display a frequency/wavelength table'}
  routine(06) = {rs, 0, '',              ''}
  routine(07) = {rs, 1, 'HELP',          'HELP'}
  routine(08) = {rs, 1, 'EXIT',          'Return to MAIN MENU'}
;
;  Create a menu of selectable operations.
;  ---------------------------------------
  menu_title = 'Spectrum operations'
  menu = STRARR(n_routines + 1)
  menu(0) = menu_title
  FOR i = 0, n_routines - 1 DO menu(i + 1) = routine(i).title
;
;  Display the menu.
;  -----------------
show_menu:
  sel = UMENU(menu, title = 0, init = sel, valid = valid)
;
;  Invoke the routine that implements the selected operation.
;  ----------------------------------------------------------
  routine_name = routine(sel - 1).name
  IF(routine_name EQ '')     THEN GOTO, show_menu
  IF(routine_name EQ 'EXIT') THEN RETURN
  IF(routine_name EQ 'HELP') THEN uimage_help, menu_title $
                             ELSE $
    IF(routine_name NE 'AREA') THEN CALL_PROCEDURE, routine_name $
                               ELSE extract_spec, /area
  GOTO, show_menu
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


