PRO datamana
;+
;  DATAMANA - a UIMAGE-specific routine.  This routine puts up the
;  menu labeled 'DATA I/O AND MANAGEMENT', and takes appropriate
;  action depending on selections made by the user.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 06 93  Add more comments.  J Ewing
;  SPR 10456  Jan 13 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10829  Apr  1993 Added 2 new fields, Assocwts and hidden. J Newmark
;  SPR 10912  May 10 93  Add new field, MKFACE. J Newmark
;  SPR 11031  Jun 08 93  Add call to REM_ALL J Newmark
;  SPR 11103  Jun 28 93  Add call to MERGE2D J. Newmark
;--------------------------------------------------------------------------
;
;  Define routines which can be called from here, as well as the
;  phrases, to appear in a menu, which identify those routines.
;  [The VT field below tells whether that routine should be accessible
;  from a VT200 (non-X-window) terminal.]
;  -------------------------------------------------------------------
  n_routines = 14
  routine = REPLICATE({rs, vt: 0, name: '', title: ''}, n_routines)
  routine(00) = {rs, 1, 'IMPORT',        'Read a data set...'}
  routine(01) = {rs, 1, 'EXPORT',        'Write a data file...'}
  routine(02) = {rs, 1, 'UPOSTSCRIPT',   'Create a PostScript file'}
  routine(03) = {rs, 1, 'ASSOCWTS',      'Associate weight with an object'}
  routine(04) = {rs, 1, 'MKFACE',        'Extract Face from skymap'}
  routine(05) = {rs, 1, 'MERGE2D',       'Combine 2D objects into 3D object'}
  routine(06) = {rs, 0, '',              ''}
  routine(07) = {rs, 1, 'DISPLAY_CHAR',  'Report object attributes'}
  routine(08) = {rs, 1, 'DELETE_OBJECT', 'Remove an object'}
  routine(09) = {rs, 1, 'HIDDEN',        'Display/hide an object'}
  routine(10) = {rs, 1, 'REM_ALL',       'Remove ALL objects'}
  routine(11) = {rs, 0, '',              ''}
  routine(12) = {rs, 1, 'HELP',          'HELP'}
  routine(13) = {rs, 1, 'EXIT',          'Return to MAIN MENU'}
;
;  Create a menu of selectable operations.
;  ---------------------------------------
  menu_title = 'DATA I/O AND MANAGEMENT'
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
  routine_name = routine(sel-1).name
  IF(routine_name EQ '')     THEN GOTO, show_menu
  IF(routine_name EQ 'EXIT') THEN RETURN
  IF(routine_name EQ 'HELP') THEN uimage_help, menu_title $
                             ELSE CALL_PROCEDURE, routine_name
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


