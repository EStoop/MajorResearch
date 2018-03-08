PRO change_labels, name
;+
;  CHANGE_LABELS - a UIMAGE-specific routine.  This routine gives the
;  user the opportunity to change that X and Y axis labels for a
;  graph.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 06 93  Add more comments.  J Ewing
;  SPR 10456  Jan 13 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
;
;  If the name of graph to be operated upon was not supplied when
;  this routine was invoked, then put up a menu in which the user
;  can identify the desired graph.
;  --------------------------------------------------------------
get_argument:
  menu_title = 'Change graph axis labels'
  IF(N_PARAMS(0) GT 0) THEN GOTO, have_name
  name = select_object(/graph, /help, /exit, title = menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, get_argument
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no graphs available.', /CONT
    RETURN
  ENDIF
;
;  Extract pertinent information from the GRAPH structure.
;  -------------------------------------------------------
have_name:
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('window = ' + name + '.window')
  j = EXECUTE('x_title = ' + name + '.x_title')
  j = EXECUTE('y_title = ' + name + '.y_title')
  IF(uimage_version LT 4) THEN BEGIN
    PRINT, ' '
    PRINT, 'Selected graph:  ' + bold(title)
    PRINT, ' '
  ENDIF
;
;  Display the X-axis label, and allow the user to enter in a new
;  value if desired.
;  --------------------------------------------------------------
  PRINT, 'X-axis label:  ' + bold(x_title)
  menu = ['Do you wish to change the X-axis label?', 'Yes', 'No']
  selx = one_line_menu(menu)
  IF(selx EQ 1) THEN BEGIN
    READ, 'Enter X-axis label:  ', x_title
    j = EXECUTE(name + '.x_title = x_title')
  ENDIF
  PRINT, ' '
;
;  Display the Y-axis label, and allow the user to enter in a new
;  value if desired.
;  --------------------------------------------------------------
  PRINT, 'Y-axis label:  ' + bold(y_title)
  menu = ['Do you wish to change the Y-axis label?', 'Yes', 'No']
  sely = one_line_menu(menu)
  IF(sely EQ 1) THEN BEGIN
    READ, 'Enter Y-axis label:  ', y_title
    j = EXECUTE(name + '.y_title = y_title')
  ENDIF
;
;  Re-display the graph.
;  ---------------------
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN WSET, window
  ERASE
  plot_graph, name, /redraw
;
;  If journaling is enabled, then send some info to the journal file.
;  ------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title
    IF(selx EQ 1) $
      THEN PRINTF, luj, '  new X-axis label:  ' + x_title $
      ELSE PRINTF, luj, '  the X-axis label was left unchanged'
    IF(sely EQ 1) $
      THEN PRINTF, luj, '  new Y-axis label:  ' + y_title $
      ELSE PRINTF, luj, '  the Y-axis label was left unchanged'
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
  ENDIF
  IF(uimage_version EQ 2) THEN GOTO, get_argument
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


