PRO change_cval, name
;+
;  CHANGE_CVAL - a UIMAGE-specific routine.  This routine allows
;  an X-window user to select a graph and then to select a new color
;  value at which the graph will be plotted.
;#
;  Written by John Ewing,  Nov 92.
;  SPR 10292  Dec 03 92  Re-worded a message.  J Ewing
;  SPR 10456  Jan 13 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON color_values,c_badpix,c_draw,c_scalemin
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    RETURN
  ENDIF
  first = 1
;  menu_title = 'Change an associated color value'
  menu_title = "Change a graph's color"
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Show a menu in which the user can select an available graph.
;  ------------------------------------------------------------
show_menu:
  name = select_object(/graph, /help, /exit, title = menu_title)
  IF(name EQ 'EXIT') THEN BEGIN
    IF(journal_on AND (first NE 1)) THEN $
      PRINTF, luj, '----------------------------------------' + $
                   '--------------------------------------'
    RETURN
  ENDIF
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no graphs currently available.', /CONT
    RETURN
  ENDIF
have_name:
  j = EXECUTE('hidden = ' + name + '.hidden')
  IF (hidden EQ 1) THEN BEGIN
    PRINT,'The selected graph has been hidden. Please display it' + $
        'before choosing this operation.'
    GOTO, show_menu
  ENDIF
;
;  Let the user select a color value to be used for that graph.
;  ------------------------------------------------------------
  cval = select_cval()
  IF(cval LT 0) THEN RETURN
;
;  Re-plot the graph (using the new color value).
;  ----------------------------------------------
  j = EXECUTE(name + '.color = cval')
  j = EXECUTE('window = ' + name + '.window')
  WSET, window
  plot_graph, name
;
;  If journaling is enabled, then report what was done to the journal file.
;  ------------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    IF(first) THEN PRINTF, luj, menu_title
    ig = FIX(STRMID(name, 6, STRPOS(name, ')') - 6))
    PRINTF, luj, '  color value ' + STRTRIM(STRING(cval),2) + $
      ' was assigned to:  ' + graph(ig).title
  ENDIF
  first = 0
  IF(uimage_version EQ 2) THEN GOTO, show_menu
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


