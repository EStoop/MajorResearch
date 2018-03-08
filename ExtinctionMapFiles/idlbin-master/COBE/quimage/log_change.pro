PRO log_change
;+
;  LOG_CHANGE a UIMAGE-specific routine.  This routine allows a user
;  to specify the type axis for a graph.
;#
; SPR 10953 Written by J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
;
;  Have the user select a graph.
;  -----------------------------
show_menu:
  menu_title = 'Change Log Scaling'
  name = select_object(/graph, /help, /exit, title = menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no graphs available.', /CONT
    RETURN
  ENDIF
;
;  Extract attributes.
;  -------------------
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('window = ' + name + '.window')
  j = EXECUTE('multi = ' + name + '.multi')
  PRINT, ' '
  PRINT, 'Selected graph:  ' + bold(title)
;
;  Let the user set the X range if he wants to.
;  --------------------------------------------
  lmenu = ['Choose Type of Scaling',$
          'Linear-Linear','Linear-Log','Log-Linear','Log-Log']
  sellog = umenu(lmenu,title=0,init=sellog)
  IF (multi NE 0) THEN BEGIN
   FOR i=0,multi-1 DO BEGIN
    j=EXECUTE('tempid=' + name + '.m_id(i)')
    tempname=get_name(tempid)
    CASE sellog OF
     1: j=EXECUTE(tempname + '.logflag = 0')
     2: j=EXECUTE(tempname + '.logflag = 1')
     3: j=EXECUTE(tempname + '.logflag = 2')
     4: j=EXECUTE(tempname + '.logflag = 3')
    ENDCASE
   ENDFOR
  ENDIF
  CASE sellog OF
   1: j=EXECUTE(name + '.logflag = 0')
   2: j=EXECUTE(name + '.logflag = 1')
   3: j=EXECUTE(name + '.logflag = 2')
   4: j=EXECUTE(name + '.logflag = 3')
  ENDCASE
;
;  Re-display the graph.
;  ---------------------
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN WSET, window
  PLOT_GRAPH, name
;
;  If journaling is enabled, then send some info to the journal file.
;  ------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title
    IF(sellog EQ 1) THEN PRINTF, luj, '  Linear-Linear Plot' 
    IF(sellog EQ 2) THEN PRINTF, luj, '  Linear-Log Plot' 
    IF(sellog EQ 3) THEN PRINTF, luj, '  Log-Linear Plot' 
    IF(sellog EQ 4) THEN PRINTF, luj, '  Log-Log Plot' 
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
  ENDIF
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


