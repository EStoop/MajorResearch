PRO resize_graph, name
;+
; NAME:	
;	RESIZE_GRAPH
; PURPOSE:
;       To allow a UIMAGE user to set the size of a graphics window.
; CATEGORY:
;       Graphics.
; CALLING SEQUENCE:
;       resize_graph, name
; INPUTS:
;       NAME      = Name of a UIMAGE structure (e.g., "MAP6(0)").
; OUTPUTS:
;       None.
;#
; COMMON BLOCKS:
;       JOURNAL, UIMAGE_DATA, XWINDOW
; RESTRICTIONS:
;       To be called from UIMAGE-routines only.
; PROCEDURE:
;       The SELECT_OBJECT function is called to allow the user to select
;       the desired data object.  If a selection is made, then the data
;       associated with that object is returned through the variable ARR.
;       If there were no data objects stored in the common statement, then
;       a message is printed out and ARR will be undefined.
; SUBROUTINES CALLED:
;       make_window, plot_graph, select_object, umenu
; MODIFICATION HISTORY:
;  Creation:  John A. Ewing, ARC, February 1992.
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  IF(NOT defined(journal_on)) THEN journal_on = 0
  first = 1
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    RETURN
  ENDIF
;
;  If the name of the operand was not supplied when this routine was invoked,
;  then put up a menu in which the user can select the desired graph.
;  --------------------------------------------------------------------------
  IF(N_PARAMS(0) EQ 0) THEN BEGIN
get_argument:
    menu_title = 'Resize a graph'
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
    j = EXECUTE('title = ' + name + '.title')
    PRINT, ' '
    PRINT, 'Selected graph:  ' + bold(title)
  ENDIF
  orig_bu = block_usage
;
;  Make a list of potential sizes, depending on the window magnification factor.
;  -----------------------------------------------------------------------------
  CASE zoom_index OF
    0: numblocks = [4,5,6,7,8,9,10,11,12]
    1: numblocks = [6,8,10,12,14]
    2: numblocks = [8,12,16]
    3: numblocks = [8,16]
  ENDCASE
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('window = ' + name + '.window')
resize:
;
;  Put up a menu in which the user can select the desired X size.
;  --------------------------------------------------------------
  m_title = 'Select X size'
  menu_size = [m_title, STRTRIM(numblocks, 2), ' ', 'HELP', $
    'Exit this menu']
  selx = umenu(menu_size, title = 0) - 1
  IF(menu_size(selx + 1) EQ ' ') THEN GOTO, resize
  IF(menu_size(selx + 1) EQ 'Exit this menu') THEN RETURN
  IF(menu_size(selx + 1) EQ 'HELP') THEN BEGIN
    uimage_help, m_title
    GOTO, resize
  ENDIF
get_y:
;
;  Put up a menu in which the user can select the desired Y size.
;  --------------------------------------------------------------
  m_title = 'Select Y size'
  menu_size = [m_title, STRTRIM(numblocks, 2), ' ', 'HELP', $
    'Exit this menu']
  sely = umenu(menu_size, init = selx + 1, title = 0) - 1
  IF(menu_size(sely + 1) EQ ' ') THEN GOTO, get_y
  IF(menu_size(sely + 1) EQ 'Exit this menu') THEN RETURN
  IF(menu_size(sely + 1) EQ 'HELP') THEN BEGIN
    uimage_help, m_title
    GOTO, get_y
  ENDIF
  nblocks_x_abs = numblocks(selx)
  nblocks_y_abs = numblocks(sely)
;
;  Free up the space in BLOCK_USAGE that was taken up previously by the graph.
;  ---------------------------------------------------------------------------
  w = WHERE(block_usage EQ window + 1)
  IF(w(0) NE -1) THEN block_usage(w) = 0
;
;  Create a new window which has the specified size.
;  -------------------------------------------------
  make_window, window, nblocks_x_abs, nblocks_y_abs, title
  j = EXECUTE(name + '.nblocks_x = numblocks(selx) / (2. ^ zoom_index)')
  j = EXECUTE(name + '.nblocks_y = numblocks(sely) / (2. ^ zoom_index)')
;
;  If the operand was a graph, then display the graph.
;  ---------------------------------------------------
  IF(STRMID(name, 0, 1) EQ 'G') THEN plot_graph, name
  sel = umenu([' ', 'Resize again', 'Exit this menu'], init = 1, title = 0)
  IF(sel EQ 1) THEN BEGIN
    block_usage = orig_bu
    GOTO, resize
  ENDIF
  IF((N_PARAMS(0) EQ 0) AND journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title
    PRINTF, luj, '  X size (blocks) = ' + STRTRIM(nblocks_x_abs, 2)
    PRINTF, luj, '  Y size (blocks) = ' + STRTRIM(nblocks_y_abs, 2)
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
  ENDIF
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


