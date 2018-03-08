PRO ushow3, name
;+
; NAME:	
;	USHOW3
; PURPOSE:
;       To call SHOW3 while UIMAGE is running.
; CATEGORY:
;	Image display.
; CALLING SEQUENCE:
;       USHOW3
; INPUTS:
;       None.
; OUTPUTS:
;       None.
;#
; COMMON BLOCKS:
;	history, journal, uimage_data, uimage_data2, xwindow
; RESTRICTIONS:
;       To be called from UIMAGE only.
; PROCEDURE:
;       Call SELECT_OBJECT so the user can identify the desired object.
;       Then call SHOW3 with the data from that object.
; SUBROUTINES CALLED:
;       defined, select_object, show3, uimage_help
; MODIFICATION HISTORY:
;  Creation:  John A. Ewing, ARC, January 1992.
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  COMMON zback,zbgr,zdsrcmap,zbsub
  IF(NOT defined(journal_on)) THEN journal_on = 0
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    RETURN
  ENDIF
show_menu:
;
;  If the name of the desired operand was not supplied when this routine
;  was called, then put up a menu in which the user can select the desired
;  image.
;  -----------------------------------------------------------------------
  menu_title = '3-D surface plot'
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
  name=select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
    /face8,/face9,/zoomed,/proj_map,/proj2_map,/help,/exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no appropriate data objects to choose from.', /CONT
    RETURN
  ENDIF
have_name:
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('window = ' + name + '.window')
  WSET, window
;
;  Define numerous system variables.
;  ---------------------------------
  xstyle = !X.STYLE
  ystyle = !Y.STYLE
  zstyle = !Z.STYLE
  !X.STYLE = 1 OR 4
  !Y.STYLE = 1 OR 4
  !Z.STYLE = 1 OR 4
  !X.TITLE = ' '
  !Y.TITLE = ' '
  !Z.TITLE = ' '
;
;  Call SHOW3 in order to display the "3-D" representation.
;  --------------------------------------------------------
  IF(STRMID(name, 0, 6) NE 'ZOOMED') THEN BEGIN
    j = EXECUTE('pmin = ' + name + '.scale_min')
    j = EXECUTE('pmax = ' + name + '.scale_max')
    j = EXECUTE('SHOW3, BYTSCL(' + name + '.data, pmin, pmax)')
  ENDIF ELSE BEGIN
    j = EXECUTE('win_orig = ' + name + '.win_orig')
    j = EXECUTE('start_x = ' + name + '.start_x')
    j = EXECUTE('stop_x = ' + name + '.stop_x')
    j = EXECUTE('start_y = ' + name + '.start_y')
    j = EXECUTE('stop_y = ' + name + '.stop_y')
    start_x = start_x > 0
    start_y = start_y > 0
    bflag=0
    IF (win_orig LT 0) THEN BEGIN
      win_orig=ABS(win_orig)
      bflag=1
    ENDIF
    j=EXECUTE('zoomflag = ' + name + '.zoomflag')
    IF (zoomflag NE 0) THEN bflag=1
    name_orig = get_name(win_orig)
    j = EXECUTE('sz = SIZE(' + name_orig + '.data)')
    dim1 = sz(1)
    dim2 = sz(2)
    stop_x = stop_x < (dim1 - 1)
    stop_y = stop_y < (dim2 - 1)
    j = EXECUTE('data = ' + name_orig + '.data(start_x:stop_x, start_y:stop_y)')
    j = EXECUTE('pmin = ' + name_orig + '.scale_min')
    j = EXECUTE('pmax = ' + name_orig + '.scale_max')
    IF (bflag EQ 1) THEN BEGIN
       IF (zoomflag EQ 1) THEN data=zbgr
       IF (zoomflag EQ 2) THEN data=zdsrcmap
       IF (zoomflag EQ 3) THEN data=zbsub
       pmin=MIN(data,MAX=pmax)
    ENDIF
    SHOW3, BYTSCL(data, pmin, pmax)
  ENDELSE
;
;  Reset system variables.
;  -----------------------
  !X.STYLE = xstyle
  !Y.STYLE = ystyle
  !Z.STYLE = zstyle
  !X.TITLE = ''
  !Y.TITLE = ''
  !Z.TITLE = ''
  !X.RANGE = 0
  !Y.RANGE = 0
  !Z.RANGE = 0
  !P.TITLE = ''
;
;  If journaling is enabled, then send some info to the journal file.
;  ------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title
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


