PRO change_title, name
;+
; NAME:	
;	CHANGE_TITLE
; PURPOSE:
;       To allow a UIMAGE user to change the title associated with
;       a selected object.
; CALLING SEQUENCE:
;       CHANGE_TITLE
; INPUTS:
;       None.
; OUTPUTS:
;       None.
; COMMON BLOCKS:
;	UIMAGE_DATA, XWINDOW
; RESTRICTIONS:
;       To be called from UIMAGE only.
; PROCEDURE:
;       Call SELECT_OBJECT so the user can identify the object whose
;       title he/she wants to change.  Then have the user type in the
;       new title.  Set the title field in the appropriate structure.
;       Redisplay the image in a ne window if an X-window terminal is
;       being used.
;#
; SUBROUTINES CALLED:
;       SELECT_OBJECT, UIMAGE_HELP, XDISPLAY
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, January 1992.
;  SPR 10383  Jan 06 93  Add more comments.  J Ewing
;  SPR 10456  Jan 13 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 11 93  Change journaling.  J Ewing
;  SPR 10829  Apr  1993  Change info_3d common block
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = 'Change a title'
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Put up a menu so the user can select an operand.
;  ------------------------------------------------
show_menu:
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7, $
           /face8,/face9,/graph,/zoomed,/proj_map,/proj2_map,/object3d,$
           /help,/exit, title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no objects available.', /CONT
    RETURN
  ENDIF
have_name:
  j = EXECUTE('oldtitle = ' + name + '.title')
  title = ''
get_title:
;
;  Have the user type in the new title.
;  ------------------------------------
  READ, 'Enter new title:  ', title
  IF (title EQ '') THEN BEGIN
    PRINT,'Entered a null string. Please enter a character string for ' + $
       'the new title.'
    GOTO, get_title
  ENDIF
;
;  Reject any user-defined titles that contain brackets (to make sure
;  that brackets in titles will always enclose an integer number).
;  ------------------------------------------------------------------
  ic1 = STRPOS(title, ']')
  ic2 = STRPOS(title, '[')
  IF((ic1 NE -1) OR (ic2 NE -1)) THEN BEGIN
    PRINT, bold('Please do not put brackets in user-defined titles.')
    PRINT, ' '
    GOTO, get_title
  ENDIF
;
;  Update the structure to contain the new title.
;  ----------------------------------------------
  j = EXECUTE(name + '.title = title')
;
;  For X-windows, display image with the new title.
;  ------------------------------------------------
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
    j = EXECUTE('window = ' + name + '.window')
    w = WHERE(block_usage EQ window + 1)
    IF(w(0) NE -1) THEN block_usage(w) = 0
    WDELETE, window
    first_letter = STRMID(name, 0, 1)
    CASE first_letter OF
      'M' : xdisplay, name
      'F' : xdisplay, name
      'O' : xdisplay, name
      'P' : xdisplay, name
      'G' : BEGIN
              zoom = 2.^zoom_index
              j = EXECUTE('nblocks_x = ' + name + '.nblocks_x')
              j = EXECUTE('nblocks_y = ' + name + '.nblocks_y')
              make_window, window, nblocks_x*zoom, nblocks_y*zoom, title
              plot_graph, name
            END
      'Z' : BEGIN
              zoom = 2.^zoom_index
              j = EXECUTE('nblocks_x = ' + name + '.nblocks_x')
              j = EXECUTE('nblocks_y = ' + name + '.nblocks_y')
              make_window, window, nblocks_x*zoom, nblocks_y*zoom, title
              show_zoomed, name
            END
      ELSE:
    ENDCASE
  ENDIF
;
;  If journaling is enabled, then send some info to the journal file.
;  ------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:    ' + oldtitle
    PRINTF, luj, '  new title:  ' + title
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


