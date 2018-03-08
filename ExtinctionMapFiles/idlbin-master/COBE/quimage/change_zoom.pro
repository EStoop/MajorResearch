PRO change_zoom
;+
; NAME:
;      CHANGE_ZOOM
; PURPOSE:
;      Change the ZOOM_INDEX variable.
; CATEGORY:
;      Image processing.
; CALLING SEQUENCE:
;      Change_zoom
; INPUTS:
;      None.
; OUTPUTS:
;      None.
;#
; COMMON BLOCKS:
;      XWINDOW
; RESTRICTIONS:
;      None.
; PROCEDURE:
;      Call UMENU to allow the user to select a zoom level.
; SUBROUTINES CALLED:
;      UIMAGE_HELP,UMENU,XDISPLAY_ALL
; MODIFICATION HISTORY: 
;      Original:  John Ewing, ARC, January 1992.
;  SPR 10442  Feb 04 93  Change journaling.  J Ewing
;----------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(journal_on)) THEN journal_on = 0
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  first = 1
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    RETURN
  ENDIF
  IF(NOT defined(zoom_index)) THEN BEGIN
    MESSAGE, 'ZOOM_INDEX is undefined.', /CONT
    RETURN
  ENDIF
;
;  Put up a menu in which the user can select the desired window magnification
;  factor [1, 2, 4, or 8].
;  ---------------------------------------------------------------------------
show_menu:
  menu_title = 'Window magnification factor'
  zmenu = [menu_title,'1','2','4','8',' ','HELP','Return to previous menu']
  zsel = umenu(zmenu, TITLE = 0, INIT = zoom_index + 1)
  CASE zsel OF
    6   : BEGIN
            uimage_help, menu_title
            GOTO, show_menu
          END
    7   : BEGIN
            IF(journal_on AND (first NE 1)) THEN $
              PRINTF, luj, '----------------------------------------' + $
                           '--------------------------------------'
            RETURN
          END
    ELSE: 
  ENDCASE
  IF((zsel - 1) GT zoom_index) THEN BEGIN
;
;  When switching to a higher window magnification factor, check to
;  see if any big images are stored in the data-environment.  If any
;  image is too big to fit on the screen, then do not use the selected
;  window magnification factor, and inform the user that it won't be
;  used because an image would be too big to fit on the screen.
;  -------------------------------------------------------------------
    m7 = SIZE(map7)
    m8 = SIZE(map8)
    m9 = SIZE(map9)
    f8 = SIZE(face8)
    f9 = SIZE(face9)
    CASE (zsel - 1) OF
      1: IF(m9(2) EQ 8) THEN GOTO, too_big_zoom
      2: BEGIN
           IF(m9(2) EQ 8) THEN GOTO, too_big_zoom
           IF(m8(2) EQ 8) THEN GOTO, too_big_zoom
           IF(f9(2) EQ 8) THEN GOTO, too_big_zoom
         END
      3: BEGIN
           IF(m9(2) EQ 8) THEN GOTO, too_big_zoom
           IF(m8(2) EQ 8) THEN GOTO, too_big_zoom
           IF(m7(2) EQ 8) THEN GOTO, too_big_zoom
           IF(f9(2) EQ 8) THEN GOTO, too_big_zoom
           IF(f8(2) EQ 8) THEN GOTO, too_big_zoom
         END
      ELSE:
    ENDCASE
    GOTO, ok_zoom
too_big_zoom:
    PRINT, 'Magnification request denied.  At least one image  ' + $
      'would not be able to'
    PRINT, 'fit on the screen.'
    GOTO, show_menu
ok_zoom:
;
;  Look at the contents of the array BLOCK_USAGE.  This array tells what
;  sections of the screen are being used.
;  ---------------------------------------------------------------------
    sz = SIZE(block_usage)
;    IF(sz(0) EQ 2) THEN nlevels = 1 ELSE nlevels = sz(3)
;    IF(nlevels NE 1) THEN GOTO, problem
    dim1 = sz(1)
    dim2 = sz(2)
;
;  Set MAXLEN_X and MAXLEN_Y to be the largest number of contiguous
;  blocks used by any object in the X and Y dimensions respectively.
;  -----------------------------------------------------------------
    maxlen_x = 0
    maxlen_y = 0
    FOR j = 0, dim2 - 1 DO BEGIN
      FOR k = 1, 50 DO BEGIN
        temp = WHERE(block_usage(*, j) EQ k)
        IF(temp(0) NE -1) THEN BEGIN
          sz = SIZE(temp)
          len_y = sz(1)
          IF(len_y GT maxlen_y) THEN maxlen_y = len_y
        ENDIF
      ENDFOR
    ENDFOR
    FOR i = 0, dim1 - 1 DO BEGIN
      temp = WHERE(block_usage(i, *) NE 0)
      sz = SIZE(temp)
      len_x = sz(1)
      IF(len_x GT maxlen_x) THEN maxlen_x = len_x
    ENDFOR
;    diff = (zsel - 1) - zoom_index
;    factor = 2^diff
;    newlen_x = maxlen_x * factor
;    newlen_y = maxlen_y * factor
;    IF((newlen_x LE dim1) AND  (newlen_y LE dim2)) THEN GOTO, doitanyway
    IF((maxlen_x LE dim1) AND  (maxlen_y LE dim2)) THEN GOTO, doitanyway
problem:
;
;  Print out a message asking the user if he still wants to use the
;  selected window magnification factor, eventhough there probably
;  will be some problems.
;  ----------------------------------------------------------------
    PRINT, 'Switching to the selected magnification factor may cause', $
      ' some images/graphs to overlap'
    PRINT, 'or perhaps be too big to fit on the screen.'
    PRINT, ' '
    PRINT, 'You may wish to remove some objects prior to', $
      ' switching to this magnification factor,'
    PRINT, 'or else remain at the present magnification factor.'
    menu = ['Switch to requested magnification factor anyway?', 'Yes', 'No']
    sel = one_line_menu(menu, init = 2)
    IF(sel EQ 2) THEN GOTO, show_menu
  ENDIF
doitanyway:
;
;  Assign the selected value to ZOOM_INDEX, and then re-display all
;  data-objects.
;  ----------------------------------------------------------------
  zoom_index = zsel - 1
  xdisplay_all
  IF(journal_on) THEN BEGIN
    IF(first EQ 1) THEN PRINTF, luj, menu_title
    PRINTF, luj, '  magnification factor was set to ' + $
      STRTRIM(STRING(2^zoom_index), 2)
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


