FUNCTION select_cval, dummy
;+
;  SELECT_CVAL - a UIMAGE-specific routine.  This routine puts up a menu
;  so a user can select a color value (ranging from 1 to C_SCALEMIN-1).
;  This routine is similar to BIGWMENU.  The color associated with each
;  color value is shown on the menu.
;#
;  Written by John Ewing, Nov 92.
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;--------------------------------------------------------------------------
  COMMON color_values,c_badpix,c_draw,c_scalemin
  COMMON colors,r_orig,g_orig,b_orig,r_curr,g_curr,b_curr
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE,'This routine can only be run on an X-window terminal.',/CONT
  ENDIF
  locolor = 7
  hicolor = 8
  arr = BYTARR(80, 16)
  old_window = !D.WINDOW
  height = 30
  csize = 1.6
  xsize = 200
  num_options = c_scalemin + 2
;
;  See to it that the menu gets placed at the center of the screen.
;  ----------------------------------------------------------------
  DEVICE, GET_SCREEN_SIZE = scrsiz
  xpos = (scrsiz(0) - 200) / 2
  ypos = (scrsiz(1) - (16 + height * num_options)) / 2
;
;  Put up the menu.
;  ----------------
  ysize = num_options * height
  menu_title = 'Select a color value'
make_menu:
  WINDOW, /FREE, XSIZE = 200, YSIZE = ysize, XPOS = xpos, $
    YPOS = ypos, TITLE = menu_title
  options = STRARR(num_options)
  FOR i = 0, c_scalemin - 2 DO options(i) = STRTRIM(STRING(i + 1), 2)
  options(c_scalemin - 1) = ''
  options(c_scalemin)     = 'HELP'
  options(c_scalemin + 1) = 'Exit this menu'
  FOR i = 0, c_scalemin - 2 DO BEGIN
    XYOUTS, 30, ysize - (i + 1) * height + 10, options(i), SIZE = csize, /DEVICE
    TV, arr + (i + 1), 60, ysize - (i + 1) * height + 8
  ENDFOR
  FOR i = c_scalemin, c_scalemin + 1 DO $
    XYOUTS, 10, ysize - (i + 1) * height + 10, options(i), SIZE = csize, /DEVICE
  xx = [2, 2, 2, xsize - 2, xsize - 2, xsize - 2, xsize - 2, 2]
  sel = 0
  osel = 0
  WHILE(options(osel) EQ '') DO osel = osel + 1
  ybox = INTARR(8, num_options)
  FOR i = 0, num_options - 1 DO BEGIN
    y1 = (num_options - i) * height - (height - 2)
    y2 = (num_options - i) * height - 2
    ybox(0, i) = [y1, y2, y2, y2, y2, y1, y1, y1]
  ENDFOR
;
;  Monitor mouse movements.
;  ------------------------
  first = 1
  !ERR = 0
  WHILE(!ERR EQ 0) DO BEGIN
    CURSOR, x, y, 0, /DEVICE
    IF((y NE -1) OR (first EQ 1)) THEN BEGIN
      sel = num_options - FIX(y / height) - 1
      IF(options(sel) EQ '') THEN sel = osel
      PLOTS, xx, ybox(*, sel), /DEVICE, THICK = 2.5, COLOR = hicolor
      FOR i = 0, num_options - 1 DO IF((i NE sel) AND (options(i) NE '')) THEN $
        PLOTS, xx, ybox(*, i), /DEVICE, THICK = 2.5, COLOR = locolor
    ENDIF
    first = 0
    osel = sel
  ENDWHILE
  WDELETE
  IF(options(sel) EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, make_menu
  ENDIF
  IF(old_window NE -1) THEN WSET, old_window
  IF(sel GE c_scalemin) THEN sel = -10
  RETURN, sel + 1
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


