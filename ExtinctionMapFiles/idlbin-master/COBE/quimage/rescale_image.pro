PRO rescale_image, name
;+
;  RESCALE_IMAGE - a UIMAGE-specific routine.  This routine allows the
;  user to scale an image within a specified range.
;#
;  Written by John A. Ewing, ARC, January 1992.
;  SPR 10398  Jan 11 93  Allow an easy set back to default values.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 04 93  Change journaling.  J Ewing
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
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    RETURN
  ENDIF
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = 'Stretch contrast'
  IF(N_PARAMS(0) GE 1) THEN GOTO, have_name
show_menu:
;
;  Put up a menu in which the user can identify the desired image.
;  ---------------------------------------------------------------
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
           /face8,/face9,/zoomed,/proj_map,/proj2_map,/help,$
           /exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no objects to rescale.', /CONT
    RETURN
  ENDIF
have_name:
;
;  Extract pertinent info from the appropriate structure.
;  ------------------------------------------------------
  j = EXECUTE('title = ' + name + '.title')
  PRINT, ' '
  PRINT, 'Selected image:  ' + bold(title)
  j = EXECUTE('scale_min = ' + name + '.scale_min')
  j = EXECUTE('scale_max = ' + name + '.scale_max')
  j = EXECUTE('window = ' + name + '.window')
  IF(STRMID(name, 0, 6) NE 'ZOOMED') THEN BEGIN
    j = EXECUTE('data_min = ' + name + '.data_min')
    j = EXECUTE('data_max = ' + name + '.data_max')
    IF(data_min GT data_max) THEN BEGIN
      PRINT, 'Image consists of all bad pixels.'
      GOTO, show_menu
    ENDIF
    IF(data_min EQ data_max) THEN BEGIN
      PRINT, 'Image has no dynamic range.'
      GOTO, show_menu
    ENDIF
  ENDIF ELSE BEGIN
    j = EXECUTE('start_x = ' + name + '.start_x')
    j = EXECUTE('start_y = ' + name + '.start_y')
    j = EXECUTE('stop_x = ' + name + '.stop_x')
    j = EXECUTE('stop_y = ' + name + '.stop_y')
    j = EXECUTE('win_orig = ' + name + '.win_orig')
    bflag=0
    IF (win_orig LT 0) THEN BEGIN
      win_orig=ABS(win_orig)
      bflag=1
    ENDIF
    j=EXECUTE('zoomflag = ' + name + '.zoomflag')
    IF (zoomflag NE 0) THEN bflag=1
    name_orig = get_name(win_orig)
    j = EXECUTE('data = ' + name_orig + '.data')
    IF (bflag EQ 1) THEN BEGIN
      IF (zoomflag EQ 1) THEN data(start_x:stop_x,start_y:stop_y)=zbgr
      IF (zoomflag EQ 2) THEN data(start_x:stop_x,start_y:stop_y)=zdsrcmap
      IF (zoomflag EQ 3) THEN data(start_x:stop_x,start_y:stop_y)=zbsub
    ENDIF
    sz = SIZE(data)
    dim1_orig = sz(1)
    dim2_orig = sz(2)
    x0 = start_x > 0
    x1 = stop_x < (dim1_orig-1)
    y0 = start_y > 0
    y1 = stop_y < (dim2_orig-1)
    data = data(x0:x1, y0:y1)
    j = EXECUTE('badpixval = ' + name_orig + '.badpixval')
    good = WHERE(data NE badpixval)
    IF(good(0) NE -1) THEN BEGIN
      data_min = MIN(data(good), MAX = data_max)
    ENDIF ELSE BEGIN
      PRINT, 'Image consists of all bad pixels.  Operation not performed.'
      GOTO, show_menu
    ENDELSE
  ENDELSE
;
;  Have the user enter the new scaling range.
;  ------------------------------------------
  PRINT, 'The valid data min and max are:  ', STRTRIM(STRING(data_min), 2), $
    ',  ', STRTRIM(STRING(data_max), 2)
  print, 'The current scaling min and max are:  ', STRTRIM(STRING(scale_min), $
    2), ',  ', STRTRIM(STRING(scale_max), 2)
  IF((scale_min NE data_min) or (scale_max NE data_max)) THEN BEGIN
    menu = ['Set the scaling min & max to the full range of the data?', $
      'Yes', 'No']
    sel = one_line_menu(menu)
    IF(sel EQ 1) THEN BEGIN
      scale_min = data_min
      scale_max = data_max
      GOTO, set_new_values
    ENDIF
  ENDIF
  str = ' '
  val = 0
  WHILE(val EQ 0) DO BEGIN
    READ, 'Enter new scaling min:  ', str
    IF(str EQ '') THEN val = 1 ELSE BEGIN $
      val = validnum(str)
      IF(val EQ 1) THEN scale_min = FLOAT(str) $
                   ELSE PRINT, 'Invalid number, please re-enter.'
    ENDELSE
  ENDWHILE
  val = 0
  WHILE(val EQ 0) DO BEGIN
    READ, 'Enter new scaling max:  ', str
    IF(str EQ '') THEN val = 1 ELSE BEGIN $
      val = validnum(str)
      IF(val EQ 1) THEN scale_max = FLOAT(str) $
                   ELSE PRINT, 'Invalid number, please re-enter.'
    ENDELSE
  ENDWHILE
  IF(scale_min GT scale_max) THEN BEGIN
    PRINT, 'Scaling min must be less than scaling max.  No changes made.'
    GOTO, show_menu
  ENDIF
;
;  Update the structure to contain the new scaling range & redisplay the image.
;  ----------------------------------------------------------------------------
set_new_values:
  j = EXECUTE(name + '.scale_min = scale_min')
  j = EXECUTE(name + '.scale_max = scale_max')
  IF(STRMID(name, 0, 6) NE 'ZOOMED') THEN xdisplay, name, /rescale $
                                     ELSE show_zoomed, name
;
;  If journaling is enabled, then send some info to the journal file.
;  ------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ', title
    PRINTF, luj, '  new scale_min = ', STRTRIM(scale_min, 2)
    PRINTF, luj, '  new scale_max = ', STRTRIM(scale_max, 2)
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


