PRO udispvt, name
;+NAME/ONE LINE DESCRIPTION:
;    Displays images on non-X window terminals.
;
;DESCRIPTION:
;    This routine calls VTGREY with inputs inorder to display
;    skymaps on non-x window terminals.
;
;CALLING SEQUENCE:
;    udispvt
;
;ARGUMENTS (I=input, o=output, []=optional): 
;    None
;
;WARNINGS:
;    None
;
;EXAMPLE:
;    UDISPVT
;    The invocation is as follows,
;    IDL> UDISPVT
;#
;COMMON BLOCKS:  journal, uimage_data, uimage_data2, xwindow
;
;PROCEDURE:
;    This routine will be called from UIMAGE.
;
;REVISION HISTORY:
;    Written by 
;    S.R.K. VIDYA SAGAR (ARC) MAR 92
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;
;.TITLE
;Routine UDISPVT
;-
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  COMMON zback,zbgr,zdsrcmap,zbsub
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = 'Tek or VT200 grey scale plot'
  IF(NOT KEYWORD_SET(name)) THEN BEGIN
;
;  Have the user identify which image he wants to work with.
;  ---------------------------------------------------------
show_menu:
    name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
      /face8,/face9,/proj_map,/proj2_map,/zoomed,/help,/exit,title=menu_title)
    IF(name EQ 'HELP') THEN BEGIN
      UIMAGE_HELP, menu_title
      GOTO, show_menu
    ENDIF
    IF(name EQ 'EXIT') THEN RETURN
    IF(name EQ 'NO_OBJECTS') THEN BEGIN
      MESSAGE, 'No images are currently available.', /CONTINUE
      RETURN
    ENDIF
  ENDIF
  j = EXECUTE('title = ' + name + '.title')
  PRINT, ' '
  PRINT, 'Selected image:  ' + bold(title)
  j = EXECUTE('scale_min = ' + name + '.scale_min')
  j = EXECUTE('scale_max = ' + name + '.scale_max')
  j = EXECUTE('window = ' + name + '.window')
  IF(STRMID(name, 0, 6) NE 'ZOOMED') THEN BEGIN
;
;  Set up for MAP6,MAP7,MAP8,MAP9,FACE6,FACE7,FACE8,FACE9,PROJ_MAP.
;  ----------------------------------------------------------------
    j = EXECUTE('data = ' + name + '.data')
    j = EXECUTE('badpixval = ' + name + '.badpixval')
    good = WHERE(data NE badpixval)
    j = EXECUTE('pos_x = ' + name + '.pos_x')
    j = EXECUTE('pos_y = ' + name + '.pos_y')
    j = EXECUTE('data_min = ' + name + '.data_min')
    j = EXECUTE('data_max = ' + name + '.data_max')
    IF(data_min GT data_max) THEN BEGIN
      PRINT, 'Image consists of all bad pixels.'
      IF(uimage_version EQ 2) THEN GOTO, show_menu ELSE RETURN
    ENDIF
    IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN zoom = 2 ^ zoom_index
  ENDIF ELSE BEGIN
;
;  Set up for a ZOOMED.
;  --------------------
    j = EXECUTE('start_x = ' + name + '.start_x')
    j = EXECUTE('start_y = ' + name + '.start_y')
    j = EXECUTE('stop_x = ' + name + '.stop_x')
    j = EXECUTE('stop_y = ' + name + '.stop_y')
    j = EXECUTE('win_orig = ' + name + '.win_orig')
    j = EXECUTE('specific_zoom = ' + name + '.specific_zoom')
    zoom = specific_zoom * 2 ^ zoom_index
    IF(start_x GE 0) THEN pos_x = 0 ELSE pos_x = (-start_x) * zoom
    IF(start_y GE 0) THEN pos_y = 0 ELSE pos_y = (-start_y) * zoom
    pos_x = pos_x + zoom / 2
    pos_y = pos_y + zoom / 2
    bflag=0
    IF (win_orig LT 0) THEN BEGIN
     win_orig=ABS(win_orig)
     bflag=1
    ENDIF
    j=EXECUTE('zoomflag = ' + name + '.zoomflag')
    IF (zoomflag NE 0) THEN bflag=1
    name_orig = get_name(win_orig)
    j = EXECUTE('data = ' + name_orig + '.data')
    sz = SIZE(data)
    dim1_orig = sz(1)
    dim2_orig = sz(2)
    x0 = start_x > 0
    x1 = stop_x < (dim1_orig - 1)
    y0 = start_y > 0
    y1 = stop_y < (dim2_orig - 1)
    data = data(x0:x1, y0:y1)
    IF (bflag EQ 1) THEN BEGIN
      IF (zoomflag EQ 1) THEN data=zbgr
      IF (zoomflag EQ 2) THEN data=zdsrcmap
      IF (zoomflag EQ 3) THEN data=zbsub
    ENDIF
    j = EXECUTE('badpixval = ' + name_orig + '.badpixval')
    good = WHERE(data NE badpixval)
    IF(good(0) NE -1) THEN BEGIN
      data_min = MIN(data(good), MAX = data_max)
    ENDIF ELSE BEGIN
      PRINT, 'Image consists of all bad pixels.  Operation not performed.'
      IF(uimage_version EQ 2) THEN GOTO, show_menu ELSE RETURN
    ENDELSE
  ENDELSE
  IF(data_min EQ data_max) THEN BEGIN
    PRINT, 'Image has no dynamic range.  Operation not performed.'
    IF(uimage_version EQ 2) THEN GOTO, show_menu ELSE RETURN
  ENDIF
  sz = SIZE(data)
  dim1 = sz(1)
  dim2 = sz(2)
;
;  Define TERMX and TERMY for the call to VTGREY.
;  ----------------------------------------------
  CASE !D.NAME OF
    'TEK'   : BEGIN & termx = 429 & termy = 326 & END
    'REGIS' : BEGIN & termx = 319 & termy = 186 & END
    'X'     : BEGIN
                termx = dim1
                termy = dim2
              END
    'WIN'   : BEGIN
                termx = dim1
                termy = dim2
              END
    ELSE    : BEGIN
                MESSAGE, 'This function is accessible only under TEK ' + $
                  'and REGIS.', /CONT
                RETURN
              END
  ENDCASE
;
;  Reduce array resolution if it's too big to fit on the screen.
;  -------------------------------------------------------------
  IF((dim1 GT termx) OR (dim2 GT termy)) THEN BEGIN
    facx = FIX(dim1 / termx)
    IF(facx * termx NE dim1) THEN facx = facx + 1
    facy = FIX(dim1 / termy)
    IF(facy * termy NE dim1) THEN facy = facy + 1
    factor = MAX([facx, facy])
    PRINT, 'The image will be temporarily reduced in size by a factor of ' + $
      bold(STRTRIM(STRING(factor), 2))
    PRINT, 'in order to be displayed.'
    data = CONGRID(data, dim1 / factor, dim2 / factor)
    good = WHERE(data NE badpixval)
  ENDIF
;
;  Give the user the opportunity to change the scaling min and max.
;  ----------------------------------------------------------------
  PRINT, 'The valid data min and max are:  ', STRTRIM(STRING(data_min), 2), $
    ',  ', STRTRIM(STRING(data_max), 2)
  PRINT, 'The current scaling min and max are:  ', $
    STRTRIM(STRING(scale_min), 2), ',  ', STRTRIM(STRING(scale_max), 2)
  menu = ['Do you want to change the scaling range?', 'Yes', 'No']
  sel = one_line_menu(menu, init = 2)
  IF(sel EQ 1) THEN BEGIN
    str = ' '
    val = 0
    WHILE(val EQ 0) DO BEGIN
      READ, 'Enter new scaling min:  ', str
      val = validnum(str)
      IF(val EQ 1) THEN scale_min = FLOAT(str) $
                   ELSE PRINT, 'Invalid number, please re-enter.'
    ENDWHILE
    val = 0
    WHILE(val EQ 0) DO BEGIN
      READ, 'Enter new scaling max:  ', str
      val = validnum(str)
      IF(val EQ 1) THEN scale_max = FLOAT(str) $
                   ELSE PRINT, 'Invalid number, please re-enter.'
    ENDWHILE
    IF(scale_min GT scale_max) THEN BEGIN
      PRINT, 'Scaling min must be less than scaling max.  No changes made.'
      IF(uimage_version EQ 2) THEN GOTO, show_menu ELSE RETURN
    ENDIF
    j = EXECUTE(name + '.scale_min = scale_min')
    j = EXECUTE(name + '.scale_max = scale_max')
  ENDIF
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
    WSET, window
    pos = [pos_x, pos_y, pos_x + zoom * dim1 + 1, pos_y + zoom * dim2 + 1]
  ENDIF
;
;  Make the call to VTGREY.
;  ------------------------
  scaled_data = FIX(data) * 0
  scaled_data(good) = FIX(BYTSCL(data(good), scale_min, scale_max))
  exstr = 'vtgrey,scaled_data,title=title,termx=termx,termy=termy,pos=pos'
  j = EXECUTE(exstr)
;
;  Special after-plot handling.
;  ----------------------------
  CASE !D.NAME OF
    'TEK'   : BEGIN
                c = GET_KBRD(1) 
                text
                PRINT, ' '
              END
    'REGIS' : BEGIN
                c = GET_KBRD(1) 
                ERASE
                PRINT, ' '
              END
    ELSE    : BEGIN
;
;  Indicate the resolution on the lower left corner of the X-window.
;  -----------------------------------------------------------------
                pos1 = STRPOS(name, '(')
                aos = STRMID(name, 0, pos1)
                res = 'RES' + STRMID(aos, STRLEN(aos) - 1, 1)
                CASE strupcase(aos) OF
                  'FACE6' : image_type = 0
                  'FACE7' : image_type = 1
                  'FACE8' : image_type = 2
                  'FACE9' : image_type = 3
                  'MAP6'  : image_type = 4
                  'MAP7'  : image_type = 5
                  'MAP8'  : image_type = 6
                  'MAP9'  : image_type = 7
                  ELSE    : GOTO, update_journal
                ENDCASE
                x = dim1 * zoom * .05
                IF((x/2) * 2 NE x) THEN x = x + 1
                y = (dim2 * zoom * .05) < 18
                IF((y/2) * 2 NE y) THEN y = y + 1
                char_size = [1.2,1.2,1.6,1.8,1.3,1.4,1.9,3.2]
                index = (image_type + zoom_index) < 7
                IF((STRUPCASE(STRMID(name, 0, 8)) NE 'PROJ_MAP') OR $
                   (STRUPCASE(STRMID(name, 0, 9)) NE 'PROJ2_MAP')) THEN $
                  XYOUTS, x, y, res, SIZE = char_size(index), /DEV, COLOR = 1
              END
  ENDCASE
update_journal:
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:    ', title
    PRINTF, luj, '  scale_min = ', STRTRIM(scale_min, 2)
    PRINTF, luj, '  scale_max = ', STRTRIM(scale_max, 2)
    PRINTF, luj, '  the following statement was executed:'
    PRINTF, luj, '  ' + exstr
    PRINTF, luj, '  where:'
    PRINTF, luj, '    scaled_data is byte-scaled data from the operand'
    PRINTF, luj, '    title = ' + title
    PRINTF, luj, '    termx = ' + STRTRIM(termx, 2)
    PRINTF, luj, '    termy = ' + STRTRIM(termy, 2)
    IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN $
      PRINTF, luj, '    pos   = ' , pos
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


