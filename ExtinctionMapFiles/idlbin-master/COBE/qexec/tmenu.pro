FUNCTION tmenu,strings,title=title,init=init,xpos=xpos,ypos=ypos,valid=valid,$
  center=center,retain=retain,ysize=ysize,exit_on_refresh=exit_on_refresh,$
  just_reg=just_reg,tababort=tababort,demowait=demowait,xsize=xsize
;+NAME
;	TMENU
; PURPOSE:
;       This function displays a menu whose choices are given by the
;       elements of a string array and then returns the index of that
;       choice which the user selects.  This routine is designed to
;       correspond to the intrinsic WMENU routine, only in this case the
;       menu is displayed in a TEKTRONIX environment.
; CATEGORY:
;	User interface, Menu.
; CALLING SEQUENCE:
;	i=TMENU(['Select one','Item_1','Item_2'],title=0,init=1,xpos=10)
; INPUTS:
;       STRINGS  = A string array with each element containing the text
;       of one menu choice.  Each menu choice can consist of up to 76
;       characters, and there can be up to 20 menu choices.
;       TITLE    = The index of the STRINGS element that is the title,
;       normally 0.  The title element is not selectable and is displayed
;       bolded and underlined.  If this parameter is omitted, all items
;       are selectable.
;       INIT     = The index of the initial selection.  If this parameter
;       is specified and within the range of STRINGS indices (excluding
;       the TITLE index), then the initial menu display is made with the
;       designated item selected.
;       XPOS     = The X position (in units of characters, there being 80
;       characters per line) at which the menu should begin.  There are
;       some limitations on the size of XPOS due to the fact that IDL will
;       not print out more that 80 characters per line.  Basically, the
;       wider the menu, the smaller the maximum allowable value of XPOS.
;       If this parameter is omitted, the menu will border on the left
;       edge of the screen.
;       VALID    = List of 0's and 1's, which designate valid choices.
;       (optional)  If supplied, then invalid items will not be
;       selectable.
;       CENTER   = 1 if options should be centered.
;       RETAIN   = 1 if menu should remain on the screen after exiting.
;       XSIZE    = Number of columns to be contained within the menu.
;       YSIZE    = Number of rows to be contained within the menu,
;       should not exceed 20.
;       EXIT_ON_REFRESH = If set to 0 then a cntl-W will refresh the menu.
;       If set to 1, then TMENU will exit when cntl-W is pressed.
;       JUST_REG = If set to 1, then TMENU will just draw the menu and exit.
;       DEMOWAIT = If non-zero, then this is interpreted as the number of
;       seconds TMENU should wait before exiting, at which point it returns
;       the value given for INIT (for demo purposes only).
; OUTPUTS:
;       The index number of the selected menu item.  A value of -1 is
;       returned if this function is called incorrectly.
;#
; COMMON BLOCKS:
;	None.
; RESTRICTIONS:
;       Menu choices must be less than 77 characters.
;       There must be less than 21 menu choices.
; MODIFICATION HISTORY:
;  Creation:  John Ewing, ARC, November 1991.
;  Allow screen refresh:  John Ewing, ARC, January 1992.
;  Incorporate VALID keyword:  John Ewing, ARC, January 1992.
;  Erase menu with display commands: K Turpie, GSC, March 1992.
;  String options can be as big as 76 characters: J Ewing, March 1992.
;  Loop from top to bottom & vice versa:  John Ewing, June 1992.
;  Skip over blank lines:  John Ewing, June 1992.
;  Changed to run on ULTRIX also:  John Ewing, July 1992.
;  Scrolling capability:  John Ewing, July 1992.
;  JUST_REG & TABABORT keywords:  John Ewing, July 1992.
;  SPR 10206  Nov    92  Repair scrolling capability.  J Ewing.
;  SPR 10509  Jan 25 93  Work in 8-bit mode, eliminate GET_KBRD(0)s.  J Ewing
;  SPR 11003  Jun 02 93  Changed !version.os to !cgis_os. J Newmark
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;----------------------------------------------------------------------------
;
;  Define escape sequences for highlighting, bolding, line drawing, etc.
;  ---------------------------------------------------------------------
;
  esc   = STRING(27b)
  csi   = STRING(155b)
  IF((!cgis_os EQ 'vms') OR (!cgis_os EQ 'windows')) THEN cr = STRING(13b) $
    ELSE cr = STRING(10b)
  cntlw = STRING(23b)
  tab   = STRING(9b)
  bel   = STRING(7b)
  wb = esc + '[1m' & wl = esc + '[4m' & wr = esc + '[7m' & wz = esc + '[m'
  cs = esc + '(0'  & cz = esc + '(B'
  cub = esc + '[D' & cud = esc + '[B' & cuf = esc + '[C' & cuu = esc + '[A'
  blank = '                                        ' + $
          '                                        '
;
;  Check the validity of what was supplied for the "STRINGS" argument.
;  -------------------------------------------------------------------
  sz = SIZE(strings)
  IF(sz(0) NE 1) OR (sz(2) NE 7) THEN BEGIN
    MESSAGE, 'An inappropriate argument was supplied.', /CONT
    RETURN, -1
  ENDIF
;
;  Get the row of the current cursor position.
;  -------------------------------------------
get_cursor_pos:
  str = ''
  PRINT, esc + '[A' + esc + '[6n'
  char = ''
  dum = GET_KBRD(1)
  WHILE(char NE 'R') DO BEGIN
    char = GET_KBRD(1)
    str = str + char
  ENDWHILE
  ic1 = STRPOS(str, '[')
  ic2 = STRPOS(str, ';')
  ic3 = STRPOS(str, 'R')
  IF(ic3 EQ -1) THEN GOTO, get_cursor_pos
  str_row = STRMID(str, ic1 + 1, ic2 - ic1 - 1)
  FOR i = 0, STRLEN(str_row) - 1 DO BEGIN
    bt = BYTE(STRMID(str_row, i, 1))
    IF((bt(0) LT 48) OR (bt(0) GT 57)) THEN GOTO, get_cursor_pos
  ENDFOR
  row_curr = FIX(str_row) + 1
;
;  Determine YSIZE.
;  ----------------
  dim1 = sz(1)
  IF(KEYWORD_SET(ysize)) THEN BEGIN
    IF(ysize LT 1) THEN BEGIN
      MESSAGE, 'The supplied value for YSIZE was too small.', /CONT
      RETURN, -1
    ENDIF
    IF(ysize GT dim1) THEN ysize = dim1
    IF(ysize GT 20) THEN ysize = 20
  ENDIF ELSE BEGIN
    IF(dim1 GT 20) THEN ysize = 20 ELSE ysize = dim1
  ENDELSE
;
;  Check what, if anything, has been passed in for "VALID"
;  -------------------------------------------------------
   IF(KEYWORD_SET(valid)) THEN BEGIN
     sz = SIZE(valid)
     IF((sz(0) NE 1) OR (sz(1) NE dim1)) THEN BEGIN
       MESSAGE, 'An invalid argument has been passed in for VALID.', /CONT
       RETURN, -1
     ENDIF
     w = WHERE(valid EQ 1)
     IF(w(0) EQ -1) THEN BEGIN
       MESSAGE, 'Some of the VALID entries must be set to 1.', /CONT
       RETURN, -1
     ENDIF
   ENDIF ELSE BEGIN
     valid = INTARR(dim1) + 1
   ENDELSE
;
;  Determine the length of the longest menu choice.
;  ------------------------------------------------
  maxlen = 0
  FOR i = 0, dim1 - 1 DO BEGIN
    len = STRLEN(strings(i))
    IF(len GT maxlen) THEN maxlen = len
    IF(maxlen GT 76) THEN BEGIN
        maxlen=76
        MESSAGE, 'The length of menu choices can not exceed 76 characters.', $
        /CONT
        strings(i)=strmid(strings(i),0,75)
    ENDIF
  ENDFOR
  IF(KEYWORD_SET(xsize)) THEN BEGIN
    IF(xsize GT maxlen) THEN maxlen = xsize - 4
  ENDIF
  IF(ysize LT dim1) THEN maxlen = maxlen > 3
;
;  Define TITLE_STR and LIST.
;  --------------------------
  sz = SIZE(title)
  IF(sz(0) + sz(1) NE 0) THEN BEGIN
    IF((title LT 0) OR (title GE dim1)) THEN BEGIN
      MESSAGE, 'TITLE has a value outside the appropriate range.', /CONT
      RETURN, -1
    ENDIF
    IF(dim1 LE 1) THEN BEGIN
      MESSAGE, 'The dimension of STRINGS is too small.', /CONT
      RETURN, -1
    ENDIF
    title_str = strings(title)
    temp = STRMID(blank, 0, (maxlen - STRLEN(title_str)) / 2) + title_str
    title_str = temp + STRMID(blank, 0, (maxlen - STRLEN(temp)))
    nlist = dim1 - 1
    CASE title OF
               0 : index = INDGEN(nlist) + 1
      (dim1 - 1) : index = INDGEN(nlist)
            ELSE : index = [INDGEN(title), INDGEN(dim1 - title - 1) + title + 1]
    ENDCASE
    list = strings(index)
  ENDIF ELSE BEGIN
    title_str = ''
    nlist = dim1
    list = strings
    index = INDGEN(dim1)
  ENDELSE
  temp = STRMID(blank, 0, (maxlen - 3) / 2) + '...'
  dotdotdot = ' ' + temp + STRMID(blank, 0, (maxlen - STRLEN(temp))) + ' '
  val = valid(index)
;
;  Make sure that there is at least one good (selectable) option.
;  --------------------------------------------------------------
  a_good_option = 0
  FOR i = 0, nlist - 1 DO $
    IF((val(i) EQ 1) AND (STRTRIM(list(i), 2) NE '')) THEN a_good_option = 1
  IF(NOT a_good_option) THEN BEGIN
    MESSAGE, 'No selectable options were supplied.', /CONT
    RETURN, -1
  ENDIF
;
;  Define increment and decrement functions.
;  -----------------------------------------
  increment = INTARR(nlist) - 1
  decrement = INTARR(nlist) - 1
  FOR i = 0, nlist - 1 DO BEGIN
    IF((val(i) EQ 1) AND (STRTRIM(list(i), 2) NE '')) THEN BEGIN
      ic = i + 1
      IF(ic EQ nlist) THEN ic = 0
      WHILE((val(ic) NE 1) OR (STRTRIM(list(ic), 2) EQ '')) DO BEGIN
        ic = ic + 1
        IF(ic EQ nlist) THEN ic = 0
      ENDWHILE
      increment(i) = ic
      ic = i - 1
      IF(ic EQ -1) THEN ic = nlist - 1
      WHILE((val(ic) NE 1) OR (STRTRIM(list(ic), 2) EQ '')) DO BEGIN
        ic = ic - 1
        IF(ic EQ -1) THEN ic = nlist - 1
      ENDWHILE
      decrement(i) = ic
    ENDIF
  ENDFOR
;
;  Pad option-strings with spaces.
;  -------------------------------
  IF(KEYWORD_SET(center)) THEN FOR i = 0, nlist - 1 DO $
    list(i) = STRMID(blank, 0, (maxlen - STRLEN(list(i))) / 2) + list(i)
  FOR i = 0, nlist - 1 DO $
    list(i) = list(i) + STRMID(blank, 0, (maxlen - STRLEN(list(i))))
;
;  Use the value of "XPOS" (if supplied) to determine the amount of
;  initial spacing in the horizontal dimension.
;  ----------------------------------------------------------------
  IF(KEYWORD_SET(xpos)) THEN BEGIN
    IF(maxlen + xpos + 4 GT 80) THEN BEGIN
      MESSAGE, 'The supplied value for XPOS is too big.', /CONT
      xpos = 80 - (maxlen + 4)
    ENDIF
    col_str = STRTRIM(STRING((xpos > 0) < (78 - maxlen)), 2)
  ENDIF ELSE col_str = '0'
;
;  Use the value of "YPOS" (if supplied) to determine the row of the
;  top bar of the menu.
;  ----------------------------------------------------------------
  IF(KEYWORD_SET(ypos)) THEN row = (ypos > 0) < 23 ELSE row = row_curr
;
;  Determine which menu item should be selected by default.
;  --------------------------------------------------------
  sz = SIZE(init)
  IF(sz(0) + sz(1) NE 0) THEN BEGIN
    sel = init
    IF(title_str NE '') THEN IF(init GT title) THEN sel = sel - 1
    sel = (sel > 0) < (nlist - 1)
  ENDIF ELSE sel = 0
  WHILE((val(sel) EQ 0) OR (STRTRIM(list(sel), 2) EQ '')) DO BEGIN
    sel = sel + 1
    IF(sel EQ nlist) THEN sel = 0
  ENDWHILE
  noptshow = ysize
  IF(title_str NE '') THEN noptshow = noptshow - 1
;
;  Determine the screen number and row number for each menu-item.
;  --------------------------------------------------------------
  IF(title_str NE '') THEN rowoffset = 1 ELSE rowoffset = 0
  scrarr = INTARR(nlist)
  rowarr = INTARR(nlist)
  layout = INTARR(noptshow, nlist) - 9
  is = 0
  ic = 0
  FOR i = 0, nlist - 1 DO BEGIN
    scrarr(i) = is
    rowarr(i) = row + rowoffset + 1 + ic
    layout(ic, is) = i
    ic = ic + 1
    IF((ic EQ noptshow) AND (i NE nlist - 1)) THEN BEGIN
      layout(ic - 1, is) = -1
      is = is + 1
      layout(0, is) = -1
      scrarr(i) = is
      rowarr(i) = row + rowoffset + 2
      layout(1, is) = i
      ic = 2
    ENDIF
  ENDFOR
  nscreens = is + 1
;
;  Print the menu out on the terminal screen.
;  ---------------------------------------
show_menu:
  screen = scrarr(sel)
  line = 'q' & WHILE(STRLEN(line) LT (maxlen + 2)) DO line = line + 'q'
  trow = row
  PRINT, esc + '[' + STRTRIM(STRING(trow), 2) + ';' + col_str + 'f' + $
    cs + 'l' + line + 'k' + cz, FORMAT = '(a)'
  IF(title_str NE '') THEN BEGIN
    trow = trow + 1
    PRINT, esc + '[' + STRTRIM(STRING(trow), 2) + ';' + col_str + 'f' + $
      cs + 'x' + cz + wb + wl + ' ' + title_str + ' ' + wz + cs + 'x' + cz, $
      FORMAT = '(a)'
  ENDIF
  FOR i = 0, noptshow - 1 DO BEGIN
    trow = row + rowoffset + 1 + i
    ic = layout(i, screen)
    CASE ic OF
       -1 : PRINT, esc + '[' + STRTRIM(STRING(trow), 2) + ';' + col_str + $
              'f' + cs + 'x' + cz + dotdotdot + cs + 'x' + cz, FORMAT = '(a)'
       -9 : PRINT, esc + '[' + STRTRIM(STRING(trow), 2) + ';' + col_str + $
              'f' + cs + 'x' + cz + STRMID(blank, 0, maxlen + 2) + cs + $
              'x' + cz, FORMAT = '(a)'
     ELSE : IF((ic EQ sel) AND NOT KEYWORD_SET(just_reg)) THEN BEGIN
              PRINT, esc + '[' + STRTRIM(STRING(trow), 2) + ';' + col_str + $
                'f' + cs + 'x' + cz + ' ' + wr + list(sel) + wz + ' ' + $
                cs + 'x' + cz, FORMAT = '(a)'
            ENDIF ELSE BEGIN
              PRINT, esc + '[' + STRTRIM(STRING(trow), 2) + ';' + col_str + $
                'f' + cs + 'x' + cz + ' ' + list(ic) + ' ' + cs + 'x' + $
                cz, FORMAT = '(a)'
            ENDELSE
    ENDCASE
  ENDFOR
  trow = row + ysize + 1
  IF(KEYWORD_SET(just_reg)) THEN str = '' ELSE str = esc + '[?25l'
  PRINT, esc + '[' + STRTRIM(STRING(trow), 2) + ';' + col_str + 'f' + $
    cs + 'm' + line + 'j' + cz + cuu + cuu + str, FORMAT = '(a)'
  IF(KEYWORD_SET(just_reg)) THEN RETURN, index(sel)
;
;  Correct for any scrolling which occured.
;  ----------------------------------------
  IF(trow GT 24) THEN BEGIN
    row = row - (trow - 24)
    rowarr = rowarr - (trow - 24)
  ENDIF
;
;  Just wait and exit if DEMOWAIT is set.
;  --------------------------------------
  IF(KEYWORD_SET(demowait)) THEN BEGIN
    WAIT, demowait
    GOTO, made_selection
  ENDIF
;
;  Monitor and respond to key presses until a carriage return is sent.
;  -------------------------------------------------------------------
watch_keyboard:
  WHILE 1 DO BEGIN
    CASE GET_KBRD(1) of
      cr   : GOTO, made_selection
      esc  : IF(GET_KBRD(1) EQ '[') THEN GOTO, handle_arrowkeys
      csi  : GOTO, handle_arrowkeys
      tab  : IF(KEYWORD_SET(tababort)) THEN BEGIN
               PRINT, esc + '[' + STRTRIM(STRING(rowarr(sel)), 2) + ';' + $
                 col_str + 'f' + cs + 'x' + cz + ' ' + list(sel) + ' ' + $
                 cs + 'x' + cz, FORMAT = '(a)'
               RETURN, index(sel) - 1000
             ENDIF ELSE GOTO, down_arrow
      cntlw: GOTO, refresh
      ELSE :
    ENDCASE
  ENDWHILE
refresh:
  IF(KEYWORD_SET(exit_on_refresh)) THEN BEGIN
    RETURN, - (index(sel) + 1)
  ENDIF ELSE BEGIN
    PRINT, esc + '[2J' + esc + '[f'
    GOTO, show_menu
  ENDELSE
handle_arrowkeys:
  CASE GET_KBRD(1) of
    'A' : BEGIN  ;up
            PRINT, esc + '[' + STRTRIM(STRING(rowarr(sel)), 2) + ';' + $
              col_str + 'f' + cs + 'x' + cz + ' ' + list(sel) + ' ' + $
              cs + 'x' + cz, FORMAT = '(a)'
            sel = decrement(sel)
            oscreen = screen
            screen = scrarr(sel)
            IF(screen NE oscreen) THEN GOTO, show_menu
            PRINT, esc + '[' + STRTRIM(STRING(rowarr(sel)), 2) + ';' + $
              col_str + 'f' + cs + 'x' + cz + ' ' + wr + list(sel) + wz + $
              ' ' + cs + 'x' + cz, FORMAT = '(a)'
          END
    'B' : BEGIN  ;down
down_arrow:
            PRINT, esc + '[' + STRTRIM(STRING(rowarr(sel)), 2) + ';' + $
              col_str + 'f' + cs + 'x' + cz + ' ' + list(sel) + ' ' + $
              cs + 'x' + cz, FORMAT = '(a)'
            sel = increment(sel)
            oscreen = screen
            screen = scrarr(sel)
            IF(screen NE oscreen) THEN GOTO, show_menu
            PRINT, esc + '[' + STRTRIM(STRING(rowarr(sel)), 2) + ';' + $
              col_str + 'f' + cs + 'x' + cz + ' ' + wr + list(sel) + wz + $
              ' ' + cs + 'x' + cz, FORMAT = '(a)'
          END
    '5' : IF(screen GT 0) THEN BEGIN  ; scroll up
            currow = rowarr(sel)
            pos = layout(*, screen - 1)
            pos = pos(WHERE(pos GE 0))
            w = WHERE(rowarr(pos) EQ currow)
            IF(w(0) NE -1) THEN sel = pos(w(0)) ELSE sel = pos(0)
            GOTO, show_menu
          ENDIF ELSE PRINT, esc + '[A' + bel
    '6' : IF(screen LT nscreens - 1) THEN BEGIN  ; scroll down
            currow = rowarr(sel)
            pos = layout(*, screen + 1)
            pos = pos(WHERE(pos GE 0))
            w = WHERE(rowarr(pos) EQ currow)
            IF(w(0) NE -1) THEN sel = pos(w(0)) ELSE sel = pos(0)
            GOTO, show_menu
          ENDIF ELSE PRINT, esc + '[A' + bel
    ELSE: GOTO, watch_keyboard
  ENDCASE
  GOTO, watch_keyboard
;
;  Erase menu, return the index of the selected choice.
;  ----------------------------------------------------
made_selection:
  PRINT, esc + '[A' + esc + '[?25h'
  IF(NOT KEYWORD_SET(retain)) THEN BEGIN
    FOR i = 0, ysize + 1 DO BEGIN
      trow = row + ysize + 1 - i
      PRINT, esc + '[' + STRTRIM(STRING(trow), 2) + ';' + col_str + 'f' + $
        esc + '[0K' + esc + '[A', FORMAT = '(a)'
    ENDFOR
    PRINT, esc + '[' + STRTRIM(STRING(row - 1), 2) + 'f'
  ENDIF ELSE BEGIN
    PRINT, esc + '[' + STRTRIM(STRING(rowarr(sel)), 2) + ';' + col_str + $
      'f' + cs + 'x' + cz + ' ' + wb + list(sel) + wz + ' ' + cs + 'x' + $
      cz, FORMAT = '(a)'
  ENDELSE
  RETURN, index(sel)
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


