FUNCTION pmenu, desc, demo = demo, ypos = ypos, noclear = noclear
;+
; NAME
;        PMENU
;
; PURPOSE:
;        This function displays a pull-down menu whose choices are given
;        by the elements of a supplied string array.  After the user makes
;        a selection (via the arrow and return keys), this function returns
;        the selected choice (a string).  This routine is designed to
;        correspond to RSI's XPDMENU, and the array-of-strings which needs
;        to be supplied for this routine is identical in format to an
;        argument for XPDMENU, with the following exception:  THIS VERSION
;        OF PMENU ALLOWS THE CREATION OF MENUS DOWN TO ONE LEVEL ONLY.
;
; CALLING SEQUENCE:
;
;        pmenu(desc)
;
; INPUTS:
;        DESC:    An array-of-strings which defines the contents of the
;                 pull-down menu.  An example is shown below:
;
;                 example = ['"Line Style" {', $
;                                '"Solid"', $
;                                '"Dashed"', $
;                                '"Dot-Dash"', $
;                                '}', $
;                            '"Color" {', $
;                                '"Red"', $
;                                '"Green"', $
;                                '"Blue"', $
;                                '}']
;
; KEYWORDS:
;        DEMO:    If set, then a default/demo pull-down menu is shown,
;                 regardless of whether or not anything was supplied
;                 for DESC.
;
;        YPOS:    A integer scalar that indicates at what row the top
;                 of the pull-down menu should be placed.
;
;        NOCLEAR: If set, then the screen will not be cleared prior to
;                 when the pull-down menu is displayed.
;
; OUTPUTS:
;        This function returns the text for the option which was selected.
;
; COMMON BLOCKS:
;        None.
;
; SIDE EFFECTS:
;        None.
;
; RESTRICTIONS:
;        Very little syntax checking is done on the array-of-strings
;        supplied for DESC.  Incorrectly formatted input can lead to
;        unexpected results.  The same is true for the keywords.
;        A screen width of 80 characters is assumed.
;
; EXAMPLE:
;       print, pmenu(/demo)
;
;       or
;
;       sample = ['"Category 1" {', '"option 1"', '"option 2"', '}', $
;         '"Category 2" {', '"option 3"', '"option 4"', '"option 5"', '}', $
;         '"Category 3" {', '"option 6"', '"option 7"', '}']
;       selected = pmenu(sample)
;
; MODIFICATION HISTORY:
;  Written by John Ewing (ARC), June 12 1992.
;  Prepared for delivery to the external community.  J Ewing, March 29 1993.
;  SPR 11003  Jun 02 93  Changed !version.os to !cgis_os. J Newmark
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;---------------------------------------------------------------------------
;
;  If the DEMO keyword was set, then use the following menu
;  --------------------------------------------------------
  IF(KEYWORD_SET(demo)) THEN BEGIN
    desc = [ $
      '"Breakfast"   {', '"Egg Biscuit"', '"Sausage & Egg Biscuit"', $
                         '"Croissant"', '"Pancakes"', '"Hash Browns"', '}', $
      '"Sandwiches"  {', '"Hamburger"', '"Cheeseburger"', '"Chicken"', $
                         '"Fish"', '"Roast Beef"', '"Turkey"', '}', $
      '"Side Orders" {', '"French Fries"', '"Onion Rings"', '}', $
      '"Desserts"    {', '"Apple Pie"', '"Brownie"', '"Cake"', '"Cookies"', $
                         '"Eclaire"', '"Ice Cream"', '}', $
      '"Drinks"      {', '"Soda"', '"Milk"', '"Orange Juice"', '"Coffee"', $
                         '"Iced Tea"', '}']
  ENDIF ELSE BEGIN
;
;  Do a partial check on the supplied array of strings.
;  ----------------------------------------------------
    sz = SIZE(desc)
    IF(sz(0) + sz(1) EQ 0) THEN $
      MESSAGE, 'One argument (an array of strings) is required.'
    IF((sz(0) NE 1) OR (sz(2) NE 7)) THEN $
      MESSAGE, 'The argument must be an array of strings.'
  ENDELSE
  IF(KEYWORD_SET(ypos)) THEN ypos = ypos > 1 ELSE ypos = 1
;
;  From the supplied array-of-strings, DESC, generate an array-of-strings
;  named WORD that contains all menu options, and an array-of-integers
;  named PARENT that, for each element of word, contains the index of the
;  parent of that element.
;  ----------------------------------------------------------------------
  list = ''
  sz = SIZE(desc)
  FOR i = 0, sz(1) - 1 DO list = list + desc(i)
  level = 0 & wordmode = 0 & cur_word = '' & cur_num = -1
  word = [''] & parent = [0]
  stack = INTARR(20) & stack(0) = -1
  FOR i = 0, STRLEN(list) - 1 DO BEGIN
    char = STRMID(list, i, 1)
    CASE char OF
      '"' : BEGIN
              wordmode = 1 - wordmode
              IF(wordmode EQ 0) THEN BEGIN
                word = [word, ' ' + cur_word + ' ']
                parent = [parent, stack(level)]
                cur_word = ''
              ENDIF ELSE cur_num = cur_num + 1
            END
      '{' : BEGIN
              level = level + 1
              stack(level) = cur_num
            END
      '}' : level = level - 1
      ELSE: IF(wordmode EQ 1) THEN cur_word = cur_word + char
    ENDCASE
  ENDFOR
  word = word(1:*)
  parent = parent(1:*)
  sz = SIZE(word)
  nword = sz(1)
;
;  Define special escape sequences
;  -------------------------------
  esc   = STRING(27b)
  csi   = STRING(155b)
  IF((!cgis_os EQ 'vms') OR (!cgis_os EQ 'windows')) THEN cr = STRING(13b) $
    ELSE cr = STRING(10b)
  cntlw = STRING(23b)
  tab   = STRING(9b)
  wb = esc + '[1m' & wl = esc + '[4m' & wr = esc + '[7m' & wz = esc + '[m'
  cs = esc + '(0'  & cz = esc + '(B'
  blank='                                                                '
  qstr = 'qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq'+$
    'qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq'
  menu = STRARR(8, 4)
  sel = 0
  top = WHERE(parent EQ -1)
  sz = SIZE(top)
  ntop = sz(1)
  line = '  '
  FOR i = 0, ntop - 1 DO line = line + word(top(i)) + '  '
  topwidth = STRLEN(line)
;  IF(topwidth GT 78) THEN PRINT, esc + '[?3h'
  IF(topwidth GT 78) THEN MESSAGE, $
    'The requested menu requires too much horizontal space.'
  PRINT, esc + '[?25l', FORMAT = '(a)'
clear_and_draw:
  IF(NOT KEYWORD_SET(noclear)) THEN PRINT, esc + '[2J'
  noclear = 0
;
;  Display the pull down menu.
;  ---------------------------
draw:
  line = '  '
  w = WHERE(top EQ sel)
  IF(w(0) NE -1) THEN BEGIN
    FOR i = 0, ntop - 1 DO $
      IF(top(i) NE sel) THEN line = line + word(top(i)) + '  ' $
                        ELSE line = line + wr + word(top(i)) + wz + '  '
    PRINT, esc+'['+STRTRIM(ypos,2)+'f'+cs+'l'+STRMID(qstr,0,topwidth)+'k'+cz, $
      FORMAT='(a)'
    PRINT, esc+'['+STRTRIM(ypos+1,2)+'f'+cs+'x'+cz+line+cs+'x'+cz, FORMAT='(a)'
    PRINT, esc+'['+STRTRIM(ypos+2,2)+'f'+cs+'m'+STRMID(qstr,0,topwidth)+ $
      'j'+cz, FORMAT='(a)'
  ENDIF ELSE BEGIN
    tophit = sel
    WHILE(parent(tophit) NE -1) DO tophit = parent(tophit)
    FOR i = 0, ntop - 1 DO $
      IF(top(i) NE tophit) THEN line = line + word(top(i)) + '  ' $
                           ELSE line = line + wb + word(top(i)) + wz + '  '
    tline = '  '
    FOR i = 0, tophit - 1 DO IF(parent(i) EQ -1) THEN tline=tline+word(i)+'  '
    left_margin = STRLEN(tline) + 1
    lm = STRTRIM(STRING(left_margin), 2)
    width = STRLEN(word(tophit))
    FOR i = 0, nword - 1 DO BEGIN
      IF(parent(i) EQ tophit) THEN BEGIN
        len = STRLEN(word(i))
        IF(len GT width) THEN width = len
      ENDIF
    ENDFOR
    PRINT,esc+'['+STRTRIM(ypos,2)+'f'+cs+'l'+STRMID(qstr,0,STRLEN(line)-7)+ $
      'k'+cz,FORMAT='(a)'
    PRINT,esc+'['+STRTRIM(ypos+1,2)+'f'+cs+'x'+cz+line+cs+'x'+cz,FORMAT='(a)'
    IF(left_margin+width LT topwidth) THEN BEGIN
      temp = 'm'+STRMID(qstr, 0, topwidth)+'j'
      STRPUT, temp, 'w', left_margin - 1
      STRPUT, temp, 'w', left_margin + width
    ENDIF ELSE BEGIN
      temp = 'm'+STRMID(qstr, 0, left_margin+width-1)+'k'
      STRPUT, temp, 'w', left_margin - 1
      STRPUT, temp, 'v', topwidth + 1
    ENDELSE
    PRINT, esc+'['+STRTRIM(ypos+2,2)+'f'+cs+temp+cz, FORMAT='(a)'
draw_below_bar:
    lno = 4
    FOR i = 0, nword - 1 DO BEGIN
      IF(parent(i) EQ tophit) THEN BEGIN
        temp = word(i)
        filler = ''
        WHILE(STRLEN(temp+filler) LT width) DO filler = filler + ' '
        IF(sel EQ i) THEN temp = wr + temp + wz
        temp = temp + filler
        ln = STRTRIM(STRING(lno+ypos-1), 2)
        PRINT,esc+'['+ln+';'+lm+'f'+cs+'x'+cz+temp+cs+'x'+cz,FORMAT='(a)'
        lno = lno + 1
      ENDIF
    ENDFOR
    ln = STRTRIM(STRING(lno+ypos-1), 2)
    PRINT,esc+'['+ln+';'+lm+'f'+cs+'m'+STRMID(qstr,0,width)+'j'+cz,FORMAT='(a)'
  ENDELSE
;
;  Monitor and respond to key presses until a carriage return is sent.
;  -------------------------------------------------------------------
watch_keyboard:
  WHILE 1 DO BEGIN
    CASE GET_KBRD(1) of
      cr   : BEGIN
carriage_return:
               w = WHERE(parent EQ sel)
               IF(w(0) NE -1) THEN BEGIN
                 sel = w(0)
                 GOTO, draw
               ENDIF ELSE BEGIN
                 PRINT, esc + '[?25h', FORMAT='(a)'
                 RETURN, STRTRIM(word(sel), 2)
               ENDELSE
             END
      esc  : IF(GET_KBRD(0) EQ '[') THEN GOTO, handle_arrowkeys
      csi  : GOTO, handle_arrowkeys
      tab  : GOTO, down_arrow
      cntlw: GOTO, clear_and_draw
      ELSE :
    ENDCASE
  ENDWHILE
handle_arrowkeys:
  CASE GET_KBRD(0) of
    'A' : BEGIN  ;up
            IF(parent(sel) NE -1) THEN BEGIN
              temp = parent(sel)
              sel = sel - 1
              IF(sel EQ -1) THEN sel = temp
              WHILE((sel NE TEMP) AND (parent(sel) NE temp)) DO BEGIN
                sel = sel - 1
                IF(sel EQ -1) THEN sel = temp
              ENDWHILE
              IF(sel EQ temp) THEN BEGIN
                PRINT,esc+'['+STRTRIM(ypos+2,2)+'f'+cs+'m'+ $
                  STRMID(qstr,0,topwidth)+'j'+cz+esc+'[0K',FORMAT='(a)'
                FOR i = 4, lno DO PRINT, esc + '[2K'
                GOTO, draw
              ENDIF ELSE GOTO, draw_below_bar
            ENDIF
          END
    'B' : BEGIN  ;down
down_arrow:
            IF(parent(sel) NE -1) THEN BEGIN
              temp = parent(sel)
              sel = sel + 1
              IF(sel EQ nword) THEN sel = temp
              WHILE((sel NE TEMP) AND (parent(sel) NE temp)) DO BEGIN
                sel = sel + 1
                IF(sel EQ nword) THEN sel = temp
              ENDWHILE
              IF(sel EQ temp) THEN BEGIN
                PRINT,esc+'['+STRTRIM(ypos+2,2)+'f'+cs+'m'+ $
                  STRMID(qstr,0,topwidth)+'j'+cz+esc+'[0K',FORMAT='(a)'
                FOR i = 4, lno DO PRINT, esc + '[2K'
                GOTO, draw
              ENDIF ELSE GOTO, draw_below_bar
            ENDIF ELSE GOTO, carriage_return
          END
    'C' : BEGIN
            IF(parent(sel) NE -1) THEN BEGIN
              PRINT,esc+'['+STRTRIM(ypos+2,2)+'f'+cs+'m'+ $
                STRMID(qstr,0,topwidth)+'j'+cz+esc+'[0K',FORMAT='(a)'
              FOR i = 4, lno DO PRINT, esc + '[2K'
            ENDIF
            sel = sel + 1
            IF(sel EQ nword) THEN sel = 0
            WHILE(parent(sel) NE -1) DO BEGIN
              sel = sel + 1
              IF(sel EQ nword) THEN sel = 0
            ENDWHILE
            GOTO, draw
          END
    'D' : BEGIN
            IF(parent(sel) NE -1) THEN BEGIN
              PRINT,esc+'['+STRTRIM(ypos+2,2)+'f'+cs+'m'+ $
                STRMID(qstr,0,topwidth)+'j'+cz+esc+'[0K',FORMAT='(a)'
              FOR i = 4, lno DO PRINT, esc + '[2K'
              sel = tophit
            ENDIF
            sel = sel - 1
            IF(sel EQ -1) THEN sel = nword - 1
            WHILE(parent(sel) NE -1) DO BEGIN
              sel = sel - 1
              IF(sel EQ -1) THEN sel = nword - 1
            ENDWHILE
            GOTO, draw
          END
    ELSE: GOTO, watch_keyboard
  ENDCASE
  GOTO, watch_keyboard
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


