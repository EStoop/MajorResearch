PRO uimage_help, menu_title
;+
;  UIMAGE_HELP - a UIMAGE-specific routine.  This routine sees to it that
;  a user can examine the help-text associated with the menu identified
;  by the supplied menu title string MENU_TITLE.  If MENU_TITLE is not
;  supplied, then the help-text for UIMAGE's main menu is displayed.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 08 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10433  Jan 21 93  Replaced first entry in 'possible_file' array
;                        with getenv('uimage_help') so it will work 
;                        under UNIX
;  SPR 10442  Feb 08 93  Change journaling.  J Ewing
;--------------------------------------------------------------------------
  COMMON uimage_help, helpfile
  COMMON journal, journal_on, luj
  IF(NOT defined(journal_on)) THEN journal_on = 0
  esc = STRING(27b)
  erase_screen = esc + '[2J'
  cursor_home = esc + '[H'
;
;  Set up the title for the scroll window.
;  ---------------------------------------
  IF(N_PARAMS(0) EQ 0) THEN BEGIN
    menu_title = 'MAIN MENU'
    title = 'Help for the UIMAGE main menu'
  ENDIF ELSE title = 'Help for "' + menu_title + '"'
  text = STRARR(150)
  num = 0
;
;  Try to identify the name of the file which contains all the help-text
;  for the UIMAGE package.
;  ---------------------------------------------------------------------
  possible_file = [getenv('uimage_help'), $
                   'cgis$uitree:uimage_help.dat', $
                   'uimage_help.dat' ]
  filename = ''
  FOR i = 0, N_ELEMENTS(possible_file) - 1 DO BEGIN
    IF(filename EQ '') THEN BEGIN
      find = FINDFILE(possible_file(i))
      IF(find(0) NE '') THEN filename = possible_file(i)
    ENDIF
  ENDFOR
  IF(filename EQ '') THEN BEGIN
    PRINT, 'The file containing help info for UIMAGE (UIMAGE_HELP.DAT) was not'
    PRINT, 'found.  Please have the person maintaining UIMAGE at your site set up'
    PRINT, 'a logical named UIMAGE_HELP which points to that file.'
    RETURN
  ENDIF
;
;  Open the help file.
;  -------------------
  OPENR, lun, filename, /GET_LUN
  str = ''
  curr_title = ''
  cap_menu_title = STRUPCASE(menu_title)
  READF, lun, curr_title
loop:
  IF(curr_title EQ cap_menu_title) THEN BEGIN
;
;  Read in the help info into TEXT.
;  --------------------------------
    READF, lun, str
    WHILE((STRMID(str, 0, 5) NE '-----') AND (NOT EOF(lun))) DO BEGIN
      text(num) = str
      num = num + 1
      READF, lun, str
    ENDWHILE
    text = text(0:(num - 1))
    FREE_LUN, lun
;
;  If run on X-windows, then eliminate special character sequences, such
;  as "{B}", which indicate to the software that special character
;  attributes are to be enabled for VT200-type terminals.
;  ---------------------------------------------------------------------
    IF(!D.NAME EQ 'X' or !D.NAME EQ 'WIN') THEN BEGIN
      FOR i = 0, num - 1 DO BEGIN
        ic = STRPOS(text(i), '{')
        WHILE(ic NE -1) DO BEGIN
          ic2 = STRPOS(text(i), '}')
          text(i) = STRMID(text(i), 0, ic) + $
                    STRMID(text(i), ic2 + 1, STRLEN(text(i)))
          ic = STRPOS(text(i), '{')
        ENDWHILE
      ENDFOR
    ENDIF
;
;  Call USCROLL to display the text.
;  ---------------------------------
    uscroll, text, title = title
    IF(journal_on) THEN BEGIN
      PRINTF, luj, 'On-line Help   ---   ' + title
      FOR i = 0, N_ELEMENTS(text) - 1 DO PRINTF, luj, '  ' + text(i)
      PRINTF, luj, '----------------------------------------' + $
                   '--------------------------------------'
    ENDIF
    RETURN
  ENDIF ELSE BEGIN
;
;  Read in a menu title.  Ignore anything in the file which is not preceded
;  by a row of hyphens.  (Menu titles are preceded by a row of hyphens).
;  ------------------------------------------------------------------------
    READF, lun, str
    WHILE((STRMID(str, 0, 5) NE '-----') AND (NOT EOF(lun))) DO READF, lun, str
    IF(NOT EOF(lun)) THEN BEGIN
      READF, lun, curr_title
      curr_title = STRUPCASE(curr_title)
      GOTO, loop
    ENDIF ELSE BEGIN
      PRINT, '<UIMAGE is still under development.  HELP information has not yet'
      PRINT, '             been incorporated for this menu.>'
    ENDELSE
  ENDELSE
  FREE_LUN, lun
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


