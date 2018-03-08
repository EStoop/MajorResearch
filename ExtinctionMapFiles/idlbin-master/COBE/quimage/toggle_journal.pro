PRO toggle_journal
;+
;  TOGGLE_JOURNAL - a UIMAGE-specific routine.  This routine gives the
;  user to enable or disable journaling.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;--------------------------------------------------------------------------
  COMMON journal, journal_on, luj
  IF(NOT defined(journal_on)) THEN journal_on = 0
  PRINT, ' '
  IF(journal_on EQ 0) THEN BEGIN
;
;  Journaling is disabled. See if the user wants to enable it.
;  -----------------------------------------------------------
    PRINT, 'Journaling is currently set to ' + bold('OFF') + '.'
    menu = ['Change it to ON?', 'Yes', 'No']
    sel = one_line_menu(menu)
    IF(sel EQ 1) THEN BEGIN
;
;  Have the user type in the name of the journal file.
;  ---------------------------------------------------
      err = 1
      filename = ''
      PRINT, 'Enter the complete name of the file to which the journaled text'
      PRINT, '  should be sent.  Example:  SYS$LOGIN:today.sav'
      WHILE(err NE 0) DO BEGIN
        READ, underline('Enter file name:') + '  ', filename
        OPENW, luj, filename, ERROR = err, /GET_LUN
        IF(err NE 0) THEN BEGIN
          PRINT, bold('An inappropriate file name was entered.')
          PRINT, ' '
        ENDIF
      ENDWHILE
      PRINTF, luj, '----------------------------------------' + $
                   '--------------------------------------'
      journal_on = 1
      PRINT, 'Journaling is now set to ' + bold('ON') + '.' 
    ENDIF ELSE PRINT, 'Journaling is still set to ' + bold('OFF') + '.' 
  ENDIF ELSE BEGIN
;
;  Journaling is enabled. See if the user wants to disable it.
;  -----------------------------------------------------------
    PRINT, 'Journaling is currently set to ' + bold('ON') + '.'
    menu = ['Change it to OFF?', 'Yes', 'No']
    sel = one_line_menu(menu)
    IF(sel EQ 1) THEN BEGIN
      journal_on = 0
      FREE_LUN, luj
      PRINT, 'Journaling is now set to ' + bold('OFF') + '.' 
    ENDIF ELSE PRINT, 'Journaling is still set to ' + bold('ON') + '.' 
  ENDELSE
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


