PRO uimage
;+NAME/ONE LINE DESCRIPTION:
;    UIMAGE is a package that allows one to manipulate COBE data.
;
;DESCRIPTION:
;    UIMAGE is an image analysis program that works with the COBE
;    initial products and project data sets from all three instruments 
;    (FIRAS, DIRBE, and DMR), as well as some other approved 
;    astronomical data bases. The program is invocable from both 
;    X-window and non-X-window terminals, although the functionality 
;    is somewhat different between the two.
;
;    The UIMAGE package has been provided so that Guest Investigators
;    can have a means to perform numerous operations on data, whether or
;    not they know IDL.  A menu-driven environment is incorporated into
;    the program, to facilitate ease of use.
;
;    Data may be inputted into UIMAGE from a disk file, or else via
;    main-level UIDL variables.  The types of disk files from which
;    input may be acquired are COBETRIEVE-type files, COBE-specific FITS
;    files, or COBE-specific IDL savesets.  Main-level UIDL variables
;    may be transferred into the UIMAGE environment by using the
;    UPUTDATA routine while at the UIDL prompt.
;
;    The forms of data which can be inputted into UIMAGE are sky maps or
;    individual sky map faces in resolution 6, 7, 8, or 9, in the
;    sky-cube projection, and in either the ecliptic, equatorial, or
;    galactic coordinate system.
;
;    Data may be outputted from UIMAGE into a disk file or else
;    (indirectly) into main-level UIDL variables.  The types of output
;    disk files are COBE-specific FITS files or COBE-specific IDL
;    savesets.  Data may be transferred into main-level UIDL variables
;    by calling the UGETDATA routine from the UIDL prompt.
;
;CALLING SEQUENCE:
;    uimage
;
;ARGUMENTS (I=input, O=output, []=optional):
;    none
;
;EXAMPLE:
;    uimage
;
;WARNINGS:
;    Trying to input data that is not pixelized in the Quadrilateralized
;    sky-cube format will cause the data input routines to crash!
;#
;COMMON BLOCKS:  colors, nonxwindow
;
;LIBRARY CALLS:  dataman, dispmani, imagenha, instr_spec, lineplot,
;    one_line_menu, specoper, toggle_journal, ualgebra, uimage_help,
;    umenu.
;
;PROCEDURE:
;    Set !QUIET to 1.  If an X-window terminal is being used, then define
;    a default color table.  Define an array of structures which tell
;    what procedures will be invocable from the main menu, under what
;    terminal types they will be allowed to be invoked, and with what
;    phrases they will be referenced.  If a non-X window terminal is
;    being used, then have the user identify the terminal type.  Call
;    UMENU to display the main menu, and then execute any selected
;    option.  Repeat the calls to UMENU until an exit request is selected.
;
;REVISION HISTORY:
;    Written by John Ewing, ARC, January 1992.
;
;  SPR 10332 Dec 21, 92   Switch help & exit options on main menu.  J Ewing
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;  SPR 11314 15 Sep 1993  Set journaling on as default. J. Newmark
;  SPR 11757 17 May 1994  Update for DEC Alpha. J. Newmark
;  SPR 11905 07-Sep-1994  Ingest DIRBE PDS's. J. Newmark
;.TITLE
;Routine UIMAGE
;-
  COMMON nonxwindow, term_type
  COMMON colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
  COMMON journal, journal_on, luj
  COMMON color_values, c_badpix, c_draw, c_scalemin
  COMMON history, uimage_version
  old_quiet = !QUIET
  !QUIET = 1
  uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 1
;
;  If a non-X-window terminal is being used, then have the user
;  identify the terminal type.
;  ------------------------------------------------------------
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    IF(NOT defined(term_type)) THEN BEGIN
      supported_terms=['Identify terminal type','REGIS','TEK','Other']
      sel = one_line_menu(supported_terms)
      term_type = STRUPCASE(supported_terms(sel))
      IF(term_type EQ 'TEK') THEN SET_PLOT, 'TEK'
      IF(term_type EQ 'REGIS') THEN SET_PLOT, 'REGIS'
    ENDIF
  ENDIF ELSE BEGIN
;
  IF (!D.NAME EQ 'WIN') THEN DEVICE, RETAIN = 2 ELSE $ 
      DEVICE, RETAIN = 2, PSEUDO_COLOR = 8
;  Define color values which will be used for special purposes.
;  ------------------------------------------------------------
    IF(NOT defined(c_scalemin)) THEN BEGIN
      c_badpix = 0  ;  never change this setting.
      c_draw = 1
      c_scalemin = 10
    ENDIF
    IF(defined(r_orig)) THEN BEGIN
;
;  If the color table is already defined, then make sure the graphics
;  color values (1..c_scalemin-1) are distinct from the background color
;  value (0).
;  ---------------------------------------------------------------------
      r_term = FLOAT(r_orig(1:c_scalemin-1) - r_orig(0)) ^ 2
      g_term = FLOAT(g_orig(1:c_scalemin-1) - g_orig(0)) ^ 2
      b_term = FLOAT(b_orig(1:c_scalemin-1) - b_orig(0)) ^ 2
      contrast = SQRT(r_term + g_term + b_term)
      low_contrast = WHERE(contrast LT 50, nlow)
      IF(nlow GT 0) THEN BEGIN
        r_spec = [0, 255, 255,   0,   0, 255, 255, 182, 255, 255]
        g_spec = [0, 255,   0, 255,   0,   0, 255,  11, 170, 255]
        b_spec = [0, 255,   0,   0, 255, 255,   0,   0,  78, 255]
        FOR i = 0, c_scalemin - 2 DO BEGIN
          IF(i EQ 0) THEN BEGIN
            r_orig(1) = 255 - r_orig(0)
            g_orig(1) = 255 - g_orig(0)
            b_orig(1) = 255 - b_orig(0)
          ENDIF ELSE BEGIN
            r_orig(i + 1) = r_spec(i + 1)
            g_orig(i + 1) = g_spec(i + 1)
            b_orig(i + 1) = b_spec(i + 1)
          ENDELSE
        ENDFOR
        TVLCT, r_orig, g_orig, b_orig
        TVLCT, r_orig, g_orig, b_orig
        PRINT, 'A change to the color table has been made.'
      ENDIF
    ENDIF ELSE uloadct, 3
    xdisplay_all
  ENDELSE
;
; Set up default journaling
;
  IF(!cgis_os EQ 'vms') THEN filename='SYS$LOGIN:uimage.jnl' $
     ELSE BEGIN
        IF(!cgis_os EQ 'unix') THEN filename='$HOME/uimage.jnl' $
        ELSE IF (!cgis_os EQ 'windows') THEN filename='c:\uimage.jnl'
     ENDELSE
   err=1
   WHILE(err NE 0) DO BEGIN
      OPENW, luj, filename, ERROR = err, /GET_LUN
      IF(err NE 0) THEN BEGIN
          PRINT, bold('An inappropriate file name was entered.')
          PRINT, ' '
          READ, prompt='Enter Log file name:  ', filename
      ENDIF
   ENDWHILE
   PRINTF, luj, '----------------------------------------' + $
                   '--------------------------------------'
;
;  Create a menu of selectable operations.
;  ---------------------------------------
  n_routines = 14
  routine = REPLICATE({rs, vt: 0, name: '', title: ''}, n_routines)
  routine(00) = {rs, 1, 'DATAMANA',       'Data I/O and Management...'}
  routine(01) = {rs, 0, '',               ''}
  routine(02) = {rs, 1, 'DISPMANI',       'Display Manipulation...'}
  routine(03) = {rs, 1, 'IMAGENHA',       'Image Enhancement...'}
  routine(04) = {rs, 1, 'UALGEBRA',       'Algebraic Operations...'}
  routine(05) = {rs, 1, 'SPECOPER',       'Spectrum Operations...'}
  routine(06) = {rs, 1, 'LINEPLOT',       'Line Plots and Statistics...'}
  routine(07) = {rs, 1, 'MODELFIT',       'Modeling and Fitting...'}
  routine(08) = {rs, 0, '',               ''}
  str = 'Journal Enable/Disable   '
  IF(journal_on EQ 1) THEN str = str + 'ON' ELSE str = str + 'OFF'
  routine(09) = {rs, 1, 'TOGGLE_JOURNAL', str}
  routine(10) = {rs, 1, 'SENDCOMMENT',    'Report Problems or Comments'}
  routine(11) = {rs, 0, '',               ''}
  routine(12) = {rs, 1, 'UIMAGE_HELP',    'HELP'}
  routine(13) = {rs, 1, 'EXIT',           'Exit UIMAGE'}
;
  menu_title = 'UIMAGE [MAIN MENU]'
  menu = STRARR(n_routines + 1)
  menu(0) = menu_title
  sel = 1
remake_menu:
  str = 'Journal Enable/Disable   '
  IF(journal_on EQ 1) THEN str = str + 'ON' ELSE str = str + 'OFF'
  routine(09) = {rs, 1, 'TOGGLE_JOURNAL', str}
  FOR i = 0, n_routines - 1 DO menu(i + 1) = routine(i).title
;
;  Display the menu.
;  -----------------
show_menu:
  sel = UMENU(menu, title = 0, init = sel, valid = valid)
;
;  Invoke the routine that implements the selected operation.
;  ----------------------------------------------------------
  routine_name = routine(sel-1).name
  IF(routine_name EQ '')     THEN GOTO, show_menu
  IF(routine_name EQ 'EXIT') THEN BEGIN
    !quiet = old_quiet
    uimage_version = -1
    RETURN
  ENDIF
  CALL_PROCEDURE, routine_name
  GOTO, remake_menu
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


