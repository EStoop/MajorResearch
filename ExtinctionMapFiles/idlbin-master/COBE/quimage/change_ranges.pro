PRO change_ranges
;+
;  CHANGE_RANGES - a UIMAGE-specific routine.  This routine allows a user
;  to change the X and Y ranges of a graph.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 06 93  Add more comments.  J Ewing
;  SPR 10456  Jan 13 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 10956  May 17 93  Properly handle composite graphs. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
;
;  Have the user select a graph.
;  -----------------------------
show_menu:
  menu_title = 'Change X and Y axis ranges'
  name = select_object(/graph, /help, /exit, title = menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no graphs available.', /CONT
    RETURN
  ENDIF
;
;  Extract attributes.
;  -------------------
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('scatter = ' + name + '.scatter')
  j = EXECUTE('x = ' + name + '.data(*, 0)')
  j = EXECUTE('y = ' + name + '.data(*, 1)')
  j = EXECUTE('xstart = ' + name + '.xstart')
  j = EXECUTE('xstop = ' + name + '.xstop')
  j = EXECUTE('ystart = ' + name + '.ystart')
  j = EXECUTE('ystop = ' + name + '.ystop')
  j = EXECUTE('window = ' + name + '.window')
  j = EXECUTE('multi = ' + name + '.multi')
  PRINT, ' '
  PRINT, 'Selected graph:  ' + bold(title)
;
;  Determine the min and max of the full X and Y ranges.
;  -----------------------------------------------------
  IF(scatter NE 1) THEN BEGIN
    xmin = MIN(x)
    xmax = MAX(x)
    ymin = MIN(y)
    ymax = MAX(y)
  ENDIF ELSE BEGIN
    j = EXECUTE('win_orig1 = ' + name + '.win_orig1')
    j = EXECUTE('win_orig2 = ' + name + '.win_orig2')
    name1 = get_name(win_orig1)
    name2 = get_name(win_orig2)
    j = EXECUTE('badpixval1 = ' + name1 + '.badpixval')
    j = EXECUTE('badpixval2 = ' + name2 + '.badpixval')
    j = EXECUTE('good = WHERE((' + name1 + '.data NE badpixval1) AND (' + $
      name2 + '.data NE badpixval2))')
    IF(good(0) EQ -1) THEN BEGIN
      MESSAGE, 'The two images have no overlapping good pixels.', /CONT
      RETURN
    ENDIF
    j = EXECUTE('xmin = MIN(' + name1 + '.data(good))')
    j = EXECUTE('xmax = MAX(' + name1 + '.data(good))')
    j = EXECUTE('ymin = MIN(' + name2 + '.data(good))')
    j = EXECUTE('ymax = MAX(' + name2 + '.data(good))')
  ENDELSE
  IF (multi NE 0) THEN BEGIN
   mxmin=FLTARR(multi)
   mxmax=FLTARR(multi)
   mymin=FLTARR(multi)
   mymax=FLTARR(multi)
   j=EXECUTE('m_id = ' +name + '.m_id')
   FOR i=0,multi-1 DO BEGIN
     mname=get_name(m_id(i))
     j=EXECUTE('mx='+mname+'.data(*,0)')
     j=EXECUTE('my='+mname+'.data(*,1)')
     mxmin(i)=MIN(mx)
     mxmax(i)=MAX(mx)
     mymin(i)=MIN(my)
     mymax(i)=MAX(my)
   ENDFOR
   xmin=MIN(mxmin)
   xmax=MAX(mxmax)
   ymin=MIN(mymin)
   ymax=MAX(mymax)
  ENDIF 
;
;  Let the user set the X range if he wants to.
;  --------------------------------------------
  PRINT, 'The X data ranges from ' + STRTRIM(STRING(xmin), 2) + ' to ' + $
    STRTRIM(STRING(xmax), 2)
  menu = ['Do you wish to set the X display range?', 'Yes', 'No']
  selx = one_line_menu(menu)
  IF(selx EQ 1) THEN BEGIN
    str = ''
    val = 0
    WHILE(val EQ 0) DO BEGIN
      READ, bold('minimum') + ':  ', str
      val = validnum(str)
      IF(val EQ 1) THEN BEGIN
        xstart = FLOAT(str)
      ENDIF ELSE BEGIN
        PRINT, 'Invalid number, please re-enter.'
      ENDELSE
    ENDWHILE
    val = 0
    WHILE(val EQ 0) DO BEGIN
      READ, bold('maximum') + ':  ', str
      val = validnum(str)
      IF(val EQ 1) THEN BEGIN
        xstop = FLOAT(str)
        IF(xstop LT xstart) THEN BEGIN
          val = 0
          PRINT, 'Invalid number, please re-enter.'
        ENDIF
      ENDIF ELSE BEGIN
        PRINT, 'Invalid number, please re-enter.'
      ENDELSE
    ENDWHILE
    j = EXECUTE(name + '.xstart = xstart')
    j = EXECUTE(name + '.xstop = xstop')
  ENDIF
;
;  Let the user set the Y range if he wants to.
;  --------------------------------------------
  PRINT, 'The Y data ranges from ' + STRTRIM(STRING(ymin), 2) + ' to ' + $
    STRTRIM(STRING(ymax), 2)
  menu = ['Do you wish to set the Y display range?', 'Yes', 'No']
  sely = one_line_menu(menu)
  IF(sely EQ 1) THEN BEGIN
    str = ''
    val = 0
    WHILE(val EQ 0) DO BEGIN
      READ, bold('minimum') + ':  ', str
      val = validnum(str)
      IF(val EQ 1) THEN BEGIN
        ystart = FLOAT(str)
      ENDIF ELSE BEGIN
        PRINT, 'Invalid number, please re-enter.'
      ENDELSE
    ENDWHILE
    val = 0
    WHILE(val EQ 0) DO BEGIN
      READ, bold('maximum') + ':  ', str
      val = validnum(str)
      IF(val EQ 1) THEN BEGIN
        ystop = FLOAT(str)
        IF(ystop LT ystart) THEN BEGIN
          val = 0
          PRINT, 'Invalid number, please re-enter.'
        ENDIF
      ENDIF ELSE BEGIN
        PRINT, 'Invalid number, please re-enter.'
      ENDELSE
    ENDWHILE
    j = EXECUTE(name + '.ystart = ystart')
    j = EXECUTE(name + '.ystop = ystop')
  ENDIF
;
;  Re-display the graph.
;  ---------------------
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN WSET, window
  PLOT_GRAPH, name
;
;  If journaling is enabled, then send some info to the journal file.
;  ------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title
    IF(selx EQ 1) THEN PRINTF, luj, '  new X range:  ' + $
      STRTRIM(xstart, 2) + ' to ' + STRTRIM(xstop,2) $
      ELSE PRINTF, luj, '  the X range was left unchanged'
    IF(sely EQ 1) THEN PRINTF, luj, '  new Y range:  ' + $
      STRTRIM(ystart, 2) + ' to ' + STRTRIM(ystop,2) $
      ELSE PRINTF, luj, '  the Y range was left unchanged'
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


