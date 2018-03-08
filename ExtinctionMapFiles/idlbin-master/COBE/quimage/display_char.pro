PRO display_char, name
;+
; NAME:	
;	DISPLAY_CHAR
; PURPOSE:
;       For UIMAGE users, display the characteristics associated
;       with a selected object.
; CATEGORY:
;	Data access.
; CALLING SEQUENCE:
;       display_char
; INPUTS:
;       None.
; OUTPUTS:
;       None.
;#
; COMMON BLOCKS:
;	uimage_data
; RESTRICTIONS:
;       For use with UIMAGE only.
; PROCEDURE:
;       Have the user identify an object via SELECT_OBJECT.  Then print
;       out the value of various fields.
; SUBROUTINES CALLED:
;       SELECT_OBJECT
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, March 1992.
; SPR 10442  Feb 04 93  Change journaling.  J Ewing
; SPR 10708  Mar 17 93  Add in code for new .LINKWEIGHT field.  J Ewing
; SPR 10740  Mar 25 93  .LINK3D = -1 means no linkage.  J Ewing
; SPR 10746  Mar 25 93  Remove obsolete code from DISPLAY_CHAR.  J Ewing
; SPR 10829  Apr 1993   Add in code for new .HIDDEN field. J. Newmark
; SPR 10843  Apr 22 93  Correct code for graphs. J Nemark
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;----------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
     freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON journal,journal_on,luj
  COMMON xscreen,magnify,block_usage,scrsiz
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  esc = STRING(27b)
  cuu = esc + '[A'
  menu_title = 'Report object attributes'
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Put up a menu so the user can select an operand.
;  ------------------------------------------------
show_menu:
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
           /face8,/face9,/graph,/zoomed,/proj_map,/proj2_map,/object3d,$
           /help,/exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no objects available.', /CONT
    RETURN
  ENDIF
have_name:
  text = [' ']
  j = EXECUTE('title     = ' + name + '.title')
  IF(uimage_version NE 4) THEN BEGIN
    PRINT, ' '
    PRINT, underline('Attributes of ' + bold(title)), FORMAT = '(a)'
  ENDIF
  first_letter = STRMID(name, 0, 1)
  IF((first_letter EQ 'M') OR (first_letter EQ 'F') OR (first_letter EQ 'P'))$
    THEN BEGIN
;
;  Print out fields for an element of MAP* or FACE* or PROJ_MAP
;  ------------------------------------------------------------
    j = EXECUTE('badpixval = ' + name + '.badpixval')
    j = EXECUTE('bandno    = ' + name + '.bandno')
    j = EXECUTE('bandwidth = ' + name + '.bandwidth')
    j = EXECUTE('coordinate_system = ' + name + '.coordinate_system')
    j = EXECUTE('data_min  = ' + name + '.data_min')
    j = EXECUTE('data_max  = ' + name + '.data_max')
    j = EXECUTE('faceno    = ' + name + '.faceno')
    j = EXECUTE('frequency = ' + name + '.frequency')
    j = EXECUTE('instrume  = ' + name + '.instrume')
    j = EXECUTE('link3d    = ' + name + '.link3d')
    j = EXECUTE('linkweight= ' + name + '.linkweight')
    j = EXECUTE('hidden    = ' + name + '.hidden')
    j = EXECUTE('orient    = ' + name + '.orient')
    j = EXECUTE('projection= ' + name + '.projection')
    j = EXECUTE('scale_min = ' + name + '.scale_min')
    j = EXECUTE('scale_max = ' + name + '.scale_max')
    j = EXECUTE('units     = ' + name + '.units')
    j = EXECUTE('version   = ' + name + '.version')
    CASE first_letter OF
      'M' : res = STRMID(name, 3, 1)
      'F' : res = STRMID(name, 4, 1)
      ELSE:
    ENDCASE
    j = EXECUTE('siz = SIZE(' + name + '.data)')
    dim1 = siz(1)
    dim2 = siz(2)
    dim1      = STRTRIM(STRING(dim1), 2)
    dim2      = STRTRIM(STRING(dim2), 2)
    badpixval = STRTRIM(STRING(badpixval), 2)
    bandno    = STRTRIM(STRING(bandno), 2)
    bandwidth = STRTRIM(STRING(bandwidth), 2)
    IF(data_min LE data_max) THEN BEGIN
      data_min  = STRTRIM(STRING(data_min), 2)
      data_max  = STRTRIM(STRING(data_max), 2)
    ENDIF ELSE BEGIN
      data_min = '<undefined>'
      data_max = '<undefined>'
    ENDELSE
    faceno    = STRTRIM(STRING(faceno), 2)
    frequency = STRTRIM(STRING(frequency), 2)
    IF(link3d NE -1) THEN link3dname = object3d(link3d).title $
                     ELSE link3dname = '<no link to a 3-D object>'
    scale_min = STRTRIM(STRING(scale_min), 2)
    scale_max = STRTRIM(STRING(scale_max), 2)
    text = [text, '  Array dimensions:      ' + dim1 + ' x ' + dim2]
    IF(first_letter EQ 'P') THEN BEGIN
      j = EXECUTE('win_orig = ' + name + '.win_orig')
      name_orig = get_name(win_orig)
      j = EXECUTE('title_orig = ' + name_orig + '.title')
      text = [text, '  Associated sky-cube:   ' + title_orig]
    ENDIF                                                      
    IF(linkweight NE -1) THEN BEGIN
      temp = get_name(linkweight)
      j = EXECUTE('linkwtname = ' + temp + '.title')
    ENDIF ELSE linkwtname = '<no link to a weight array>'
    IF (hidden NE -1) THEN hidekey='<object is not displayed>' $
      ELSE hidekey='<object is displayed>'
    text = [text, '  Associated 3-D object: ' + link3dname]
    text = [text, '  Associated weights:    ' + linkwtname]
    text = [text, '  Hidden Object:         ' + hidekey]
    text = [text, '  Bad pixel value:       ' + badpixval]
    text = [text, '  Band number:           ' + bandno]
    text = [text, '  Band width:            ' + bandwidth]
    text = [text, '  Coordinate-system:     ' + coordinate_system]
    text = [text, '  Data minimum:          ' + data_min]
    text = [text, '  Data maximum:          ' + data_max]
    IF(first_letter EQ 'F') THEN $
      text = [text, '  Face number:           ' + faceno]
    text = [text, '  Frequency:             ' + frequency]
    text = [text, '  Instrument:            ' + instrume]
    text = [text, '  Orientation:           ' + orient]
    text = [text, '  Projection:            ' + projection]
    IF(first_letter NE 'P') THEN $
      text = [text, '  Resolution Index:    ' + res]
    text = [text, '  Scaling minimum:       ' + scale_min]
    text = [text, '  Scaling maximum:       ' + scale_max]
    text = [text, '  Units:                 ' + units]
  ENDIF
  IF(first_letter EQ 'G') THEN BEGIN
;
;  Print out fields for an element of GRAPH
;  ----------------------------------------
    j = EXECUTE('linkweight= ' + name + '.linkweight')
    j = EXECUTE('hidden    = ' + name + '.hidden')
    IF(linkweight NE -1) THEN BEGIN
      temp = get_name(linkweight)
      j = EXECUTE('linkwtname = ' + temp + '.title')
    ENDIF ELSE linkwtname = '<no link to a weight array>'
    IF (hidden NE -1) THEN hidekey='<object is not displayed>' $
      ELSE hidekey='<object is displayed>'
    j = EXECUTE('badpixval = ' + name + '.badpixval')
    j = EXECUTE('logflag   = ' + name + '.logflag')
    j = EXECUTE('num       = ' + name + '.num')
    j = EXECUTE('psymbol   = ' + name + '.psymbol')
    j = EXECUTE('special   = ' + name + '.special')
    j = EXECUTE('x_title   = ' + name + '.x_title')
    j = EXECUTE('y_title   = ' + name + '.y_title')
    IF(logflag NE 0) THEN ptype='LOG' ELSE ptype='LINEAR'
    badpixval = STRTRIM(STRING(badpixval), 2)
    IF(num LE 1024) THEN num = FIX(num)
    num       = STRTRIM(STRING(num), 2)
    psymbol   = STRTRIM(STRING(psymbol), 2)
    IF(special EQ 1) THEN text = [text, '  Bad pixel value:   ' + badpixval]
    text = [text, '  Number of points:  ' + num]
    text = [text, '  Plot type:         ' + ptype]
    text = [text, '  Psymbol:           ' + psymbol]
    text = [text, '  X-axis title:      ' + x_title]
    text = [text, '  Y-axis title:      ' + y_title]
    text = [text, '  Associated weights:    ' + linkwtname]
    text = [text, '  Hidden Object:         ' + hidekey]
  ENDIF
  IF(first_letter EQ 'Z') THEN BEGIN
;
;  Print out fields for an element of ZOOMED
;  -----------------------------------------
    j = EXECUTE('specific_zoom = ' + name + '.specific_zoom')
    j = EXECUTE('start_x   = ' + name + '.start_x')
    j = EXECUTE('start_y   = ' + name + '.start_y')
    j = EXECUTE('stop_x    = ' + name + '.stop_x')
    j = EXECUTE('stop_y    = ' + name + '.stop_y')
    j = EXECUTE('win_orig  = ' + name + '.win_orig')
    name_orig = get_name(ABS(win_orig))
    j = EXECUTE('title_orig= ' + name_orig + '.title')
    specific_zoom = STRTRIM(STRING(specific_zoom), 2)
    start_x   = STRTRIM(STRING(start_x), 2)
    start_y   = STRTRIM(STRING(start_y), 2)
    stop_x    = STRTRIM(STRING(stop_x), 2)
    stop_y    = STRTRIM(STRING(stop_y), 2)
    text = [text, '  Data origin:       ' + title_orig]
    text = [text, '  X range:           ' + start_x + ' to ' + stop_x]
    text = [text, '  Y range:           ' + start_y + ' to ' + stop_y]
    text = [text, '  Zoom factor:       ' + specific_zoom]
  ENDIF
;
;  Print out fields for an OBJECT3D.
;  ---------------------------------
  IF(first_letter EQ 'O') THEN BEGIN
    str_index3d=STRMID(name,9,1)
    j = EXECUTE('havewt = wt3d_' + str_index3d)
    wsz=SIZE(havewt)
    j = EXECUTE('badpixval = ' + name + '.badpixval')
    j = EXECUTE('coordinate_system = ' + name + '.coordinate_system')
    j = EXECUTE('dim1      = ' + name + '.dim1')
    j = EXECUTE('dim2      = ' + name + '.dim2')
    j = EXECUTE('dim3      = ' + name + '.dim3')
    j = EXECUTE('faceno    = ' + name + '.faceno')
    j = EXECUTE('instrume  = ' + name + '.instrume')
    j = EXECUTE('linkweight = ' + name + '.linkweight')
    j = EXECUTE('hidden    = ' + name + '.hidden')
    j = EXECUTE('orient    = ' + name + '.orient')
    j = EXECUTE('projection= ' + name + '.projection')
    j = EXECUTE('units     = ' + name + '.units')
    j = EXECUTE('version   = ' + name + '.version')
    IF(linkweight NE -1) THEN BEGIN
      temp = get_name(linkweight)
      j = EXECUTE('linkwtname = ' + temp + '.title')
    ENDIF ELSE IF (wsz(0) EQ 3) THEN linkwtname='<3D weight array>' $
               ELSE linkwtname = '<no link to a weight array>'
    dim1      = STRTRIM(STRING(dim1), 2)
    dim2      = STRTRIM(STRING(dim2), 2)
    dim3      = STRTRIM(STRING(dim3), 2)
    badpixval = STRTRIM(STRING(badpixval), 2)
    faceno    = STRTRIM(STRING(faceno), 2)
    text = [text, '  Array dimensions:  ' + dim1 + ' x ' + dim2 + ' x ' + dim3]
    text = [text, '  Bad pixel value:   ' + badpixval]
    text = [text, '  Associated weights:    ' + linkwtname]
    text = [text, '  Coordinate-system: ' + coordinate_system]
    IF(dim1 EQ dim2) THEN $
      text = [text, '  Face number:       ' + faceno]
    text = [text, '  Instrument:        ' + instrume]
    IF (hidden NE -1) THEN hidekey='<object is not displayed>' $
      ELSE hidekey='<object is displayed>'
    text = [text, '  Hidden Object:     ' + hidekey]
    text = [text, '  Orientation:       ' + orient]
    text = [text, '  Projection:        ' + projection]
    text = [text, '  Units:             ' + units]
  ENDIF
  sz = SIZE(text)
  dim1 = sz(1)
;
;  Display the text as appropriate.
;  --------------------------------
  IF((!D.NAME EQ 'X' OR !D.NAME EQ 'WIN') AND (uimage_version EQ 4)) THEN BEGIN
    offset = get_location(nblocks_x = 14, nblocks_y = 14)
    base = WIDGET_BASE(TITLE = 'Attributes of ' + title, $
      XOFFSET = offset(0), YOFFSET = offset(1), XSIZE = 34*12)
    w1 = WIDGET_BUTTON(base, VALUE = "Done")
    w2 = WIDGET_TEXT(base, XSIZE = 34, YSIZE = 15, VALUE = text, /SCROLL, $
      XOFFSET = 0, YOFFSET = 26)
    WIDGET_CONTROL, base, /REALIZE
    XMANAGER, 'UIMAGE4', base
    w = WHERE(block_usage EQ -1)
    block_usage(w) = base
  ENDIF ELSE BEGIN
    FOR i = 1, dim1 - 1 DO PRINT, text(i)
  ENDELSE
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    PRINT, ' '
    PRINT, '<hit any key to continue>'
    c = GET_KBRD(1)
    PRINT, cuu, '                             ', cuu
  ENDIF
;
;  If journaling is enabled, then send some info to the journal file.
;  ------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  Attributes of ' + title
    FOR i = 1, dim1 - 1 DO PRINTF, luj, '  ' + text(i)
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


