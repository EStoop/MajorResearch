PRO uscatter
;+
;  USCATTER - a UIMAGE-specific routine.  This routine allows a UIMAGE user
;  to create a scatter plot by selecting two images of the same size.  A
;  GRAPH-type data-object is generated.
;#
;  Written by John Ewing.
;  SPR 10385  Dec 22 92  Set XSTART, XSTOP, etc for graph-object.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  nblocks_x = 5
  nblocks_y = 5
  first = 1
;
;  Have the user select the image for the X-axis.
;  ----------------------------------------------
  menu_title = 'Scatter plot - select image for X-axis'
show_menu:
  name1 = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
           /face8,/face9,/proj_map,/proj2_map,/help,/exit,title=menu_title)
  IF(name1 EQ 'EXIT') THEN RETURN
  IF(name1 EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name1 EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no objects available.', /CONT
    RETURN
  ENDIF
  j = EXECUTE('title1 = ' + name1 + '.title')
  j = EXECUTE('window1 = ' + name1 + '.window')
;
;  Have the user select the image for the Y-axis.
;  ----------------------------------------------
  menu2_title = 'Scatter plot - select image for Y-axis'
show_menu2:
  type = STRMID(name1, 0, STRPOS(name1, ')') - 2)
  j = EXECUTE('name2 = select_object(/' + type + ', /help, /exit, ' + $
    'title = menu2_title)')
  IF(name2 EQ 'EXIT') THEN RETURN
  IF(name2 EQ 'HELP') THEN BEGIN
    uimage_help, menu2_title
    GOTO, show_menu2
  ENDIF
  j = EXECUTE('title2 = ' + name2 + '.title')
  j = EXECUTE('window2 = ' + name2 + '.window')
;
;  Set-up and plot the scatter graph.
;  ----------------------------------
  j = EXECUTE('xstart = ' + name1 + '.data_min')
  j = EXECUTE('xstop  = ' + name1 + '.data_max')
  j = EXECUTE('ystart = ' + name2 + '.data_min')
  j = EXECUTE('ystop  = ' + name2 + '.data_max')
  j = EXECUTE('sz = SIZE(' + name1 + '.data)')
  dim1 = sz(1)
  dim2 = sz(2)
  num = FLOAT(dim1) * FLOAT(dim2)
  title = append_number('Scatter plot')
  name = setup_graph([0.], [0.], num, title = title, psymbol = 1, $
    x_title = title1, y_title = title2, /scatter, win_orig1 = window1, $
    win_orig2 = window2, nblocks_x = nblocks_x, nblocks_y = nblocks_y, $
    xstart = xstart, xstop = xstop, ystart = ystart, ystop = ystop)
  plot_graph, name
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
    IF(first) THEN BEGIN
      sel = umenu(['Resize graph?', $
        'Yes / Resize the graph', 'No / Exit this menu'], title = 0)
      IF(sel EQ 1) THEN BEGIN
        resize_graph, name
        j = EXECUTE('nblocks_x = ' + name + '.nblocks_x')
        j = EXECUTE('nblocks_y = ' + name + '.nblocks_y')
      ENDIF
    ENDIF
  ENDIF
  IF(journal_on) THEN BEGIN
    PRINTF, luj, 'Scatter plot'
    PRINTF, luj, '  1st image:  ' + title1
    PRINTF, luj, '  2nd image:  ' + title2
    PRINTF, luj, '  resulting graph:  ' + title
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
  ENDIF
  first = 0
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


