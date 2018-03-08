PRO combine_graphs
;+
;  COMBINE_GRAPHS - a UIMAGE-specific routine.  This routine allows users
;  to select various graphs which will be overlaid together.
;#
;  Written by John Ewing, Nov 92.
;  SPR 10291  Dec 01 92  Added an abort option to the menu.  J Ewing
;  SPR 10385  Dec 22 92  Changed handling of XSTART, XSTOP, etc.  J Ewing
;  SPR 10456  Jan 13 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 10958  May 18 93  Added ability to choose plotting linestyle. J Newmark
;-----------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON control,stay_in_menus
  COMMON color_values,c_badpix,c_draw,c_scalemin
  COMMON colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
  IF(NOT defined(journal_on)) THEN journal_on = 0
;
;  Catch the situation in which there are less than 2 graphs hanging around.
;  -------------------------------------------------------------------------
  sz = SIZE(graph)
  IF(sz(2) NE 8) THEN GOTO, toofewgraphs
  IF(sz(3) LT 2) THEN GOTO, toofewgraphs
  GOTO, enoughgraphs
toofewgraphs:
  PRINT, 'There must be at least 2 graphs present in the UIMAGE data ' + $
   'environment'
  PRINT, 'prior to being able to access this operation.'
  RETURN
enoughgraphs:
  nblocks_x = 5
  nblocks_y = 5
  PRINT, underline('List of graphs to be overlaid:')
;
;  Set up a menu in which a user can select graphs to be overlaid.
;  ---------------------------------------------------------------
  menu_title = 'Overlay graphs'
  m_id = INTARR(10)
  list = STRARR(50)
  name_a = STRARR(50)
  list(0) = menu_title
  ic = 1
  sz = SIZE(graph)
  IF(sz(2) EQ 8) THEN BEGIN
    num_graph = sz(3)
    FOR i = 0, num_graph - 1 DO BEGIN
      list(ic) = graph(i).title
      name_a(ic) = 'GRAPH(' + strtrim(string(i), 2) + ')'
      ic = ic + 1
    ENDFOR
  ENDIF
  IF(ic EQ 0) THEN BEGIN
    MESSAGE, 'No graphs are available.', /CONT
    RETURN
  ENDIF
  list(ic) = ''
  name_a(ic) = ''
  ic = ic + 1
  list(ic) = 'HELP'
  name_a(ic) = 'HELP'
  ic = ic + 1
  list(ic) = 'Finished selecting - overlay the graphs'
  name_a(ic) = 'FINISH'
  ic = ic + 1
  list(ic) = 'Abort - return to previous menu'
  name_a(ic) = 'ABORT'
  ic = ic + 1
  list = list(0:(ic - 1))
  name_a = name_a(0:(ic - 1))
  multi = 0
;
;  Show that menu.
;  ---------------
show_menu:
  sel = umenu(list, title = 0)
  IF(name_a(sel) EQ '') THEN GOTO, show_menu
  name = name_a(sel)
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'ABORT') THEN RETURN
  IF(name NE 'FINISH') THEN BEGIN
;
;  Add this graph to the list of graphs to be overlaid.
;  ----------------------------------------------------
    ig = FIX(STRMID(name, 6, STRPOS(name, ')') - 6))
    PRINT, '  ' + graph(ig).title
    IF(graph(ig).multi LE 0) THEN BEGIN
      m_id(multi) = graph(ig).window
      multi = multi + 1
    ENDIF ELSE BEGIN
      FOR j = 0, graph(ig).multi - 1 DO m_id(multi + j) = graph(ig).m_id(j)
      multi = multi + graph(ig).multi
    ENDELSE
    GOTO, show_menu
  ENDIF
  IF(multi LT 1) THEN BEGIN
    PRINT, 'You must select at least 2 graphs prior to selecting ' + $
      '"Make composite graph".'
    multi = 0
    GOTO, show_menu
  ENDIF
  PRINT, '<end of list>'
  title = append_number('Composite graph')
  xstart = 1.e30
  xstop = -1.e30
  ystart = 1.e30
  ystop = -1.e30
  m_color = INTARR(multi)
  m_linetype = INTARR(multi)
  FOR i = 0, multi - 1 DO BEGIN
    tempname = get_name(m_id(i))
    ig = FIX(STRMID(tempname, 6, STRPOS(tempname, ')') - 6))
    m_color(i) = graph(ig).color
    PRINT,'Enter below the plotting linestyle (0-5) for graph',i
    READ,tmpline
    m_linetype(i)=tmpline
;
;  For each graph, determine values for X1, X2, Y1, Y2, which denote
;  the range of the graph in the X and Y dimensions.
;  -----------------------------------------------------------------
    x1 = graph(ig).xstart
    x2 = graph(ig).xstop
    y1 = graph(ig).ystart
    y2 = graph(ig).ystop
;
;  Set the range for the composite graph (XSTART, XSTOP, YSTART, YSTOP)
;  to be the extreme values from the ranges from the individual graphs.
;  --------------------------------------------------------------------
    IF(x1 LT xstart) THEN xstart = x1
    IF(x2 GT xstop)  THEN xstop  = x2
    IF(y1 LT ystart) THEN ystart = y1
    IF(y2 GT ystop)  THEN ystop  = y2
;
;  If all X labels are the same, then set X_TITLE_PROB equal to 0, else
;  set it equal to 1.  Do the same thing for Y_TITLE_PROB.
;  --------------------------------------------------------------------
    IF(i EQ 0) THEN BEGIN
      x_title = graph(ig).x_title
      y_title = graph(ig).y_title
      x_title_prob = 0
      y_title_prob = 0
    ENDIF ELSE BEGIN
      IF(graph(ig).x_title NE x_title) THEN x_title_prob = 1
      IF(graph(ig).y_title NE y_title) THEN y_title_prob = 1
    ENDELSE
  ENDFOR
;
;  If an axis label doesn't match for each graph to be overlaid, then
;  prompt the user for a new axis label.
;  ------------------------------------------------------------------
  IF(x_title_prob EQ 1) THEN BEGIN
    PRINT, 'The selected graphs have different X-axis labels.'
    PRINT, 'Please enter below the X-axis label for the composite graph.'
    READ, x_title
  ENDIF
  IF(y_title_prob EQ 1) THEN BEGIN
    PRINT, 'The selected graphs have different Y-axis labels.'
    PRINT, 'Please enter below the Y-axis label for the composite graph.'
    READ, y_title
  ENDIF
;
;  Instantiate the combined-graph.
;  -------------------------------
  name = setup_graph([0.],[0.],0,title=title,multi=multi,m_color=m_color,$
    m_id=m_id,nblocks_x=nblocks_x,nblocks_y=nblocks_y,xstart=xstart,$
    xstop=xstop,ystart=ystart,ystop=ystop,x_title=x_title,y_title=y_title,$
    m_linetype=m_linetype)
;
;  Plot the graph.
;  ---------------
  plot_graph, name
;
;  Allow X-window users the opportunity to re-size the graph.
;  ----------------------------------------------------------
  IF(!D.NAME EQ 'X' OR !D.NAME EQ 'WIN') THEN BEGIN
    sel = UMENU(['Resize graph?', $
      'Yes / Resize the graph', 'No / Exit this menu'], title = 0)
    IF(sel EQ 1) THEN BEGIN
      resize_graph, name
      j = EXECUTE('nblocks_x = ' + name + '.nblocks_x')
      j = EXECUTE('nblocks_y = ' + name + '.nblocks_y')
    ENDIF
  ENDIF
;
;  If journaling is enabled, then report what occurred to the journal file.
;  ------------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  the following graphs were combined together:'
    FOR i = 0, multi - 1 DO BEGIN
      tempname = get_name(m_id(i))
      ig = FIX(STRMID(tempname, 6, STRPOS(tempname, ')') - 6))
      PRINTF, luj, '    ', graph(ig).title
    ENDFOR
    ig = FIX(STRMID(name, 6, STRPOS(name, ')') - 6))
    PRINTF, luj, '  output:  ' + graph(ig).title
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
  ENDIF
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


