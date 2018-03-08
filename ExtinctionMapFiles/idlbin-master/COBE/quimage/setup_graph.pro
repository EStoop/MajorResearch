FUNCTION setup_graph, x, y, num, $
  title = title, $
  x_title = x_title, $
  y_title = y_title, $
  nblocks_x = nblocks_x, $
  nblocks_y = nblocks_y, $
  logflag = logflag, $
  psymbol = psymbol, $
  topmin = topmin, $
  topmax = topmax, $
  label_title = label_title, $
  label_size = label_size, $
  badpixval = badpixval, $
  linkweight = linkweight, $
  hidden = hidden, $
  special = special, $
  scatter = scatter, $
  win_orig1 = win_orig1, $
  win_orig2 = win_orig2, $
  xstart = xstart, $
  xstop = xstop, $
  ystart = ystart, $
  ystop = ystop, $
  multi = multi, $
  m_color = m_color, $
  m_id = m_id, $
  color = color, $
  m_linetype = m_linetype
;+
; NAME:	
;	setup_graph
; PURPOSE:
;       To initialize a UIMAGE graphics structure and, if an X-window
;       terminal is being used, to open a window of a certain size and
;       to position it on an X-window screen in accordance with the UIMAGE
;       X-window screen management system.
; CALLING SEQUENCE:
;       name = setup_graph( x, y, num, <optional keywords>)
; INPUTS:
;       X        = a 1-D array of X values (up to size 1024)
;       Y        = a 1-D array of Y values (up to size 1024)
;       NUM      = the number of elements in X and Y  (0 < NUM < 1025)
; OUTPUTS:
;       NAME     = A string which references where a structure
;       associated with a graph is stored (e.g., "GRAPH(0)").
;#
; COMMON BLOCKS:
;       uimage_data,xwindow
; RESTRICTIONS:
;       To be called from UIMAGE only.
; PROCEDURE:
;       Insert a new structure/element into the array-of-structures
;       named GRAPH.  Open a new window on the X-window screen for
;       the new graph.
; SUBROUTINES CALLED:
;       make_window, w_available
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, February 1992.
;  SPR 10385  Dec 22 92  Special handling for XSTART, XSTOP, etc.  J Ewing
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10829  Apr  1993  Change to add .HIDDEN field
;  SPR 10843  Apr 22 93  Add .LINKWEIGHT field
;  SPR 10958  May 18 93  Add M_LINETYPE to call. J Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON color_values,c_badpix,c_draw,c_scalemin
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(uimage_version EQ 4) THEN BEGIN
    nblocks_x = 4
    nblocks_y = 3
    c_draw = 1
  ENDIF
  IF(N_PARAMS(0) NE 3) THEN BEGIN
    MESSAGE, 'Three arguments must be supplied.', /CONT
    RETURN, 'UNDEFINED'
  ENDIF
  IF(NOT defined(c_draw)) THEN c_draw = 1
  w_available, window
  IF(NOT KEYWORD_SET(title)) THEN title='Object #'+strtrim(string(window+1),2)
;
;  Create a new instance of the graph-structure and store it in GRAPH. 
;  If GRAPH does not exist, then set it equal to this new instance,
;  else append the new instance onto the end of GRAPH.
;  -------------------------------------------------------------------
  sz = SIZE(graph)
  IF(sz(2) NE 8) THEN BEGIN
    ic = 0
    graph = {gr_struct, data: FLTARR(1024,2), num: 0., window: 0, $
             title: '', x_title: '', y_title: '', logflag: 0, psymbol: 0, $
             nblocks_x: 0, nblocks_y: 0, topmin: 0., topmax: 0., $
             label_title: '', label_size: 0., badpixval: 0., special: 0, $
             hidden: 0, scatter: 0, win_orig1: 0, win_orig2: 0, xstart: 0., $
             xstop: 0., ystart: 0., ystop: 0., multi: 0, m_id: INTARR(10), $
             m_color: INTARR(10), m_linetype: INTARR(10), color: 0, $
             linkweight: 0, widget : 0L}
  ENDIF ELSE BEGIN
    ic = sz(3)
    graph = [graph, {gr_struct}]
  ENDELSE
;
;  Set various fields within the new instance of the structure.
;  ------------------------------------------------------------
  graph(ic).window = window
  graph(ic).title = title
  graph(ic).num = num
  data = FLTARR(1024, 2)
  IF((NOT KEYWORD_SET(scatter)) AND (NOT KEYWORD_SET(multi))) THEN BEGIN
    data(0, 0) = x(0 : (num - 1))
    data(0, 1) = y(0 : (num - 1))
  ENDIF
  graph(ic).data = data
  arr_elem = 'GRAPH(' + STRTRIM(STRING(ic), 2) + ')'
  IF(KEYWORD_SET(x_title)) THEN j = EXECUTE(arr_elem + '.x_title = x_title')
  IF(KEYWORD_SET(y_title)) THEN j = EXECUTE(arr_elem + '.y_title = y_title')
  IF(NOT KEYWORD_SET(nblocks_x)) THEN BEGIN
    nblocks_x = 5
  ENDIF
  IF(NOT KEYWORD_SET(nblocks_y)) THEN BEGIN
    nblocks_y = 5
  ENDIF
;
;  Determine suggested axis ranges.
;  --------------------------------
  IF(KEYWORD_SET(special)) THEN BEGIN
    good = WHERE(y NE badpixval, ngood)
    IF(ngood GT 0) THEN BEGIN
      sug_x1 = MIN(x(good), MAX = sug_x2)
      sug_y1 = MIN(y(good), MAX = sug_y2)
    ENDIF
  ENDIF ELSE BEGIN
    sug_x1 = MIN(x, MAX = sug_x2)
    sug_y1 = MIN(y, MAX = sug_y2)
  ENDELSE
;
;  Continue to various fields within the new instance of the structure.
;  --------------------------------------------------------------------
  j=EXECUTE(arr_elem+'.nblocks_x=nblocks_x')
  j=EXECUTE(arr_elem+'.nblocks_y=nblocks_y')
  IF(KEYWORD_SET(logflag))     THEN j=EXECUTE(arr_elem+'.logflag=logflag')
  IF(KEYWORD_SET(psymbol))     THEN j=EXECUTE(arr_elem+'.psymbol=psymbol')
  IF(KEYWORD_SET(topmin))      THEN j=EXECUTE(arr_elem+'.topmin=topmin')
  IF(KEYWORD_SET(topmax))      THEN j=EXECUTE(arr_elem+'.topmax=topmax')
  IF(KEYWORD_SET(label_title)) THEN j=EXECUTE(arr_elem+'.label_title=label_title')
  IF(KEYWORD_SET(label_size))  THEN j=EXECUTE(arr_elem+'.label_size=label_size')
  IF(KEYWORD_SET(badpixval))   THEN j=EXECUTE(arr_elem+'.badpixval=badpixval')
  IF(KEYWORD_SET(hidden))      THEN j=EXECUTE(arr_elem+'.hidden=hidden') $
                               ELSE j=EXECUTE(arr_elem+'.hidden=-1') 
  IF(KEYWORD_SET(linkweight))THEN j=EXECUTE(arr_elem+'.linkweight=linkweight') $
                               ELSE j=EXECUTE(arr_elem+'.linkweight=-1') 
  IF(KEYWORD_SET(special))     THEN j=EXECUTE(arr_elem+'.special=special')
  IF(KEYWORD_SET(scatter))     THEN j=EXECUTE(arr_elem+'.scatter=scatter')
  IF(KEYWORD_SET(win_orig1))   THEN j=EXECUTE(arr_elem+'.win_orig1=win_orig1')
  IF(KEYWORD_SET(win_orig2))   THEN j=EXECUTE(arr_elem+'.win_orig2=win_orig2')
  IF(KEYWORD_SET(xstart))      THEN j=EXECUTE(arr_elem+'.xstart=xstart') $
                               ELSE j=EXECUTE(arr_elem+'.xstart=sug_x1')
  IF(KEYWORD_SET(xstop))       THEN j=EXECUTE(arr_elem+'.xstop=xstop') $
                               ELSE j=EXECUTE(arr_elem+'.xstop =sug_x2')
  IF(KEYWORD_SET(ystart))      THEN j=EXECUTE(arr_elem+'.ystart=ystart') $
                               ELSE j=EXECUTE(arr_elem+'.ystart=sug_y1')
  IF(KEYWORD_SET(ystop))       THEN j=EXECUTE(arr_elem+'.ystop=ystop') $
                               ELSE j=EXECUTE(arr_elem+'.ystop =sug_y2')
  IF(KEYWORD_SET(multi))       THEN j=EXECUTE(arr_elem+'.multi=multi')
  IF(KEYWORD_SET(m_color))     THEN j=EXECUTE(arr_elem+'.m_color=m_color')
  IF(KEYWORD_SET(m_id))        THEN j=EXECUTE(arr_elem+'.m_id=m_id')
  IF(KEYWORD_SET(m_linetype))  THEN j=EXECUTE(arr_elem+'.m_linetype=m_linetype')
  IF(KEYWORD_SET(color))       THEN j=EXECUTE(arr_elem+'.color=color') $
                       ELSE IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) $
                                 THEN j=EXECUTE(arr_elem+'.color=c_draw') $
                                 ELSE j=EXECUTE(arr_elem+'.color=255')
;
;  If running on X-windows and depending on the value of UIMAGE_VERSION,
;  make a window for the graph to be plotted in.
;  ---------------------------------------------------------------------
  IF((uimage_version EQ 2) OR (uimage_version LT 0)) THEN BEGIN
    IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
      j=EXECUTE('hidden = ' +arr_elem + '.hidden')
      IF (hidden EQ 1) THEN RETURN,arr_elem
      IF(NOT defined(zoom_index)) THEN BEGIN
        MESSAGE, 'ZOOM_INDEX is undefined.', /CONT
        RETURN, -1
      ENDIF
      zoom = 2. ^ zoom_index
      make_window, window, nblocks_x * zoom, nblocks_y * zoom, title
    ENDIF
  ENDIF
;
;  Return the string which identifies the new instance of the graph-structure,
;  e.g., "GRAPH(1)".
;  ---------------------------------------------------------------------------
  RETURN, arr_elem
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


