PRO plot_graph, name, color = color, xstart = xstart, xstop = xstop, $
  xstyle = xstyle, ystart = ystart, ystop = ystop, ystyle = ystyle, $
  mplot = mplot, redraw = redraw, linestyle = linestyle
;+
; NAME:	
;	plot_graph
; PURPOSE:
;       (UIMAGE-specific) To plot a graph in an X-window.
; CATEGORY:
;       Graphics.
; CALLING SEQUENCE:
;       plot_graph,name
; INPUTS:
;       NAME     = Name of a graphics structure (e.g., "GRAPH(1)").
; OUTPUTS:
;       None.
; COMMON BLOCKS:
;	uimage_data,xwindow
; RESTRICTIONS:
;       For use with UIMAGE only.
;#
; PROCEDURE:
;       Extract plot-related info from the structure and then make
;       the plot.
; SUBROUTINES CALLED:
;       None.
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, February 1992.
;  SPR 10290  Dec 03, 92  Use /YNOZERO in calls to PLOT.  J Ewing
;  SPR 10384  Dec 21, 92  Remove XSTYLE keyword in calls to OPLOT.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10953  May 17 93  Add check for different types of LOGFLAG. J Newmark
;  SPR 10958  May 18 93  Add ability to use different plotting linestyles.
;                        J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON history,uimage_version
  COMMON color_values,c_badpix,c_draw,c_scalemin
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    dummy,zoom_index
  COMMON xscreen,magnify,block_usage,scrsiz
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
;
;  Set IG to be the index of the appropriate structure-element in GRAPH.
;  ---------------------------------------------------------------------
  ic = STRPOS(name, '(')
  ig = FIX(STRMID(name, ic + 1, STRPOS(name, ')') - ic - 1))
  IF (graph(ig).hidden EQ 1) THEN RETURN
;
;  Special handling for different versions of UIMAGE.
;  --------------------------------------------------
  IF(uimage_version EQ 4) THEN BEGIN
    IF(NOT KEYWORD_SET(redraw)) THEN BEGIN
      offset = get_location(name)
      base = WIDGET_BASE(TITLE = graph(ig).title, UVALUE = name, $
        XOFFSET = offset(0), YOFFSET = offset(1))
      XPDMENU, get_menu('GRAPH'), base
      draw = WIDGET_DRAW(base, XSIZE = graph(ig).nblocks_x*32*magnify, $
        YSIZE = graph(ig).nblocks_y*32*magnify, YOFFSET= 40)
      WIDGET_CONTROL, /REALIZE, base
      WIDGET_CONTROL, GET_VALUE = window, draw
      graph(ig).window = window
      graph(ig).widget = base
      w = WHERE(block_usage EQ -1)
      block_usage(w) = base
    ENDIF
    WSET, graph(ig).window
  ENDIF ELSE BEGIN
    IF(NOT defined(zoom_index)) THEN zoom_index = 1
    magnify = 2.^zoom_index
  ENDELSE
;
;  Extract various pieces of information from the structure which will
;  be supplied in the call to PLOT.
;  -------------------------------------------------------------------
  IF(NOT KEYWORD_SET(color)) THEN color = graph(ig).color
  x = [0., 0.] & y = [0., 0.]
  cutoff = 0 > (FIX(graph(ig).num) - 1) < 1023
  x = graph(ig).data(0:cutoff, 0)
  y = graph(ig).data(0:cutoff, 1)
  IF(NOT KEYWORD_SET(xstart)) THEN xstart = graph(ig).xstart
  IF(NOT KEYWORD_SET(xstop))  THEN xstop  = graph(ig).xstop
  IF(NOT KEYWORD_SET(ystart)) THEN ystart = graph(ig).ystart
  IF(NOT KEYWORD_SET(ystop))  THEN ystop  = graph(ig).ystop
  psymbol = graph(ig).psymbol
  IF(NOT KEYWORD_SET(xstyle)) THEN $
    IF((graph(ig).topmin NE 0.) OR (graph(ig).topmax NE 0.)) THEN xstyle = 8 $
                                                             ELSE xstyle = 0
  IF(NOT KEYWORD_SET(ystyle)) THEN ystyle = 0
  IF(NOT KEYWORD_SET(linestyle)) THEN linestyle = 0
  IF(NOT KEYWORD_SET(mplot)) THEN label_title = graph(ig).label_title $
                             ELSE label_title = ''
  IF(NOT KEYWORD_SET(mplot)) THEN $
    graph_region, graph(ig).nblocks_x, graph(ig).nblocks_y, magnify, $
      label_title, ypos
    !X.TITLE = graph(ig).x_title
    !Y.TITLE = graph(ig).y_title
;
;  Handle multiple plots recursively.
;  ----------------------------------
  IF(graph(ig).multi NE 0) THEN BEGIN
    IF (graph(ig).logflag EQ 0) THEN $
      PLOT, [-1.e30, -1.e-30], [1.e-30, 1.e-30], /YNOZERO, $
        XRANGE = [xstart, xstop], YRANGE = [ystart, ystop], COLOR = color
    IF (graph(ig).logflag EQ 1) THEN $
      PLOT_IO, [-1.e30, -1.e-30], [1.e-30, 1.e-30], /YNOZERO, $
        XRANGE = [xstart, xstop], YRANGE = [ystart, ystop], COLOR = color
    IF (graph(ig).logflag EQ 2) THEN $
      PLOT_OI, [-1.e30, -1.e-30], [1.e-30, 1.e-30], /YNOZERO, $
        XRANGE = [xstart, xstop], YRANGE = [ystart, ystop], COLOR = color
    IF (graph(ig).logflag EQ 3) THEN $
      PLOT_OO, [-1.e30, -1.e-30], [1.e-30, 1.e-30], /YNOZERO, $
        XRANGE = [xstart, xstop], YRANGE = [ystart, ystop], COLOR = color
    !NOERAS = 1
    FOR i = 0, graph(ig).multi - 1 DO BEGIN
      tempname = get_name(graph(ig).m_id(i))
      plot_graph, tempname, color = graph(ig).m_color(i),  $
        xstart = graph(ig).xstart, xstop = graph(ig).xstop, /mplot, $
        ystart = graph(ig).ystart, ystop = graph(ig).ystop, xstyle = 4, $
        ystyle = 4, linestyle = graph(ig).m_linetype(i)
    ENDFOR
    !NOERAS = 0
    GOTO, made_plot
  ENDIF
;
;  If SCATTER=1, then make a scatter plot.
;  ---------------------------------------
  IF(graph(ig).scatter EQ 1) THEN BEGIN
    name1 = get_name(graph(ig).win_orig1)
    name2 = get_name(graph(ig).win_orig2)
    j = EXECUTE('badpixval1 = ' + name1 + '.badpixval')
    j = EXECUTE('badpixval2 = ' + name2 + '.badpixval')
    j = EXECUTE('good = WHERE((' + name1 + '.data NE badpixval1) AND (' + $
      name2 + '.data NE badpixval2))')
    IF(good(0) EQ -1) THEN BEGIN
      MESSAGE, 'The two images have no overlapping good pixels.', /CONT
      RETURN
    ENDIF
    IF (graph(ig).logflag EQ 0) THEN $
      j = EXECUTE('PLOT, ' + name1 + '.data(good), ' + name2 + $
       '.data(good), PSYM = psymbol, /YNOZERO, XRANGE = [xstart, xstop], ' + $
       'YRANGE = [ystart, ystop], COLOR = color, XSTYLE = xstyle, ' + $
       'YSTYLE = ystyle')
    IF (graph(ig).logflag EQ 1) THEN $
      j = EXECUTE('PLOT_IO, ' + name1 + '.data(good), ' + name2 + $
       '.data(good), PSYM = psymbol, /YNOZERO, XRANGE = [xstart, xstop], ' + $
       'YRANGE = [ystart, ystop], COLOR = color, XSTYLE = xstyle, ' + $
       'YSTYLE = ystyle')
    IF (graph(ig).logflag EQ 2) THEN $
      j = EXECUTE('PLOT_OI, ' + name1 + '.data(good), ' + name2 + $
       '.data(good), PSYM = psymbol, /YNOZERO, XRANGE = [xstart, xstop], ' + $
       'YRANGE = [ystart, ystop], COLOR = color, XSTYLE = xstyle, ' + $
       'YSTYLE = ystyle')
    IF (graph(ig).logflag EQ 3) THEN $
      j = EXECUTE('PLOT_OO, ' + name1 + '.data(good), ' + name2 + $
       '.data(good), PSYM = psymbol, /YNOZERO, XRANGE = [xstart, xstop], ' + $
       'YRANGE = [ystart, ystop], COLOR = color, XSTYLE = xstyle, ' + $
       'YSTYLE = ystyle')
    GOTO, made_plot
  ENDIF
;
;  Determine the X-range and Y-range.
;  ----------------------------------
  IF((xstart EQ 0.) AND (xstop EQ 0.)) THEN BEGIN
    xstart = MIN(x)
    xstop = MAX(x)
  ENDIF
  index_range = WHERE((x GE xstart) AND (x LE xstop))
  IF(index_range(0) EQ -1) THEN BEGIN
    MESSAGE, 'The are no points within the set X range.', /CONT
    RETURN
  ENDIF
  IF((ystart EQ 0.) AND (ystop EQ 0.)) THEN BEGIN
    ystart = MIN(y(index_range))
    ystop = MAX(y(index_range))
    set_yrange = 0
  ENDIF ELSE set_yrange = 1
;
;  Handle a "SPECIAL" plot - one which has gaps in it.
;  ---------------------------------------------------
  IF(graph(ig).special EQ 1) THEN BEGIN
    badpixval = graph(ig).badpixval
    num = graph(ig).num
    good = WHERE(y NE badpixval)
    bad = WHERE(y EQ badpixval)
    IF(bad(0) EQ -1) THEN BEGIN
      IF(graph(ig).logflag EQ 0) $
        THEN PLOT, x, y, PSYM = psymbol, /YNOZERO, XSTYLE = xstyle, $
               YSTYLE = ystyle, XRANGE = [xstart, xstop], $
               YRANGE = [ystart, ystop], COLOR = color, linestyle=linestyle 
      IF(graph(ig).logflag EQ 1) $
        THEN PLOT_IO, x, y ,PSYM = psymbol, XSTYLE = xstyle, YSTYLE = ystyle, $
               XRANGE = [xstart, xstop], YRANGE = [ystart, ystop], $
               COLOR = color, linestyle=linestyle
      IF(graph(ig).logflag EQ 2) $
        THEN PLOT_OI, x, y ,PSYM = psymbol, XSTYLE = xstyle, YSTYLE = ystyle, $
               XRANGE = [xstart, xstop], YRANGE = [ystart, ystop], $
               COLOR = color, linestyle=linestyle
      IF(graph(ig).logflag EQ 3) $
        THEN PLOT_OO, x, y ,PSYM = psymbol, XSTYLE = xstyle, YSTYLE = ystyle, $
               XRANGE = [xstart, xstop], YRANGE = [ystart, ystop], $
               COLOR = color, linestyle=linestyle
    ENDIF ELSE BEGIN
      IF(good(0) NE -1) THEN BEGIN
      minimum = MIN(y(good))
      maximum = MAX(y(good))
      p0 = 0
      WHILE(y(p0) EQ badpixval) DO BEGIN
        p0 = p0 + 1
        IF(p0 GE (num - 1)) THEN BEGIN
          MESSAGE, 'Graphics data consists of all bad values.', /CONT
          RETURN
        ENDIF
      ENDWHILE
      p1 = p0 + 1
      WHILE((y(p1) NE badpixval) AND (p1 LT (num - 1))) DO p1 = p1 + 1
      IF(y(p1) EQ badpixval) THEN p1 = p1 - 1
      IF(set_yrange NE 1) THEN BEGIN
        ystart = minimum
        ystop = maximum
      ENDIF
      IF (graph(ig).logflag EQ 0) THEN $
        PLOT, x(p0:p1), y(p0:p1), XRANGE = [xstart, xstop], $
          YRANGE = [ystart, ystop], XSTYLE = xstyle, YSTYLE = ystyle, $
          COLOR = color, PSYM = psymbol, /YNOZERO, linestyle=linestyle
      IF (graph(ig).logflag EQ 1) THEN $
        PLOT_IO, x(p0:p1), y(p0:p1), XRANGE = [xstart, xstop], $
          YRANGE = [ystart, ystop], XSTYLE = xstyle, YSTYLE = ystyle, $
          COLOR = color, PSYM = psymbol, /YNOZERO, linestyle=linestyle
      IF (graph(ig).logflag EQ 2) THEN $
        PLOT_OI, x(p0:p1), y(p0:p1), XRANGE = [xstart, xstop], $
          YRANGE = [ystart, ystop], XSTYLE = xstyle, YSTYLE = ystyle, $
          COLOR = color, PSYM = psymbol, /YNOZERO, linestyle=linestyle
      IF (graph(ig).logflag EQ 3) THEN $
        PLOT_OO, x(p0:p1), y(p0:p1), XRANGE = [xstart, xstop], $
          YRANGE = [ystart, ystop], XSTYLE = xstyle, YSTYLE = ystyle, $
          COLOR = color, PSYM = psymbol, /YNOZERO, linestyle=linestyle
loop:
      p0 = p1 + 1
      IF(p0 GT (num - 1)) THEN GOTO, made_plot
      WHILE((y(p0) EQ badpixval) AND (p0 LT (num - 1))) DO p0 = p0 + 1
      IF((p0 EQ (num - 1)) AND (y(p0) EQ badpixval)) THEN GOTO, made_plot
      p1 = p0 + 1
      IF(p1 GT (num - 1)) THEN p1 = (num - 1) ELSE $
        WHILE((y(p1) NE badpixval) AND (p1 LT (num - 1))) DO p1 = p1 + 1
      IF(y(p1) EQ badpixval) THEN p1 = p1 - 1
      OPLOT, x(p0:p1), y(p0:p1), COLOR = color, PSYM = psymbol
      IF(p1 EQ (num - 1)) THEN GOTO, made_plot
      GOTO, loop
      ENDIF
    ENDELSE
  ENDIF ELSE BEGIN
;
;  Handle a normal (non-special) plot.
;  -----------------------------------
      IF(graph(ig).logflag EQ 0) $
        THEN PLOT, x, y, PSYM = psymbol, /YNOZERO, XSTYLE = xstyle, $
               YSTYLE = ystyle, XRANGE = [xstart, xstop], $
               YRANGE = [ystart, ystop], COLOR = color, linestyle=linestyle 
      IF(graph(ig).logflag EQ 1) $
        THEN PLOT_IO, x, y ,PSYM = psymbol, XSTYLE = xstyle, YSTYLE = ystyle, $
               XRANGE = [xstart, xstop], YRANGE = [ystart, ystop], $
               COLOR = color, linestyle=linestyle
      IF(graph(ig).logflag EQ 2) $
        THEN PLOT_OI, x, y ,PSYM = psymbol, XSTYLE = xstyle, YSTYLE = ystyle, $
               XRANGE = [xstart, xstop], YRANGE = [ystart, ystop], $
               COLOR = color, linestyle=linestyle
      IF(graph(ig).logflag EQ 3) $
        THEN PLOT_OO, x, y ,PSYM = psymbol, XSTYLE = xstyle, YSTYLE = ystyle, $
               XRANGE = [xstart, xstop], YRANGE = [ystart, ystop], $
               COLOR = color, linestyle=linestyle
  ENDELSE
made_plot:
  !X.TITLE = ''
  !Y.TITLE = ''
;
;  For histogram plots, put a special X-axis at the top of the plot,
;  along with a special label.
;  -----------------------------------------------------------------
  IF((label_title NE '') AND (NOT KEYWORD_SET(mplot))) THEN BEGIN
    IF((graph(ig).topmin NE 0.) OR (graph(ig).topmax NE 0.)) THEN BEGIN
      topmin = graph(ig).topmin
      topmax = topmin + (graph(ig).topmax - topmin)*!X.CRANGE(1) / graph(ig).num
      AXIS, XAXIS = 1, XRANGE = [topmin, topmax], XSTYLE = 1, COLOR = color
    ENDIF
    XYOUTS, /NORM, !X.WINDOW(0) + (!X.WINDOW(1) - !X.WINDOW(0)) / 2., ypos, $
      label_title, COLOR = color, $
      SIZE = graph(ig).label_size, ALIGNMENT = .5
  ENDIF
;new
  IF(NOT KEYWORD_SET(mplot)) THEN $
    set_viewport, 0.150005, 0.955005, 0.110005, 0.945005
;
;  For non-X-window terminals, wait for a key-press before continuing.
;  -------------------------------------------------------------------
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    CASE !D.NAME OF
       'TEK'   : BEGIN
		   IF (NOT KEYWORD_SET(mplot)) THEN i = GET_KBRD(1)
                   text
                   PRINT, ' '
		 END      
       'REGIS' : BEGIN
		   IF (NOT KEYWORD_SET(mplot)) THEN i = GET_KBRD(1)
                   IF (NOT KEYWORD_SET(mplot)) THEN ERASE
                   PRINT, ' '
		 END
       ELSE    :
    ENDCASE
  ENDIF ELSE BEGIN
    IF((uimage_version EQ 4) AND (NOT KEYWORD_SET(redraw))) THEN $
      XMANAGER, 'GRAPH', base, /JUST_REG
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


