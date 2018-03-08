PRO uskycut, name_sel
;+
;  USKYCUT - a UIMAGE-specific routine.  This is UIMAGE's driver for the
;  SKYCUT procedure.
;#
;  Written by J Ewing & J Gales
;  SPR 10385  Dec 22 92  Set XSTART, XSTOP, etc for graph-object.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;  SPR 10645  May 02 94  Remove window message for lon-lat input. J. Newmark
;--------------------------------------------------------------------------
  COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
    cube_side, proj, coord, sz_proj
  COMMON arc, x_plot,y_plot,win_title,fij_input,n_pix_arc,i_arc,j_arc, $
    vec_arc,end_pixs
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON nonxwindow,term_type
  COMMON journal,journal_on,luj
  COMMON color_values,c_badpix,c_draw,c_scalemin
  COMMON history,uimage_version
  COMMON zback,zbgr,zdsrcmap,zbsub
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  first = 1
  arc_type = 'S'
  noinpt = 0
  n_plot = 0
  nblocks_x = 5
  nblocks_y = 5
  color = 1
;
;  If the name of the desired operand was not supplied when this routine
;  was invoked, then put up a menu in which the user can identify the
;  desired image.
;  ---------------------------------------------------------------------
  menu_title = 'Cross sections, Sky cuts'
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
show_topmenu:
  mask = 0
  name_sel = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
           /face8,/face9,/proj_map,/proj2_map,/zoomed,/help,$
           /exit,title=menu_title)
  IF(name_sel EQ 'EXIT') THEN BEGIN
    IF(journal_on AND (first NE 1)) THEN $
      PRINTF, luj, '----------------------------------------' + $
                   '--------------------------------------'
    RETURN
  ENDIF
  IF(name_sel EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_topmenu
  ENDIF
  IF(name_sel EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no objects available.', /CONT
    RETURN
  ENDIF
have_name:
;
;  At this point, NAME_SEL identifies the image on which the user will
;  mark the two arc-endpoints.
;  -------------------------------------------------------------------
  j = EXECUTE('title_sel = ' + name_sel + '.title')
  IF(journal_on) THEN BEGIN
    IF(first EQ 1) THEN PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title_sel
  ENDIF
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN zoom = 2.^zoom_index $
    ELSE zoom = 1
;
;  Define NAME_BASE to be equal to NAME_SEL if NAME_SEL does not identify
;  a zoomed image, else set NAME_BASE so it identifies the image which
;  had been zoomed in on.
;  ----------------------------------------------------------------------
  type_sel = STRMID(name_sel, 0, STRPOS(name_sel, '('))
  IF(type_sel EQ 'ZOOMED') THEN BEGIN
    j = EXECUTE('win_orig = ' + name_sel + '.win_orig')
    bflag=0
    IF (win_orig LT 0) THEN BEGIN
      win_orig=ABS(win_orig)
      bflag=1
    ENDIF
    j=EXECUTE('zoomflag = ' + name_sel + '.zoomflag')
    IF (zoomflag NE 0) THEN bflag=1
    name_base = get_name(win_orig)
  ENDIF ELSE name_base = name_sel

  type_base = STRMID(name_base, 0, STRPOS(name_base, '('))
  j = EXECUTE('coordinate = ' + name_base + '.coordinate_system')
  coord = STRMID(coordinate, 0, 1)
  IF(coordinate EQ 'EQUATORIAL') THEN coor = 'Q'
;
;  If NAME_BASE identifies a projected image, then set NAME_ORIG so it
;  identifies the image which had been reprojected.
;  -------------------------------------------------------------------
  IF((type_base EQ 'PROJ_MAP') OR (type_base EQ 'PROJ2_MAP')) THEN BEGIN
    j = EXECUTE('projection = ' + name_base + '.projection')
    CASE projection OF
      'AITOFF'            : proj = 'A'
      'GLOBAL SINUSOIDAL' : proj = 'S'
      'MOLLWEIDE'         : proj = 'M'
       ELSE               : BEGIN
                              PRINT, 'Unsupported projection.'
                              GOTO, show_topmenu
                            END
    ENDCASE
    p_coord = coord
    p_type = proj
    j = EXECUTE('win_orig =' + name_base + '.win_orig')
    name_orig = get_name(win_orig)
    type_orig = STRMID(name_orig, 0, STRPOS(name_orig, '('))
    IF(type_orig EQ 'FACE') THEN $
      j = EXECUTE('face_num = ' + name_orig + '.faceno')
    j = EXECUTE('nat_cube = ' + name_orig + '.data')
    re_proj = 1
  ENDIF ELSE re_proj = 0
  j = EXECUTE('img_win = ' + name_sel + '.window')
  j = EXECUTE('orient = ' + name_base + '.orient')
  j = EXECUTE('pos_x = ' + name_base + '.pos_x')
  j = EXECUTE('pos_y = ' + name_base + '.pos_y')
  j = EXECUTE('badpixval = ' + name_base + '.badpixval')
  j = EXECUTE('scale_min = ' + name_base + '.scale_min')
  j = EXECUTE('scale_max = ' + name_base + '.scale_max')
  j = EXECUTE('y_title = ' + name_base + '.units')
  j = EXECUTE('input = ' + name_base + '.data')
  IF(type_sel EQ 'ZOOMED') THEN BEGIN
    j = EXECUTE('start_x = ' + name_sel + '.start_x')
    j = EXECUTE('start_y = ' + name_sel + '.start_y')
    j = EXECUTE('stop_x = ' + name_sel + '.stop_x')
    j = EXECUTE('stop_y = ' + name_sel + '.stop_y')
    j = EXECUTE('specific_zoom = ' + name_sel + '.specific_zoom')
    pos_x = -start_x * zoom
    pos_y = -start_y * zoom
    IF (bflag EQ 1) THEN BEGIN
      IF (zoomflag EQ 1) THEN input(start_x:stop_x,start_y:stop_y)=zbgr
      IF (zoomflag EQ 2) THEN input(start_x:stop_x,start_y:stop_y)=zdsrcmap
      IF (zoomflag EQ 3) THEN input(start_x:stop_x,start_y:stop_y)=zbsub
    ENDIF
  ENDIF ELSE specific_zoom = 1.
  w_pos = [pos_x, pos_y]
  sz = SIZE(input)
  dim1 = sz(1)
  dim2 = sz(2)
  input_l = zoom * dim1
  input_h = zoom * dim2
  IF((type_base EQ 'PROJ_MAP') OR (type_base EQ 'PROJ2_MAP'))THEN BEGIN
    re_proj = 1
    sz_proj = SIZE(input)
    sz_proj(1) = sz_proj(1) * zoom
    sz_proj(2) = sz_proj(2) * zoom
    res_input = 0
  ENDIF ELSE BEGIN
    re_proj = 0
    ic = STRPOS(name_base, '(')
    res_input = FIX(STRMID(name_base, ic-1, 1))
    IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) $
      THEN res_input = res_input + zoom_index
    cube_side = zoom * dim2 / 3.
  ENDELSE
  face_num = -1
  IF(STRMID(type_base, 0, 4) EQ 'FACE') THEN BEGIN
    j = EXECUTE('face_num = ' + name_base + '.faceno')
    offx = [0,0,0,0,0,0]
    offy = [0,0,0,0,0,0]
    cube_side = zoom * dim2
  ENDIF ELSE BEGIN
    offx = [0,0,1,2,3,0]
    offy = [2,1,1,1,1,0]
  ENDELSE
  input = CONGRID(input, sz(1)*zoom, sz(2)*zoom)
  IF(type_base EQ 'FACE') THEN j = EXECUTE('face_num = '+name_base+'.faceno')
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN noshow = 1 ELSE noshow = 0
  nostore = 1
;
;  Preparatory set-up is completed, now call Joel's primitives.
;  ------------------------------------------------------------
  ret_stat = getinpar(input, nat_cube, res, re_proj, psize, mask, face_num)
e_mode:
  arc_type = 'S'
  getemode, entry_mode, noshow, coord, infmt, re_proj
loop:
  IF(noshow EQ 0) THEN BEGIN
    IF (entry_mode eq 1) THEN $
      PRINT, 'Please mark points on the window titled:  ' + bold(title_sel)
    WSET, img_win
    WSHOW, img_win
  ENDIF
  ret_stat = getarc(entry_mode, input, nat_cube, mask, res, face_num, w_pos, $
    re_proj, arc_type, lon_flag, lat_flag, noinpt, color, zoom = specific_zoom)
  IF (ret_stat EQ 1) THEN GOTO, loop
  IF(noshow EQ 0) THEN BEGIN
    IF(first NE 1) THEN refresh_image, name_sel
    drawcirc, fij_input, n_pix_arc, i_arc, j_arc, vec_arc, $
      re_proj, w_pos, face_num, color, zoom = specific_zoom
  ENDIF
;
;  Catch the case in which the Y values for the plot are all zeros.
;  ----------------------------------------------------------------
  good = WHERE(y_plot NE badpixval, ngood)
  IF(ngood LE 0) THEN BEGIN
    PRINT, 'A plot could not be made for the selected points.'
    GOTO, wish
  ENDIF
  xstart = MIN(x_plot, MAX = xstop)
  ystart = MIN(y_plot(good), MAX = ystop)
;
;  Instantiate the graph-object, plot the graph, & allow it to be resized.
;  -----------------------------------------------------------------------
  win_title=append_number(win_title)
  gname = setup_graph(x_plot, y_plot, N_ELEMENTS(y_plot), $
    title = win_title, bad = badpixval, x_title = 'Degrees along arc', $
    y_title = y_title, psym = -3, nblocks_x = nblocks_x, nblocks_y = $
    nblocks_y, xstart = xstart, xstop = xstop, ystart = ystart, ystop = $
    ystop, /special)
  plot_graph, gname
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
    WSET, img_win
    IF(first EQ 1) THEN BEGIN
      sel = UMENU(['Resize graph?', $
        'Yes / Resize the graph', 'No / Exit this menu'], title = 0)
      IF(sel EQ 1) THEN BEGIN
        resize_graph, gname
        j = EXECUTE('nblocks_x = ' + gname + '.nblocks_x')
        j = EXECUTE('nblocks_y = ' + gname + '.nblocks_y')
      ENDIF
    ENDIF
  ENDIF
  IF(journal_on) THEN PRINTF, luj, '  output graph:  ' + win_title
  first = 0
;
;  A cross-cut was performed.  See what the user wants to do next.
;  ---------------------------------------------------------------
wish:
  wish_title = 'Next action for C.S. plot'
  wish = umenu([wish_title, $
                'Specify new cross-section', $
                'View "complementary" cross-section', $
                'Change entry mode', $
                '', $
                'HELP',$
                'Return to previous menu'], title = 0, init = 1)
  CASE wish OF
    1 : arc_type = 'S'
    2 : BEGIN
          IF((lat_flag EQ 0) AND (lon_flag EQ 0) AND $
            (face_num EQ -1) AND (type_sel NE 'ZOOMED')) THEN BEGIN
            IF (arc_type EQ 'S') THEN arc_type = 'L' $
                                 ELSE arc_type = 'S'
            noinpt = 1
          ENDIF ELSE BEGIN
            PRINT, $
              'This operation not available for single faces, zoomed '+$
              'images or latitude cuts.'
            GOTO, wish
          ENDELSE
        END
    3 : BEGIN
          IF(noshow EQ 1) THEN BEGIN 
            PRINT, $
              'Cursor entry mode not allowed for this invocation or terminal.'
            GOTO, wish
          ENDIF ELSE GOTO, e_mode
        END
    5 : BEGIN
          uimage_help, wish_title
          GOTO, wish
        END
    6 : IF(uimage_version EQ 2) THEN GOTO, show_topmenu
  ENDCASE
  GOTO, loop
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


