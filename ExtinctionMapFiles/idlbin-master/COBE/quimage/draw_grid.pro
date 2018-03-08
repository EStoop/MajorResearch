PRO draw_grid, name
;+
;  DRAW_GRID - a UIMAGE-specific routine.  This routine allows a grid
;  to be drawn on a window associated with a skycube map, skycube face,
;  or reprojected map.  This routine is UIMAGE's driver to the DRAW_PM
;  routine.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 06 93  Add more comments.  J Ewing
;  SPR 10456  Jan 13 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 10325  Mar 05 93  Coordinate grids on zoomed images.  J Ewing
;  SPR 11105  Jun 28 93  Add choice of entering delta coordinates rather
;                        than absolute coordinates. J. Newmark
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  COMMON image_parms,input_l,input_h,cur_2_face,offx,offy, $
    cube_side, proj, coord, sz_proj
  COMMON zback,zbgr,zdsrcmap,zbsub
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    RETURN
  ENDIF
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = 'Draw coordinate grid'
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Have the user identify the image on which the grid will be drawn.
;  -----------------------------------------------------------------
show_menu:
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
           /face8,/face9,/proj_map,/proj2_map,/zoomed,/help,/exit,$
           title=menu_title)
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
;
;  Extract a few fields out of the operand's structure.
;  ----------------------------------------------------
  face_num = -1
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('window = ' + name + '.window')
  type = STRMID(name, 0, STRPOS(name, '('))
  IF(type EQ 'ZOOMED') THEN BEGIN
    j = EXECUTE('win_orig = ' + name + '.win_orig')
    bflag=0
    IF (win_orig LT 0) THEN BEGIN
      win_orig=ABS(win_orig)
      bflag=1
    ENDIF
    j=EXECUTE('zoomflag = ' + name + '.zoomflag')
    IF (zoomflag NE 0) THEN bflag=1
    name_orig = get_name(win_orig)
    j = EXECUTE('start_x = ' + name + '.start_x')
    j = EXECUTE('start_y = ' + name + '.start_y')
    j = EXECUTE('stop_x = ' + name + '.stop_x')
    j = EXECUTE('stop_y = ' + name + '.stop_y')
    subset = [start_x, stop_x, start_y, stop_y]
    subset = subset * 2.^zoom_index
    j = EXECUTE('specific_zoom = ' + name + '.specific_zoom')
    pos_x = -start_x * 2.^zoom_index
    pos_y = -start_y * 2.^zoom_index
  ENDIF ELSE BEGIN
    j = EXECUTE('pos_x = ' + name + '.pos_x')
    j = EXECUTE('pos_y = ' + name + '.pos_y')
    name_orig = name
  ENDELSE
  j = EXECUTE('data = ' + name_orig + '.data')
  IF(type EQ 'ZOOMED') THEN BEGIN
    IF (bflag EQ 1) THEN BEGIN
      IF (zoomflag EQ 1) THEN data(start_x:stop_x,start_y:stop_y)=zbgr
      IF (zoomflag EQ 2) THEN data(start_x:stop_x,start_y:stop_y)=zdsrcmap
      IF (zoomflag EQ 3) THEN data(start_x:stop_x,start_y:stop_y)=zbsub
    ENDIF
  ENDIF
  j = EXECUTE('projection = ' + name_orig + '.projection')
  j = EXECUTE('coordinate_system = ' + name_orig + '.coordinate_system')
  w_pos = [pos_x,  pos_y]
  sz = SIZE(data)
  dim1 = sz(1)
  dim2 = sz(2)
  input_l = 2.^zoom_index * dim1
  input_h = 2.^zoom_index * dim2
  CASE projection OF
    'AITOFF'            : proj = 'A'
    'GLOBAL SINUSOIDAL' : proj = 'S'
    'MOLLWEIDE'         : proj = 'M'
     ELSE :
  ENDCASE
  CASE coordinate_system OF
    'ECLIPTIC'   : coord = 'E'
    'EQUATORIAL' : coord = 'Q'
    'GALACTIC'   : coord = 'G'
    ELSE:
  ENDCASE
;
;  Have the user select the grid coordinate type.
;  ----------------------------------------------
  menu = ['Grid coordinate type:', 'Ecliptic', 'Equatorial', 'Galactic']
  sel = one_line_menu(menu)
  gcoord = STRMID('EQG', sel - 1, 1)
;
;  Determine what parallels will be drawn.
;  ---------------------------------------
  para = [0, 30, 60]
  para_inc = '30'
  PRINT, 'The default parallels are at ' + bold('0, 30, & 60') + ' degrees.'
get_parallels:
  menu = ['Use default parallels or set new values?', 'Use default', $
    'Set new values']
  sel = one_line_menu(menu)
  IF(sel EQ 2) THEN BEGIN
    inmenu=['Enter specific values or increments','Enter specifics', $
    'Enter increments']
    sel=one_line_menu(inmenu)
    IF (sel EQ 1) THEN BEGIN
       para_str = '0,30,60'
       PRINT, 'The default parallels are at ' + bold('0, 30, & 60')+' degrees.'
         PRINT, 'Enter parallels separated by commas (e.g. ' + $
           bold('0,15,30,45,60,75') + ')'
         READ, para_str
         IF(para_str EQ '') THEN BEGIN
           para = -1
         ENDIF ELSE BEGIN
           para = [-1]
           str = para_str + ','
           ic = STRPOS(str, ',')
           WHILE(ic NE -1) DO BEGIN
             sub = STRMID(str, 0, ic)
             str = STRMID(str, ic + 1, STRLEN(str))
             IF(validnum(sub) EQ 0) THEN BEGIN
               PRINT, 'Invalid numerical value.'
               GOTO, get_parallels
             ENDIF
             value = FLOAT(sub)
             IF((value LT 0.) OR (value GT 90.)) THEN BEGIN
               PRINT, 'A value is outside the allowable range [0..90].'
               GOTO, get_parallels
             ENDIF
             para = [para, value]
             ic = STRPOS(str, ',')
           ENDWHILE
           para = para(1:*)
         ENDELSE
    ENDIF ELSE BEGIN
     READ,'Enter new increment (degrees) for parallels: ',para_inc
     IF(para_inc EQ '') THEN BEGIN
      para = -1
     ENDIF ELSE BEGIN
      para_flt= FLOAT(para_inc)
      step= FIX(90 / para_flt)
      para=FLTARR(step+1)
      index=0
      FOR i=0,step DO para(i)=i*para_flt
     ENDELSE
    ENDELSE
  ENDIF
;
;  Determine what medidians will be drawn.
;  ---------------------------------------
  merd = [0, 60, 120, 180]
  merd_inc = '60'
  PRINT, 'The default meridians are at ' + bold('0, 60, 120, & 180') + $
    ' degrees.'
get_meridians:
  menu = ['Use default meridians or set new values?', 'Use default', $
    'Set new values']
  sel = one_line_menu(menu)
  IF(sel EQ 2) THEN BEGIN
    inmenu=['Enter specific values or increments','Enter specifics', $
    'Enter increments']
    sel=one_line_menu(inmenu)
    IF (sel EQ 1) THEN BEGIN
       merd_str = '0,60,120,180'
       PRINT, 'The default meridians are at ' + bold('0, 60, 120, & 180') + $
           ' degrees.'
         PRINT, 'Enter meridians separated by commas (e.g. ' + $
           bold('0,90,180') + ')'
         READ, merd_str
         IF(merd_str EQ '') THEN BEGIN
           merd = -1
         ENDIF ELSE BEGIN
           merd = [-1]
           str = merd_str + ','
           ic = STRPOS(str, ',')
           WHILE(ic NE -1) DO BEGIN
             sub = STRMID(str, 0, ic)
             str = STRMID(str, ic + 1, STRLEN(str))
             IF(validnum(sub) EQ 0) THEN BEGIN
               PRINT, 'Invalid numerical value.'
               GOTO, get_meridians
             ENDIF
             value = FLOAT(sub)
             IF((value LT 0.) OR (value GT 180.)) THEN BEGIN
               PRINT, 'A value is outside the allowable range [0..180].'
               GOTO, get_meridians
             ENDIF
             merd = [merd, value]
             ic = STRPOS(str, ',')
           ENDWHILE
           merd = merd(1:*)
         ENDELSE
    ENDIF ELSE BEGIN
     READ,'Enter new increment (degrees) for meridians: ',merd_inc
     IF(merd_inc EQ '') THEN BEGIN
      merd = -1
     ENDIF ELSE BEGIN
      merd_flt= FLOAT(merd_inc)
      step= FIX(180 / merd_flt)
      merd=FLTARR(step+1)
      index=0
      FOR i=0,step DO merd(i)=i*merd_flt
     ENDELSE
    ENDELSE
  ENDIF
;
;  Get ready for the call to DRAW_PM by defining its arguments and by
;  seting up the IMAGE_PARMS common block.
;  ------------------------------------------------------------------
  IF(STRMID(name_orig, 0, 1) EQ 'P') THEN BEGIN
    re_proj = 1
    sz_proj = SIZE(data)
    sz_proj(1) = sz_proj(1) * 2.^zoom_index
    sz_proj(2) = sz_proj(2) * 2.^zoom_index
    res_input = 0
  ENDIF ELSE BEGIN
    re_proj = 0
    ic = STRPOS(name_orig, '(')
    res_input = FIX(STRMID(name_orig, ic - 1, 1)) + zoom_index 
    cube_side = 2.^zoom_index * dim2/3.
  ENDELSE
  IF(STRMID(name_orig, 0, 1) EQ 'F') THEN BEGIN
    j = EXECUTE('face_num = ' + name_orig + '.faceno')
    offx = [0,0,0,0,0,0]
    offy = [0,0,0,0,0,0]
    cube_side = 2.^zoom_index * dim2
  ENDIF ELSE BEGIN
    offx = [0,0,1,2,3,0]
    offy = [2,1,1,1,1,0]
  ENDELSE
  WSET, window
  WSHOW, window
;
;  Call DRAW_PM so that the grid will be drawn.
;  --------------------------------------------
  exstr = 'draw_pm, para, merd, re_proj, face_num, w_pos, ' + $
    'res_input, gcoord, color = 255, subset = subset, zoom = specific_zoom'
  j = EXECUTE(exstr)
;
;  Handle journaling if appropriate.
;  ---------------------------------
  IF(journal_on EQ 1) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:   ' + title
    PRINTF, luj, '  parallels: ' 
    PRINTF, luj, STRTRIM(para,2)
    PRINTF, luj, '  meridians: ' 
    PRINTF, luj, STRTRIM(merd,2)
    PRINTF, luj, '  the following statement was executed:'
    PRINTF, luj, '  ' + exstr
    PRINTF, luj, '  where:'
    PRINTF, luj, '    para is a 1-D array containing the parallels'
    PRINTF, luj, '    merd is a 1-D array containing the meridians'
    PRINTF, luj, '    re_proj   = ' + STRTRIM(re_proj, 2)
    PRINTF, luj, '    face_num  = ' + STRTRIM(face_num, 2)
    PRINTF, luj, '    w_pos     = ' , w_pos
    PRINTF, luj, '    res_input = ' + STRTRIM(res_input, 2)
    PRINTF, luj, '    gcoord    = ' + gcoord
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
  ENDIF
  subset = 0
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


