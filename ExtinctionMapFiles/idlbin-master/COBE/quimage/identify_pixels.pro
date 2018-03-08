PRO identify_pixels, name
;+
;  IDENTIFY_PIXELS - a UIMAGE-specific routine.  This routine allows an
;  X-window user to mark points on a selected image, and to have the
;  coordinates and data values for the selected points be printed out
;  in a table.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10405  Jan 11 93  Handle points beyond the image.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 11104  Jun 28 93  Add ability to accept keyboard selection of pixels.
;                         J. Newmark
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  COMMON zback,zbgr,zdsrcmap,zbsub
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  keymenu=['Select input format','Mouse','Keyboard']
  sel=one_line_menu(keymenu,init=1)
  IF (sel EQ 2) THEN GOTO, keyboard
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on X-windows.', /CONT
    RETURN
  ENDIF
  menu_title = 'Single pixel information'
  IF(N_PARAMS(0) GE 1) THEN GOTO, have_name
show_menu:
;
;  Put up a menu in which the user can select the image he wants to
;  work with.
;  ----------------------------------------------------------------
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
  first_letter = STRMID(name, 0, 1)
show_menu2:
;
;  Put up a menu in which the user can select the corrdinate system to
;  be used for the outputted pairs of {lon,lat} which are derived for
;  the selected pixels.
;  -------------------------------------------------------------------
  menu_title2 = 'Select output coordinate system'
  menu = [menu_title2, 'Ecliptic', 'Equatorial', 'Galactic', $
    ' ', 'HELP', 'Return to previous menu']
  sel = umenu(menu, title = 0)
  CASE sel OF
    1: outco = 'E'
    2: outco = 'Q'
    3: outco = 'G'
    5: BEGIN
         uimage_help, menu_title2
         GOTO, show_menu2
       END
    6: GOTO, show_menu
    ELSE: GOTO, show_menu2
  ENDCASE
  outcoor = STRUPCASE(menu(sel))
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('window = ' + name + '.window')
  IF(first_letter EQ 'F') THEN j = EXECUTE('faceno = ' + name + '.faceno') $
                          ELSE faceno = -1
  ic = STRPOS(name, '(')
;
;  Handle complexity for different types of data-objects.
;  ------------------------------------------------------
  stype = STRUPCASE(STRMID(name, 0, ic))
  CASE stype OF
    'ZOOMED' : BEGIN
       j = EXECUTE('win_orig = ' + name + '.win_orig')
       bflag=0
       IF (win_orig LT 0) THEN BEGIN
         win_orig=ABS(win_orig)
         bflag=1
       ENDIF
       j=EXECUTE('zoomflag = ' + name + '.zoomflag')
       IF (zoomflag NE 0) THEN bflag=1
       name_orig = get_name(win_orig)
       IF(name_orig EQ 'UNDEFINED') THEN BEGIN
         MESSAGE, $
         'The user has deleted the data associated with a zoom window.', /CONT
         GOTO, show_menu
       ENDIF
       IF(STRMID(name_orig, 0, 1) EQ 'F') THEN $
         j = EXECUTE('faceno = ' + name_orig + '.faceno')
       j = EXECUTE('start_x = ' + name + '.start_x')
       j = EXECUTE('start_y = ' + name + '.start_y')
       j = EXECUTE('stop_x = ' + name + '.stop_x')
       j = EXECUTE('stop_y = ' + name + '.stop_y')
       j = EXECUTE('specific_zoom = ' + name + '.specific_zoom')
       j = EXECUTE('orient = ' + name_orig + '.orient')
       j = EXECUTE('data = ' + name_orig + '.data')
       IF (bflag EQ 1) THEN BEGIN
          IF (zoomflag EQ 1) THEN data(start_x:stop_x,start_y:stop_y)=zbgr
          IF (zoomflag EQ 2) THEN data(start_x:stop_x,start_y:stop_y)=zdsrcmap
          IF (zoomflag EQ 3) THEN data(start_x:stop_x,start_y:stop_y)=zbsub
       ENDIF
       sz = SIZE(data)
       dim1_orig = sz(1)
       dim2_orig = sz(2)
       x0 = start_x > 0
       x1 = stop_x < (dim1_orig-1)
       y0 = start_y > 0
       y1 = stop_y < (dim2_orig-1)
       zoom = FIX(2.^zoom_index)*specific_zoom
       IF(start_x GE 0) THEN pos_x = 0 ELSE pos_x = -start_x*zoom
       IF(start_y GE 0) THEN pos_y = 0 ELSE pos_y = -start_y*zoom
       ic = STRPOS(name_orig, '(')
       dtype = STRUPCASE(STRMID(name_orig, 0, ic))
       IF((dtype EQ 'PROJ_MAP') OR (dtype EQ 'PROJ2_MAP')) THEN BEGIN
         PRINT,'WARNING...information from zoomed images of reprojections'+$
            ' is not photometrically accurate.'
         re_proj = 1 
         j = EXECUTE('pcs = ' + name_orig + '.coordinate_system')
         CASE pcs OF
           'ECLIPTIC'   : p_coord = 'E'
           'EQUATORIAL' : p_coord = 'Q'
           'GALACTIC'   : p_coord = 'G'
           ELSE :
         ENDCASE
         j = EXECUTE('projection = ' + name_orig + '.projection')
         CASE projection OF
           'AITOFF'            : proj = 'A'
           'GLOBAL SINUSOIDAL' : proj = 'S'
           'MOLLWEIDE'         : proj = 'M'
           ELSE                :
         ENDCASE
         j = EXECUTE('win_orig_orig = ' + name_orig + '.win_orig')
         name_orig_orig = get_name(win_orig_orig)
         IF (bflag EQ 1) THEN BEGIN
           IF (zoomflag EQ 1) THEN data(start_x:stop_x,start_y:stop_y)=zbgr
           IF (zoomflag EQ 2) THEN data(start_x:stop_x,start_y:stop_y)=zdsrcmap
           IF (zoomflag EQ 3) THEN data(start_x:stop_x,start_y:stop_y)=zbsub
         ENDIF
         ic = STRPOS(name_orig_orig, '(')
         dtype = STRUPCASE(STRMID(name_orig_orig, 0, ic))
       ENDIF ELSE re_proj = 0
     END
    'PROJ_MAP' :BEGIN
       re_proj = 1
       j = EXECUTE('win_orig = ' + name + '.win_orig')
       name_orig = get_name(win_orig)
       IF(name_orig EQ 'UNDEFINED') THEN BEGIN
         MESSAGE, 'The user has deleted the data associated with a ' + $
           'reprojected image.', /CONT
         GOTO, show_menu
       ENDIF
       j = EXECUTE('pos_x = ' + name + '.pos_x')
       j = EXECUTE('pos_y = ' + name + '.pos_y')
       j = EXECUTE('data = ' + name + '.data')
       j = EXECUTE('orient = ' + name + '.orient')
       j = EXECUTE('pcs = ' + name + '.coordinate_system')
       CASE pcs OF
         'ECLIPTIC'   : p_coord = 'E'
         'EQUATORIAL' : p_coord = 'Q'
         'GALACTIC'   : p_coord = 'G'
         ELSE :
       ENDCASE
       j = EXECUTE('projection = ' + name + '.projection')
       CASE projection OF
         'AITOFF'            : proj = 'A'
         'GLOBAL SINUSOIDAL' : proj = 'S'
         'MOLLWEIDE'         : proj = 'M'
         ELSE                :
       ENDCASE
       zoom = FIX(2.^zoom_index)
       ic = STRPOS(name_orig, '(')
       dtype = STRUPCASE(STRMID(name_orig, 0, ic))
     END
    'PROJ2_MAP': BEGIN
       re_proj = 1
       j = EXECUTE('win_orig = ' + name + '.win_orig')
       name_orig = get_name(win_orig)
       IF(name_orig EQ 'UNDEFINED') THEN BEGIN
         MESSAGE, 'The user has deleted the data associated with a ' + $
           'reprojected image.', /CONT
         GOTO, show_menu
       ENDIF
       j = EXECUTE('pos_x = ' + name + '.pos_x')
       j = EXECUTE('pos_y = ' + name + '.pos_y')
       j = EXECUTE('data = ' + name + '.data')
       j = EXECUTE('orient = ' + name + '.orient')
       j = EXECUTE('pcs = ' + name + '.coordinate_system')
       CASE pcs OF
         'ECLIPTIC'   : p_coord = 'E'
         'EQUATORIAL' : p_coord = 'Q'
         'GALACTIC'   : p_coord = 'G'
         ELSE :
       ENDCASE
       j = EXECUTE('projection = ' + name + '.projection')
       CASE projection OF
         'AITOFF'            : proj = 'A'
         'GLOBAL SINUSOIDAL' : proj = 'S'
         'MOLLWEIDE'         : proj = 'M'
         ELSE                :
       ENDCASE
       zoom = FIX(2.^zoom_index)
       ic = STRPOS(name_orig, '(')
       dtype = STRUPCASE(STRMID(name_orig, 0, ic))
     END
     ELSE : BEGIN
       re_proj = 0
       j = EXECUTE('pos_x = ' + name + '.pos_x')
       j = EXECUTE('pos_y = ' + name + '.pos_y')
       j = EXECUTE('orient = ' + name + '.orient')
       j = EXECUTE('data = ' + name + '.data')
       zoom = FIX(2.^zoom_index)
       ic = STRPOS(name, '(')
       dtype = STRUPCASE(STRMID(name, 0, ic))
     END
  ENDCASE
  sz = SIZE(data)
  dim1_orig = sz(1)
  dim2_orig = sz(2)
;
;  Store the resolution number in RES, for the call to XY2PIX.
;  -----------------------------------------------------------
  CASE STRUPCASE(dtype) OF
    'MAP6'  : res = 6
    'MAP7'  : res = 7
    'MAP8'  : res = 8
    'MAP9'  : res = 9
    'FACE6' : res = 6
    'FACE7' : res = 7
    'FACE8' : res = 8
    'FACE9' : res = 9
     ELSE   : BEGIN
                MESSAGE, 'Unsupported resolution.', /CONT
                GOTO, show_menu
              END
  ENDCASE
  inco = orient + STRING(res, '(i1)')
  WSET, window
  WSHOW, window
  PRINT, ' '
  PRINT, bold(outcoor) + ' coordinates have been selected.'
  PRINT, 'Select points on the image labeled:  ' + bold(title)
  PRINT, ' '
  PRINT, 'Hit ' + bold('LEFT') + ' mouse button to select points, ' + $
    bold('RIGHT') + ' mouse button to exit.'
  oldx = -1
  oldy = -1
  point_number = 0
read_xy:
;
;  Monitor presses of mouse buttons.
;  ---------------------------------
  TVRDC, screenx, screeny, /DEV
  x = (screenx - pos_x) / zoom
  y = (screeny - pos_y) / zoom
  IF(STRMID(name, 0, 6) EQ 'ZOOMED') THEN BEGIN
    x = (start_x > 0) + x
    y = (start_y > 0) + y
  ENDIF
  CASE !ERR OF
    1: BEGIN
;
;  A point was selected by a press of the leftmost mouse button.  Calculate
;  the pixel number and {lon,lat} for that point, but first reject any
;  points which are the same as the last selected point, or any which are
;  beyond the bounds of the image.
;  ------------------------------------------------------------------------
         IF((screenx EQ oldx) AND (screeny EQ oldy)) THEN GOTO, read_xy
         oldx = screenx
         oldy = screeny
         IF((x LT 0) OR (x GE dim1_orig)) THEN GOTO, read_xy
         IF((y LT 0) OR (y GE dim2_orig)) THEN GOTO, read_xy
         IF(re_proj EQ 0) THEN BEGIN
           IF(faceno NE -1) $
             THEN pix = xy2pix([x], [y], res = res, face = faceno) $
             ELSE pix = xy2pix([x], [y], res = res)
           IF(pix(0) EQ -1) THEN GOTO, read_xy
           lonlat = coorconv(pix, infmt='P', outfmt='L', inco=inco, outco=outco)
         ENDIF ELSE BEGIN
           pnt = proj2uv([[x], [y]], proj, SIZE(data))
           pnt_tran = TRANSPOSE(pnt)
           pix = coorconv(pnt_tran, infmt = 'U', outfmt = 'P', inco = p_coord, $
             outco = inco)
           lonlat = coorconv(pix, infmt = 'P', outfmt = 'L', inco = inco, $
             outco = outco)
         ENDELSE
;
;  Put a dot on the image to identify the selected point.
;  ------------------------------------------------------
         TV, [1], screenx, screeny
;
;  If this if the first point to be selected, then print out the
;  top part of the point-information table.
;  -------------------------------------------------------------
         IF(point_number EQ 0) THEN BEGIN
           PRINT, ' '
           str = 'Point #   Pixel #   Lat.(deg)   Lon.(deg)    Data value '
           PRINT, underline(str)
           IF(journal_on) THEN BEGIN
             PRINTF, luj, menu_title
             PRINTF, luj, '  operand:  ' + title
             PRINTF, luj, '  ' + outcoor + ' coordinates have been selected.'
             PRINTF, luj, '  ' + str
           ENDIF
         ENDIF
         point_number = point_number + 1
;
;  Print out the information for this point.
;  -----------------------------------------
         str = STRING(point_number,'(i4)')+STRING(pix,'(i11)')+'    '+$
           STRING(lonlat(1),'(f8.2)')+'    '+STRING(lonlat(0),'(f8.2)')+$
           '    '+STRING(data(x,y),'(g12.4)')
         PRINT, str
         IF(journal_on) THEN PRINTF, luj, '  ' + str
       END
    4: BEGIN
;
;  The user indicated that he is finished selection points by hitting
;  the rightmost mouse button.  Give him the opportunity to clean up
;  the image if desired (i.e., remove the marked points).
;  ------------------------------------------------------------------
         menu = ['Remove marked points from displayed image?', 'Yes', 'No']
         sel = one_line_menu(menu)
         IF(sel EQ 1) THEN IF(STRMID(name, 0, 6) NE 'ZOOMED') $
           THEN xdisplay, name, /rescale ELSE show_zoomed, name
         IF(journal_on) THEN $
           PRINTF, luj, '----------------------------------------' + $
                        '--------------------------------------'
         IF(uimage_version EQ 2) THEN GOTO, show_menu ELSE RETURN
       END
    ELSE:
  ENDCASE
  GOTO, read_xy
keyboard:
 datpix
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


