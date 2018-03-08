PRO extract_spec, area = area
;+
;  EXTRACT_SPEC - a UIMAGE-specific routine.  This routine allows a user
;  to make a spectral plot from a 3-D data-object.  If the AREA keyword
;  is not set, then this routine assumes that the spectra will be generated
;  for a single pixel (i.e., for a single spatial coordinate), else if it
;  is set, then this routine will allow the user to mark off a polygon, and
;  the spectral plot will be created using data-values which are averaged
;  within the defined spatial area.  In either case, the point or polygon
;  is defined on a 2-D object which was generated from (linked to) the
;  selected 3-D object.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 08 93  Change journaling.  J Ewing
;  SPR 10829  Apr  1993  Add code to use any associated weights when
;                        performing an area average. J. Newmark
;  SPR 10906  May 07 93  Add coed to extract a weight spectrum. J Newmark
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11174  Jun 14 94  Upgrade for use of PXINPOLY, J. Newmark
;--------------------------------------------------------------------------
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  first = 1
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    IF(KEYWORD_SET(area)) THEN BEGIN
      MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
      RETURN
    ENDIF
  ENDIF
  IF(NOT KEYWORD_SET(area)) THEN menu_title = 'Extract spectrum for a pixel' $
                            ELSE menu_title = 'Extract spectrum for an area'
;
;  Have the user identify the 3-D object he wants to work with.
;  ------------------------------------------------------------
  select_object3d, index3d
  IF(index3d EQ -1) THEN RETURN
  keyflag=1
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
   IF(NOT KEYWORD_SET(area)) THEN BEGIN
    menu = ['Choose a pixel-selection mode', 'Mouse', 'Keyboard']
    sel = one_line_menu(menu)
    IF(journal_on) THEN PRINTF, luj, '  Pixel selection mode:  ' + menu(sel)
    IF(sel EQ 2) THEN GOTO, input_data ELSE keyflag=0
   ENDIF
;
;  Check if any 2-D objects were derived from the 3-D object.
;  ----------------------------------------------------------
    trace_3d_links, index3d, nlink3d, names
;
;  If there are no linked 2-D objects, then print out a message and exit.
;  ----------------------------------------------------------------------
    IF(nlink3d EQ 0) THEN BEGIN
      PRINT, 'You need to make a 2-D representation of the desired 3-D ' + $
        'data prior to'
      PRINT, 'selecting this operation.  Please make a 2-D representation ' + $
        'by selecting'
      PRINT, 'either '+bold('Extract a frequency/wavelength slice')+ ' or ' + $
        bold('Integrate over frequency/wavelength') + ' from'
      PRINT, 'the ' + bold('Spectrum operations') + ' menu.'
      PRINT,'This will be used in selecting the appropriate pixel or area' + $
         ' with a mouse.'
      RETURN
    ENDIF
;
;  See if any zoomed images exist which are derived from the linked 2-D objects.
;  -----------------------------------------------------------------------------
    nzlink = 0
    sz = SIZE(zoomed)
    IF(sz(2) EQ 8) THEN BEGIN
      dim1_zoomed = sz(1)
      linklist = INTARR(dim1_zoomed)
      FOR i = 0, dim1_zoomed - 1 DO BEGIN
        win_orig = ABS(zoomed(i).win_orig)
        name_orig = get_name(win_orig)
        include_zoom = 0
        FOR k = 0, nlink3d - 1 DO $
          IF(name_orig EQ names(k)) THEN include_zoom = 1
        IF(include_zoom EQ 1) THEN BEGIN
          linklist(nzlink) = i
          nzlink = nzlink + 1
        ENDIF
      ENDFOR
    ENDIF
;
;  Show the menu of all (unzoomed & zoomed) linked 2-D objects.
;  ------------------------------------------------------------
show_menu:
    menu = STRARR(nlink3d + nzlink + 4)
    exit_option = 'Return to previous menu'
    objects = STRARR(nlink3d + nzlink + 4)
    menu(0) = menu_title
    FOR i = 0, nlink3d - 1 DO BEGIN
      objects(i + 1) = names(i)
      j = EXECUTE('title = ' + names(i) + '.title')
      menu(i + 1) = title
    ENDFOR
    FOR i = 0, nzlink - 1 DO BEGIN
      name = 'ZOOMED(' + STRTRIM(STRING(linklist(i)), 2) + ')'
      objects(i + 1 + nlink3d) = name
      j = EXECUTE('title = ' + name + '.title')
      menu(i + 1 + nlink3d) = title
    ENDFOR
    menu(nlink3d + nzlink + 2) = 'HELP'
    menu(nlink3d + nzlink + 3) = exit_option
    sel = umenu(menu, title = 0)
    IF(menu(sel) EQ exit_option) THEN BEGIN
      IF(journal_on AND (first NE 1)) THEN $
        PRINTF, luj, '----------------------------------------' + $
                     '--------------------------------------'
      RETURN
    ENDIF
    IF(menu(sel) EQ 'HELP') THEN BEGIN
      uimage_help, menu_title
      GOTO, show_menu
    ENDIF
    name = objects(sel)
    title2d = menu(sel)
    IF(STRMID(name, 0, 1) EQ 'Z') THEN BEGIN
      using_zoomed = 1
      j = EXECUTE('win_orig = ' + name + '.win_orig')
      name_orig = get_name(win_orig)
    ENDIF ELSE BEGIN
      using_zoomed = 0
      name_orig = name
    ENDELSE
  ENDIF
;
;  Determine various attributes of the selected 3-D object.
;  --------------------------------------------------------
input_data:
  dim1 = object3d(index3d).dim1
  dim2 = object3d(index3d).dim2
  dim3 = object3d(index3d).dim3
  frequnits = object3d(index3d).frequnits
  str_index3d = STRING(index3d, '(i1)')
  j = EXECUTE('frequency = freq3D_' + str_index3d)
  j = EXECUTE('havewt = wt3D_' + str_index3d )
  wsz=SIZE(havewt)
  linkweight=-1
  IF((dim2 / 3) * 4 EQ dim1) THEN BEGIN
;
;  MAP6, MAP7, MAP8, or MAP9
;  -------------------------
    type = 'MAP'
    CASE dim1 OF
       128 : res = 6
       256 : res = 7
       512 : res = 8
      1024 : res = 9
      ELSE : 
    ENDCASE
  ENDIF
  IF(dim1 EQ dim2) THEN BEGIN
;
;  FACE6, FACE7, FACE8, or FACE9
;  -----------------------------
    type = 'FACE'
    faceno = object3d(index3d).faceno
    CASE dim1 OF
        32 : res = 6
        64 : res = 7
       128 : res = 8
       256 : res = 9
      ELSE : 
    ENDCASE
  ENDIF
  badpixval = object3d(index3d).badpixval
  y_title = object3d(index3d).units
  ocoord = 'R' + STRING(res, '(i1)')
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN') AND (keyflag NE 1)) THEN BEGIN
    j = EXECUTE('title = ' + name + '.title')
    j = EXECUTE('data_window = ' + name + '.window')
    j = EXECUTE('pos_x = ' + name_orig + '.pos_x')
    j = EXECUTE('pos_y = ' + name_orig + '.pos_y')
    IF(using_zoomed EQ 1) THEN BEGIN
      j = EXECUTE('start_x = ' + name + '.start_x')
      j = EXECUTE('start_y = ' + name + '.start_y')
      j = EXECUTE('stop_x = ' + name + '.stop_x')
      j = EXECUTE('stop_y = ' + name + '.stop_y')
      j = EXECUTE('specific_zoom = ' + name + '.specific_zoom')
      zoom = (2 ^ zoom_index) * specific_zoom
      pos_x = 0.
      pos_y = 0.
    ENDIF ELSE zoom = 2. ^ zoom_index
    nblocks_x = 5
    nblocks_y = 5
  ENDIF
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  3-D operand:  ' + object3d(index3d).title
    IF (keyflag EQ 0) THEN PRINTF, luj, '  2-D operand:  ' + title2d
  ENDIF
;
;  Handle "singe pixel" case.
;  --------------------------
  IF(NOT KEYWORD_SET(area)) THEN BEGIN
    IF (keyflag NE 1) THEN BEGIN
;
;  Monitor mouse activity.
;  -----------------------
      j=EXECUTE('hidden =' + name + '.hidden')
      IF (hidden EQ 1) THEN BEGIN
         PRINT,' The chosen image is hidden, please redisplay it.'
         RETURN
      ENDIF
      PRINT, ' '
      PRINT, 'Please select points on the image with the following label:'
      PRINT, bold(title)
      PRINT, ' '
      PRINT, 'Use ' + bold('LEFT') + ' button for point selection, ' + $
        bold('RIGHT') + ' button to quit.', FORMAT = '(a)'
monitor_mouse:
      WSET, data_window
      WSHOW, data_window
      TVRDC, cur_x, cur_y, /DEV
      CASE !err OF
        4 : BEGIN
              IF(journal_on AND (first NE 1)) THEN $
                PRINTF, luj, '----------------------------------------' + $
                             '--------------------------------------'
              RETURN
            END
        1 : BEGIN
;
;  A point was selected - see if it's a valid point.
;  -------------------------------------------------
              ix = FIX((cur_x - pos_x) / zoom)
              iy = FIX((cur_y - pos_y) / zoom)
              IF(using_zoomed) THEN BEGIN
                ix = start_x + ix
                iy = start_y + iy
              ENDIF
              IF((ix LT 0) OR (ix GT (dim1 - 1)) OR (iy LT 0) OR $
                (iy GT (dim2 - 1))) THEN BEGIN
                PRINT, 'An inappropriate point was selected.  ', $
                  'Please pick another point.'
                GOTO, monitor_mouse
              ENDIF
              j = EXECUTE('data_y = data3D_' + str_index3d + '(ix, iy, *)')
              allbad = 1
              FOR i = 0, dim3 - 1 DO IF(data_y(i) NE badpixval) THEN allbad = 0
              IF(allbad EQ 1) THEN BEGIN
                PRINT, 'A bad pixel was selected.  Please pick another pixel.'
                GOTO, monitor_mouse
              ENDIF
;
;  It is a valid point - get the pixel number and make the graph.
;  --------------------------------------------------------------
              TV, [[1]], cur_x, cur_y
              IF(type EQ 'FACE') $
                THEN pix = xy2pix([ix], [iy], res = res, face = faceno) $
                ELSE pix = xy2pix([ix], [iy], res = res)
              title = 'Pixel spectrum [' + STRTRIM(STRING(pix(0)), 2) + ']'
;   First set up weight object if weights exist.
              IF (wsz(0) EQ 3) THEN BEGIN
               j = EXECUTE('wt_y = wt3D_' + str_index3d + '(ix, iy, *)')
               wtitle='Weights for ' +title
               wname = setup_graph(frequency, wt_y, dim3,$
                 title = wtitle, x_title = frequnits, y_title = y_title, $
                 nblocks_x = nblocks_x, nblocks_y = nblocks_y,hidden=1)
               j=EXECUTE('linkweight=' + wname +'.window')
              ENDIF
              gname = setup_graph(frequency, data_y, dim3,$
                title = title, x_title = frequnits, y_title = y_title, $
                nblocks_x = nblocks_x, nblocks_y = nblocks_y,$
                linkweight=linkweight)
              j=EXECUTE('winnum =' + gname + '.window')
              IF ((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN WSET,winnum
              plot_graph, gname
              IF(journal_on) THEN PRINTF, luj, '  outputted graph:  ' + title
              IF(first EQ 1) THEN BEGIN
                sel = umenu(['Resize graph?', 'Yes / Resize the graph', $
                  'No / Exit this menu'], title = 0)
                IF(sel EQ 1) THEN BEGIN
                  resize_graph, gname
                  j = EXECUTE('nblocks_x = ' + gname + '.nblocks_x')
                  j = EXECUTE('nblocks_y = ' + gname + '.nblocks_y')
                ENDIF
                first = 0
              ENDIF
;
; Print out spectrum to a file if requested.
        menu=['Do you wish to the spectrum to be sent to a file?','Yes','No']
        ans=one_line_menu(menu,init=2)
        IF (ans EQ 1) THEN BEGIN
         err = 1
         filename = ''
         PRINT, 'Enter the complete name of the file to which the spectrum'
         PRINT, '  should be sent.  Example:  SYS$LOGIN:spectrum.firas'
         WHILE(err NE 0) DO BEGIN
          READ, underline('Enter file name:') + '  ', filename
          OPENW, lus, filename, ERROR = err, /GET_LUN
          IF(err NE 0) THEN BEGIN
            PRINT, bold('An inappropriate file name was entered.')
            PRINT, ' '
          ENDIF
         ENDWHILE
         PRINTF,lus,title    
         PRINTF,lus,' '
         IF (wsz(0) EQ 3) THEN BEGIN
          PRINTF,lus,'Frequency/wavelength    Value   Weight'
          FOR i=0,dim3-1 DO PRINTF,lus,frequency(i),data_y(i),wt_y(i)
         ENDIF ELSE BEGIN
          PRINTF,lus,'Frequency/wavelength     Value'
          FOR i=0,dim3-1 DO PRINTF,lus,frequency(i),data_y(i)
         ENDELSE
      FREE_LUN,lus
        ENDIF
       END
      ENDCASE
      GOTO, monitor_mouse
    ENDIF ELSE BEGIN
keyboard_entry:
;
;  Have the user type in the {lon,lat} coordinates of the desired point.
;  ---------------------------------------------------------------------
      PRINT, 'You will be asked to enter a longitude/latitude coordinate pair'
      PRINT, 'which identifies a location for which a spectrum will be plotted.'
      PRINT, 'Please select the coordinate system you will use.'
      menu = ['Coordinate system:', 'Ecliptic', 'Equatorial', 'Galactic']
      sel = one_line_menu(menu)
      IF(journal_on) THEN PRINTF, luj, '  coordinate system:  ' + menu(sel)
      icoord = STRMID('EQG', sel - 1, 1)
      PRINT, 'Next enter the lon/lat pair in degrees (example: ' + $
        bold('20. -30.') + ').  [<CR> to exit]'
      str = ''
read_lonlat:
      READ, 'lon lat:  ', str
      str = STRTRIM(str, 2)
      IF(str EQ '') THEN RETURN
      ic = STRPOS(str, ' ')
      IF(ic EQ -1) THEN ic = STRPOS(str, ',')
      IF(ic EQ -1) THEN BEGIN
        PRINT, 'Invalid delimiter (should be " " or ",").  Please enter again.'
        GOTO, read_lonlat
      ENDIF
      str1 = STRMID(str, 0, ic)
      str2 = STRMID(str, ic + 1, STRLEN(str))
      IF(validnum(str1) + validnum(str2) NE 2) THEN BEGIN
        PRINT, 'Invalid numerical value(s).  Please enter again.'
        GOTO, read_lonlat
      ENDIF
      lon = FLOAT(str1)
      lat = FLOAT(str2)
      IF((lat LT -90.) OR (lat GT 90.)) THEN BEGIN
        PRINT, 'Latitude must be between -90 and 90.  Please enter again.'
        GOTO, read_lonlat
      ENDIF
      IF(journal_on) THEN PRINTF, luj, '  lon, lat = ' + str
;
;  Convert {lon,lat} to array indices.
;  -----------------------------------
      pix = coorconv([[lon], [lat]], infmt = 'L', outfmt = 'P', $
        inco = icoord, outco = ocoord)
      IF(type EQ 'FACE') THEN BEGIN
        onface = FIX(pix(0) / (4. ^ (res - 1)))
        IF(onface NE faceno) THEN BEGIN
          PRINT, 'Selected point does not lie on the appropriate face.' + $
            '  Please enter again.'
          GOTO, read_lonlat
        ENDIF
        pix2xy, pix, ix, iy, res = res, /face
      ENDIF ELSE BEGIN
        pix2xy, pix, ix, iy, res = res
      ENDELSE
;
;  Make the graph of the spectrum.
;  -------------------------------
      j = EXECUTE('data_y = data3D_' + str_index3d + '(ix, iy, *)')
      allbad = 1
      FOR i = 0, dim3 - 1 DO IF(data_y(i) NE badpixval) THEN allbad = 0
      IF(allbad EQ 1) THEN BEGIN
        PRINT, 'A bad pixel was selected.  Please pick another pixel.'
        GOTO, read_lonlat
      ENDIF
      title = 'Pixel spectrum [' + STRTRIM(STRING(pix(0)), 2) + ']'
; First set up weight object if weights exist.
      IF (wsz(0) EQ 3) THEN BEGIN
         j = EXECUTE('wt_y = wt3D_' + str_index3d + '(ix, iy, *)')
         wtitle='Weights for ' +title
         wname = setup_graph(frequency, wt_y, dim3,$
           title = wtitle, x_title = frequnits, y_title = y_title, $
           nblocks_x = nblocks_x, nblocks_y = nblocks_y,hidden=1)
         j=EXECUTE('linkweight=' + wname +'.window')
      ENDIF
      gname = setup_graph(frequency,data_y,dim3,linkweight=linkweight,$
        title = title, x_title = frequnits, y_title = y_title)
      j=EXECUTE('winnum =' + gname + '.window')
      IF ((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN WSET,winnum
      plot_graph, gname
      IF(journal_on) THEN PRINTF, luj, '  outputted graph:  ' + title
      IF(first AND ((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN'))) THEN BEGIN
        sel = umenu(['Resize graph?', $
          'Yes / Resize the graph', 'No / Exit this menu'], title = 0)
        IF(sel EQ 1) THEN BEGIN
          resize_graph, gname
          j = EXECUTE('nblocks_x = ' + gname + '.nblocks_x')
          j = EXECUTE('nblocks_y = ' + gname + '.nblocks_y')
        ENDIF
        first = 0
      ENDIF
;
; Print out spectrum to a file if requested.
      menu=['Do you wish to the spectrum to be sent to a file?','Yes','No']
      ans=one_line_menu(menu,init=2)
      IF (ans EQ 1) THEN BEGIN
         err = 1
         filename = ''
         PRINT, 'Enter the complete name of the file to which the spectrum'
         PRINT, '  should be sent.  Example:  SYS$LOGIN:spectrum.firas'
         WHILE(err NE 0) DO BEGIN
          READ, underline('Enter file name:') + '  ', filename
          OPENW, lus, filename, ERROR = err, /GET_LUN
          IF(err NE 0) THEN BEGIN
            PRINT, bold('An inappropriate file name was entered.')
            PRINT, ' '
          ENDIF
         ENDWHILE
         PRINTF,lus,title    
         PRINTF,lus,' '
         IF (wsz(0) EQ 3) THEN BEGIN
          PRINTF,lus,'Frequency/wavelength    Value   Weight'
          FOR i=0,dim3-1 DO PRINTF,lus,frequency(i),data_y(i),wt_y(i)
         ENDIF ELSE BEGIN
          PRINTF,lus,'Frequency/wavelength     Value'
          FOR i=0,dim3-1 DO PRINTF,lus,frequency(i),data_y(i)
         ENDELSE
         FREE_LUN,lus
       ENDIF
    ENDELSE
  ENDIF ELSE BEGIN
;
;  Handle "area" case.
;  -------------------
mark_an_area:
;
;  The DEFINE_AREA routine will let them identify the pixels in the area.
;  ------------------------------------------------------------------------
    j=EXECUTE('hidden =' + name + '.hidden')
    IF (hidden EQ 1) THEN BEGIN
       PRINT,' The chosen image is hidden, please redisplay it.'
       RETURN
    ENDIF
    dum = define_area(name, poly_pix)
    yval = FLTARR(dim3)
    wtval= FLTARR(dim3)
    numpix=INTARR(dim3)
    allbad = 1
    links=object3d(index3d).linkweight
    IF (links NE -1) THEN BEGIN
      wtname=get_name(links)
      j=EXECUTE('weight= ' + wtname + '.data')
    ENDIF
;
;  For each frequency, find the mean value of good pixels in the area.
;  -------------------------------------------------------------------
;
;  If there are any weights associated with the data, use them in
;   performing the area average.
    IF ((wsz(0) EQ 3) OR (links NE -1)) THEN BEGIN
     FOR i = 0, dim3 - 1 DO BEGIN
      j = EXECUTE('data = data3D_' + str_index3d + '(*, *, i)')
      IF (wsz(0) EQ 3) THEN j = EXECUTE('weight = wt3D_' + str_index3d +$
       '(*, *, i)')
      area = pix2dat(pixel=poly_pix,raster=data)
      good = WHERE(area NE badpixval)
      wt = pix2dat(pixel=poly_pix,raster=weight)
      scale=TOTAL(wt(good))
      IF((good(0) EQ -1) OR (scale EQ 0)) THEN mean = badpixval ELSE BEGIN
         yval(i)=TOTAL(area(good)*wt(good))/scale
         allbad = 0
         wsd=STDEV(wt(good), mean)
         wtval(i)=mean
         sz=SIZE(good)
         numpix(i)=sz(1)
      ENDELSE
     ENDFOR
    ENDIF ELSE BEGIN
     FOR i = 0, dim3 - 1 DO BEGIN
      j = EXECUTE('data = data3D_' + str_index3d + '(*, *, i)')
      area = pix2dat(pixel=poly_pix,raster=data)
      good = WHERE(area NE badpixval)
      IF(good(0) EQ -1) THEN mean = badpixval ELSE BEGIN
        sd = STDEV(area(good), mean)
        allbad = 0
        sz=SIZE(good)
        numpix(i)=sz(1)
     ENDELSE
      yval(i) = mean
     ENDFOR
    ENDELSE
    IF(allbad EQ 1) THEN BEGIN
      PRINT, ' '
      PRINT, bold('An area full of bad pixels was selected.  ' + $
        'Please select another area.')
      GOTO, mark_an_area
    ENDIF
;
;  Plot the averaged spectrum for the area.
;  ----------------------------------------
    title = append_number('Area spectrum')
    IF (wsz(0) EQ 3) THEN BEGIN
      wtitle='Weights for ' +title
      wname = setup_graph(frequency, wtval, dim3,$
        title = wtitle, x_title = frequnits, y_title = y_title, $
        nblocks_x = nblocks_x, nblocks_y = nblocks_y,hidden=1)
      j=EXECUTE('linkweight=' + wname +'.window')
    ENDIF
    gname = setup_graph(frequency, yval, dim3, $
      title = title, x_title = frequnits, linkweight=linkweight,$
      y_title = y_title, nblocks_x = nblocks_x, nblocks_y = nblocks_y)
    j=EXECUTE('winnum =' + gname + '.window')
    IF ((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN WSET,winnum
    plot_graph, gname
    sdnum=STDEV(numpix,mean)
    pixels=FIX(mean)
    PRINT,'# of good pixels in area = ',STRTRIM(STRING(pixels),2)
    IF(journal_on) THEN PRINTF, luj, '  outputted graph:  ' + title
    IF(first EQ 1) THEN BEGIN
      sel = umenu(['Resize graph?', $
        'Yes / Resize the graph', 'No / Exit this menu'], title = 0)
      IF(sel EQ 1) THEN BEGIN
        resize_graph, gname
        j = EXECUTE('nblocks_x = ' + gname + '.nblocks_x')
        j = EXECUTE('nblocks_y = ' + gname + '.nblocks_y')
      ENDIF
      first = 0
    ENDIF
;
;  If desired, remove the polygon border from the image.
;  -----------------------------------------------------
    menu = ['Remove area boundary from image?', 'Yes', 'No']
    sel = one_line_menu(menu, init = 1)
    IF(sel EQ 1) THEN IF(using_zoomed NE 1) THEN xdisplay, name, /rescale $
                                            ELSE show_zoomed, name
;
; Print out spectrum to a file if requested.
      menu=['Do you wish to the spectrum to be sent to a file?','Yes','No']
      ans=one_line_menu(menu,init=2)
      IF (ans EQ 1) THEN BEGIN
         err = 1
         filename = ''
         PRINT, 'Enter the complete name of the file to which the spectrum'
         PRINT, '  should be sent.  Example:  SYS$LOGIN:spectrum.firas'
         WHILE(err NE 0) DO BEGIN
          READ, underline('Enter file name:') + '  ', filename
          OPENW, lus, filename, ERROR = err, /GET_LUN
          IF(err NE 0) THEN BEGIN
            PRINT, bold('An inappropriate file name was entered.')
            PRINT, ' '
          ENDIF
         ENDWHILE
         PRINTF,lus,title    
         PRINTF,lus,' '
         PRINTF,lus,'# of good pixels in area = ',STRTRIM(STRING(pixels),2)
         IF ((wsz(0) EQ 3) OR (links NE -1)) THEN BEGIN
          PRINTF,lus,'Frequency/wavelength    Value   Weight'
          FOR i=0,dim3-1 DO PRINTF,lus,frequency(i),yval(i),wtval(i)
         ENDIF ELSE BEGIN
          PRINTF,lus,'Frequency/wavelength     Value'
          FOR i=0,dim3-1 DO PRINTF,lus,frequency(i),yval(i)
         ENDELSE
         FREE_LUN,lus
       ENDIF
;
;  If desired, allow another area to be marked off.
;  ------------------------------------------------
    menu = ['Make a plot for another area?', 'Yes', 'No']
    sel = one_line_menu(menu, init = 2)
    IF(sel EQ 1) THEN GOTO, mark_an_area
  ENDELSE
  IF(journal_on EQ 1) THEN $
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
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


