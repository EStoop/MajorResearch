PRO ualgebra, operation, name1
;+
;  UALGEBRA - a UIMAGE-specific routine.  This routine allows a user to
;  perform algebraic operations.
;#
;  Written by Vidya Sagar & John Ewing.
;  SPR 10332  Jan 05 93  Put "init=sel" into call to UMENU.  J Ewing
;  SPR 10383  Jan 08 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 10829  Apr  1993  If data has weights associated with it, properly
;                        propagate the weights when combining objects.
;                        Modified to have weighted addition, subtraction and 
;                        averages. J Newmark
;  SPR 10923  May 10 93  Add reciprocal function. J. Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  nblocks_x = 5
  nblocks_y = 5
  resized_graph = 0
  first = 1
;
;  Have the user select which algebraic operation he'd like to employ.
;  -------------------------------------------------------------------
get_argument:
  main_menu_title = 'ALGEBRAIC OPERATIONS'
  IF(N_PARAMS(0) GE 1) THEN GOTO, have_op
  main_menu = [main_menu_title, $
              'X1 + X2', $
              'X1 - X2', $
              'X1 * X2', $
              'X1 / X2', $
              '', $
              '(X1*wt1+X2*wt2)/(wt1+wt2)', $
              '(X1*wt1-X2*wt2)/(wt1+wt2)', $
              '', $
              'SQRT(X)', $
              'LOG10(X)', $
              'ABS(X)', $
              '1./X', $
              '', $
              'C0 + C1*X1 + C2*X2 +...+ Cn*Xn', $
              'Average(X1, X2, X3, ..., Xn)', $
              'Weighted Avg((X1*wt1,..., Xn*Wtn)/Sum(wt1-n))', $
              '', $
              'HELP', $
              'Return to MAIN MENU']
  selop = UMENU(main_menu, title = 0, init = selop)
;
;  Set the variable OPERATION so it indicates which operation is to
;  be performed.
;  ----------------------------------------------------------------
  CASE selop OF
     1: operation = 'ADD'
     2: operation = 'SUB'
     3: operation = 'MUL'
     4: operation = 'DIV'
     6: operation = 'WADD'
     7: operation = 'WSUB'
     9: operation = 'SQR'
    10: operation = 'LOG'
    11: operation = 'ABS'
    12: operation = 'REC'
    14: operation = 'WSM'
    15: operation = 'AVG'
    16: operation = 'WAVG'
    18: BEGIN
          uimage_help, main_menu_title
          GOTO, get_argument
        END
    19: RETURN
    ELSE : GOTO, get_argument
  ENDCASE
  PRINT, ' '
  PRINT, 'Selected operation:  ' + bold(main_menu(selop))
have_op:
;
;  If the weighted-sum or average operation was selected, then branch
;  elsewhere.
;  ------------------------------------------------------------------
  IF((operation EQ 'WSM') OR (operation EQ 'AVG')) THEN GOTO, wsum_or_avg
;
; If weighted addition, weighted subtraction or weighted average was
;   selected, then branch elsewhere.
;
  IF((operation EQ 'WADD') OR (operation EQ 'WSUB') OR (operation EQ 'WAVG')) $
   THEN GOTO, w_avg
;
;  Have the user select an operand.
;  --------------------------------
identify_X1:
  IF((operation EQ 'SQR') OR (operation EQ 'LOG') OR (operation EQ 'ABS') $
    OR (operation EQ 'REC')) THEN id = 'X' ELSE id = 'X1'
  menu_title = 'Select ' + id
  IF(N_PARAMS(0) EQ 2) THEN GOTO, have_x1
  name1 = select_object(/map6,/map7,/map8,/map9,/face6,/face7,/face8,$
          /face9,/proj_map,/proj2_map,/simple_graph,/help,$
          /EXIT,title=menu_title)
  IF(name1 EQ 'HELP') THEN BEGIN
    UIMAGE_HELP, menu_title
    GOTO, identify_X1
  ENDIF
  IF(name1 EQ 'EXIT') THEN BEGIN
    IF(uimage_version EQ 2) THEN GOTO, get_argument ELSE RETURN
  ENDIF
  IF(name1 EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'No objects are currently available.', /CONT
    RETURN
  ENDIF
  j = EXECUTE('title1 = ' + name1 + '.title')
  PRINT, bold(id) + ':  ' + title1
have_x1:
;
;  Extract some fields out of the selected operand's structure.
;  ------------------------------------------------------------
  parent = name1
  IF((STRMID(name1, 0, 8) EQ 'PROJ_MAP') OR $
     (STRMID(name1, 0, 9) EQ 'PROJ2_MAP'))THEN BEGIN
    name_proj = name1
    j = EXECUTE('win_orig = ' + name1 + '.win_orig')
    name1 = get_name(win_orig)
    j = EXECUTE('projection = ' + name_proj + '.projection')
    j = EXECUTE('coordinate_system = ' + name_proj + '.coordinate_system')
    j = EXECUTE('sz = SIZE(' + name_proj + '.data)')
    IF (sz(1) EQ 1024) THEN psize='L' ELSE psize='S'
    proj = STRUPCASE(STRMID(projection, 0, 1))
    IF(projection EQ 'GLOBAL SINUSOIDAL') THEN proj = 'S'
    coor = STRUPCASE(STRMID(coordinate_system, 0, 1))
    IF(STRUPCASE(coordinate_system) EQ 'EQUATORIAL') THEN coor = 'Q'
    IF(STRMID(name1, 0, 4) EQ 'FACE') THEN $
      j = EXECUTE('faceno = ' + name1 + '.faceno')
  ENDIF ELSE name_proj = ''
  IF(STRMID(name1, 0, 5) EQ 'GRAPH') THEN graphs = 1 ELSE graphs = 0
;
;  Handle the situation for the square root, log, and absolute value
;  operations.
;  -----------------------------------------------------------------
  IF((operation EQ 'SQR') OR (operation EQ 'LOG') OR (operation EQ 'ABS') $
    OR (operation EQ 'REC')) THEN BEGIN
    j = EXECUTE('badpixval = ' + name1 + '.badpixval')
    j = EXECUTE('title1 = ' + parent + '.title')
    IF(NOT graphs) THEN BEGIN
      j = EXECUTE('scale_min = ' + name1 + '.scale_min')
      j = EXECUTE('scale_max = ' + name1 + '.scale_max')
    ENDIF ELSE j = EXECUTE('n1 = FIX(' + name1 + '.num) - 1')
    CASE operation OF
      'SQR': BEGIN
           title = 'SQRT("' + title1 + '")'
           IF(STRLEN(title) GT 40) THEN $
             title = append_number('Square Root Result')
            IF(NOT graphs) THEN BEGIN
             exstr = 'result = SQUAREROOT(' + name1 + $
               '.data, badval_in = badpixval, outbadval = outbadval)'
            ENDIF ELSE BEGIN
             exstr = 'result = SQUAREROOT(' + name1 + $
               '.data(0:n1, *), badval_in = badpixval, outbadval = outbadval)'
            ENDELSE
           j = EXECUTE(exstr)
         END
      'LOG': BEGIN
           title = 'LOG("' + title1 + '")'
           IF(STRLEN(title) GT 40) THEN title = append_number('Log Result')
            IF(NOT graphs) THEN BEGIN
             exstr = 'result = LOGBASE10(' + name1 + $
               '.data, badval_in = badpixval, outbadval = outbadval)'
            ENDIF ELSE BEGIN
             exstr = 'result = LOGBASE10(' + name1 + $
               '.data(0:n1, *), badval_in = badpixval, outbadval = outbadval)'
            ENDELSE
           j = EXECUTE(exstr)
         END
      'ABS': BEGIN
           title = 'ABS("' + title1 + '")'
           IF(STRLEN(title) GT 40) THEN title = append_number('Abs Result')
           j = EXECUTE('neg = WHERE(' + name1 + '.data LT 0)')
            IF(NOT graphs) THEN BEGIN
             exstr = 'result = ABSVALUE(' + name1 + $
               '.data, badval_in = badpixval, outbadval = outbadval)'
            ENDIF ELSE BEGIN
             exstr = 'result = ABSVALUE(' + name1 + $
               '.data(0:n1, *), badval_in = badpixval, outbadval = outbadval)'
            ENDELSE
           j = EXECUTE(exstr)
         END
      'REC': BEGIN
           title = 'RECIPROCAL("' + title1 + '")'
           IF(STRLEN(title) GT 40) THEN title = append_number('Rec Result')
           j=EXECUTE('result=' +name1 + '.data')
           IF(NOT graphs) THEN BEGIN
             good=WHERE((result NE badpixval) AND (result NE 0))
             exstr='result(good) = 1./result(good)'
           ENDIF ELSE BEGIN
             good=WHERE((result(0:n1,1) NE badpixval) AND (result(0:n1,1) NE 0))
             exstr='result(good,1) = 1./result(good,1)'
           ENDELSE
           j = EXECUTE(exstr)
           outbadval=badpixval
         END
    ENDCASE
;
;  If the result is an image, then store that image in the UIMAGE data
;  environment, and then display it (if on an X-window terminal).
;  -------------------------------------------------------------------
    IF(NOT graphs) THEN BEGIN
      good = WHERE(result NE badpixval)
      IF(good(0) NE -1) THEN BEGIN
        scale_min = MIN(result(good), MAX = scale_max)
      ENDIF ELSE BEGIN
        scale_min = badpixval
        scale_max = badpixval
      ENDELSE
      newname=setup_image(data=result, title=title, badpixval=outbadval, $
        scale_min=scale_min, scale_max=scale_max, parent=parent)
      IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, newname
      IF(name_proj NE '') THEN BEGIN
        reproj,result,result,proj=proj,coor=coor,psize=psize,face=faceno,/noshow
        title = 'Reprojected ' + title
        IF(STRLEN(title) GT 40) THEN title = append_number('Reprojected result')
        newnameproj=setup_image(data=result,title=title,badpixval=outbadval,$
          projection=projection,coordinate_system=coordinate_system, $
          scale_min=scale_min,scale_max=scale_max,parent=newname)
        IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, newnameproj
      ENDIF
    ENDIF ELSE BEGIN
;
;  If the result is an graph, then store that graph in the UIMAGE data
;  environment, and then display it.
;  -------------------------------------------------------------------
      good = WHERE(result(*, 1) NE outbadval)
      IF(good(0) EQ -1) THEN BEGIN
        PRINT, 'Resulting graph consists of all bad values - operation ' + $
          'aborted.'
        GOTO, done
      ENDIF
      j = EXECUTE('x_title = ' + name1 + '.x_title')
      j = EXECUTE('y_title = ' + name1 + '.y_title')
      j = EXECUTE('psymbol = ' + name1 + '.psymbol')
      sz = SIZE(result)
      newname=setup_graph(result(*,0),result(*,1),sz(1),title=title, $
        /special,badpixval=outbadval,x_title=x_title,y_title=y_title,$
        nblocks_x=nblocks_x,nblocks_y=nblocks_y,psymbol=psymbol)
      plot_graph, newname
      IF(((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) $
        AND (uimage_version EQ 2)) THEN BEGIN
        IF(resized_graph EQ 0) THEN BEGIN
          sel = umenu(['Resize graph?', $
            'Yes / Resize the graph', 'No / Exit this menu'], title = 0)
          IF(sel EQ 1) THEN BEGIN
            resize_graph, newname
            j = EXECUTE('nblocks_x = ' + newname + '.nblocks_x')
            j = EXECUTE('nblocks_y = ' + newname + '.nblocks_y')
          ENDIF
          resized_graph = 1
        ENDIF
      ENDIF
    ENDELSE
    IF(N_PARAMS(0) NE 2) THEN PRINT, bold('result') + ':  ' + title
;
;  If journaling is enabled, then write some info to the journal file.
;  -------------------------------------------------------------------
    IF(journal_on) THEN BEGIN
      PRINTF, luj, main_menu(selop)
      j = EXECUTE('oldtitle = ' + name1 + '.title')
      j = EXECUTE('newtitle = ' + newname + '.title')
      PRINTF, luj, '  the following statement was executed:'
      PRINTF, luj, '  ' + STRMID(exstr, 0, 76)
      exstr = STRMID(exstr, 76, STRLEN(exstr))
      WHILE(STRLEN(exstr) GE 1) DO BEGIN
        PRINTF, luj, '    ' + STRMID(exstr, 0, 74)
        exstr = STRMID(exstr, 74, STRLEN(exstr))
      ENDWHILE
      PRINTF, luj, '  where:'
      PRINTF, luj, '    ' + name1 + ' --> ' + oldtitle
      PRINTF, luj, '    badpixval = ', STRTRIM(badpixval, 2)
      PRINTF, luj, '    result --> ' + newtitle
      PRINTF, luj, '----------------------------------------' + $
                   '--------------------------------------'
    ENDIF
    IF(uimage_version EQ 2) THEN GOTO, get_argument ELSE RETURN
  ENDIF ELSE BEGIN
;
;  The selected operation requires the selection of a second operand.
;  Provide a mechanism below for that second operand to be selected.
;  ------------------------------------------------------------------
    types = STRMID(parent, 0, STRPOS(parent, '('))
    badval_in = FLTARR(2)
    titlen = 'Select X2'
identify_X2:
    j = EXECUTE('name2=select_object(/'+types+',/help,/EXIT,title=' + $
      'titlen)')
    IF(name2 EQ 'HELP') THEN BEGIN
      uimage_help, 'Select Xn'
      GOTO, identify_X2
    ENDIF
    IF(name2 EQ 'EXIT') THEN BEGIN
      IF(uimage_version EQ 2) THEN GOTO, get_argument ELSE RETURN
    ENDIF
    parent2 = name2
    IF((STRMID(name2, 0, 8) EQ 'PROJ_MAP') OR $ 
      (STRMID(name2, 0, 9) EQ 'PROJ2_MAP'))THEN BEGIN
      j = EXECUTE('win_orig = ' + name2 + '.win_orig')
      name2 = get_name(win_orig)
    ENDIF
    j = EXECUTE('title2 = ' + parent2 + '.title')
    PRINT, bold('X2') + ':  ' + title2
;
;  Get bad pixel values and scaling parameters for both operands.
;  --------------------------------------------------------------
    j = EXECUTE('badval_in(0) = ' + name1 + '.badpixval')
    j = EXECUTE('badval_in(1) = ' + name2 + '.badpixval')
    IF(NOT graphs) THEN BEGIN
      j = EXECUTE('scale_min1 = ' + name1 + '.scale_min')
      j = EXECUTE('scale_min2 = ' + name2 + '.scale_min')
      j = EXECUTE('scale_max1 = ' + name1 + '.scale_max')
      j = EXECUTE('scale_max2 = ' + name2 + '.scale_max')
    ENDIF
    IF(graphs) THEN BEGIN
      j = EXECUTE('n1 = FIX(' + name1 + '.num) - 1')
      j = EXECUTE('n2 = FIX(' + name2 + '.num) - 1')
    ENDIF
;
; Check if there are any weights associated with the objects.
j=EXECUTE('link1= ' + name1 + '.linkweight')
j=EXECUTE('link2= ' + name2 + '.linkweight')
newlink=-1
;
;  Handle the situation for the binary operations (add, subtract, multiply,
;  and divide). In each case. if both objects have weights associated with
;  them, then create a new weight object for the output.
;  ------------------------------------------------------------------------
    CASE operation OF
      'ADD': BEGIN
            IF(NOT graphs) THEN BEGIN
              exstr = 'result = WEIGHTEDSUM(' + name1 + '.data, ' + name2 + $
               '.data, coef=[0.,1.,1.], badval_in = badval_in, ' + $
               'outbadval = outbadval)'
            ENDIF ELSE BEGIN
              exstr = 'result = WEIGHTEDSUM(' + name1 + '.data(0:n1,*), ' + $
               name2 + '.data(0:n2,*), coef = [0.,1.,1.], badval_in = ' + $
               'badval_in, outbadval = outbadval)'
            ENDELSE
            IF ((link1 NE -1) AND (link2 NE -1)) then begin
              newlink=mkweight(link1,link2)
            ENDIF
            j = EXECUTE(exstr)
            sz = SIZE(result)
            IF((sz(0) NE 2) OR (sz(1) LT 2)) THEN BEGIN
              IF(uimage_version EQ 2) THEN GOTO, get_argument ELSE RETURN
            ENDIF
            title = '"' + title1 + '" + "' + title2 + '"'
            IF(STRLEN(title) GT 40) THEN title=append_number('Addition Result')
            good = WHERE(result NE outbadval)
            IF(good(0) NE -1) THEN BEGIN
              scale_min = MIN(result(good), MAX = scale_max)
            ENDIF ELSE BEGIN
              scale_min = outbadval
              scale_max = outbadval
            ENDELSE
          END
      'SUB' : BEGIN
            IF(NOT graphs) THEN BEGIN
              exstr = 'result = WEIGHTEDSUM(' + name1 + '.data, ' + name2 + $
               '.data, coef = [0.,1.,-1.], badval_in = badval_in, ' + $
               'outbadval = outbadval)'
            ENDIF ELSE BEGIN
              exstr = 'result = WEIGHTEDSUM(' + name1 + '.data(0:n1,*),' + $
               name2 + '.data(0:n2,*), coef=[0.,1.,-1.], badval_in = ' + $
               'badval_in, outbadval = outbadval)'
            ENDELSE
            IF ((link1 NE -1) AND (link2 NE -1)) then begin
                 newlink=mkweight(link1,link2)
            ENDIF
            j = EXECUTE(exstr)
            sz = SIZE(result)
            IF((sz(0) NE 2) OR (sz(1) LT 2)) THEN BEGIN
              IF(uimage_version EQ 2) THEN GOTO, get_argument ELSE RETURN
            ENDIF
            title = '"' + title1 + '" - "' + title2 + '"'
            IF(STRLEN(title) GT 40) THEN $
              title = append_number('Subtraction Result')
            good = WHERE(result NE outbadval)
            IF(good(0) NE -1) THEN BEGIN
              scale_min = MIN(result(good), MAX = scale_max)
            ENDIF ELSE BEGIN
              scale_min = outbadval
              scale_max = outbadval
            ENDELSE
          END
      'MUL' : BEGIN
            IF(NOT graphs) THEN BEGIN
              exstr = 'result = MULTIPLY(' + name1 + '.data, ' + name2 + $
                '.data, badval_in = badval_in, outbadval = outbadval)'
            ENDIF ELSE BEGIN
              exstr = 'result = MULTIPLY(' + name1 + '.data(0:n1,*), ' + $
                name2 + '.data(0:n2,*), badval_in = badval_in, ' + $
                'outbadval = outbadval)'
            ENDELSE
            IF ((link1 NE -1) AND (link2 NE -1)) then begin
                 newlink=mkweight(link1,link2)
            ENDIF
            j = EXECUTE(exstr)
            sz = SIZE(result)
            IF((sz(0) NE 2) OR (sz(1) LT 2)) THEN BEGIN
              IF(uimage_version EQ 2) THEN GOTO, get_argument ELSE RETURN
            ENDIF
            title = '"' + title1 + '" * "' + title2 + '"'
            IF(STRLEN(title) GT 40) THEN $
              title = append_number('Multiplication Result')
            IF(NOT graphs) THEN BEGIN
              good = WHERE(result NE outbadval)
              IF(good(0) NE -1) THEN BEGIN
                sig = STDEV(result(good), mean)
                scale_min = mean - 2.*sig
                scale_max = mean + 2.*sig
              ENDIF ELSE BEGIN
                scale_min = outbadval
                scale_max = outbadval
              ENDELSE
            ENDIF
          END
      'DIV' : BEGIN
            IF(NOT graphs) THEN BEGIN
              exstr = 'result = DIVIDE(' + name1 + '.data, ' + name2 + $
                '.data, badval_in = badval_in, outbadval = outbadval)'
            ENDIF ELSE BEGIN
              exstr = 'result = DIVIDE(' + name1 + '.data(0:n1,*), ' + $
                name2 + '.data(0:n2,*), badval_in = badval_in, ' + $
                'outbadval = outbadval)'
            ENDELSE
            IF ((link1 NE -1) AND (link2 NE -1)) then begin
                 newlink=mkweight(link1,link2)
            ENDIF
            j = EXECUTE(exstr)
            sz = SIZE(result)
            IF((sz(0) NE 2) OR (sz(1) LT 2)) THEN BEGIN
              IF(uimage_version EQ 2) THEN GOTO, get_argument ELSE RETURN
            ENDIF
            title = '"' + title1 + '" / "' + title2 + '"'
            IF(STRLEN(title) GT 40) THEN title=append_number('Division Result')
            IF(NOT graphs) THEN BEGIN
              good = WHERE(result NE outbadval)
              IF(good(0) NE -1) THEN BEGIN
                sig = STDEV(result(good), mean)
                scale_min = mean - 2.*sig
                scale_max = mean + 2.*sig
              ENDIF ELSE BEGIN
                scale_min = outbadval
                scale_max = outbadval
              ENDELSE
            ENDIF
          END
    ENDCASE
;
;  If the result is an image, then store that image in the UIMAGE data
;  environment, and then display it (if on an X-window terminal).
;  -------------------------------------------------------------------
    IF(NOT graphs) THEN BEGIN
      newname=setup_image(data=result,title=title,projection='SKY-CUBE', $
        badpixval=outbadval,scale_min=scale_min,scale_max=scale_max, $
        linkweight=newlink,parent=name1)
      IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, newname
      IF(name_proj NE '') THEN BEGIN
        reproj,result,result,proj=proj,coor=coor,psize=psize,face=faceno,/noshow
        title = 'Reprojected ' + title
        IF(STRLEN(title) GT 40) THEN title = append_number('Reprojected result')
        newnameproj=setup_image(data=result,title=title, $
          projection=projection,coordinate_system=coordinate_system, $
          badpixval=outbadval,scale_min=scale_min,scale_max=scale_max, $
          parent=newname)
        IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, newnameproj
      ENDIF
    ENDIF ELSE BEGIN
;
;  If the result is an graph, then store that graph in the UIMAGE data
;  environment, and then display it.
;  -------------------------------------------------------------------
      good = WHERE(result(*, 1) NE outbadval)
      IF(good(0) EQ -1) THEN BEGIN
        PRINT, 'Resulting graph consists of all bad values - operation ' + $
          'aborted.'
        GOTO, done
      ENDIF
      j = EXECUTE('x_title1 = ' + name1 + '.x_title')
      j = EXECUTE('x_title2 = ' + name1 + '.x_title')
      j = EXECUTE('y_title1 = ' + name1 + '.y_title')
      j = EXECUTE('y_title2 = ' + name1 + '.y_title')
      j = EXECUTE('psymbol1 = ' + name1 + '.psymbol')
      j = EXECUTE('psymbol2 = ' + name2 + '.psymbol')
      IF(x_title1 EQ x_title2) THEN x_title = x_title1 ELSE x_title = ''
      IF(y_title1 EQ y_title2) THEN y_title = y_title1 ELSE x_title = ''
      IF(psymbol1 EQ psymbol2) THEN psymbol = psymbol1 ELSE psymbol = 0
      sz = SIZE(result)
      newname=setup_graph(result(*,0),result(*,1),sz(1),title=title, $
        /special,badpixval=outbadval,x_title=x_title,y_title=y_title,$
        nblocks_x=nblocks_x,nblocks_y=nblocks_y,psymbol=psymbol,$
        linkweight=newlink)
      plot_graph,newname
      IF(((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) $
        AND (uimage_version EQ 2)) THEN BEGIN
        IF(resized_graph EQ 0) THEN BEGIN
          sel = umenu(['Resize graph?', $
            'Yes / Resize the graph', 'No / Exit this menu'], title = 0)
          IF(sel EQ 1) THEN BEGIN
            resize_graph, newname
            j = EXECUTE('nblocks_x = ' + newname + '.nblocks_x')
            j = EXECUTE('nblocks_y = ' + newname + '.nblocks_y')
          ENDIF
          resized_graph = 1
        ENDIF
      ENDIF
    ENDELSE
    PRINT, bold('result') + ':  ' + title
;
;  If journaling is enabled, then write some info to the journal file.
;  -------------------------------------------------------------------
    IF(journal_on) THEN BEGIN
      PRINTF, luj, main_menu(selop)
      j = EXECUTE('newtitle = ' + newname + '.title')
      PRINTF, luj, '  the following statement was executed:'
      PRINTF, luj, '  ' + STRMID(exstr, 0, 76)
      exstr = STRMID(exstr, 76, STRLEN(exstr))
      WHILE(STRLEN(exstr) GE 1) DO BEGIN
        PRINTF, luj, '    ' + STRMID(exstr, 0, 74)
        exstr = STRMID(exstr, 74, STRLEN(exstr))
      ENDWHILE
      PRINTF, luj, '  where:'
      PRINTF, luj, '    ' + name1 + ' --> ' + title1
      PRINTF, luj, '    ' + name2 + ' --> ' + title2
      PRINTF, luj, '    badval_in = ', badval_in
      PRINTF, luj, '    result --> ' + newtitle
      PRINTF, luj, '----------------------------------------' + $
                   '--------------------------------------'
    ENDIF
    first = 0
    IF(uimage_version EQ 2) THEN GOTO, get_argument ELSE RETURN
  ENDELSE
;
;  For the weighted sum or average operation, call the UWEIGHTEDSUM procedure.
;  ---------------------------------------------------------------------------
wsum_or_avg:
  uweightedsum, operation, resized_graph
  GOTO, done
;
;  For the weighted addition, weighted subtraction or weighted average, call
;   the UWEIGHTAVG procedure.
;
w_avg:
 uweightavg, operation
;
;  The operation has been completed.
;  ---------------------------------
done:
  IF(uimage_version EQ 2) THEN GOTO, get_argument
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


