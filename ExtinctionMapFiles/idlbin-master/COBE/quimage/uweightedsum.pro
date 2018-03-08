PRO uweightedsum, operation, resized_graph
;+
;  UWEIGHTEDSUM - a UIMAGE-specific routine.  The code for this routine
;  was once part of UALGEBRA, but was split off separately into this
;  routine to keep UALGEBRA from getting too big.  This routine handles
;  the weighted-sum and averaging operations.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Changed journaling.  J Ewing
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
  IF(operation EQ 'WSM') THEN title = append_number('Weighted-sum Result') $
                         ELSE title = append_number('Averaging Result')
;
;  Have the user select the number of operands.
;  --------------------------------------------
  num_images = one_line_menu(['Select number of objects','1','2','3','4',$
    '5','6','7','8','9','10'])
  badval_in = FLTARR(num_images)
  coef = FLTARR(num_images + 1)
  names = STRARR(num_images)
  IF(operation EQ 'WSM') THEN BEGIN
    str = ' '
    val = 0
    PRINT, 'Please enter values for coefficients when prompted for below.'
    WHILE(val EQ 0) DO BEGIN
      READ, bold('C0') + ':  ', str
      val = validnum(str)
      IF(val EQ 1) THEN coef(0) = FLOAT(str) $
                   ELSE PRINT, 'Invalid number, please re-enter.'
    ENDWHILE
  ENDIF ELSE coef(0) = 0.
;
;  Have the user select the first operand.
;  ---------------------------------------
identify_1:
  names(0) = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
    /face8,/face9,/proj_map,/proj2_map,/simple_graph,/help,$
    /EXIT,title='Select X1')
  IF(names(0) EQ 'HELP') THEN BEGIN
    UIMAGE_HELP, menu_title
    GOTO, identify_1
  ENDIF
  IF(names(0) EQ 'EXIT') THEN RETURN
  parent = names(0)
  IF((STRMID(names(0), 0, 8) EQ 'PROJ_MAP') OR $
     (STRMID(names(0), 0, 9) EQ 'PROJ2_MAP')) THEN BEGIN
    name_proj = names(0)
    j = EXECUTE('win_orig = ' + names(0) + '.win_orig')
    names(0) = get_name(win_orig)
    j = EXECUTE('projection = ' + name_proj + '.projection')
    j = EXECUTE('coordinate_system = ' + name_proj + '.coordinate_system')
    j = EXECUTE('sz = SIZE(' + name_proj + '.data)')
    IF (sz(1) EQ 1024) THEN psize='L' ELSE psize='S'
    proj = STRUPCASE(STRMID(projection, 0, 1))
    IF(projection EQ 'GLOBAL SINUSOIDAL') THEN proj = 'S'
    coor = STRUPCASE(STRMID(coordinate_system, 0, 1))
    IF(STRUPCASE(coordinate_system) EQ 'EQUATORIAL') THEN coor = 'Q'
    IF(STRMID(names(0), 0, 4) EQ 'FACE') THEN $
      j = EXECUTE('faceno = ' + names(0) + '.faceno')
  ENDIF ELSE name_proj = ''
  IF(STRMID(names(0), 0, 5) EQ 'GRAPH') THEN graphs = 1 ELSE graphs = 0
  j = EXECUTE('badval_in(0) = ' + names(0) + '.badpixval')
  j = EXECUTE('temptitle = ' + names(0) + '.title')
  PRINT, bold('X1') + ':  ' + temptitle
  IF(operation EQ 'WSM') THEN BEGIN
    str = ' '
    val = 0
    WHILE(val EQ 0) DO BEGIN
      READ, bold('C1') + ':  ', str
      val = validnum(str)
      IF(val EQ 1) THEN coef(1) = FLOAT(str) $
                   ELSE PRINT, 'Invalid number, please re-enter.'
    ENDWHILE
  ENDIF ELSE coef(1) = 1. / num_images
  types = STRMID(parent, 0, STRPOS(parent, '('))
;
;  Have the user select all other operands such that their data arrays
;  match the size of the data array for the first operand.
;  -------------------------------------------------------------------
  IF(num_images GT 1) THEN BEGIN
    FOR i = 2, num_images DO BEGIN
identify_n:
      str_i = STRMID(STRCOMPRESS(STRING(i)), 1, 1)
      titlen = 'Select X' + str_i
      j = EXECUTE('names(i-1)=select_object(/' + types + ',/help,/EXIT,' + $
        'title=titlen)') 
      IF(names(i - 1) EQ 'HELP') THEN BEGIN
        UIMAGE_HELP, 'Select Xn'
        GOTO, identify_n
      ENDIF
      IF(names(i - 1) EQ 'EXIT') THEN RETURN
      parentn = names(i - 1)
      IF((STRMID(names(i - 1), 0, 8) EQ 'PROJ_MAP') OR $
         (STRMID(names(i - 1), 0, 9) EQ 'PROJ2_MAP')) THEN BEGIN
        j = EXECUTE('win_orig = ' + names(i - 1) + '.win_orig')
        names(i - 1) = get_name(win_orig)
      ENDIF
      j = EXECUTE('badval_in(i - 1) = ' + names(i - 1) + '.badpixval')
      j = EXECUTE('temptitle = ' + parentn + '.title')
      PRINT, bold('X' + str_i) + ':  ' + temptitle
      IF(operation EQ 'WSM') THEN BEGIN
        str = ' '
        val = 0
        WHILE(val EQ 0) DO BEGIN
          READ, bold('C' + str_i) + ':  ', str
          val = validnum(str)
          IF(val EQ 1) THEN coef(i) = FLOAT(str) $
                       ELSE PRINT, 'Invalid number, please re-enter.'
        ENDWHILE
      ENDIF ELSE coef(i) = 1./num_images
    ENDFOR
  ENDIF
;
;  Call WEIGHTEDSUM, which will return the desired summed result.
;  --------------------------------------------------------------
  exstr = 'result = WEIGHTEDSUM('
  IF(NOT graphs) THEN FOR i=0,num_images-1 DO exstr=exstr+names(i)+'.data, ' $
                 ELSE FOR i=0,num_images-1 DO $
           exstr = exstr + names(i) + '.data(0:' + names(i) + '.num-1,*), '
  exstr = exstr + 'coef = coef, badval_in = badval_in, outbadval = outbadval)'
  j = EXECUTE(exstr)
  sz = SIZE(result)
  IF (sz(0) NE 2) THEN RETURN
  IF ((sz(1) lt 2) OR (sz(2) lt 2)) THEN RETURN
;
;  If the result is an image, then store it in an image data-object,
;  and display that image if run on an X-window terminal.
;  -----------------------------------------------------------------
  IF(NOT graphs) THEN BEGIN
    good = WHERE(result NE outbadval)
    IF(good(0) NE -1) THEN BEGIN
      scale_min = MIN(result(good), MAX = scale_max)
    ENDIF ELSE BEGIN
      scale_min = outbadval
      scale_max = outbadval
    ENDELSE
    newname=setup_image(data=result,title=title,projection='SKY-CUBE', $
      badpixval=outbadval,scale_min=scale_min,scale_max=scale_max,$
      parent=names(0))
    IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, newname
;
;  If the operands were reprojected maps, then reproject the result.
;  -----------------------------------------------------------------
    IF(name_proj NE '') THEN BEGIN
      reproj,result,result,proj=proj,coor=coor,psize=psize,face=faceno,/noshow
      newnameproj=setup_image(data=result,title=title, $
        projection=projection,coordinate_system=coordinate_system, $
        badpixval=outbadval,scale_min=scale_min,scale_max=scale_max,$
        parent=newname)
      IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, newnameproj
    ENDIF
  ENDIF ELSE BEGIN
;
;  If the result is a graph, then store it in a graph data-object,
;  and then plot the graph.
;  ---------------------------------------------------------------
    good = WHERE(result(*, 1) NE outbadval)
    IF(good(0) EQ -1) THEN BEGIN
      PRINT, 'Resulting graph consists of all bad values - operation ' + $
        'aborted.'
      RETURN
    ENDIF
    j = EXECUTE('x_title = ' + names(0) + '.x_title')
    j = EXECUTE('y_title = ' + names(0) + '.y_title')
    j = EXECUTE('psymbol = ' + names(0) + '.psymbol')
    FOR i = 1, num_images - 1 DO BEGIN
      j = EXECUTE('temp = ' + names(i) + '.x_title')
      IF(temp NE x_title) THEN x_title = ''
      j = EXECUTE('temp = ' + names(i) + '.y_title')
      IF(temp NE x_title) THEN y_title = ''
      j = EXECUTE('temp = ' + names(i) + '.psymbol')
      IF(temp NE psymbol) THEN psymbol = 0
    ENDFOR
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
  PRINT, bold('result') + ':  ' + title
;
;  If journaling is enabled, then send out info to the journal file.
;  -----------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    IF(operation EQ 'WSM') THEN PRINTF, luj, 'C0 + C1*X1 + C2*X2 +...+ Cn*Xn' $
      ELSE PRINTF, luj, 'Average(X1, X2, X3, ..., Xn)'
    j = EXECUTE('newtitle = ' + newname + '.title')
    PRINTF, luj, '  the following statement was executed:'
    PRINTF, luj, '  ' + STRMID(exstr, 0, 76)
    exstr = STRMID(exstr, 76, STRLEN(exstr))
    WHILE(STRLEN(exstr) GE 1) DO BEGIN
      PRINTF, luj, '    ' + STRMID(exstr, 0, 74)
      exstr = STRMID(exstr, 74, STRLEN(exstr))
    ENDWHILE
    PRINTF, luj, '  where:'
    FOR i = 0, num_images - 1 DO BEGIN
      j = EXECUTE('tit = ' + names(i) + '.title')
      PRINTF, luj, '    ' + names(i) + ' --> ' + tit
    ENDFOR
    PRINTF, luj, '    coef = ', coef
    PRINTF, luj, '    badval_in = ', badval_in
    PRINTF, luj, '    result --> ' + newtitle
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


