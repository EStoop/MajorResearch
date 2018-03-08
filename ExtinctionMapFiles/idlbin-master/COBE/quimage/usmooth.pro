PRO usmooth, name
;+
;  USMOOTH - a UIMAGE-specific routine.  This routine serves as UIMAGE's
;  driver to the SMOOTHCUBE routine,  allowing the creation of a data-
;  object containing a smoothed image derived from a sky cube, a sky
;  cube face, or a reprojected map.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 04 93  Change journaling.  J Ewing
;  SPR 10829  Apr 04 93  Add weighting to calculations, if data has an
;                        associated weight object. J. Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON smooth,kernel
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = 'Smooth an image'
  IF(N_PARAMS(0) GE 1) THEN GOTO, have_name
show_menu:
;
;  Put up a menu in which the user can identify the desired image.
;  ---------------------------------------------------------------
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,/face8,$
    /face9,/proj_map,/proj2_map,/help,/exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no appropriate maps available.', /CONT
    RETURN
  ENDIF
;
;  Extract some fields from the operand's structure.
;  -------------------------------------------------
have_name:
  parent = name
  IF((STRMID(name, 0, 8) EQ 'PROJ_MAP') OR $
     (STRMID(name, 0, 9) EQ 'PROJ2_MAP')) THEN BEGIN
    name_proj = name
    j = EXECUTE('win_orig = ' + name + '.win_orig')
    name = get_name(win_orig)
    j = EXECUTE('projection = ' + name_proj + '.projection')
    j = EXECUTE('coordinate_system = ' + name_proj + '.coordinate_system')
    j = EXECUTE('sz = SIZE(' + name_proj + '.data)')
    IF (sz(1) EQ 1024) THEN psize='L' ELSE psize='S'
    proj = STRUPCASE(STRMID(projection, 0, 1))
    IF(projection EQ 'GLOBAL SINUSOIDAL') THEN proj = 'S'
    coor = STRUPCASE(STRMID(coordinate_system, 0, 1))
    IF(STRUPCASE(coordinate_system) EQ 'EQUATORIAL') THEN coor = 'Q'
    IF(STRMID(name,0,4) EQ 'FACE') THEN j=EXECUTE('faceno = '+name+'.faceno')
  ENDIF ELSE name_proj = ''
  j = EXECUTE('orient = ' + name + '.orient')
  j = EXECUTE('title = ' + parent + '.title')
  IF(orient NE 'R') THEN BEGIN
    MESSAGE, 'Sky-cube is not in the RIGHT orientation.', /CONT
    RETURN
  ENDIF
  PRINT, ' '
  PRINT, 'Selected image:  ' + bold(title)
  filename = ''
;
;  If a kernel was defined previously within this IDL session, retained
;  via the "smooth" COMMON block, then ask the user if he'd like to
;  re-use that same kernel.
;  --------------------------------------------------------------------
  IF(defined(kernel)) THEN BEGIN
    menu = ['Re-use the previously-defined kernel?', 'Yes', 'No']
    reuse = one_line_menu(menu, init = 1)
    IF(reuse EQ 1) THEN GOTO, do_smooth
  ENDIF
;
;  Ask the user if he wants to type in the kernel, or if he'd like to
;  read one in from an IDL save set.
;  ------------------------------------------------------------------
re_enter:
  menu=['How should kernel be defined?','Keyboard entry','Read from a file']
  mode = one_line_menu(menu, init = 1)
  IF(mode EQ 2) THEN GOTO, get_from_file
;
;  Have the user select the size of the kernel.
;  --------------------------------------------
  menu = ['Select dimension for a square kernel', '3', '5', '7']
  kdim1 = FIX(menu(one_line_menu(menu, init = 2)))
;
;  Show the user how the kernel will look with regards to symmetry.
;  ----------------------------------------------------------------
  CASE kdim1 OF
    3: BEGIN & x1 = 2 & x2 = 4 & list = 'A,B,C' & END
    5: BEGIN & x1 = 1 & x2 = 5 & list = 'A,B,C,D,E,F' & END
    7: BEGIN & x1 = 0 & x2 = 6 & list = 'A,B,C,D,E,F,G,H,I,J' & END
  ENDCASE
  y1 = x1
  y2 = x2
  r0 = ['       J','       I','       H','       G','       H','       I',$
        '       J']
  r1 = ['       I','       F','       E','       D','       E','       F',$
        '       I']
  r2 = ['       H','       E','       C','       B','       C','       E',$
        '       H']
  r3 = ['       G','       D','       B','       A','       B','       D',$
        '       G']
  r4 = ['       H','       E','       C','       B','       C','       E',$
        '       H']
  r5 = ['       I','       F','       E','       D','       E','       F',$
        '       I']
  r6 = ['       J','       I','       H','       G','       H','       I',$
        '       J']
  matrix = [[r0], [r1], [r2], [r3], [r4], [r5], [r6]]
  PRINT, 'Shape of kernel:'
  PRINT, ' '
  PRINT, '', matrix(x1:x2, y1:y2)
  PRINT, ' '
  PRINT, 'Enter values for ' + list + ' below.'
  a = 0 & b = 0 & c = 0 & d = 0 & e = 0 & f = 0 & g = 0 & h = 0 & i = 0 & j = 0
  str = ' '
;
;  Have the user enter in values for elements of the kernel.
;  ---------------------------------------------------------
  WHILE(list NE '') DO BEGIN
    varname = STRMID(list, 0, 1)
    list = STRMID(list, 2, STRLEN(list))
    val = 0
    WHILE(val EQ 0) DO BEGIN
      READ, varname + ':  ', str
      val = validnum(str)
      IF(val) THEN dummy = EXECUTE(varname + '= FLOAT(str)') $
                   ELSE PRINT, 'Invalid number, please re-enter.'
    ENDWHILE
  ENDWHILE
;
;  Create the kernel with the specified values.
;  --------------------------------------------
  array=[[j,i,h,g,h,i,j],[i,f,e,d,e,f,i],[h,e,c,b,c,e,h],[g,d,b,a,b,d,g],$
    [h,e,c,b,c,e,h],[i,f,e,d,e,f,i],[j,i,h,g,h,i,j]]
  kernel = array(x1:x2, y1:y2)
  PRINT, ' '
  PRINT, 'kernel:'
  PRINT, ' '
  PRINT, kernel
  PRINT, ' '
  menu = ['Do you want to re-enter those values?', 'Yes', 'No']
  redo = one_line_menu(menu, init = 2)
  IF(redo EQ 1) THEN BEGIN
    PRINT, ' '
    GOTO, re_enter
  ENDIF ELSE BEGIN
    GOTO, do_smooth
  ENDELSE
get_from_file:
;
;  Have the user enter in a file name of an IDL saveset which is supposed to
;  contain a 2-D floating point array named KERNEL.
;  -------------------------------------------------------------------------
  err = 1
  WHILE(err NE 0) DO BEGIN
stay_in_loop:
    PRINT, 'Enter file name of an IDL saveset containing a 2-D floating point array'
    PRINT, '  named KERNEL, or else enter "E" to exit.  (Note - program blowup will'
    PRINT, '  occur if the indicated file is not an IDL saveset.)'
    READ, underline('File name:') + '  ', filename
    kernel = 0
    filename = STRUPCASE(filename)
    IF(filename EQ 'E') THEN GOTO, show_menu
    OPENR, lun, filename, ERROR = err, /GET_LUN
    IF(err EQ 0) THEN BEGIN
      FREE_LUN, lun
      RESTORE, filename
      sz = SIZE(kernel)
      IF(sz(0) NE 2) THEN BEGIN
        PRINT, bold('KERNEL is not a 2-D array.')
        PRINT, ' '
        GOTO, stay_in_loop
      ENDIF
      IF(sz(3) NE 4) THEN BEGIN
        PRINT, bold('KERNEL is not of type FLOAT.')
        PRINT, ' '
        GOTO, stay_in_loop
      ENDIF
      kdim1 = sz(1)
      kdim2 = sz(2)
      IF(kdim1 NE kdim2) THEN BEGIN
        PRINT, bold('KERNEL is not a square array.')
        PRINT, ' '
        GOTO, stay_in_loop
      ENDIF
    ENDIF ELSE BEGIN
      PRINT, bold('That file could not be opened.')
      PRINT, ' '
    ENDELSE
  ENDWHILE
;
;  Make a call to SMOOTHCUBE, which will return a smoothed image.
;  --------------------------------------------------------------
do_smooth:
  j = EXECUTE('badpixval = ' + name + '.badpixval')
  j = EXECUTE('good = WHERE(' + name + '.data NE badpixval)')
  IF(good(0) EQ -1) THEN BEGIN
    MESSAGE, 'Image consists of all bad pixels.', /CONT
    RETURN
  ENDIF
  j=execute('linkweight= '+name+'.linkweight')
  IF (linkweight ne -1 ) THEN BEGIN
;
; If there are weights then Result= SMOOTH(data*weight)/SMOOTH(weight)
;  
    namew=get_name(linkweight)
    parent=namew
    j = EXECUTE('titlew = ' + parent + '.title')
    print, 'Weights will be taken from ' + bold(titlew)
    badval_in = FLTARR(2)
    j = EXECUTE('badval_in(0) = ' + name + '.badpixval')
    j = EXECUTE('badval_in(1) = ' + namew + '.badpixval')
; multiply data by weights
    exstr = 'result1 = MULTIPLY(' + name + '.data, ' + namew + $
        '.data, badval_in = badval_in, outbadval = badval1)'
    j = EXECUTE(exstr)
; create smoothed array=data*weight
    exstr = 'result2 = SMOOTHCUBE(result1, kernel, bad = badval1)'
    j = EXECUTE(exstr)
; create smoothed weight array and setup new weight structure
    j = EXECUTE('badval2 = ' + namew + '.badpixval')
    exstr = 'result3 = SMOOTHCUBE(' + namew + '.data,kernel,bad=badval2)'
    j = EXECUTE(exstr)
    newtitle = 'Smoothed version of "' + titlew + '"'
    IF(STRLEN(newtitle) GT 40) THEN newtitle=append_number('Smoothed Result')
    newname=setup_image(data=result3,title=newtitle,bad=badval2, $
      hidden=1,parent=parent)
    j=EXECUTE('newlink= '+ newname + '.window')
; create output array= smoothed (data*wt)/ smoothed wt
    badval_in(0)=badval1
    badval_in(1)=badval2
    exstr = 'result = DIVIDE(result2,result3, badval_in = badval_in,' + $
        'outbadval = outbadval)'
    j = EXECUTE(exstr)
;
;  Store the result in a UIMAGE data-object, and display it if an X-window
;  terminal is being used.
;  -----------------------------------------------------------------------
    parent=name
    newtitle = 'Smoothed version of "' + title + '"'
    IF(STRLEN(newtitle) GT 40) THEN newtitle=append_number('Smoothed Result')
    newname=setup_image(data=result,title=newtitle,bad=outbadval,$
         linkweight=newlink,parent=parent)
    IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, newname
;
  ENDIF ELSE BEGIN
; no weighting
    exstr = 'result = SMOOTHCUBE(' + name + '.data, kernel, bad = badpixval)'
    j = EXECUTE(exstr)
;
;  Store the result in a UIMAGE data-object, and display it if an X-window
;  terminal is being used.
;  -----------------------------------------------------------------------
    newtitle = 'Smoothed version of "' + title + '"'
    IF(STRLEN(newtitle) GT 40) THEN newtitle=append_number('Smoothed Result')
    newname=setup_image(data=result,title=newtitle,bad=badpixval,parent=parent)
    IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, newname
  ENDELSE
;
;  If the operand was a reprojected map, then create a smoothed reprojected
;  map as output, in addition to the smoothed sky cube.
;  ------------------------------------------------------------------------
  IF(name_proj NE '') THEN BEGIN
    reproj,result,result,proj=proj,coor=coor,psize=psize,face=faceno,/noshow
    newtitle = 'Reprojected ' + newtitle
    IF(STRLEN(newtitle) GT 40) THEN newtitle = $
     append_number('Reprojected smoothed result')
    newnameproj = setup_image(data=result,title=newtitle,bad=badpixval,$
      parent=newname,proj=projection)
    IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, newnameproj
  ENDIF
;
;  If journaling is enabled, then send out some info to the journal file.
;  ----------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title
    PRINTF, luj, '  the following statement was executed:'
    PRINTF, luj, '  ' + exstr
    PRINTF, luj, '  where:'
    PRINTF, luj, '    ' + name + " is the name of the operand's structure"
    IF(filename EQ '') THEN BEGIN
      PRINTF, luj, '    kernel: '
      PRINTF, luj, kernel
    ENDIF ELSE BEGIN
      PRINTF, luj, '    kernel was read in from ' + filename
    ENDELSE
    PRINTF, luj, '    badpixval = ', STRTRIM(badpixval, 2)
    PRINTF, luj, '    result --> ', newtitle
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
  ENDIF
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


