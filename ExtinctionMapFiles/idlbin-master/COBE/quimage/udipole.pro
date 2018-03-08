PRO udipole, name
;+
;  UDIPOLE - a UIMAGE-specific routine.  This is UIMAGE's driver routine
;  for the DIPOLE procedure.
;#
;  Written by John Ewing.
;  SPR 10294  Dec 03 92  Won't work on a reprojected face.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 10829  Apr  1993  If data has an associated weight object then
;                        weighting is automatically used. J. Newmark
;  SPR 10835  Apr 21 93  Fix call to DIPOLE for new parameters. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON journal,journal_on,luj
  COMMON xscreen,magnify,block_usage,scrsiz
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = 'Subtract a dipole'
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Put up a menu of potential operands.
;  ------------------------------------
show_menu:
  name = select_object(/map6,/map7,/map8,/map9,/pmap_noface, $
           /pmap2_noface,/help,/exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no sky-cubes available.', /CONT
    RETURN
  ENDIF
;
;  An operand has been selected - make sure it is valid.
;  -----------------------------------------------------
have_name:
  j = EXECUTE('data_min = ' + name + '.data_min')
  j = EXECUTE('data_max = ' + name + '.data_max')
  IF(data_min EQ data_max) THEN BEGIN
    PRINT, 'Image has no dynamic range.'
    IF(uimage_version EQ 2) THEN GOTO, show_menu ELSE RETURN
  ENDIF
  IF(data_min GT data_max) THEN BEGIN
    PRINT, 'Image consists of all bad pixels.'
    IF(uimage_version EQ 2) THEN GOTO, show_menu ELSE RETURN
  ENDIF
  parent = name
  IF ((STRMID(name, 0, 8) EQ 'PROJ_MAP') OR $
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
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('badpixval = ' + name + '.badpixval')
;
;  Find out what will be the galactic exclusion angle.
;  ---------------------------------------------------
  PRINT, 'Enter a galactic exclusion angle below in degrees (default = 20).'
  str = ' '
  val = 0
  WHILE(val EQ 0) DO BEGIN
    READ, bold('angle') + ':  ', str
    IF(str EQ '') THEN BEGIN
      val = 1
      galexc = 20.
      PRINT, 'Using default value of 20 degrees.'
    ENDIF ELSE BEGIN
      val = validnum(str)
      IF(val EQ 1) THEN BEGIN
        galexc = FLOAT(str)
        IF((galexc LT 0.) OR (galexc GT 90.)) THEN BEGIN
          val = 0
          PRINT, 'Value must be between 0 and 90, please re-enter.'
        ENDIF
      ENDIF ELSE BEGIN
        PRINT, 'Invalid number, please re-enter.'
      ENDELSE
    ENDELSE
  ENDWHILE
;
;  Find out if an array of weights will be used.
;  ---------------------------------------------
  wsel=0
  j=execute('linkweight= ' +name+ '.linkweight')
  IF( linkweight ne -1) THEN BEGIN
    wsel=1
    namew=get_name(linkweight)
    j = EXECUTE('titlew = ' + namew + '.title')
    PRINT, 'Weights will be taken from ' + bold(titlew)
    j = EXECUTE('badvalw = ' + namew + '.badpixval')
    j = EXECUTE('weights = ' + namew + '.data')
    goodw = WHERE(weights NE badvalw)
    badw = WHERE(weights EQ badvalw)
    IF(goodw(0) NE -1) THEN BEGIN
      minw = MIN(weights(goodw))
      IF(badw(0) NE -1) THEN weights(badw) = minw
    ENDIF
    exstr = 'dipole, ' + name + '.data, monopole, mono_sig, diamp1,' + $
      'diamp1_sig, diglon, diglat, digl_sig,' + $
      'galexc = galexc, residmap = residmap, badval = badpixval, ' + $
      'weights = weights'
  ENDIF ELSE BEGIN
    exstr = 'dipole, ' + name + '.data, monopole, mono_sig, diamp1,' + $
      'diamp1_sig, diglon, diglat, digl_sig,' + $
      'galexc = galexc, residmap = residmap, badval = badpixval'
  ENDELSE
;
;  Call DIPOLE.
;  ------------
  residmap=-1
  j = EXECUTE(exstr)
  bsz=SIZE(residmap)
  IF (bsz(0) EQ 0) THEN RETURN
;
;  Store the result in a UIMAGE data-object.
;  -----------------------------------------
  j = EXECUTE('title = ' + name + '.title')
  titleout = '(' + title + ') - DIPOLE'
  IF(STRLEN(titleout) GT 40) THEN $
    titleout = append_number('Dipole subtraction result')
  good = WHERE(residmap NE 0.)
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(residmap(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = 0.
    scale_max = 0.
  ENDELSE
  nameout = setup_image(data = residmap, title = titleout, badpixval = 0., $
    scale_min = scale_min, scale_max = scale_max, parent = parent)
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, nameout
;
;  If a reprojected data-object was the operand, then reproject the result.
;  ------------------------------------------------------------------------
  IF(name_proj NE '') THEN BEGIN
    reproj, residmap, residmap, proj = proj, coor = coor, psize = psize, $
      face = faceno, /noshow
    titleout = 'Reprojected ' + titleout
    IF(STRLEN(titleout) GT 40) THEN $
      titleout = append_number('Reprojected dipole subtr. result')
    good = WHERE(residmap NE 0.)
    IF(good(0) NE -1) THEN BEGIN
      scale_min = MIN(residmap(good), MAX = scale_max)
    ENDIF ELSE BEGIN
      scale_min = 0.
      scale_max = 0.
    ENDELSE
    nameoutproj = setup_image(data = residmap, title = titleout, $
      badpixval = 0., scale_min = scale_min, scale_max = scale_max, $
      parent = nameout)
    IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, nameoutproj
  ENDIF
;
;  Print out details of the calculations.
;  --------------------------------------
  text = STRARR(7)
  text(0) = ' '
  text(1) = 'Mean CMBR intensity:  ' + STRTRIM(STRING(monopole, '(g18.3)'), 2)
  text(2) = 'Sigma of CMBR inten.:  ' + STRTRIM(STRING(mono_sig, '(g18.3)'), 2)
  text(3) = 'Dipole amplitude:     ' + STRTRIM(STRING(diamp1, '(g18.6)'), 2)
  text(4) = 'Sigma of Dipole amp.: ' + STRTRIM(STRING(diamp1_sig, '(g18.6)'), 2)
  text(5) = 'Galactic lon & lat of dipole pole (degrees):  ' + $
    STRTRIM(STRING(diglon, '(f9.2)'), 2) + ', ' + $
    STRTRIM(STRING(diglat, '(f9.2)'), 2)
  text(6) = 'Sigma of Dipole dir.: ' + STRTRIM(STRING(digl_sig, '(f9.2)'), 2)
 FOR i = 0, 6 DO PRINT, text(i)
;
;  If journaling is enabled, then write info to the journal file.
;  --------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title
    PRINTF, luj, '  the following statement was executed:'
    PRINTF, luj, '  ' + STRMID(exstr, 0, 76)
    exstr = STRMID(exstr, 76, STRLEN(exstr))
    WHILE(STRLEN(exstr) GE 1) DO BEGIN
      PRINTF, luj, '    ' + STRMID(exstr, 0, 74)
      exstr = STRMID(exstr, 74, STRLEN(exstr))
    ENDWHILE
    PRINTF, luj, '  where:'
    PRINTF, luj, '    galexc    = ' + STRTRIM(STRING(galexc), 2)
    PRINTF, luj, '    badpixval = ' + STRTRIM(STRING(badpixval), 2)
    IF(wsel EQ 1) THEN BEGIN
      PRINTF, luj, '    weights came from:  ' + titlew
    ENDIF
    FOR i = 1, 6 DO PRINTF, luj, '  ' + text(i)
    PRINTF, luj, '  residmap --> ' + titleout
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


