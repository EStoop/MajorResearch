PRO integ_power
;+
;  INTEG_POWER - a UIMAGE-specific routine.  This routine allows the user
;  to create a 2-D object which is formed by summing together constant
;  frequency images from a selected 3-D object which are within a specified
;  range of indices.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 08 93  Change journaling.  J Ewing
;  SPR 10740  Mar 25 93  .LINK3D = 0 means no linkage.  J Ewing
;  SPR 10829  Apr  1993  Change info_3d common block
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
   freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON journal,journal_on,luj
;
;  Have the user select the desired 3-D object.
;  --------------------------------------------
  select_object3d, index3d
  IF(index3d EQ -1) THEN RETURN
  PRINT, ' '
  PRINT, 'Title of 3-D object:  ' + bold(object3d(index3d).title)
;
;  Have the user enter the desired range of indices.
;  -------------------------------------------------
get_range:
  PRINT, 'Please enter numbers to indicate the desired range of indices.'
  dim3 = object3d(index3d).dim3
  PRINT, 'The start index should be >= ' + bold('1') + $
    ' and the stop index should be <= ' + bold(STRTRIM(STRING(dim3), 2))
  str = ''
  val = 0
  WHILE(val EQ 0) DO BEGIN
    READ, bold('start index') + ': ', str
    val = validnum(str)
    IF(val EQ 1) THEN BEGIN
      istart = FIX(str)
      IF((istart LT 1) OR (istart GT dim3)) THEN BEGIN
        val = 0
        PRINT, 'Number not in proper range, please re-enter.'
      ENDIF
    ENDIF ELSE BEGIN
      PRINT, 'Invalid number, please re-enter.'
    ENDELSE
  ENDWHILE
  val = 0
  WHILE(val EQ 0) DO BEGIN
    READ, bold('stop index') + ':  ', str
    val = validnum(str)
    IF(val EQ 1) THEN BEGIN
      istop = FIX(str)
      IF((istop LT istart) OR (istop GT dim3)) THEN BEGIN
        val = 0
        PRINT, 'Number not in proper range, please re-enter.'
      ENDIF
    ENDIF ELSE BEGIN
      PRINT, 'Invalid number, please re-enter.'
    ENDELSE
  ENDWHILE
;
;  Sum the data over the specified range of indices.
;  -------------------------------------------------
  i0 = istart - 1
  i1 = istop - 1
  dataname = 'DATA3D_' + STRING(index3d, '(i1)')
  j = EXECUTE('result = ' + dataname + '(*, *, i0)')
  IF(i1 GT i1) THEN FOR i = 1 + i0, i1 DO $
    j = EXECUTE('result = result + ' + dataname + '(*, *, i)')
  badpixval = object3d(index3d).badpixval
  FOR i = 0, object3d(index3d).dim3 - 1 DO BEGIN
    j = EXECUTE('bad = WHERE(' + dataname + '(*, *, i) EQ badpixval)')
    IF(bad(0) NE -1) THEN result(bad) = badpixval
  ENDFOR
;
;  Store the resulting 2-D array in the UIMAGE data environment.
;  -------------------------------------------------------------
  title = 'Integrated power (indices ' + STRTRIM(STRING(istart), 2) + $
    ' through ' + STRTRIM(STRING(istop), 2) + ')'
  title=append_number(title)
  name = setup_image( $
    badpixval = object3d(index3d).badpixval, $
    coordinate_system = object3d(index3d).coordinate_system, $
    data = result, $
    faceno = object3d(index3d).faceno, $
    instrume = object3d(index3d).instrume, $
    link3d = index3d, $
    orient = object3d(index3d).orient, $
    projection = object3d(index3d).projection, $
    title = title, $
    units = object3d(index3d).units)
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, name
  IF(journal_on EQ 1) THEN BEGIN
    PRINTF, luj, 'Integrate over frequency/wavelength'
    PRINTF, luj, '  operand (3-D):  ' + object3d(index3d).title
    PRINTF, luj, '  <note:  indices start at 1>'
    PRINTF, luj, '  start index =   ' + STRTRIM(istart, 2)
    PRINTF, luj, '  stop index  =   ' + STRTRIM(istop, 2)
    PRINTF, luj, '  output (2-D):   ' + title
  ENDIF
  menu = ['Integrate another range from the same 3-D object?', 'Yes', 'No']
  sel = one_line_menu(menu, init = 2)
  IF(sel EQ 1) THEN GOTO, get_range
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


