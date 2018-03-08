PRO extract_slice
;+
;  EXTRACT_SLICE - a UIMAGE-specific routine.  This routine allows a user
;  to 2-D data-object from a selected 3-D data object after indicating
;  a desired index number.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 08 93  Change journaling.  J Ewing
;  SPR 10740  Mar 25 93  .LINK3D = -1 means no linkage.  J Ewing
;--------------------------------------------------------------------------
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,freq3D_2
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map
  COMMON journal,journal_on,luj
  first = 1
;
;  Have the user select a 3-D object.
;  ----------------------------------
  select_object3d, index3d
  IF(index3d EQ -1) THEN RETURN
  IF(journal_on EQ 1) THEN PRINTF, luj, 'Extract slice of 3-D object'
  dim1 = object3d(index3d).dim1
  dim2 = object3d(index3d).dim2
  dim3 = object3d(index3d).dim3
  maxslice = STRTRIM(STRING(dim3), 2)
  PRINT, ' '
  PRINT, 'Title of 3-D object:  ' + bold(object3d(index3d).title)
;
;  Have the user enter the index of the desired slice.
;  ---------------------------------------------------
get_sliceindex:
  PRINT, 'Please indicate the index number of the desired slice.'
  PRINT, '(This number should be between ' + bold('1') + ' and ' + $
    bold(maxslice) + ').'
  str = ''
  val = 0
  WHILE(val EQ 0) DO BEGIN
    READ, bold('index') + ':  ', str
    val = validnum(str)
    IF(val EQ 1) THEN BEGIN
      index = FIX(str)
      IF((index LT 0) OR (index GT maxslice)) THEN BEGIN
        val = 0
        PRINT, 'Index value is not within the proper range.'
      ENDIF
    ENDIF ELSE BEGIN
      PRINT, 'Invalid number, please re-enter.'
    ENDELSE
  ENDWHILE
;
;  Extract the slice and put it in the UIMAGE data environment.
;  ------------------------------------------------------------
  title = 'Slice #' + STRTRIM(STRING(index), 2) + ' from "' + $
    object3d(index3d).title + '"'
  str = STRING(index3d, '(i1)')
  j = EXECUTE('data = data3D_' + str + '(*, *, index - 1)')
  j = EXECUTE('freq = freq3D_' + str + '(index - 1)')
  name = setup_image( $
    badpixval = object3d(index3d).badpixval, $
    coordinate_system = object3d(index3d).coordinate_system, $
    data = data, $
    faceno = object3d(index3d).faceno, $
    frequency = freq, $
    instrume = object3d(index3d).instrume, $
    link3d = index3d, $
    orient = object3d(index3d).orient, $
    projection = object3d(index3d).projection, $
    title = title, $
    units = object3d(index3d).units)
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, name
  IF(journal_on EQ 1) THEN BEGIN
    PRINTF, luj, '  operand (3-D):  ' + object3d(index3d).title
    PRINTF, luj, '  output (2-D):   ' + title
  ENDIF
  menu = ['Extract another slice from the same 3-D object?', 'Yes', 'No']
  sel = one_line_menu(menu, init = 2)
  IF(sel EQ 1) THEN GOTO, get_sliceindex
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


