PRO show_object3d, name
;+
;  SHOW_OBJECT3D - a UIMAGE-specific routine.  This routine will create a
;  window on a X-window screen to indicate the presence of a specified
;  3-D object (specified via the NAME argument, which would have a value
;  such as "OBJECT3D(0)").  The window will have a title which matches
;  the .TITLE field of the 3-D object's structure, and will show small
;  images of the first and last spectral slices of the 3-D data array.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10829  Apr  1993  Change info_3d common block
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;--------------------------------------------------------------------------
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON color_values,c_badpix,c_draw,c_scalemin
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    RETURN
  ENDIF
  pos1 = STRPOS(name, '(')
  pos2 = STRPOS(name, ')')
  IF((pos1 EQ -1) OR (pos2 EQ -1)) THEN BEGIN
    MESSAGE, 'Invalid argument.', /CONT
    RETURN
  ENDIF
  index = FIX(STRMID(name, pos1 + 1, pos2 - pos1 - 1))
;
;  Set up the display window.
;  --------------------------
  j = EXECUTE('window = ' + name + '.window') 
  j = EXECUTE('title = ' + name + '.title') 
  nblocks_x = 6
  nblocks_y = 2
  make_window, window, nblocks_x, nblocks_y, title, error_status
  IF(error_status EQ -1) THEN RETURN
;
;  Put out the "3-D Object" text.
;  ------------------------------
  XYOUTS, 12, 30, '!53-D!3', COLOR = 1, SIZE = 2.2, /DEV
  XYOUTS, 04, 08, '!5Object!3', COLOR = 1, SIZE = 2.2, /DEV
;
;  Determine size of small images to be displayed.
;  -----------------------------------------------
  j = EXECUTE('dim1 = ' + name + '.dim1') 
  j = EXECUTE('dim2 = ' + name + '.dim2')
  j = EXECUTE('dim3 = ' + name + '.dim3')
  j = EXECUTE('badpixval = '+name+'.badpixval')
  d2 = 46
  IF(dim1 EQ dim2) THEN d1 = 46 ELSE d1 = 62
  dataname = 'DATA3D_' + STRING(index, '(i1)')
;
;  Show a small image of the first slice.
;  --------------------------------------
  j = EXECUTE('slice = CONGRID(' + dataname + '(*,*,0),d1,d2)')
  good = WHERE(slice NE badpixval)
  bad = WHERE(slice EQ badpixval)
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(slice(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = badpixval
    scale_max = badpixval
  ENDELSE
  IF(scale_min LT scale_max) THEN BEGIN
    vimage = BYTSCL(slice, min = scale_min, max = scale_max, $
      top = 255 - c_scalemin) + c_scalemin
  ENDIF ELSE BEGIN
    vimage = BYTE(slice * 0.) + 255
  ENDELSE
  IF(bad(0) NE -1) THEN vimage(bad) = 0
  XYOUTS, 84, 36, 'First', COLOR=1, SIZE=1.2, /DEV
  XYOUTS, 84, 22, 'Slice:', COLOR=1, SIZE=1.2, /DEV
  TV, vimage, 120, 0
;
;  Show a small image of the second slice.
;  ---------------------------------------
  j = EXECUTE('slice = CONGRID(' + dataname + '(*,*,dim3-1),d1,d2)')
  good = WHERE(slice NE badpixval)
  bad = WHERE(slice EQ badpixval)
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(slice(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = badpixval
    scale_max = badpixval
  ENDELSE
  IF(scale_min LT scale_max) THEN BEGIN
    vimage = BYTSCL(slice, min = scale_min, max = scale_max, $
      top = 255 - c_scalemin) + c_scalemin
  ENDIF ELSE BEGIN
    vimage = BYTE(slice * 0.) + 255
  ENDELSE
  IF(bad(0) NE -1) THEN vimage(bad) = 0
  XYOUTS, 200, 36, 'Last', COLOR=1, SIZE=1.2, /DEV
  XYOUTS, 200, 22, 'Slice:', COLOR=1, SIZE=1.2, /DEV
  TV, vimage, 236, 0
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


