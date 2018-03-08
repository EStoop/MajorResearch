PRO xdisplay,name,rescale=rescale,xpos=xpos,ypos=ypos,noerase=noerase
;+
;  XDISPLAY - a UIMAGE-specific routine.  This routine displays a UIMAGE
;  image-object (map, face, or reprojected map) on an X-windows screen.
;#
;  Written by John Ewing.
;  SPR 10294  Dec 03 92  Adjusted to use new BLOCKSIZE_X.  J Ewing
;  SPR 10383  Jan 11 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10829  Apr  1993  Check if object should be hidden or not. J. Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON color_values,c_badpix,c_draw,c_scalemin
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    RETURN
  ENDIF
  j = EXECUTE('sz = SIZE(' + name + ')')
  IF((sz(2) NE 8) OR (STRMID(name, 0, 5) EQ 'GRAPH')) THEN BEGIN
    MESSAGE, 'An invalid object name was supplied.', /CONT
    RETURN
  ENDIF
  IF(KEYWORD_SET(rescale)) THEN rescale = 1 ELSE rescale = 0
  DEVICE, RETAIN = 2
;
;  Set the variable TYPE to contain a string which is the name of the
;  array-of-structures that contains the structure associated with the
;  operand.  Example: "MAP6".  Set the variable INDEX to contain an
;  integer number which tells the index of that element of the array-of-
;  structures which contains the operand's structure-instance.
;  ---------------------------------------------------------------------
  pos1 = STRPOS(name, '(')
  pos2 = STRPOS(name, ')')
  IF((pos1 EQ -1) OR (pos2 EQ -1)) THEN BEGIN
    MESSAGE, 'Invalid argument.', /CONT
    RETURN
  ENDIF
  type = STRMID(name, 0, pos1)
  res = 'RES' + STRMID(type, STRLEN(type) - 1, 1)
  index = FIX(STRMID(name, pos1 + 1, pos2 - pos1 - 1))
  j = EXECUTE('sz = SIZE(' + type + ')')
  IF(sz(2) NE 8) THEN BEGIN
    MESSAGE, 'An invalid or empty array-of-structures was designated.', /CONT
    RETURN
  ENDIF
;
; Check to see if object should be hidden
j=EXECUTE('hidden = ' + name + '.hidden')
IF (hidden EQ 1) THEN RETURN
;
;  Enforce an overall zoom factor of 1 if any res 9 maps are present.
;  ------------------------------------------------------------------
  size_map9 = SIZE(map9)
  IF(size_map9(2) EQ 8) THEN zoom_index = 0
;
  num_entries = sz(3)
  IF(num_entries LE index) THEN BEGIN
    MESSAGE, 'An invalid index number was supplied.'
    RETURN
  ENDIF
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('window = ' + name + '.window')
;
;  Handle the situation for 3-D objects.
;  -------------------------------------
  IF(STRMID(name, 0, 8) EQ 'OBJECT3D') THEN BEGIN
    show_object3d, name
    RETURN
  ENDIF
  j = EXECUTE('instrume = ' + name + '.instrume')
  j = EXECUTE('projection = ' + name + '.projection')
  IF(NOT defined(zoom_index)) THEN BEGIN
    MESSAGE, 'ZOOM_INDEX is undefined.', /CONT
    RETURN
  ENDIF
;
;  NUMBLOCKS_X and NUMBLOCKS_Y are defined below to indicate the size (in
;  blocks) for the windows associated with 9 different types of data-objects,
;  at four possible values of the overall window magnification factor.
;  Each sub-array of 9 values corresponds to a certain window magnification
;  factor.  The are four of these sub-arrays, for values of ZOOM_INDEX
;  ranging from 0 to 3.  The nine values within each sub-array are for
;  the following types of data-objects:  FACE6, FACE7, FACE8, FACE9, MAP6,
;  MAP7, MAP8, MAP9, and PROJ_MAP.
;  -------------------------------------------------------------------------
 numblocks_x=[[02,02,03,06,03,06,10,19,10,19],[02,03,06,10,06,10,19,99,19,99],$
              [03,06,10,99,10,19,99,99,99,99],[06,10,99,99,19,99,99,99,99,99]]
 numblocks_y=[[03,03,05,08,04,06,11,20,07,20],[03,05,08,14,06,11,20,99,13,99],$
              [05,08,14,99,11,20,99,99,99,99],[08,14,99,99,20,99,99,99,99,99]]
;
;  X0 and Y0 below tell the positions, within a window, that images of
;  various types, at various window mgnification factors, should be placed.
;  These have the same arrangement as described above for NBLOCKS_X and
;  NBLOCKS_Y.  These values were heuristically determined so that images
;  will basically be centered within their windows.
;  ------------------------------------------------------------------------
  x0 = [[26,10,06,28,06,28,14,18,14,18],[10,06,28,16,28,14,18,-1,20,-1], $
        [06,28,16,-1,14,18,-1,-1,-1,-1],[28,16,-1,-1,18,-1,-1,-1,-1,-1]]
  y0 = [[38,22,38,50,32,33,55,75,10,75],[22,38,50,65,33,55,75,-1,10,-1], $
        [38,50,65,-1,55,75,-1,-1,-1,-1],[50,65,-1,-1,75,-1,-1,-1,-1,-1]]

  CASE strupcase(type) OF
    'FACE6' : image_type = 0
    'FACE7' : image_type = 1
    'FACE8' : image_type = 2
    'FACE9' : image_type = 3
    'MAP6'  : image_type = 4
    'MAP7'  : image_type = 5
    'MAP8'  : image_type = 6
    'MAP9'  : image_type = 7
    'PROJ_MAP' : image_type = 8
    'PROJ2_MAP' : image_type = 9
  ENDCASE
  nblocks_x = numblocks_x(image_type, zoom_index)
  nblocks_y = numblocks_y(image_type, zoom_index)
;
;  If the RESCALE keyword was not set, then create the window in which
;  the image will be displayed.
;  -------------------------------------------------------------------
  IF(NOT rescale) THEN BEGIN
    make_window, window, nblocks_x, nblocks_y, title, error_status
    IF(error_status EQ -1) THEN RETURN
  ENDIF ELSE WSET, window
  xspace = blocksize_x * nblocks_x
  yspace = blocksize_y * nblocks_y
  windowsize_x = xspace - 20
  windowsize_y = yspace - 22 - 16
  j = EXECUTE('scale_min = ' + name + '.scale_min')
  j = EXECUTE('scale_max = ' + name + '.scale_max')
  j = EXECUTE('bad = WHERE(' + name + '.data EQ ' + name + '.badpixval)')
;
; Set border of image to white (255).
;
  IF(NOT KEYWORD_SET(noerase)) THEN BEGIN
    ERASE
    tmp=TVRD(0,0)
    tmp(*,*)=255
    tv,tmp
  ENDIF
;
;  Set the variable IMAGE to contain a BYTSCLed version of the data array.
;  -----------------------------------------------------------------------
  IF(scale_min LT scale_max) THEN BEGIN
    j = EXECUTE('vimage = BYTSCL(' + name + $
      '.data,min=scale_min,max=scale_max,top=255-c_scalemin)+c_scalemin')
  ENDIF ELSE BEGIN
    j = EXECUTE('vimage = BYTE(' + name + '.data*0.) + 255')
  ENDELSE
  IF(bad(0) NE -1) THEN vimage(bad) = c_badpix
  sz=SIZE(vimage)
  ndim = sz(0)
  dim1 = sz(1)
  dim2 = sz(2)
  IF((ndim EQ 2) AND ((dim2/3)*4 EQ dim1)) THEN BEGIN
     CASE dim1 OF
       128 : ires = 6
       256 : ires = 7
       512 : ires = 8
       1024 : ires = 9
     ENDCASE
     pixel=FINDGEN((dim1/2.)*dim2)
     pixel=LONG(pixel)
     tmp=pixel*0.+255
     pix2xy,pixel,res=ires,data=tmp,raster=out
     badp=WHERE(out eq 0)
     vimage(badp)=255
  endif
;
;  Call CONGRID so that the array contained within IMAGE is made to be
;  at the appropriate size.
;  -------------------------------------------------------------------
  zoom = 2.^zoom_index
  sz = SIZE(vimage)
  dim1 = sz(1)
  dim2 = sz(2)
  IF(zoom_index NE 0) THEN vimage = CONGRID(vimage, dim1*zoom, dim2*zoom)
;
;  Call TV to display the image data stored within the IMAGE variable.
;  -------------------------------------------------------------------
  xpos = x0(image_type, zoom_index)
  ypos = y0(image_type, zoom_index)
  TV, vimage, xpos, ypos

  j = EXECUTE(name + '.pos_x = xpos')
  j = EXECUTE(name + '.pos_y = ypos')
;
;  For data in the sky cube format, display the resolution size in the lower
;  left corner of the window.
;  -------------------------------------------------------------------------
  x = dim1 * zoom * .05
  IF((x / 2) * 2 NE x) THEN x = x + 1
  y = (dim2 * zoom * .05)<18
  IF((y / 2) * 2 NE y) THEN y = y + 1
  char_size = [1.2,1.2,1.6,1.8,1.3,1.4,1.9,3.2]
  index = (image_type + zoom_index) < 7
  IF(STRUPCASE(STRMID(name, 0, 8)) NE 'PROJ_MAP') AND $ 
    (STRUPCASE(STRMID(name, 0, 9)) NE 'PROJ2_MAP') THEN $
    XYOUTS, x, y, res, SIZE = char_size(index), /DEV, COLOR = 7
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


