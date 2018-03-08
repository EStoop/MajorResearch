FUNCTION setup_object3d, $
  badpixval = badpixval, $
  coordinate_system = coordinate_system, $
  data = data, $
  faceno = faceno, $
  frequency = frequency, $
  frequnits = frequnits, $
  weight3d = weight3d, $
  instrume = instrume, $
  linkweight=linkweight, $
  hidden = hidden, $
  orient = orient, $
  projection = projection, $
  title = title, $
  units = units, $
  wgtunits = wgtunits, $
  version = version
;+
;  SETUP_OBJECT3D - a UIMAGE-specific routine.  This routine supports the
;  instantiation of a 3-D data-object.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10829  Apr  1993  Add in .LINKWEIGHT and .HIDDEN fields. J. Newmark
;  SPR 10829  Apr  1993  Add a new structure, WT3D_0,1,2, which contains
;                        any 3D weights (e.g. SERROR for FIRAS) associated
;                        with particular data. J. Newmark
;  SPR 11491 10 Dec 1993  Add weight units field. J. Newmark
;--------------------------------------------------------------------------
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
;
;  Check if valid arguments were supplied.
;  ---------------------------------------
  IF(NOT KEYWORD_SET(data)) THEN MESSAGE, 'DATA must be supplied.', /CONT
  IF(NOT KEYWORD_SET(frequency)) THEN $
    MESSAGE, 'FREQUENCY/WAVELENGTH must be supplied.', /CONT
  siz = SIZE(data)
  IF(siz(0) NE 3) THEN BEGIN
    MESSAGE, 'DATA must be a 3-D array.', /CONT
    RETURN, 'ERROR'
  ENDIF
  dim1 = siz(1)
  dim2 = siz(2)
  dim3 = siz(3)
  siz = SIZE(frequency)
  IF(siz(0) NE 1) THEN BEGIN
    MESSAGE, 'FREQUENCY/WAVELENGTH must be a 1-D array.', /CONT
    RETURN, 'ERROR'
  ENDIF
  IF(siz(1) NE dim3) THEN BEGIN
    MESSAGE, 'The 1st dimension of FREQUENCY/WAVELENGTH must equal the 3rd ' + $
      'dimension of DATA.', /CONT
    RETURN, 'ERROR'
  ENDIF
;Check if weights are passed in.
  IF(NOT KEYWORD_SET(weight3d)) THEN weight3d=0
;
;  If OBJECT3D is undefined, then define it.
;  -----------------------------------------
  siz = SIZE(object3d)
  IF(siz(0) EQ 0) THEN object3d = REPLICATE({tres_D, badpixval: 0., $
    coordinate_system: '', dim1: 0, dim2: 0, dim3: 0, faceno: 0, $
    frequnits: '', instrume: '', inuse: 0, orient: '', projection: '', $
    linkweight: 0, hidden: 0, title: '', units: '', wgtunits: '', $
    version: '', window: 0},3)
;
;  Check to see how many 3-D objects are already in existence.  If there
;  are already 3 (UIMAGE can only support up to 3), then give the user
;  the option to delete one of the existing ones.  If he agrees to do this
;  then supply a menu of the 3 objects, and then delete the one he
;  selected.  If he declines to do this, then abort the instantiation
;  process.
;  -----------------------------------------------------------------------
  index3d = 0
  WHILE(object3d(index3d).inuse EQ 1) DO BEGIN
    index3d = index3d + 1
    IF(index3d EQ 3) THEN BEGIN
      PRINT, 'There currently are three 3-D objects in memory.'
      PRINT, 'No more than three 3-D objects can be manipulated at one time.'
      PRINT, 'You may designate that one of the existing 3-D objects will be'
      PRINT, 'removed in order to make room for the new 3-D object.'
      menu = ['Are you willing to remove an existing 3-D object?', 'Yes', 'No']
      sel = one_line_menu(menu)
      IF(sel EQ 2) THEN RETURN, 'ERROR'
remove3d_menu:
      menu_title = 'Remove a 3-D object'
      exit_control = 'Return to previous menu'
      menu = [menu_title, object3d(0).title, object3d(1).title, $
        object3d(2).title, 'HELP', exit_control]
      sel = umenu(menu, title = 0)
      IF(menu(sel) EQ 'HELP') THEN BEGIN
        uimage_help, menu_title
        GOTO, remove3d_menu
      ENDIF
      IF(menu(sel) EQ exit_control) THEN RETURN, 'ERROR'
      index3d = sel - 1
      delete_object, 'OBJECT3D(' + STRING(index3d, '(i1)') + ')'
    ENDIF
  ENDWHILE
;
;  Set the various structure fields.
;  ---------------------------------
  IF(KEYWORD_SET(badpixval)) THEN object3d(index3d).badpixval = badpixval $
                             ELSE object3d(index3d).badpixval = 0.
  IF(KEYWORD_SET(coordinate_system)) $
    THEN object3d(index3d).coordinate_system = coordinate_system $
    ELSE object3d(index3d).coordinate_system = 'ECLIPTIC'
  object3d(index3d).dim1 = dim1
  object3d(index3d).dim2 = dim2
  object3d(index3d).dim3 = dim3
  siz = SIZE(faceno)
  IF(siz(0) + siz(1) EQ 0) THEN object3d(index3d).faceno = -1 $
                           ELSE object3d(index3d).faceno = faceno
  object3d(index3d).inuse = 1
  IF(KEYWORD_SET(frequnits)) THEN object3d(index3d).frequnits = frequnits $
    ELSE object3d(index3d).frequnits = 'Frequency (GHz)'
  IF(KEYWORD_SET(instrume)) THEN object3d(index3d).instrume = instrume $
                            ELSE object3d(index3d).instrume = 'UNKNOWN'
  IF(KEYWORD_SET(orient)) THEN object3d(index3d).orient = orient $
                          ELSE object3d(index3d).orient = 'R'
  IF(KEYWORD_SET(projection)) THEN object3d(index3d).projection = projection $
                              ELSE object3d(index3d).projection = 'SKY-CUBE'
  IF(KEYWORD_SET(title)) THEN object3d(index3d).title = title $
    ELSE object3d(index3d).title = '3-D object #' + STRING(index3d + 1, '(i1)')
  IF(KEYWORD_SET(units)) THEN object3d(index3d).units = units $
    ELSE object3d(index3d).units = 'Pixel value (unknown units)'
  IF(KEYWORD_SET(wgtunits)) THEN object3d(index3d).wgtunits = wgtunits $
    ELSE object3d(index3d).wgtunits = ''
  IF(KEYWORD_SET(version)) THEN object3d(index3d).version = version
  w_available, window
 object3d(index3d).window = window
IF(KEYWORD_SET(hidden)) THEN object3d(index3d).hidden = hidden $
   ELSE object3d(index3d).hidden = -1
IF(N_ELEMENTS(linkweight) EQ 0) THEN object3d(index3d).linkweight = -1 $
   ELSE object3d(index3d).linkweight = linkweight 
;
;  Set the associated DATA3D_ and FREQ3D arrays.
;  ---------------------------------------------
  CASE index3d OF
    0 : BEGIN & data3D_0 = data & freq3D_0 = frequency & wt3d_0=weight3d &END
    1 : BEGIN & data3D_1 = data & freq3D_1 = frequency & wt3d_1=weight3d &END
    2 : BEGIN & data3D_2 = data & freq3D_2 = frequency & wt3d_2=weight3d &END
  ENDCASE
  RETURN, 'OBJECT3D(' + STRING(index3d, '(i1)') + ')'
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


