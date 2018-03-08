FUNCTION setup_image, $
  badpixval = badpixval, $
  bandno = bandno, $
  bandwidth = bandwidth, $
  coordinate_system = coordinate_system, $
  data = data, $
  faceno = faceno, $
  frequency = frequency, $
  instrume = instrume, $
  link3d = link3d, $
  linkweight = linkweight, $
  hidden = hidden, $
  orient = orient, $
  parent = parent, $
  projection = projection, $
  rescale = rescale, $
  scale_max = scale_max, $
  scale_min = scale_min, $
  title = title, $
  units = units, $
  version = version
;+
; NAME:	
;	SETUP_IMAGE
; PURPOSE:
;       To add an instance of an image-type structure to one of the
;       arrays-of-structures supported by UIMAGE.
; CALLING SEQUENCE:
;       name = setup_image(data=array,title=string,faceno=integer_number,$
;         instrume=string,projection=string,coordinate_system=string)
; INPUTS:
;       DATA     = (optional)  A 2-D array of floating-point data.
;       TITLE    = (optional)  A title string.
;       FACENO   = (optional)  A face number [0..5].
;       INSTRUME = (optional)  A string which identifies the instrument.
;                  (e.g., 'DIRBE', 'DMR', 'FIRAS')
;       PROJECTION=(optional)  A string which identifies the projection.
;                  (e.g., 'SKY-CUBE', 'AITOFF', 'GLOBAL SINUSOIDAL')
;       COORDINATE_SYSTEM=(optional) A string which identifies the coordinate
;                  system (e.g., 'ECLIPTIC', 'EQUATORIAL', 'GALACTIC')
;       FREQUENCY= (optional)  A float scalar representing frequency.
;       BANDWITH = (optional)  A float scalar representing bandwidth.
;       VERSION  = (optional)  A string which identifies the data version.
;       BADPIXVAL= (optional)  A float scalar identifying the bad pixel value.
;       BANDNO   = (optional)  An integer scalar indicating band number.
;       ORIENT   = (optional)  A string indicating orientation ('L' or 'R').
; OUTPUTS:
;       NAME     = A variable which will, upon program completion, contain
;                  a string which identifies where the object was stored.
;#
; COMMON BLOCKS:
;	UIMAGE_DATA
; RESTRICTIONS:
;       None.
; PROCEDURE:
;       W_AVAILABLE is called to obtain the next available window number.
;       An inspection of DATA (whether or not supplied) is then done in
;       order to identify which array-of-structures is to be the recipient
;       of this new "item".  If the appropriate array-of-structures does not
;       exist then it is created and given the supplied values, else the
;       supplied values are placed in a new instance of a structure which
;       is appended to the appropriate existing array-of-structures.
; SUBROUTINES CALLED:
;       W_AVAILABLE
; MODIFICATION HISTORY:
;  Creation:  John A. Ewing, ARC, December 1991.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10708  Mar 17 93  Add in .LINKWEIGHT field for certain objects.  J Ewing
;  SPR 10740  Mar 25 93  Have .LINK3D be set to -1 for no linkage.  J Ewing
;  SPR 10748  Mar 26 93  Define .LINKWEIGHT for PROJ_MAP structures.  J Ewing
;  SPR 10829  Apr  1993  Add in .HIDDEN field for certain objects. J Newmark
;  SPR 10904  May 06 93  Change default value for FACENO to -1. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;  SPR 11315  Sep 15 93  Change default scaling for DIRBE images. J Newmark
; SPR 11905 07-Sep-1994   Ingest DIRBE PDS's. J. Newmark
;------------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  w_available,window
  IF(NOT KEYWORD_SET(title)) THEN title = append_number('Object')
  sz = SIZE(data)
  IF(sz(0) + sz(1) EQ 0) THEN BEGIN
    MESSAGE, 'No data was supplied.', /CONT
    RETURN, 'ERROR'
  ENDIF
  ndim = sz(0)
  dim1 = sz(1)
  dim2 = sz(2)
;  dirbemin=[0.,0.,0.,0.,10.,15.,0.,0.,-.5,-.8]
;  dirbemax=[2.,1.5,1.,3.,20.,25.,5.,1.,.5,.8]
  valid_object = 0
;
;  See if the array dimensions indicate that the data is a reprojected map.
;  ------------------------------------------------------------------------
  IF((ndim EQ 2) AND (dim1 EQ 512) AND (dim2 EQ 256)) THEN BEGIN
;
;  Create a new instance of a small PROJ_MAP structure.
;  ----------------------------------------------
    arr_struc = 'PROJ_MAP'
    j = EXECUTE('sz = SIZE(' + arr_struc + ')')
    IF(sz(2) NE 8) THEN BEGIN
      ic = 0
      proj_map = {pm_struct, data: FLTARR(512,256),$
        badpixval: 0., bandno: 0, bandwidth: 0.,coordinate_system: '',$
        data_min: 0., data_max: 0., $
        faceno:0, frequency: 0.,instrume: '', link3d: 0, $
        linkweight: 0, hidden: 0, orient: '',$
        pos_x: 0, pos_y: 0, projection: '', scale_min: 0., scale_max: 0.,$
        title: '', units: '', version: '', window: 0, win_orig: 0.}
    ENDIF ELSE BEGIN
      ic = sz(3)
      proj_map = [proj_map, {pm_struct}]
    ENDELSE
    valid_object = 1
    IF(KEYWORD_SET(parent)) THEN BEGIN
      IF(STRMID(parent,0,8) NE 'PROJ_MAP') $
        THEN j = EXECUTE('proj_map(ic).win_orig = ' + parent + '.window') $
        ELSE j = EXECUTE('proj_map(ic).win_orig = ' + parent + '.win_orig')
    ENDIF
  ENDIF
  IF((ndim EQ 2) AND (dim1 EQ 1024) AND (dim2 EQ 512)) THEN BEGIN
;
;  Create a new instance of a large PROJ_MAP structure.
;  ----------------------------------------------
    arr_struc = 'PROJ2_MAP'
    j = EXECUTE('sz = SIZE(' + arr_struc + ')')
    IF(sz(2) NE 8) THEN BEGIN
      ic = 0
      proj2_map = {pm2_struct, data: FLTARR(1024,512),$
        badpixval: 0., bandno: 0, bandwidth: 0.,coordinate_system: '',$
        data_min: 0., data_max: 0., $
        faceno:0, frequency: 0.,instrume: '', link3d: 0, $
        linkweight: 0, hidden: 0, orient: '',$
        pos_x: 0, pos_y: 0, projection: '', scale_min: 0., scale_max: 0.,$
        title: '', units: '', version: '', window: 0, win_orig: 0.}
    ENDIF ELSE BEGIN
      ic = sz(3)
      proj2_map = [proj2_map, {pm2_struct}]
    ENDELSE
    valid_object = 1
    IF(KEYWORD_SET(parent)) THEN BEGIN
      IF(STRMID(parent,0,9) NE 'PROJ2_MAP') $
        THEN j = EXECUTE('proj2_map(ic).win_orig = ' + parent + '.window') $
        ELSE j = EXECUTE('proj2_map(ic).win_orig = ' + parent + '.win_orig')
    ENDIF
  ENDIF
;
;  See if the array dimensions indicate that the data is a sky cube.
;  -----------------------------------------------------------------
  IF((ndim EQ 2) AND ((dim2/3)*4 EQ dim1)) THEN BEGIN
;
;  Create a new instance of either a MAP6, MAP7, MAP8, or MAP9 structure.
;  ----------------------------------------------------------------------
    CASE dim1 OF
       128 : res = 6
       256 : res = 7
       512 : res = 8
      1024 : res = 9
      ELSE : GOTO, bad_data
    ENDCASE
    arr_struc = 'MAP' + STRING(res, '(i1)')
    j = EXECUTE('sz = SIZE(' + arr_struc + ')')
    IF(sz(2) NE 8) THEN BEGIN
      ic = 0
      CASE res OF
          6 : map6 = {m6_struct, data: FLTARR(128,96),$
                badpixval: 0., bandno: 0, bandwidth: 0.,coordinate_system: '',$
                data_min: 0., data_max: 0., $
                faceno:0, frequency: 0.,instrume: '', link3d: 0, $
                linkweight: 0, hidden: 0, orient: '', pos_x: 0, pos_y: 0, $
                projection: '', scale_min: 0., scale_max: 0., $
                title: '', units: '', version: '', window: 0}
          7 : map7 = {m7_struct, data: FLTARR(256,192),$
                badpixval: 0., bandno: 0, bandwidth: 0., coordinate_system: '',$
                data_min: 0., data_max: 0., $
                faceno:0, frequency: 0.,instrume: '', link3d: 0, $
                linkweight: 0, hidden: 0, orient: '', pos_x: 0, pos_y: 0, $
                projection: '', scale_min: 0., scale_max: 0., $
                title: '', units: '', version: '', window: 0}
          8 : map8 = {m8_struct, data: FLTARR(512,384),$
                badpixval: 0., bandno: 0, bandwidth: 0., coordinate_system: '',$
                data_min: 0., data_max: 0., $
                faceno:0, frequency: 0.,instrume: '', link3d: 0, $
                linkweight: 0, hidden: 0, orient: '', pos_x: 0, pos_y: 0, $
                projection: '', scale_min: 0., scale_max: 0., $
                title: '', units: '', version: '', window: 0}
          9 : map9 = {m9_struct, data: FLTARR(1024,768),$
                badpixval: 0., bandno: 0, bandwidth: 0., coordinate_system: '',$
                data_min: 0., data_max: 0., $
                faceno:0, frequency: 0.,instrume: '', link3d: 0, $
                linkweight: 0, hidden: 0, orient: '', pos_x: 0, pos_y: 0, $
                projection: '', scale_min: 0., scale_max: 0., $
                title: '', units: '', version: '', window: 0}
        ELSE: BEGIN
                MESSAGE, 'Unsupported resolution.', /CONT
                RETURN, 'ERROR'
              ENDIF
      ENDCASE
    ENDIF ELSE BEGIN
      ic = sz(3)
      CASE res OF
          6 : map6 = [map6, {m6_struct}]
          7 : map7 = [map7, {m7_struct}]
          8 : map8 = [map8, {m8_struct}]
          9 : map9 = [map9, {m9_struct}]
        ELSE: BEGIN
                MESSAGE, 'Unsupported resolution.', /CONT
                RETURN, 'ERROR'
              ENDIF
      ENDCASE
    ENDELSE
    valid_object = 1
  ENDIF
;
;  See if the array dimensions indicate that the data is a sky cube face.
;  ----------------------------------------------------------------------
  IF(dim1 EQ dim2) THEN BEGIN
;
;  Create a new instance of either a FACE6, FACE7, FACE8, or FACE9 structure.
;  --------------------------------------------------------------------------
    CASE dim1 OF
        32 : res = 6
        64 : res = 7
       128 : res = 8
       256 : res = 9
      ELSE : GOTO, bad_data
    ENDCASE
    arr_struc = 'face' + STRING(res, '(i1)')
    j = EXECUTE('sz = SIZE(' + arr_struc + ')')
    IF(sz(2) NE 8) THEN BEGIN
      ic = 0
      CASE res OF
          6 : face6 = {f6_struct, data: FLTARR(32,32),$
                badpixval: 0., bandno: 0, bandwidth: 0., coordinate_system: '',$
                data_min: 0., data_max: 0., $
                faceno: 0, frequency: 0.,instrume: '', link3d: 0, $
                linkweight: 0, hidden: 0, orient: '', pos_x: 0, pos_y: 0, $
                projection: '', scale_min: 0., scale_max: 0., $
                title: '', units: '', version: '', window: 0}
          7 : face7 = {f7_struct, data: FLTARR(64,64),$
                badpixval: 0., bandno: 0, bandwidth: 0., coordinate_system: '',$
                data_min: 0., data_max: 0., $
                faceno:0, frequency: 0.,instrume: '', link3d: 0, $
                linkweight: 0, hidden: 0, orient: '',  pos_x: 0, pos_y: 0, $
                projection: '', scale_min: 0., scale_max: 0., $
                title: '', units: '', version: '', window: 0}
          8 : face8 = {f8_struct, data: FLTARR(128,128),$
                badpixval: 0., bandno: 0, bandwidth: 0., coordinate_system: '',$
                data_min: 0., data_max: 0., $
                faceno:0, frequency: 0.,instrume: '', link3d: 0, $
                linkweight: 0, hidden: 0, orient: '', pos_x: 0, pos_y: 0, $
                projection: '', scale_min: 0., scale_max: 0., $
                title: '', units: '', version: '', window: 0}
          9 : face9 = {f9_struct, data: FLTARR(256,256),$
                badpixval: 0., bandno: 0, bandwidth: 0., coordinate_system: '',$
                data_min: 0., data_max: 0., $
                faceno: 0, frequency: 0.,instrume: '', link3d: 0, $
                linkweight: 0, hidden: 0, orient: '', pos_x: 0, pos_y: 0, $
                projection: '', scale_min: 0., scale_max: 0., $
                title: '', units: '', version: '', window: 0}
        ELSE: BEGIN
                MESSAGE, 'Unsupported resolution.', /CONT
                RETURN, 'ERROR'
              ENDIF
      ENDCASE
    ENDIF ELSE BEGIN
      ic = sz(3)
      CASE res OF
          6 : face6 = [face6, {f6_struct}]
          7 : face7 = [face7, {f7_struct}]
          8 : face8 = [face8, {f8_struct}]
          9 : face9 = [face9, {f9_struct}]
        ELSE: BEGIN
                MESSAGE, 'Unsupported resolution.', /CONT
                RETURN, 'ERROR'
              ENDIF
      ENDCASE
    ENDELSE
    valid_object = 1
  ENDIF
  IF(valid_object NE 1) THEN BEGIN
bad_data:
    MESSAGE, 'An invalid data object was supplied.', /CONT
    RETURN, 'ERROR'
  ENDIF
  name = arr_struc + '(' + STRTRIM(STRING(ic), 2) + ')'
  sz = SIZE(badpixval)
  IF(sz(0) + sz(1) EQ 0) $
    THEN IF(KEYWORD_SET(parent)) THEN $
       j = EXECUTE('badpixval = ' + parent + '.badpixval') $
    ELSE badpixval = 0.
;
;  Determine the range of valid data.
;  ----------------------------------
  data_min = 0.
  data_max = 0.
  w = WHERE(data NE badpixval)
  IF(w(0) NE -1) THEN BEGIN
    data_min = MIN(data(w), MAX = data_max)
    IF(KEYWORD_SET(instrume)) THEN BEGIN
      IF(STRMID(instrume,0,5) EQ 'DIRBE') THEN $
         good=where(data ne badpixval and data gt -16375)
      IF((STRMID(instrume,0,5) EQ 'DIRBE') AND $
          NOT KEYWORD_SET(scale_max)) THEN BEGIN
         scale_max=2*median(data(good))
;        tsize=n_elements(good)
;        line=reform(data(good),tsize)
;        sorted=line(sort(line))
;        scale_min=sorted(tsize*0.1)
;        scale_max=sorted(tsize*0.9)
      ENDIF
      IF((STRMID(instrume,0,5) EQ 'DIRBE') AND $
         NOT KEYWORD_SET(scale_min)) THEN scale_min=MIN(data(good))
    ENDIF
  ENDIF ELSE BEGIN
    data_min = 9999
    data_max = -9999
  ENDELSE
;
; If image is a full skymap then fill in empty parts.
;
;
;  Set various fields within the new instance of the structure.
;  ------------------------------------------------------------
  j = EXECUTE(name + '.badpixval = badpixval')

  sz = SIZE(bandno)
  IF(sz(0) + sz(1) NE 0) $
    THEN j = EXECUTE(name + '.bandno = bandno') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name + '.bandno = ' + parent + '.bandno')

  sz = SIZE(bandwidth)
  IF(sz(0) + sz(1) NE 0) $
    THEN j = EXECUTE(name + '.bandwidth=bandwidth') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name + '.bandwidth = ' + parent + '.bandwidth')

  IF(KEYWORD_SET(coordinate_system)) $
    THEN j = EXECUTE(name + '.coordinate_system = coordinate_system') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name+'.coordinate_system = '+parent+'.coordinate_system') $
      ELSE j = EXECUTE(name + '.coordinate_system = "ECLIPTIC"')

  j = EXECUTE(name + '.data = data')
  j = EXECUTE(name + '.data_min = data_min')
  j = EXECUTE(name + '.data_max = data_max')

  sz = SIZE(faceno)
  IF(sz(0) + sz(1) NE 0) $
    THEN j = EXECUTE(name + '.faceno = faceno') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name + '.faceno = ' + parent + '.faceno') $
      ELSE j = EXECUTE(name + '.faceno = -1')

  sz = SIZE(frequency)
  IF(sz(0) + sz(1) NE 0) $
    THEN j = EXECUTE(name + '.frequency=frequency') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name + '.frequency = ' + parent + '.frequency') $
      ELSE j = EXECUTE(name + '.frequency = 1')

  IF(KEYWORD_SET(instrume)) $
    THEN j = EXECUTE(name + '.instrume = instrume') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name + '.instrume = ' + parent + '.instrume') $
      ELSE j = EXECUTE(name + '.instrume = "UNKNOWN"')

  IF(KEYWORD_SET(orient)) $
    THEN j = EXECUTE(name + '.orient = orient') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name + '.orient = ' + parent + '.orient') $
      ELSE j = EXECUTE(name + '.orient = "R"')

  IF(KEYWORD_SET(parent)) THEN j = EXECUTE(name+'.pos_x = '+parent+'.pos_x')
  IF(KEYWORD_SET(parent)) THEN j = EXECUTE(name+'.pos_y = '+parent+'.pos_y')

  IF(KEYWORD_SET(projection)) $
    THEN j = EXECUTE(name + '.projection = projection') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name + '.projection = ' + parent + '.projection') $
      ELSE j = EXECUTE(name + '.projection = "SKY-CUBE"')

  sz = SIZE(scale_min)
  IF(sz(0) + sz(1) NE 0) $
    THEN j = EXECUTE(name + '.scale_min = scale_min') $
    ELSE IF(KEYWORD_SET(parent) AND (NOT KEYWORD_SET(rescale))) $
      THEN j = EXECUTE(name + '.scale_min = ' + parent + '.scale_min') $
      ELSE j = EXECUTE(name + '.scale_min = data_min')

  sz = SIZE(scale_max)
  IF(sz(0) + sz(1) NE 0) $
    THEN j = EXECUTE(name + '.scale_max = scale_max') $
    ELSE IF(KEYWORD_SET(parent) AND (NOT KEYWORD_SET(rescale))) $
      THEN j = EXECUTE(name + '.scale_max = ' + parent + '.scale_max') $
      ELSE j = EXECUTE(name + '.scale_max = data_max')

  j = EXECUTE(name + '.link3d = -1')
  sz = SIZE(link3d)
  IF(sz(0) + sz(1) NE 0) $
    THEN j = EXECUTE(name + '.link3d = link3d') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name + '.link3d = ' + parent + '.link3d')

  j = EXECUTE(name + '.linkweight = -1')
  sz = SIZE(linkweight)
  IF(sz(0) + sz(1) NE 0) $
    THEN j = EXECUTE(name + '.linkweight = linkweight') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name + '.linkweight = ' + parent + '.linkweight')

  j = EXECUTE(name + '.hidden = -1')
  sz = SIZE(hidden)
  IF(sz(0) + sz(1) NE 0) $
    THEN j = EXECUTE(name + '.hidden = hidden') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name + '.hidden = ' + parent + '.hidden')

  j = EXECUTE(name + '.title = title')

  IF(KEYWORD_SET(version)) $
    THEN j = EXECUTE(name + '.version = version') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name + '.version = ' + parent + '.version')

  IF(KEYWORD_SET(units)) $
    THEN j = EXECUTE(name + '.units=units') $
    ELSE IF(KEYWORD_SET(parent)) THEN $
      j = EXECUTE(name + '.units = ' + parent + '.units') $
      ELSE j = EXECUTE(name + '.units = "Pixel value (unknown units)"')

  j = EXECUTE(name + '.window = window')
  RETURN, name
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


