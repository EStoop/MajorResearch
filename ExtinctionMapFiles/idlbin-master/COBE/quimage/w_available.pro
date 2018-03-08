PRO w_available, window
;+
; NAME:	
;	W_AVAILABLE
; PURPOSE:
;       To determine the smallest unused window number >= 0.
;       (Used window numbers are those referenced in appropriate elements
;       within the arrays-of-structures listed in the UIMAGE_DATA
;       commmon statement).
; CATEGORY:
;	Low-level.
; CALLING SEQUENCE:
;       W_AVAILABLE,window
; INPUTS:
;       None.
; OUTPUTS:
;       WINDOW   = A variable which will, upon program completion,
;                  contain a window number.
;#
; COMMON BLOCKS:
;	Uimage_data
; RESTRICTIONS:
;       None.
; PROCEDURE:
;       The existing arrays-of-structures are searched through and a
;       list is compiled of the used window numbers which are referenced
;       within the structures.  The smallest number greater than or
;       equal to 0 which is not within that list is then calculated and
;       returned.
; SUBROUTINES CALLED:
;       None.
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, December, 1991.
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10829  Apr  1993  Change info_3d common block
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
;
;  Compile a list (named USED) of used window numbers.
;  ---------------------------------------------------
  used = intarr(100)
  ic = 0
  type = ['map', 'face']
  FOR it = 0, 1 DO BEGIN
    FOR res = 6, 9 DO BEGIN
      arr_struc = type(it) + STRING(res, '(i1)')
      j = EXECUTE('sz = SIZE(' + arr_struc + ')')
      IF(sz(2) EQ 8) THEN BEGIN
        num_entries = sz(3)
        FOR i = 0, num_entries - 1 DO BEGIN
          arr_elem = arr_struc + '(' + STRTRIM(STRING(i), 2) + ')'
          j = EXECUTE('used(ic) = ' + arr_elem + '.window')
          ic = ic + 1
        ENDFOR
      ENDIF
    ENDFOR
  ENDFOR
  sz = SIZE(proj_map)
  IF(sz(2) EQ 8) THEN BEGIN
    num_graph = sz(3)
    FOR i = 0, num_graph - 1 DO BEGIN
      used(ic) = proj_map(i).window
      ic = ic + 1
    ENDFOR
  ENDIF
  sz = SIZE(proj2_map)
  IF(sz(2) EQ 8) THEN BEGIN
    num_graph = sz(3)
    FOR i = 0, num_graph - 1 DO BEGIN
      used(ic) = proj2_map(i).window
      ic = ic + 1
    ENDFOR
  ENDIF
  sz = SIZE(graph)
  IF(sz(2) EQ 8) THEN BEGIN
    num_graph = sz(3)
    FOR i = 0, num_graph - 1 DO BEGIN
      used(ic) = graph(i).window
      ic = ic + 1
    ENDFOR
  ENDIF
  sz = SIZE(zoomed)
  IF(sz(2) EQ 8) THEN BEGIN
    num_zoomed = sz(3)
    FOR i = 0, num_zoomed - 1 DO BEGIN
      used(ic) = zoomed(i).window
      ic = ic + 1
    ENDFOR
  ENDIF
  sz = SIZE(object3d)
  IF(sz(2) EQ 8) THEN BEGIN
    num_object3d = sz(3)
    FOR i = 0, num_object3d - 1 DO BEGIN
      used(ic) = object3d(i).window
      ic = ic + 1
    ENDFOR
  ENDIF
  IF(ic EQ 0) THEN window = 0 ELSE BEGIN
;
;  Determine the smallest number greater than or equal to 0 which
;  is not contained in the USED array.
;  --------------------------------------------------------------
    ic = ic - 1
    window = 0
loop:
    flag = 0
    FOR i = 0, ic DO IF(window EQ used(i)) THEN flag = 1
    IF(flag EQ 0) THEN RETURN ELSE BEGIN
      window = window + 1
      GOTO, loop
    ENDELSE
  ENDELSE
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


