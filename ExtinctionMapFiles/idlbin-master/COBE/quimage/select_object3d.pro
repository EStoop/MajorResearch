PRO select_object3d, index
;+
;  SELECT_OBJECT3D - a UIMAGE-specific routine.  This routine can up a
;  menu so that a user can select one of the available 3-D objects.  A
;  value of -1 is returned if there are no 3-D objects available (in
;  which case the menu is not displayed), if there is only one 3-D object
;  around, then the number that identifies that object is returned (in
;  which case again the menu is not displayed), else a menu is displayed
;  and a value of either 0, 1, or 2 is returned to indicate which 3-D
;  object was selected.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10829  Apr  1993  Change info_3d common block
;--------------------------------------------------------------------------
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
   freq3D_2,wt3d_0,wt3d_1,wt3d_2
  menu_title = 'Select a 3-D object'
  index = -1
  sz = SIZE(object3d)
;
;  Set NUMINUSE to be the number of 3-D objects in memory (0..3).
;  Handle cases when NUMINUSE is 0, 1, or greater than 1.
;  --------------------------------------------------------------
  numinuse = 0
  IF(sz(0) NE 0) THEN FOR i = 0, 2 DO numinuse = numinuse + object3d(i).inuse
  CASE numinuse OF
     0  : BEGIN
            PRINT, 'There is no 3-D object in memory.'
            RETURN
          END
     1  : FOR i = 0, 2 DO IF(object3d(i).inuse EQ 1) THEN index = i
    ELSE: BEGIN
;
;  There are more than one 3-D objects available, so put up a menu
;  so the user can select which one he wants to work with.
;  ---------------------------------------------------------------
            exit_option = 'Return to previous menu'
            menu = [menu_title]
            indexlist = [-1]
            FOR i = 0, 2 DO IF(object3d(i).inuse EQ 1) THEN BEGIN
                              menu = [menu, object3d(i).title]
                              indexlist = [indexlist, i]
                            ENDIF
            menu = [menu, '', 'HELP', exit_option]
identify_3d:
            sel = umenu(menu, title = 0)
            IF(menu(sel) EQ 'HELP') THEN BEGIN
              uimage_help, menu_title
              GOTO, identify_3d
            ENDIF
            IF(menu(sel) EQ exit_option) THEN RETURN
            IF(menu(sel) EQ '') THEN GOTO, identify_3d
            index = indexlist(sel)
          END
  ENDCASE
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


