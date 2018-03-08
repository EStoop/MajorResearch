pro rem_all
;+
; NAME:	
;	REM_ALL
; PURPOSE:
;       To allow the user to delete all UIMAGE objects.
; CATEGORY:
;	User interface.
; CALLING SEQUENCE:
;         rem_all
; INPUTS: None.
;
; OUTPUTS: None.
;#
; COMMON BLOCKS:
;	UIMAGE_DATA, UIMAGE_DATA2,INFO_3D
; RESTRICTIONS:
;       None.
; PROCEDURE:
;       A list is compiled of titles and associated names (such as "GRAPH(1)").
;       Delete_object is then called for each member of the list.
; SUBROUTINES CALLED:
;       Delete_object, select_object, one_line_menu
; MODIFICATION HISTORY:
;       SPR 11031  Creation: Jun 08 93 J Newmark 
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON journal,journal_on,luj
 menu=['This will delete ALL objects. Continue?','Yes','No']
 verify=one_line_menu(menu, init=2)
 IF (verify EQ 2) THEN RETURN ELSE BEGIN
   name=select_object(/map6,/map7,/map8,/map9,/face6,/face7,/face8,/face9,$
          /graph,/zoomed,/proj_map,/proj2_map,/object3d,/all)
   IF (name(0) EQ 'NO_OBJECTS') THEN BEGIN
     MESSAGE, 'There are no objects to remove.',/CONT
     RETURN
   ENDIF
   sz=SIZE(name)
   IF (journal_on) THEN BEGIN
     PRINTF, luj, 'Remove ALL objects'
     PRINTF, luj, '----------------------------------------' + $
                     '--------------------------------------'
   ENDIF
   FOR i=0,sz(1)-1 DO delete_object,name(sz(1)-1-i)
 ENDELSE
 RETURN 
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


