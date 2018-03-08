PRO assocwts
;+
; NAME:	
;	ASSOCWTS
; PURPOSE:
;       For UIMAGE users, associate chosen weight object 
;       with a selected data object.
; CATEGORY:
;	Data access.
; CALLING SEQUENCE:
;       assocwts
; INPUTS:
;       None.
; OUTPUTS:
;       None.
;#
; COMMON BLOCKS:
;	uimage_data
; RESTRICTIONS:
;       For use with UIMAGE only.
; PROCEDURE:
;       Have the user identify an object via SELECT_OBJECT.  Then 
;       associate an appropriate weight object with it.  The weight 
;       object must be a 2 dimensional object. 3d weights are handled
;       separately.
; SUBROUTINES CALLED:
;       SELECT_OBJECT
; MODIFICATION HISTORY:
;       SPR 10829: Creation:  J Newmark 6-APR-93
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;----------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
     freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON journal,journal_on,luj
  COMMON xscreen,magnify,block_usage,scrsiz
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  esc = STRING(27b)
  cuu = esc + '[A'
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Put up a menu so the user can select an operand.
;  ------------------------------------------------
obj_menu:
  menu_title = 'Associate weights with an object'
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
           /face8,/face9,/graph,/zoomed,/proj_map,/proj2_map,/object3d,$
           /help,/exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, obj_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no objects available.', /CONT
    RETURN
  ENDIF
  wmenu_title = 'Choose Weight Object'
;
;  Put up a menu so the user can select an operand.
;  ------------------------------------------------
wgt_menu:
  wname = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
           /face8,/face9,/graph,/zoomed,/proj_map,/proj2_map,/object3d,$
           /none,/help,/exit,title=wmenu_title)
  IF(wname EQ 'EXIT') THEN RETURN
  IF(wname EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, wgt_menu
  ENDIF
  IF(wname EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no objects available.', /CONT
    RETURN
  ENDIF
have_name:
IF (wname EQ 'NONE') THEN BEGIN
   IF (STRMID(name,0,1) EQ 'O') THEN BEGIN
     str_index=STRMID(name,9,1)
     j=EXECUTE('wt3d_' + str_index + '= 0')
   ENDIF 
   j=EXECUTE(name + '.linkweight = -1')
   GOTO, done
ENDIF
first_letter=STRMID(wname,0,1)
IF (first_letter EQ 'O') THEN BEGIN
  str_index3d=STRMID(name,9,1)
  wstr=STRMID(wname,9,1)
  j=EXECUTE('wt3d_' + str_index3d + ' = data3d_' + wstr)
  delname='OBJECT3D(' + wstr + ')'
  delete_object,delname
ENDIF ELSE BEGIN
  j=EXECUTE('link= ' + wname + '.window')
  j=EXECUTE(name +'.linkweight = link')
ENDELSE
done:
;
;  If journaling is enabled, then send some info to the journal file.
;  ------------------------------------------------------------------
    IF(journal_on) THEN BEGIN
      j=EXECUTE('title=' + name + '.title')
      IF (wname NE 'NONE') THEN j=EXECUTE('wtitle=' + wname + '.title') $
         ELSE wtitle='None'
      PRINTF, luj, menu_title
      PRINTF, luj, '  operand:  ' + title
      PRINTF, luj, '  weight:   ' + wtitle
      PRINTF, luj, '----------------------------------------' + $
                   '--------------------------------------'
    ENDIF
IF (uimage_version EQ 2) THEN GOTO, obj_menu
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


