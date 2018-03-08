PRO ucbar, name
;+
; NAME:	
;	UCBAR
; PURPOSE:
;       For UIMAGE users, put a color bar in chosen window.
; CATEGORY:
;	Data access.
; CALLING SEQUENCE:
;       ucbar
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
;       either use defaults for BAR or choose options for COLORBAR 
;       and call it. The default minimum is set to 10 since UIMAGE
;       reserves the first 10 colors for use with graphics.
; SUBROUTINES CALLED:
;       SELECT_OBJECT, BAR, COLORBAR
; MODIFICATION HISTORY:
;       SPR 10908 Creation:  J Newmark 05-May-93
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;----------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
     block_usage,zoom_index
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  esc = STRING(27b)
  cuu = esc + '[A'
  sentname=1
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
      MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
      RETURN
  ENDIF
  IF(NOT KEYWORD_SET(name)) THEN sentname=0
   menu_title = 'Add Color Bar'
   IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Put up a menu so the user can select an operand.
;  ------------------------------------------------
obj_menu:
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
have_name:
  min=10
  title=''
  j=EXECUTE('winnum = ' + name + '.window')
  menu=[' Do you want default bar (else custom bar)? ', 'Yes', 'No']
  bars=one_line_menu(menu,init=1)
  IF (bars EQ 1) THEN BEGIN
    wset,winnum
    menu=[' Do you want a title on the bar? ', 'Yes', 'No']
    tswitch=one_line_menu(menu,init=2)
    IF (tswitch EQ 1) THEN READ,'Enter title:',title
    bar,type='xb',min=min,title=title
    GOTO,done
  ENDIF
  menu=[' Do you want to add text to window? ', 'Yes', 'No']
  lswitch=one_line_menu(menu,init=2)
  menu=[' Do you want a reverse the color bar? ', 'Yes', 'No']
  rswitch=one_line_menu(menu,init=2)
  menu=[' Do you want to use default color MIN? ', 'Yes', 'No']
  minsw=one_line_menu(menu,init=1)
  menu=[' Do you want to use default color MAX? ', 'Yes', 'No']
  maxsw=one_line_menu(menu,init=1)
  exstr='colorbar,winnum'
  IF (lswitch EQ 1) THEN exstr=exstr+',/label'
  IF (rswitch EQ 1) THEN exstr=exstr+',/reverse'
  IF (minsw EQ 2) THEN BEGIN
   READ,'Enter the new minimum color value (0-254): ',min
  ENDIF
   exstr=exstr+',min=min'
  IF (maxsw EQ 2) THEN BEGIN
   READ,'Enter the new maximum color value (1-255): ',max
   exstr=exstr+',max=max'
  ENDIF
  j=EXECUTE(exstr)
done:
IF (sentname EQ 1) THEN RETURN
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


