PRO hidden, name
;+
; NAME:	
;	HIDDEN
; PURPOSE:
;       For UIMAGE users, choose to display or hide a data-object.
; CATEGORY:
;	Data access.
; CALLING SEQUENCE:
;       hidden
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
;       choose to either display or hide the object. If object is
;       hidden free up the block usage.
; SUBROUTINES CALLED:
;       SELECT_OBJECT, XDISPLAY, WDELETE
; MODIFICATION HISTORY:
;       SPR 10829 Creation:  J Newmark 6-APR-93
;       SPR 10843 Apr 22 93 Corrected redisplay for graphs. J Newmark
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
  IF(NOT KEYWORD_SET(name)) THEN sentname=0
   menu_title = 'Display/hide an object'
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
menu=['Do you want data-object displayed or hidden?', 'Display', 'Hide']
value=one_line_menu(menu, init=1)
;
; Display object
;
IF(NOT defined(zoom_index)) THEN zoom_index = 1
zoom = 2.^zoom_index
IF (value EQ 1) THEN BEGIN
  j=EXECUTE(name + '.hidden = -1')
  IF (STRMID(name,0,1) EQ 'G') THEN BEGIN
      str_i=STRMID(name,6,1)
      i=FIX(str_i)
      make_window, graph(i).window, graph(i).nblocks_x, $
        graph(i).nblocks_y, graph(i).title, error_status
      IF(error_status NE -1) THEN plot_graph, name
  ENDIF ELSE IF (STRMID(name,0,1) EQ 'Z') THEN BEGIN
      str_i=STRMID(name,7,1)
      i=FIX(str_i)
      make_window,zoomed(i).window,zoomed(i).nblocks_x*zoom,$
        zoomed(i).nblocks_y*zoom,zoomed(i).title, error_status
      IF (error_status NE -1) THEN show_zoomed,name 
  ENDIF ELSE xdisplay,name 
;
; Hide object, but structure is still intact.
;
ENDIF ELSE BEGIN
  j=EXECUTE(name + '.hidden = 1')
  j=EXECUTE('window = ' + name + '.window')
;
;  If an X-window terminal is being used, then delete the window and
;  update the SCREEN_USE array to show that the section of the screen
;  that was occupied by it is available to be used for other windows now.
;  ----------------------------------------------------------------------
    IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
      WDELETE, window
      w = WHERE(ABS(block_usage) EQ window + 1)
      IF(w(0) NE -1) THEN block_usage(w) = 0
    ENDIF
ENDELSE
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


