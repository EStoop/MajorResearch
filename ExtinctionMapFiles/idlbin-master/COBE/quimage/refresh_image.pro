PRO refresh_image, name
;+
; NAME:	
;	REFRESH_IMAGE
; PURPOSE:
;       (A UIMAGE-specific routine)
;       To redisplay a selected image in its default form.
; CATEGORY:
;	Image display.
; CALLING SEQUENCE:
;       REFERESH_IMAGE
; INPUTS:
;       None.
; OUTPUTS:
;       None.
;#
; COMMON BLOCKS:
;	UIMAGE_DATA,XWINDOW
; RESTRICTIONS:
;       To be called from UIMAGE only.
; PROCEDURE:
;       Use SELECT_OBJECT so the user can select the desired image,
;       the call XDISPLAY to display that image.
; SUBROUTINES CALLED:
;       SELECT_OBJECT,SHOW_ZOOMED,UIMAGE_HELP,XDISPLAY
; MODIFICATION HISTORY:
;  Creation:  John A. Ewing, ARC, January 1992.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 04 93  Change journaling.  J Ewing
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON history,uimage_version
  COMMON journal,journal_on,luj
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  first = 1
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    RETURN
  ENDIF
;
;  Put up a menu so the user can select which image he wants to refresh.
;  ---------------------------------------------------------------------
show_menu:
  menu_title = 'Refresh an image'
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
    /face8,/face9,/zoomed,/proj_map,/proj2_map,/help,/exit,title=menu_title)
  IF(name EQ 'EXIT') THEN BEGIN
    IF((journal_on EQ 1) AND (first NE 1)) THEN $
      PRINTF, luj, '----------------------------------------'+$
                   '--------------------------------------'
    RETURN
  ENDIF
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no appropriate data objects to choose from.', /CONT
    RETURN
  ENDIF
have_name:
;
;  Call either XDISPLAY or SHOW_ZOOMED in order to effect the refresh.
;  -------------------------------------------------------------------
  IF(STRMID(name, 0, 6) NE 'ZOOMED') $
    THEN xdisplay, name, /rescale ELSE show_zoomed, name
  IF(journal_on AND (N_PARAMS(0) EQ 0)) THEN BEGIN
    IF(first EQ 1) THEN PRINTF, luj, menu_title
    j = EXECUTE('title = ' + name + '.title')
    PRINTF, luj, '  operand:  ' + title
  ENDIF
  first = 0
  IF((uimage_version EQ 2) AND (N_PARAMS(0) EQ 0)) THEN GOTO, show_menu
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


