PRO uedgenhance, name
;+
; NAME:	
;	UEDGENHANCE
; PURPOSE:
;       To allow a UIMAGE user to perform edge enhancement on
;       a selected object.
; CATEGORY:
;	Egde enhancement.
; CALLING SEQUENCE:
;       UEDGENHANCE
; INPUTS:
;       None.
; OUTPUTS:
;       None.
;#
; COMMON BLOCKS:
;	UIMAGE_DATA, XWINDOW
; RESTRICTIONS:
;       To be called from UIMAGE only.
; PROCEDURE:
;       Call SELECT_OBJECT so the user can identify the image on which
;       the edge enhancement is to be performed.  Then have the user
;       select either the ROBERTS or SOBEL functions.  Then call
;       EE_SKYMAP.
; SUBROUTINES CALLED:
;       ee_skymap, select_object, setup_image, uimage_help, xdisplay
; MODIFICATION HISTORY:
;  Creation:  John A. Ewing, ARC, February 1992.
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = 'Edge enhancement'
  IF(N_PARAMS(0) GE 1) THEN GOTO, show_menu2
show_menu:
;
;  Put up a menu in which the user can identify the desired image.
;  ---------------------------------------------------------------
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
	/face8,/face9,/proj_map,/proj2_map,/help,/exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no objects available.', /CONT
    RETURN
  ENDIF
show_menu2:
;
;  Put up a menu in which the user can identify the desired edge-enhancement
;  type (ROBERTS or SOBEL).
;  -------------------------------------------------------------------------
  menu_title_2 = 'Edge En. type'
  emenu = [menu_title_2, 'ROBERTS', 'SOBEL', ' ', $
    'HELP', 'Return to previous menu']
  sel = UMENU(emenu, title = 0)
  CASE sel OF
    3: GOTO, show_menu2
    4: BEGIN
         uimage_help, menu_title_2
         GOTO, show_menu2
       END
    5: GOTO, show_menu
    ELSE:
  ENDCASE
  j = EXECUTE('badpixval = ' + name + '.badpixval')
  j = EXECUTE('data_min = ' + name + '.data_min')
  j = EXECUTE('data_max = ' + name + '.data_max')
  IF(data_min EQ data_max) THEN BEGIN
    PRINT, 'Image has no dynamic range.'
    IF(uimage_version EQ 2) THEN GOTO, show_menu ELSE RETURN
  ENDIF
  IF(data_min GT data_max) THEN BEGIN
    PRINT, 'Image consists of all bad pixels.'
    IF(uimage_version EQ 2) THEN GOTO, show_menu ELSE RETURN
  ENDIF
  type = emenu(sel)
;
;  Call ENHANMAP to generate the edge-enhanced image.
;  --------------------------------------------------
  exstr = 'result = ENHANMAP(' + name + '.data, bad = ' + $
    'badpixval, type = type)'
  j = EXECUTE(exstr)
;
;  Store the resulting edge-enhanced image into the UIMAGE data-environment.
;  -------------------------------------------------------------------------
  j = EXECUTE('title = ' + name + '.title')
  IF(sel eq 1) THEN newtitle = 'ROBERTS("' + title + '")' $
               ELSE newtitle = 'SOBEL("' + title + '")'
  IF(STRLEN(title) GT 40) THEN $
    IF(sel eq 1) THEN newtitle=append_number('ROBERTS Edge Enhancement Result') $
                 ELSE newtitle=append_number('SOBEL Edge Enhancement Result')
  badpixal = 0.
  good = WHERE(result NE badpixval)
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(result(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = badpixval
    scale_max = badpixval
  ENDELSE
  newname = setup_image(data = result, title = newtitle, badpixval = $
    badpixval, scale_min = scale_min, scale_max = scale_max, parent = name)
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, newname
;
;  If journaling is enabled, then send some info to the journal file.
;  ------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title
    PRINTF, luj, '  the following statement was executed:'
    PRINTF, luj, '  ' + exstr
    PRINTF, luj, '  where:'
    PRINTF, luj, '    ' + name + " is the name of the operand's structure"
    PRINTF, luj, '    bad  = ', STRTRIM(badpixval, 2)
    PRINTF, luj, '    type = ' + type
    PRINTF, luj, '    result --> ' + newtitle
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
  ENDIF
  IF(uimage_version EQ 2) THEN GOTO, show_menu
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


