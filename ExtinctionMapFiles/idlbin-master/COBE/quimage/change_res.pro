PRO change_res
;+
; NAME:	
;	CHANGE_RES
; PURPOSE:
;       (A UIMAGE-specific routine)
;       To generate a new skymap (or face) from another skymap (or face)
;       that's at a different resolution.  Only resolutions 6 through 9
;       are handled.
; CATEGORY:
;	Resolution change.
; CALLING SEQUENCE:
;       CHANGE_RES
; INPUTS:
;       None.
; OUTPUTS:
;       None.
;#
; COMMON BLOCKS:
;	UIMAGE_DATA
; RESTRICTIONS:
;       None.
; PROCEDURE:
;       Call SELECT_OBJECT in order to identify the input image.
;       Then have the user indicate what the new resolution will be.
;       Either DOWNRES or CONGRID is used to generate the new array,
;       depending on whether a contraction or expansion is desired.
; SUBROUTINES CALLED:
;       downres,select_object,setup_image,uimage_help,umenu,xdisplay
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, January 1992.
; SPR 10442  Feb 04 93  Change journaling.  J Ewing
;----------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON journal,journal_on,luj
  IF(NOT defined(journal_on)) THEN journal_on = 0
;
;  Put up a menu in which the user can select the desired operand.
;  ---------------------------------------------------------------
show_menu:
  menu_title = 'Change array resolution'
  name=select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
    /face8,/face9,/help,/exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no appropriate data objects to choose from.', /CONT
    RETURN
  ENDIF
;
;  Extract pertinent info from the appropiate structure & make sure
;  that the image is a suitable canditate for the operation.
;  ----------------------------------------------------------------
  j = EXECUTE('badpixval = ' + name + '.badpixval')
  j = EXECUTE('good = WHERE(' + name + '.data NE badpixval)')
  IF(good(0) EQ -1) THEN BEGIN
    PRINT, 'Image consists of all bad pixels.  Operation not performed.'
    GOTO, show_menu
  ENDIF
  ic = STRPOS(name, '(')
  oldres = FIX(STRMID(name, ic - 1, 1))
;
;  Have the user identify the desired resolution of the output image.
;  (The new resolution must be between 6 and 9, and must be different
;  from the original resolution).
;  ------------------------------------------------------------------
show_menu2:
  menu_title_2 = 'Select new resolution'
  CASE oldres OF
    6 : rmenu=[menu_title_2,'7','8','9','HELP','Return to previous menu']
    7 : rmenu=[menu_title_2,'6','8','9','HELP','Return to previous menu']
    8 : rmenu=[menu_title_2,'6','7','9','HELP','Return to previous menu']
    9 : rmenu=[menu_title_2,'6','7','8','HELP','Return to previous menu']
  ENDCASE
  sel = umenu(rmenu, title = 0)
  IF(sel eq 5) THEN GOTO, show_menu
  IF(sel EQ 4) THEN BEGIN
    uimage_help, menu_title_2
    GOTO, show_menu2
  ENDIF
  j = EXECUTE('oldtitle = ' + name + '.title')
  newres = FIX(rmenu(sel))
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + oldtitle
    PRINTF, luj, '  the following statement was executed:'
  ENDIF
;
;  If there will be a decrease in resolution, then call DOWNRES, else
;  if there will be an increase in resolution, then call CONGRID.
;  ------------------------------------------------------------------
  IF(newres LT oldres) THEN BEGIN
    exstr = 'result = DOWNRES(' + name + $
      '.data, oldres - newres, bad = badpixval)'
    j = EXECUTE(exstr)
    IF(journal_on) THEN BEGIN
      PRINTF, luj, '  ' + exstr
      PRINTF, luj, '  where:'
      PRINTF, luj, '    ' + name + ' --> ' + oldtitle
      PRINTF, luj, '    oldres - newres = ', STRTRIM(oldres - newres, 2)
      PRINTF, luj, '    bad = ', STRTRIM(badpixval, 2)
    ENDIF
  ENDIF ELSE BEGIN
    j = EXECUTE('sz = SIZE(' + name + '.data)')
    dim1 = sz(1)
    dim2 = sz(2)
    CASE newres-oldres OF
      1: factor = 2
      2: factor = 4
      3: factor = 8
    ENDCASE
    newdim1 = dim1 * factor
    newdim2 = dim2 * factor
    exstr = 'result = CONGRID(' + name + '.data, newdim1, newdim2)'
    j = EXECUTE(exstr)
    IF(journal_on) THEN BEGIN
      PRINTF, luj, '  ' + exstr
      PRINTF, luj, '  where:'
      PRINTF, luj, '    ' + name + ' --> ' + oldtitle
      PRINTF, luj, '    newdim1 = ', STRTRIM(newdim1, 2)
      PRINTF, luj, '    newdim2 = ', STRTRIM(newdim2, 2)
    ENDIF
  ENDELSE
;
;  Store the result in the UIMAGE data-environment.
;  ------------------------------------------------
  newtitle = 'RES' + STRING(newres, '(i1)') + ' version of "' + oldtitle + '"'
  IF(STRLEN(newtitle) GT 40) THEN $
    newtitle = append_number('Change-resolution result')
  newname = setup_image(data = result, title = newtitle, badp = 0., $
    parent = name, link3d = -1)
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, newname
  IF(journal_on) THEN BEGIN
    PRINTF, luj, '    result --> ' + newtitle
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
  ENDIF
  GOTO, show_menu
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


