PRO uloadct, ct_number
;+
;NAME:
;      ULOADCT
;PURPOSE:
;      Change the color table used on the image display to one of 15
;      pre-defined selections.
;CATEGORY:
;      Image Analysis.
;CALLING SEQUENCE:
;      UCOLORS
;INPUTS:
;      None.
;OUTPUTS:
;      None.
;#
;MODIFICATION HISTORY:
;  Written by:  K. Galuk, COBE/DMR, STX, April 1991.
;  Reformatted, renamed:  J. Ewing, ARC, November 1991.
;  Can now accept CT_NUMBER on the line-of-invocation, J Ewing, Nov, 1992.
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 04 93  Change journaling.  J Ewing
;--------------------------------------------------------------------------
  COMMON colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
  COMMON color_values, c_badpix, c_draw, c_scalemin
  COMMON journal, journal_on, luj
  COMMON history, uimage_version
  IF(NOT defined(journal_on)) THEN journal_on = 0
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  first = 1
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    RETURN
  ENDIF
  IF(NOT defined(ct_number)) THEN given_ctnum = 0 ELSE given_ctnum = 1
  IF(NOT given_ctnum) THEN BEGIN
;
;  If a value for CT_NUMBER was not given on the line-of-invocation, then
;  put up a menu so the user can select a color-table.
;  ----------------------------------------------------------------------
    menu_title = 'Standard IDL color table'
    menu_items = [menu_title, $
                  'B-W Linear',           $  ; 0
                  'Blue/White',           $  ; 1
                  'Grn-Red-Blu-Wht',      $  ; 2
                  'Red Temperature',      $  ; 3
                  'Blue/Green/Red/Yellow',$  ; 4
                  'Std Gamma-II',         $  ; 5
                  'Prism',                $  ; 6
                  'Red-Purple',           $  ; 7
                  'Green/White Linear',   $  ; 8
                  'Green/Wht Exponential',$  ; 9
                  'Green-Pink',           $  ;10
                  'Blue-Red',             $  ;11
                  '16 Level',             $  ;12
                  'Rainbow',              $  ;13
                  'Steps',                $  ;14
                  'Stern Special',        $  ;15
                  '', $
                  'HELP',   $
                  'Return to previous menu']
show_menu:
    sel = umenu(menu_items, title = 0, init = 1)
    IF(menu_items(sel) EQ 'HELP') THEN BEGIN
      uimage_help, menu_title
      GOTO, show_menu
    ENDIF
    IF(menu_items(sel) EQ 'Return to previous menu') THEN BEGIN
      IF((journal_on EQ 1) AND (first NE 1)) THEN $
        PRINTF, luj, '----------------------------------------'+$
                     '--------------------------------------'
      RETURN
    ENDIF
    IF(journal_on) THEN BEGIN
      IF(first EQ 1) THEN PRINTF, luj, menu_title
      PRINTF, luj, '  selected color table:  ' + menu_items(sel)
    ENDIF
    ct_number = sel - 1
  ENDIF
;
;  Retain the original color values from 0 to (C_SCALEMIN-1).
;  ----------------------------------------------------------
  IF(defined(r_orig)) THEN BEGIN
    r_spec = r_orig(0 : c_scalemin - 1)
    g_spec = g_orig(0 : c_scalemin - 1)
    b_spec = b_orig(0 : c_scalemin - 1)
  ENDIF ELSE BEGIN
    r_spec = [0, 255, 255,   0,   0, 255, 255, 182, 255, 255]
    g_spec = [0, 255,   0, 255,   0,   0, 255,  11, 170, 255]
    b_spec = [0, 255,   0,   0, 255, 255,   0,   0,  78, 255]
  ENDELSE
  LOADCT, ct_number ; note - On Motif, it's neccessary to do this twice
  LOADCT, ct_number ; in order to get proper CT arrays in the common block.
  temp = [[r_orig], [g_orig], [b_orig]]
  temp = CONGRID(temp, (N_ELEMENTS(r_orig) - c_scalemin), 3)
  r_orig = [r_spec, temp(*, 0)]
  g_orig = [g_spec, temp(*, 1)]
  b_orig = [b_spec, temp(*, 2)]
  r_curr = r_orig
  g_curr = g_orig
  b_curr = b_orig
  TVLCT, r_orig, g_orig, b_orig
  TVLCT, r_orig, g_orig, b_orig
  first = 0
  IF((uimage_version EQ 2) AND (NOT given_ctnum)) THEN GOTO, show_menu
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


