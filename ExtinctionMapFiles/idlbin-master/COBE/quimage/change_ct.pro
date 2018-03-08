PRO change_ct
;+
; NAME:	
;	CHANGE_CT
; PURPOSE:
;       To allow an X-window user to change the color table. 
; CATEGORY:
;       Color tables.
; CALLING SEQUENCE:
;       CHANGE_CT
; INPUTS:
;       none.
; OUTPUTS:
;       none.
; COMMON BLOCKS:
;	Colors
; RESTRICTIONS:
;       X-window terminals only.
; PROCEDURE:
;       A menu is put up which allows the user to select which routine
;       will be called to change the color table.
; SUBROUTINES CALLED:
;       uimage_help, umenu, uloadct, xpalette, xloadct, adjct
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, January 1992.
;  SPR 10143  Dec 11, 92  Call XPALETTE instead of PALETTE.  J Ewing
;  SPR 10959  May 18 93   Add calls to XLOADCT,ADJCT. J Newmark
;-
  COMMON colors,r_orig,g_orig,b_orig,r_curr,g_curr,b_curr
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE,'This routine can only be run on an X-window terminal.',/CONT
    RETURN
  ENDIF
  menu_title='Change color table'
  options=[menu_title,$
           'Standard IDL color table',$
           'XPALETTE','XLOADCT','ADJCT',$
           ' ',$
           'HELP',$
           'Return to previous menu']
get_argument:
  sel=umenu(options,title=0)
  CASE sel OF
    1 : uloadct
    2 : xpalette
    3 : xloadct
    4 : adjct
    6 : BEGIN
          uimage_help,menu_title
          GOTO,get_argument
        END
    7 : RETURN
    ELSE:
  ENDCASE
  GOTO,get_argument
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


