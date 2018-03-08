FUNCTION get_name,window
;+
; NAME:	
;	get_name
; PURPOSE:
;       While UIMAGE is running, return the "name" associated with
;       the supplied window number.
; CATEGORY:
;	low-level
; CALLING SEQUENCE:
;       name = get_name(window)
; INPUTS:
;       WINDOW   = window number (integer)
; OUTPUTS:
;       A string that identifies a structure (e.g., "MAP6(0)")
; COMMON BLOCKS:
;	uimage_data
; RESTRICTIONS:
;       For use with UIMAGE only.
; PROCEDURE:
;       Check all existing structures to see if any has a window value
;       which matches the supplied value.
; SUBROUTINES CALLED:
;       None.
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, March 1992.
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;-
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  type=['MAP6','MAP7','MAP8','MAP9','FACE6','FACE7','FACE8','FACE9',$
    'GRAPH','ZOOMED','PROJ_MAP','PROJ2_MAP']
  FOR i=0,11 DO BEGIN
    z = EXECUTE('sz = SIZE('+type(i)+')')
    IF(sz(2) EQ 8) THEN BEGIN
      FOR j=0,sz(3)-1 DO BEGIN
        name = type(i)+'('+STRTRIM(STRING(j),2)+')'
        z = EXECUTE('curr_win = '+name+'.window')
        IF(curr_win EQ window) THEN RETURN,name
      ENDFOR
    ENDIF
  ENDFOR
  RETURN,'UNDEFINED'
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


