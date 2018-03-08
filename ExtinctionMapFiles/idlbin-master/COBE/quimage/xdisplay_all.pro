PRO xdisplay_all
;+
; NAME:	
;	XDISPLAY_ALL
; PURPOSE:
;       To display all appropriate image data stored in the "uimage_data"
;       common area on an X-window screen.
; CATEGORY:
;	X-window image display.
; CALLING SEQUENCE:
;       xdisplay_all
; INPUTS
;       None.
; OUTPUTS:
;       None.
;#
; COMMON BLOCKS:
;	UIMAGE_DATA,XWINDOW
; RESTRICTIONS:
;       None.
; PROCEDURE:
;       Repeated calls to XDISPLAY are made.
; SUBROUTINES CALLED:
;       None.
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, January 1992.
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10829  Apr  1993  Change info_3d common block
;  SPR 11127  Jul 06 93  IDL for Windows compatibility. J. Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    RETURN
  ENDIF
  IF (!D.NAME EQ 'X') THEN FOR i = 0, 31 do WDELETE, i
  block_usage = 0
;
;  Enforce an overall zoom factor of 1 if any res 9 maps are present.
;  ------------------------------------------------------------------
  sz = SIZE(map9)
  IF(sz(2) EQ 8) THEN zoom_index = 0
;
;  If the overall zoom factor is undefined, then set it to 2.
;  ----------------------------------------------------------
  IF(NOT defined(zoom_index)) THEN zoom_index = 1
  zoom = 2.^zoom_index
;
;  Display all images.
;  -------------------
  type = ['MAP', 'FACE']
  FOR it = 0, 1 DO BEGIN
    FOR res = 6, 9 DO BEGIN
      arr_struc = type(it) + STRING(res, '(i1)')
      i = EXECUTE('sz = SIZE(' + arr_struc + ')')
      IF(sz(2) EQ 8) THEN BEGIN
        num_entries = sz(3)
        FOR i = 0, num_entries - 1 DO BEGIN
          name = arr_struc + '(' + STRTRIM(STRING(i), 2) + ')'
          xdisplay, name
        ENDFOR
      ENDIF
    ENDFOR
  ENDFOR
;
;  Display all projected maps.
;  ---------------------------
  sz = SIZE(proj_map)
  IF(sz(2) EQ 8) THEN BEGIN
    num_entries = sz(3)
    FOR i = 0, num_entries - 1 DO BEGIN
      name = 'PROJ_MAP(' + STRTRIM(STRING(i), 2) + ')'
      xdisplay, name
    ENDFOR
  ENDIF
  sz = SIZE(proj2_map)
  IF(sz(2) EQ 8) THEN BEGIN
    num_entries = sz(3)
    FOR i = 0, num_entries - 1 DO BEGIN
      name = 'PROJ2_MAP(' + STRTRIM(STRING(i), 2) + ')'
      xdisplay, name
    ENDFOR
  ENDIF
;
;  Display all zoomed images.
;  --------------------------
  sz = SIZE(zoomed)
  IF(sz(2) EQ 8) THEN BEGIN
    num_entries = sz(3)
    FOR i = 0, num_entries - 1 DO BEGIN
     IF (zoomed.hidden NE 1) THEN BEGIN
      make_window, zoomed(i).window, zoomed(i).nblocks_x*zoom, $
        zoomed(i).nblocks_y*zoom, zoomed(i).title, error_status
      IF(error_status NE -1) THEN BEGIN
        name = 'ZOOMED(' + STRTRIM(STRING(i), 2) + ')'
        show_zoomed, name
      ENDIF ELSE !ERR = 0
     ENDIF
    ENDFOR
  ENDIF
;
;  Display all 3-D objects.
;  ------------------------
  sz = SIZE(object3d)
  IF(sz(2) EQ 8) THEN BEGIN
    num_entries = sz(3)
    FOR i = 0, 2 DO BEGIN
      name = 'OBJECT3D(' + STRTRIM(STRING(i), 2) + ')'
      j = EXECUTE('inuse = ' + name + '.inuse')
      IF(inuse EQ 1) THEN xdisplay, name
    ENDFOR
  ENDIF
;
;  Display all graphs.
;  -------------------
  sz = SIZE(graph)
  IF(sz(2) EQ 8) THEN BEGIN
    num_entries = sz(3)
    FOR i = 0, num_entries - 1 DO BEGIN
      IF (graph(i).hidden NE 1) THEN BEGIN
       make_window, graph(i).window, graph(i).nblocks_x*zoom, $
         graph(i).nblocks_y*zoom, graph(i).title, error_status
       IF(error_status NE -1) THEN BEGIN
         name = 'GRAPH(' + STRTRIM(STRING(i), 2) + ')'
         plot_graph, name
       ENDIF ELSE !ERR = 0
      ENDIF
    ENDFOR
  ENDIF
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


