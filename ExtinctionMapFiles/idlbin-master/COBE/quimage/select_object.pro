FUNCTION select_object, title = title, $
  map6 = m6, map7 = m7, map8 = m8, map9 = m9, $
  face6 = f6, face7 = f7, face8 = f8, face9 = f9, proj_map = pm, $
  proj2_map=pm2,graph = gr, zoomed = zm, skycube = skycube, object3d = o3, $
  simple_graph = simple_gr, pmap_noface = pmap_noface, all=all, $
  pmap2_noface=pmap2_noface, none = none, help = help, exit = exit
;+
; NAME:	
;	SELECT_OBJECT
; PURPOSE:
;       To allow the user to select an "object" from a menu of titles.
; CATEGORY:
;	User interface.
; CALLING SEQUENCE:
;       SELECT_ITEM,selection,title=string,/map6,/map7,/map8,/map9,$
;         /face6,/face7,/face8,/face9,/graph,skycube=string
; INPUTS:
;       TITLE    = If present, a title string for the menu.
;       MAP6     = If set, then resolution 6 maps may be included in the menu.
;       MAP7     = If set, then resolution 7 maps may be included in the menu.
;       MAP8     = If set, then resolution 8 maps may be included in the menu.
;       MAP9     = If set, then resolution 9 maps may be included in the menu.
;       FACE6    = If set, then resolution 6 faces may be included in the menu.
;       FACE7    = If set, then resolution 7 faces may be included in the menu.
;       FACE8    = If set, then resolution 8 faces may be included in the menu.
;       FACE9    = If set, then resolution 9 faces may be included in the menu.
;       GRAPH    = If set, then "graphs" may be included in the menu.
;       SKYCUBE  = If present, then only those skymaps which have a
;                  projection of 'SKY-CUBE' may be included in the menu.
;       ALL      = If set, then a list of all of the chosen objects types
;                  are returned. Note, no menu is displayed.
;       HELP     = If set, then put a "Help information" option in this menu.
;       EXIT     = If set, then put an "Exit this menu" option in this menu.
; OUTPUTS:
;       SELECT_ITEM= A string which identifies the array-of-structures and
;                    the element which define the selected item, e.g.,
;                    "GRAPH(1)".
;#
; COMMON BLOCKS:
;	UIMAGE_DATA
; RESTRICTIONS:
;       None.
; PROCEDURE:
;       A list is compiled of titles and associated names (such as "GRAPH(1)").
;       What is included in the list depends on what items are stored in the
;       UIMAGE_DATA common statement, and what keywords are set when this
;       function was called.  UMENU is called to allow the user to select
;       an item from the list.
; SUBROUTINES CALLED:
;       UMENU
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing, ARC, December 1991.
;  SPR 10294  Dec 03,92  Added PMAP_NOFACE keyword.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10829  Apr  1993  Change info_3d common block
;  SPR 11031  Jun 08 93  Add /all keyword.  J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
  list = STRARR(50)
  name = STRARR(50)
  ic = 0
  IF(NOT KEYWORD_SET(skycube)) THEN skycube = 0
  IF(KEYWORD_SET(title)) THEN BEGIN
    list(0) = title
    ic = 1
  ENDIF
  IF(KEYWORD_SET(m6)) THEN BEGIN
    sz = SIZE(map6)
    IF(sz(2) EQ 8) THEN BEGIN
      num_entries = sz(3)
      FOR i = 0, num_entries - 1 DO BEGIN
        IF((NOT skycube) OR map6(i).projection EQ 'SKY-CUBE') THEN BEGIN
          list(ic) = map6(i).title
          name(ic) = 'MAP6(' + STRTRIM(STRING(i), 2) + ')'
          ic = ic + 1
        ENDIF
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(m7)) THEN BEGIN
    sz = SIZE(map7)
    IF(sz(2) EQ 8) THEN BEGIN
      num_entries = sz(3)
      FOR i = 0, num_entries - 1 DO BEGIN
        IF((NOT skycube) OR map7(i).projection EQ 'SKY-CUBE') THEN BEGIN
          list(ic) = map7(i).title
          name(ic) = 'MAP7(' + strtrim(string(i), 2) + ')'
          ic = ic + 1
        ENDIF
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(m8)) THEN BEGIN
    sz = SIZE(map8)
    IF(sz(2) EQ 8) THEN BEGIN
      num_entries = sz(3)
      FOR i = 0, num_entries - 1 DO BEGIN
        IF((NOT skycube) OR map8(i).projection EQ 'SKY-CUBE') THEN BEGIN
          list(ic) = map8(i).title
          name(ic) = 'MAP8(' + strtrim(string(i), 2) + ')'
          ic = ic + 1
        ENDIF
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(m9)) THEN BEGIN
    sz = SIZE(map9)
    IF(sz(2) EQ 8) THEN BEGIN
      num_entries = sz(3)
      FOR i = 0, num_entries - 1 DO BEGIN
        IF((NOT skycube) OR map9(i).projection EQ 'SKY-CUBE') THEN BEGIN
          list(ic) = map9(i).title
          name(ic) = 'MAP9(' + STRTRIM(STRING(i), 2) + ')'
          ic = ic + 1
        ENDIF
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(f6)) THEN BEGIN
    sz = SIZE(face6)
    IF(sz(2) EQ 8) THEN BEGIN
      num_face = sz(3)
      FOR i = 0, num_face - 1 DO BEGIN
        list(ic) = face6(i).title
        name(ic) = 'FACE6(' + STRTRIM(STRING(i), 2) + ')'
        ic = ic + 1
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(f7)) THEN BEGIN
    sz = SIZE(face7)
    IF(sz(2) EQ 8) THEN BEGIN
      num_face = sz(3)
      FOR i = 0, num_face - 1 DO BEGIN
        list(ic) = face7(i).title
        name(ic) = 'FACE7(' + STRTRIM(STRING(i), 2) + ')'
        ic = ic + 1
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(f8)) THEN BEGIN
    sz = SIZE(face8)
    IF(sz(2) EQ 8) THEN BEGIN
      num_face = sz(3)
      FOR i = 0, num_face - 1 DO BEGIN
        list(ic) = face8(i).title
        name(ic) = 'FACE8(' + STRTRIM(STRING(i), 2) + ')'
        ic = ic + 1
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(f9)) THEN BEGIN
    sz = SIZE(face9)
    IF(sz(2) EQ 8) THEN BEGIN
      num_face = sz(3)
      FOR i = 0, num_face - 1 DO BEGIN
        list(ic) = face9(i).title
        name(ic) = 'FACE9(' + STRTRIM(STRING(i), 2) + ')'
        ic = ic + 1
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(gr)) THEN BEGIN
    sz = SIZE(graph)
    IF(sz(2) EQ 8) THEN BEGIN
      num_graph = sz(3)
      FOR i = 0, num_graph - 1 DO BEGIN
        list(ic) = graph(i).title
        name(ic) = 'GRAPH(' + STRTRIM(STRING(i), 2) + ')'
        ic = ic + 1
      ENDFOR
    ENDIF
  ENDIF ELSE BEGIN
    IF(KEYWORD_SET(simple_gr)) THEN BEGIN
      sz = SIZE(graph)
      IF(sz(2) EQ 8) THEN BEGIN
        num_graph = sz(3)
        FOR i = 0, num_graph - 1 DO BEGIN
          IF((graph(i).scatter EQ 0) AND (graph(i).multi EQ 0)) THEN BEGIN
            list(ic) = graph(i).title
            name(ic) = 'GRAPH(' + STRTRIM(STRING(i), 2) + ')'
            ic = ic + 1
          ENDIF
        ENDFOR
      ENDIF
    ENDIF
  ENDELSE
  IF(KEYWORD_SET(pm)) THEN BEGIN
    sz = SIZE(proj_map)
    IF(sz(2) EQ 8) THEN BEGIN
      num_entries = sz(3)
      FOR i = 0, num_entries - 1 DO BEGIN
        list(ic) = proj_map(i).title
        name(ic) = 'PROJ_MAP(' + STRTRIM(STRING(i), 2) + ')'
        ic = ic + 1
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(pmap_noface)) THEN BEGIN
    sz = SIZE(proj_map)
    IF(sz(2) EQ 8) THEN BEGIN
      num_entries = sz(3)
      FOR i = 0, num_entries - 1 DO BEGIN
        tempname = get_name(proj_map(i).win_orig)
        IF(STRMID(tempname, 0, 3) EQ 'MAP') THEN BEGIN
          list(ic) = proj_map(i).title
          name(ic) = 'PROJ_MAP(' + STRTRIM(STRING(i), 2) + ')'
          ic = ic + 1
        ENDIF
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(pm2)) THEN BEGIN
    sz = SIZE(proj2_map)
    IF(sz(2) EQ 8) THEN BEGIN
      num_entries = sz(3)
      FOR i = 0, num_entries - 1 DO BEGIN
        list(ic) = proj2_map(i).title
        name(ic) = 'PROJ2_MAP(' + STRTRIM(STRING(i), 2) + ')'
        ic = ic + 1
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(pmap2_noface)) THEN BEGIN
    sz = SIZE(proj2_map)
    IF(sz(2) EQ 8) THEN BEGIN
      num_entries = sz(3)
      FOR i = 0, num_entries - 1 DO BEGIN
        tempname = get_name(proj2_map(i).win_orig)
        IF(STRMID(tempname, 0, 3) EQ 'MAP') THEN BEGIN
          list(ic) = proj2_map(i).title
          name(ic) = 'PROJ2_MAP(' + STRTRIM(STRING(i), 2) + ')'
          ic = ic + 1
        ENDIF
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(zm)) THEN BEGIN
    sz = SIZE(zoomed)
    IF(sz(2) EQ 8) THEN BEGIN
      num_zoomed = sz(3)
      FOR i = 0, num_zoomed - 1 DO BEGIN
        list(ic) = zoomed(i).title
        name(ic) = 'ZOOMED(' + STRTRIM(STRING(i), 2) + ')'
        ic = ic + 1
      ENDFOR
    ENDIF
  ENDIF
  IF(KEYWORD_SET(o3)) THEN BEGIN
    sz = SIZE(object3d)
    IF(sz(2) EQ 8) THEN BEGIN
      FOR i = 0, 2 DO BEGIN
        IF(object3d(i).inuse EQ 1) THEN BEGIN
          list(ic) = object3d(i).title
          name(ic) = 'OBJECT3D(' + STRTRIM(STRING(i), 2) + ')'
          ic = ic + 1
        ENDIF
      ENDFOR
    ENDIF
  ENDIF
  IF(ic EQ 0) THEN BEGIN
    RETURN, 'NO_OBJECTS'
  ENDIF
  IF(KEYWORD_SET(title) AND (ic EQ 1)) THEN RETURN, 'NO_OBJECTS'
  IF(KEYWORD_SET(none)) THEN BEGIN
    list(ic) = 'No weight'
    name(ic) = 'NONE'
    ic = ic + 1
  ENDIF
  IF(KEYWORD_SET(help)) THEN BEGIN
    list(ic) = ''
    name(ic) = ''
    ic = ic + 1
    list(ic) = 'HELP'
    name(ic) = 'HELP'
    ic = ic + 1
  ENDIF
  IF(KEYWORD_SET(exit)) THEN BEGIN
    list(ic) = 'Return to previous menu'
    name(ic) = 'EXIT'
    ic = ic + 1
  ENDIF
  list = list(0:(ic - 1))
  name = name(0:(ic - 1))
;
; The keyword /all is used in the REM_ALL procedure which removes all
; UIMAGE objects simultaneously. The whole name list is returned.
;
  IF (KEYWORD_SET(all)) THEN RETURN,name
;
show_menu:
  IF(KEYWORD_SET(title)) THEN sel = umenu(list, title = 0) $
                         ELSE sel = umenu(list)
  IF(name(sel) EQ '') THEN GOTO, show_menu
  IF(sel LT 0) THEN RETURN, sel
  RETURN, name(sel)
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


