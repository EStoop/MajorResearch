PRO delete_object, name
;+
; NAME:	
;	DELETE_OBJECT
; PURPOSE:
;       To free up some IDL memory by removing a structure that is
;       stored in the "uimage_data" common area.
; CATEGORY:
;	Memory management.
; CALLING SEQUENCE:
;       DELETE_OBJECT
; INPUTS:
;       NAME    = (optional) Name of object to be deleted.
; OUTPUTS:
;       None.
;#
; COMMON BLOCKS:
;	info_3d,uimage_data,xwindow
; RESTRICTIONS:
;       None.
; PROCEDURE:
;       The SELECT_OBJECT routine is called to allow the user to select
;       an object for deletion.  If there are no objects, then a message
;       is printed out an the routine will exit.  Otherwise the selected
;       object is eliminated.  This is done in the following manner:
;       associated with the specified window number.  If there are other
;       elements in the associated  array-of-structures besides just the
;       one to be deleted, then the array-of-structures is collapsed so as
;       to exclude that element, otherwise the array-of-structures is
;       replaced with a value of 0.
; SUBROUTINES CALLED:
;       SELECT_OBJECT,UIMAGE_HELP
; MODIFICATION HISTORY:
;       Creation:  John A. Ewing (ARC) & Sarada Chintala (STX), Dec 91-Jan 92.
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 10829  Apr 16 93  Add weight structure in 3d objects. J. Newmark
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;-----------------------------------------------------------------------------
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
  del_list = [' ']
  sentname = 1
  IF(NOT KEYWORD_SET(name)) THEN BEGIN
    sentname = 0
show_menu:
    menu_title = 'Remove an object'
    name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
  	/face8,/face9,/proj_map,/proj2_map,/graph,/zoomed,/object3d,$
        /help,/exit,title=menu_title)
    IF(name EQ 'EXIT') THEN BEGIN
      sz = SIZE(del_list)
      ndel = sz(1)-1
      IF(journal_on AND (ndel GT 0)) THEN BEGIN
        PRINTF, luj, menu_title
        PRINTF, luj, '  the following object(s) were removed:'
        FOR i = 1, ndel DO PRINTF, luj, '    ' + del_list(i)
        PRINTF, luj, '----------------------------------------' + $
                     '--------------------------------------'
      ENDIF
      RETURN
    ENDIF
    IF(name EQ 'HELP') THEN BEGIN
      uimage_help, menu_title
      GOTO, show_menu
    ENDIF
    IF(name EQ 'NO_OBJECTS') THEN BEGIN
      MESSAGE, 'There are no objects to remove.', /CONT
      RETURN
    ENDIF
  ENDIF
  type = STRMID(name, 0, STRPOS(name, '('))
  j = EXECUTE('window = ' + name + '.window')
  first_letter = STRMID(name, 0, 1)
  IF(first_letter EQ 'O') THEN BEGIN
    index3d = FIX(STRMID(name, 9, 1))
    trace_3d_links, index3d, /unlink
    object3d(index3d).inuse = 0
    str = STRING(index3d, '(i1)')
    j = EXECUTE('data3d_' + str + ' = 0')
    j = EXECUTE('freq3d_' + str + ' = 0')
    j = EXECUTE('wt3d_' + str + ' = 0')
    GOTO, remove_window
  ENDIF
;
;  If a graph, then make sure any combined-graphs don't use it.
;  ------------------------------------------------------------
  IF(type EQ 'GRAPH') THEN BEGIN
    linkid = INTARR(20)
    linkgn = INTARR(20)
    nlink = 0
    sz = SIZE(graph)
    num_entries = sz(3)
    FOR i = 0, num_entries - 1 DO BEGIN
      IF(graph(i).multi GT 0) THEN BEGIN
        flag = 0
        FOR j = 0, graph(i).multi - 1 DO BEGIN
          IF(graph(i).m_id(j) EQ window AND NOT flag) THEN BEGIN
             linkid(nlink) = graph(i).window
             linkgn(nlink) = i
             nlink = nlink + 1
             flag = 1
          ENDIF
        ENDFOR
      ENDIF
    ENDFOR
    IF(nlink NE 0) THEN BEGIN
      PRINT, 'The following composite graph(s) are derived from the graph that was'
      PRINT, '  selected for deletion:'
      FOR i = 0, nlink - 1 DO PRINT, bold(graph(linkgn(i)).title)
      PRINT, ' '
      PRINT, 'If you delete the selected graph, then the composite graph(s) ', $
        'in the'
      PRINT, '  above list will automatically be affected.'
      menu = ['Should the deletion still be performed?', 'Yes', 'No']
      sel = one_line_menu(menu, init = 2)
      IF(sel EQ 2) THEN RETURN
      FOR i = 0, nlink - 1 DO BEGIN
        tempname = get_name(linkid(i))
        ig = FIX(STRMID(tempname, 6, STRPOS(tempname, ')') - 6))
        IF(graph(ig).multi EQ 1) THEN BEGIN
          del_list = [del_list, graph(ig).title]
          delete_object, tempname
        ENDIF ELSE BEGIN
          w = WHERE(graph(ig).m_id NE window)
          graph(ig).m_id = graph(ig).m_id(w)
          graph(ig).m_color = graph(ig).m_color(w)
          graph(ig).multi = graph(ig).multi - 1
          IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
            WSET, graph(ig).window
            plot_graph, tempname
          ENDIF
        ENDELSE
      ENDFOR
    ENDIF
  ENDIF
;
;  See if any zoomed images exist which are derived from the selected object.
;  --------------------------------------------------------------------------
  sz = SIZE(zoomed)
  IF(sz(2) EQ 8) THEN BEGIN
    dim1 = sz(1)
    linklist = INTARR(dim1)
    looselist = INTARR(dim1)
    nlink = 0
    nloose = 0
    FOR i = 0, dim1 - 1 DO BEGIN
      win_orig=ABS(zoomed(i).win_orig)
      IF(get_name(win_orig) EQ name) THEN BEGIN
        linklist(nlink) = i
        nlink = nlink + 1
      ENDIF ELSE BEGIN
        looselist(nloose) = i
        nloose = nloose + 1
      ENDELSE
    ENDFOR
;
;  If there are linked zoomed images, then ask if the deletion should proceed.
;  ---------------------------------------------------------------------------
    IF(nlink NE 0) THEN BEGIN
      PRINT, 'The following zoomed image(s) are derived from the image that was'
      PRINT, '  selected for deletion:'
      FOR i = 0, nlink - 1 DO PRINT, bold(zoomed(linklist(i)).title)
      PRINT, ' '
      PRINT, 'If you delete the selected image, then the zoomed image(s) ', $
        'in the'
      PRINT, '  above list will automatically be deleted too.'
      menu = ['Should the deletion still be performed?', 'Yes', 'No']
      sel = one_line_menu(menu, init = 2)
      IF(sel EQ 2) THEN GOTO, show_menu
;
;  Delete the linked zoomed images.
;  --------------------------------
      IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
        FOR i = 0, nlink - 1 DO BEGIN
          zwindow = zoomed(linklist(i)).window
          WDELETE, zwindow
          w = WHERE(ABS(block_usage) EQ zwindow + 1)
          IF(w(0) NE -1) THEN block_usage(w) = 0
          del_list = [del_list, zoomed(linklist(i)).title]
        ENDFOR
      ENDIF
      IF(nloose EQ 0) THEN zoomed = 0 ELSE BEGIN
        str = 'zoomed=['
        FOR i = 0, nloose - 1 DO str=str + 'zoomed(' + $
          STRTRIM(STRING(looselist(i)), 2) + '),'
        str = STRMID(str, 0, STRLEN(str)-1) + ']'
        z = EXECUTE(str)
      ENDELSE
    ENDIF
  ENDIF
;
;  See if any projected images exist which are derived from the selected object.
;  -----------------------------------------------------------------------------
  sz = SIZE(proj_map)
  IF(sz(2) EQ 8) THEN BEGIN
    dim1 = sz(1)
    linklist = INTARR(dim1)
    looselist = INTARR(dim1)
    nlink = 0
    nloose = 0
    FOR i = 0, dim1 - 1 DO BEGIN
      IF(get_name(proj_map(i).win_orig) EQ name) THEN BEGIN
        linklist(nlink) = i
        nlink = nlink + 1
      ENDIF ELSE BEGIN
        looselist(nloose) = i
        nloose = nloose + 1
      ENDELSE
    ENDFOR
;
;  If there are linked PROJ_MAP images, then ask if the deletion should proceed.
;  -----------------------------------------------------------------------------
    IF(nlink NE 0) THEN BEGIN
      PRINT, 'The following re-projected image(s) are derived from the image that was'
      PRINT, '  selected for deletion:'
      FOR i = 0, nlink - 1 DO PRINT, bold(proj_map(linklist(i)).title)
      PRINT, ' '
      PRINT, 'If you delete the selected image, then the re-projected image(s) ',$
        'in the'
      PRINT, '  above list will automatically be deleted too.'
      menu = ['Should the deletion still be performed?', 'Yes', 'No']
      sel = one_line_menu(menu, init = 2)
      IF(sel EQ 2) THEN GOTO, show_menu
;
;  Delete the linked re-projected images.
;  --------------------------------------
      IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
        FOR i = 0, nlink - 1 DO BEGIN
          pwindow = proj_map(linklist(i)).window
          WDELETE, pwindow
          w = WHERE(ABS(block_usage) EQ pwindow + 1)
          IF(w(0) NE -1) THEN block_usage(w) = 0
          del_list = [del_list, proj_map(linklist(i)).title]
        ENDFOR
      ENDIF
      IF(nloose EQ 0) THEN proj_map = 0 ELSE BEGIN
        str = 'proj_map=['
        FOR i = 0, nloose - 1 DO str = str + 'proj_map(' + $
          STRTRIM(STRING(looselist(i)), 2) + '),'
        str = STRMID(str, 0, STRLEN(str)-1) + ']'
        z = EXECUTE(str)
      ENDELSE
    ENDIF
  ENDIF
  sz = SIZE(proj2_map)
  IF(sz(2) EQ 8) THEN BEGIN
    dim1 = sz(1)
    linklist = INTARR(dim1)
    looselist = INTARR(dim1)
    nlink = 0
    nloose = 0
    FOR i = 0, dim1 - 1 DO BEGIN
      IF(get_name(proj2_map(i).win_orig) EQ name) THEN BEGIN
        linklist(nlink) = i
        nlink = nlink + 1
      ENDIF ELSE BEGIN
        looselist(nloose) = i
        nloose = nloose + 1
      ENDELSE
    ENDFOR
;
;  If there are linked PROJ2_MAP images, then ask if the deletion should proceed.
;  -----------------------------------------------------------------------------
    IF(nlink NE 0) THEN BEGIN
      PRINT, 'The following re-projected image(s) are derived from the image that was'
      PRINT, '  selected for deletion:'
      FOR i = 0, nlink - 1 DO PRINT, bold(proj2_map(linklist(i)).title)
      PRINT, ' '
      PRINT, 'If you delete the selected image, then the re-projected image(s) ',$
        'in the'
      PRINT, '  above list will automatically be deleted too.'
      menu = ['Should the deletion still be performed?', 'Yes', 'No']
      sel = one_line_menu(menu, init = 2)
      IF(sel EQ 2) THEN GOTO, show_menu
;
;  Delete the linked re-projected images.
;  --------------------------------------
      IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
        FOR i = 0, nlink - 1 DO BEGIN
          pwindow = proj2_map(linklist(i)).window
          WDELETE, pwindow
          w = WHERE(ABS(block_usage) EQ pwindow + 1)
          IF(w(0) NE -1) THEN block_usage(w) = 0
          del_list = [del_list, proj2_map(linklist(i)).title]
        ENDFOR
      ENDIF
      IF(nloose EQ 0) THEN proj2_map = 0 ELSE BEGIN
        str = 'proj2_map=['
        FOR i = 0, nloose - 1 DO str = str + 'proj2_map(' + $
          STRTRIM(STRING(looselist(i)), 2) + '),'
        str = STRMID(str, 0, STRLEN(str)-1) + ']'
        z = EXECUTE(str)
      ENDELSE
    ENDIF
  ENDIF
;
;  See if any scatter plots exist which are derived from the selected object.
;  --------------------------------------------------------------------------
  sz = SIZE(graph)
  IF(sz(2) EQ 8) THEN BEGIN
    dim1 = sz(1)
    linklist = INTARR(dim1)
    looselist = INTARR(dim1)
    nlink = 0
    nloose = 0
    FOR i = 0, dim1 - 1 DO BEGIN
      IF(graph(i).scatter EQ 1) THEN BEGIN
        win_orig1 = graph(i).win_orig1
        win_orig2 = graph(i).win_orig2
        name_orig1 = get_name(win_orig1)
        name_orig2 = get_name(win_orig2)
        IF((name_orig1 EQ name) OR (name_orig2 EQ name)) THEN BEGIN
          linklist(nlink) = i
          nlink = nlink + 1
        ENDIF ELSE BEGIN
          looselist(nloose) = i
          nloose = nloose + 1
        ENDELSE
      ENDIF ELSE BEGIN
        looselist(nloose) = i
        nloose = nloose + 1
      ENDELSE
    ENDFOR
;
;  If there are linked scatter plots, then ask if the deletion should proceed.
;  ---------------------------------------------------------------------------
    IF(nlink NE 0) THEN BEGIN
      PRINT, 'The following scatter plots(s) are derived from the image that was'
      PRINT, '  selected for deletion:'
      FOR i = 0, nlink - 1 DO PRINT, bold(graph(linklist(i)).title)
      PRINT, ' '
      PRINT, 'If you delete the selected image, then the scatter plot(s) ', $
        'in the'
      PRINT, '  above list will automatically be deleted too.'
      menu = ['Should the deletion still be performed?', 'Yes', 'No']
      sel = one_line_menu(menu, init = 2)
      IF(sel EQ 2) THEN GOTO, show_menu
;
;  Delete the linked scatter plots.
;  --------------------------------
      IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
        FOR i = 0, nlink - 1 DO BEGIN
          gwindow = graph(linklist(i)).window
          WDELETE, gwindow
          w = WHERE(ABS(block_usage) EQ gwindow + 1)
          IF(w(0) NE -1) THEN block_usage(w) = 0
          del_list = [del_list, graph(linklist(i)).title]
        ENDFOR
      ENDIF
      IF(nloose EQ 0) THEN graph = 0 ELSE BEGIN
        str = 'graph=['
        FOR i = 0, nloose - 1 DO str = str + 'graph(' + $
          STRTRIM(STRING(looselist(i)), 2) + '),'
        str = STRMID(str, 0, STRLEN(str)-1) + ']'
        z = EXECUTE(str)
      ENDELSE
    ENDIF
  ENDIF
  j = EXECUTE('title = ' + name + '.title')
  del_list = [del_list, title]
  ic1 = STRPOS(name, '(')
  ic2 = STRPOS(name, ')')
  arr_struc = STRMID(name, 0, ic1)
  id = FIX(STRMID(name, ic1 + 1, ic2 - ic1 - 1))
  j = EXECUTE('sz = SIZE(' + arr_struc + ')')
  IF(sz(2) EQ 8) THEN BEGIN
    num_entries = sz(3)
    IF(num_entries eq 1) THEN i = EXECUTE(arr_struc + '=0') $
    ELSE BEGIN
      str = arr_struc + '=['
      FOR i = 0, num_entries-1 DO IF(i NE id) THEN $
        str = str + arr_struc + '(' + STRTRIM(STRING(i), 2) + '),'
      str = STRMID(str, 0, STRLEN(str)-1) + ']'
      z = EXECUTE(str)
    ENDELSE
remove_window:
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
  ENDIF
  IF(sentname EQ 1) THEN RETURN
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


