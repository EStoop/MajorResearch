PRO project_image
;+
;  PROJECT_IMAGE - a UIMAGE-specific routine.  This routine serves as
;  UIMAGE's driver to the REPROJ routine, allowing reprojected maps
;  to be generated.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 10442  Feb 05 93  Change journaling.  J Ewing
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;  SPR 11412  28 Oct 93  Add big_gs_ecl,gal reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
    cube_side, proj, coord, sz_proj
  IF(NOT defined(journal_on)) THEN journal_on = 0
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
show_menu:
;
;  Put up a menu in which the user can identify the image which will
;  be reprojected.
;  -----------------------------------------------------------------
  menu_title = 'Reprojection'
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,/face8,/face9,$
    /help,/exit,title=menu_title,/skycube)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no sky-cubes available.', /CONT
    RETURN
  ENDIF
  j = EXECUTE('oldtitle = ' + name + '.title')
  j = EXECUTE('coordinate_system = ' + name + '.coordinate_system')
  IF(coordinate_system NE 'ECLIPTIC') THEN BEGIN
    MESSAGE, $
      'Currently, this only works for data in ECLIPTIC coordinates.', /CONT
    RETURN
  ENDIF
show_menu2:
;
;  Put up a menu in which the user can identify the new projection type.
;  ---------------------------------------------------------------------
  menu_title2 = 'New projection'
  pmenu = [menu_title2, 'Aitoff', 'Global Sinusoidal', 'Mollweide', $
   'HELP', 'Return to previous menu']
  psel = UMENU(pmenu, title = 0)
  CASE psel OF
    1: proj = 'A'
    2: proj = 'S'
    3: proj = 'M'
    4: BEGIN
         uimage_help, menu_title2
         GOTO, show_menu2
       END
    5: RETURN
  ENDCASE
show_menu3:
;
;  Put up a menu in which the user can identify the new coordinate-system.
;  -----------------------------------------------------------------------
  menu_title3 = 'New coordinate system'
  cmenu = [menu_title3, 'Ecliptic', 'Equatorial', 'Galactic', $
    'HELP', 'Return to previous menu']
  csel = UMENU(cmenu, title = 0)
  CASE csel of
    1: coor = 'E'
    2: coor = 'Q'
    3: coor = 'G'
    4: BEGIN
         uimage_help, menu_title3
         GOTO, show_menu3
       END
    5: RETURN
  ENDCASE
  IF (coor EQ 'Q') THEN BEGIN
    psize='S'
    GOTO, skip
  ENDIF
  size_table = ['S','L']
  psize = ' '
  i_size = umenu(['Enter Projection Size', $
		  'Small (512 x 256)', $
                  'Large (1024 x 512)'],title=0,init=1) 
  psize = size_table(i_size-1)
skip:
  j = EXECUTE('data = ' + name + '.data')
  faceno = -1
  j = EXECUTE('orient = ' + name + '.orient')
  IF(orient NE 'R') THEN BEGIN
    PRINT, 'Improper orientation.  Image should be right-oriented.'
    GOTO, show_menu
  ENDIF
  j = EXECUTE('instrume = ' + name + '.instrume')
  j = EXECUTE('badpixval = ' + name + '.badpixval')
  j = EXECUTE('scale_min = ' + name + '.scale_min')
  j = EXECUTE('scale_max = ' + name + '.scale_max')
  IF(STRMID(name, 0, 4) EQ 'FACE') THEN $
    j = EXECUTE('faceno = ' + name + '.faceno')
;
;  Perform the reprojection by calling REPROJ.
;  -------------------------------------------
  reproj, data, result, proj = proj, coord = coor, psize = psize, $
    face = faceno, /noshow
;
; stuff to set border to 0 rather than min
;
coord=coor
ret_stat=getinpar(result,data,res,re_proj,psize,mask,faceno)
border=where(mask EQ 255)
result(border)=0
;
;  Store the re-projected image into the UIMAGE data-environment.
;  --------------------------------------------------------------
  CASE psel OF
    1: proj = 'AITOFF'
    2: proj = 'GLOBAL SINUSOIDAL'
    3: proj = 'MOLLWEIDE'
  ENDCASE
  CASE csel of
    1: coor = 'ECLIPTIC'
    2: coor = 'EQUATORIAL'
    3: coor = 'GALACTIC'
  ENDCASE
  title = append_number(proj + '/' + coor)
  newname = setup_image(data = result, title = title, instrume = instrume, $
   proj = proj, coor = coor, orient = orient, bad = badpixval, parent = name)
  j = EXECUTE(name + '.scale_min = scale_min')
  j = EXECUTE(name + '.scale_max = scale_max')
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, newname
;
;  If journaling is on, then print out some info to the journal file.
;  ------------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + oldtitle
    PRINTF, luj, '  the following statement was executed:'
    exstr = 'reproj,data,result,proj=proj,coord=coor,' + $
      'psize=psize,face=faceno,/noshow'
    PRINTF, luj, '  ' + exstr
    PRINTF, luj, '  where:'
    PRINTF, luj, "    data came from the operand's structure"
    PRINTF, luj, '    proj = ', proj
    PRINTF, luj, '    coor = ', coor
    PRINTF, luj, '    face = ', STRTRIM(faceno, 2)
    PRINTF, luj, '    result --> ' + title
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


