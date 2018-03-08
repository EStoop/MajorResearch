PRO ubackmod, name
;+
;  UBACKMOD - a UIMAGE-specific routine.  This is UIMAGE's driver routine
;  for the MKBGRMOD procedure. This generates a background model from a
;  rasterized map.
;#
;  SPR 10890  May 04 93  Creation. J Newmark
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON journal,journal_on,luj
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON history,uimage_version
  COMMON color_values,c_badpix,c_draw,c_scalemin
  COMMON zback,zbgr,zdsrcmap,zbsub
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = '2-D Background Modeling Tool'
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Put up a menu of potential operands.
;  ------------------------------------
show_menu:
  name_orig = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
           /face8,/face9,/zoomed,/help,/exit,title=menu_title)
  IF(name_orig EQ 'EXIT') THEN RETURN
  IF(name_orig EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name_orig EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no sky-cubes available.', /CONT
    RETURN
  ENDIF
;
;  An operand has been selected - make sure it is valid.
;  -----------------------------------------------------
have_name:
  using_zoomed=0
  IF (STRMID(name_orig,0,1) EQ 'Z') THEN BEGIN
     using_zoomed=1
     str_ic=STRMID(name_orig,7,1)
     ic=FIX(str_ic)
     win_orig=zoomed(ic).win_orig
     start_x =zoomed(ic).start_x
     start_y =zoomed(ic).start_y
     stop_x = zoomed(ic).stop_x
     stop_y = zoomed(ic).stop_y
     IF (start_x LT 0) THEN start_x=0
     IF (start_y LT 0) THEN start_y=0
     IF (stop_x LT 0) THEN stop_x=0
     IF (stop_y LT 0) THEN stop_y=0
     nblocks_x =zoomed(ic).nblocks_x
     nblocks_y =zoomed(ic).nblocks_y
     specific_zoom = zoomed(ic).specific_zoom
     ztitle = zoomed(ic).title
     scale_min =zoomed(ic).scale_min
     scale_max =zoomed(ic).scale_max
     name=get_name(win_orig)
     j=EXECUTE('sz=SIZE('+name+'.data)')
     stop_x=stop_x < (sz(1)-1)
     stop_y=stop_y < (sz(2)-1)
     PRINT,'NOTE: operations for zoomed images are still under development.'
     PRINT,'      This routine will work for the zoomed image HOWEVER,'
     PRINT,'      the only one zoomed image can be worked on at a time.'
     PRINT,'      Please delete any previously created outputs from'
     PRINT,'      this routine before running it a second time.'
  ENDIF ELSE name=name_orig
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
  parent = name
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('badval = ' + name + '.badpixval')
;
; Enter in various input choices:
;
pix_parm=BYTARR(9)
PRINT,'The user has a choice of fitting using a quintic interpolation'
PRINT,'to pixels specified interactively using a mouse or more '
PRINT,'specialized fits.'
nomenu=['Choose type of fit','Interactive quintic interpolation','Specialized']
sel=one_line_menu(nomenu,init=2)
IF (sel EQ 1) THEN GOTO,backfit
get_parms:
parmenu=['Define Input Parameters ',$
            'Order of polynomial surface fit (per dimension)', $
            'Size of Spatial Filter in Pixels', $
            'Number of standard deviations used to ID a source (cut)',$
            'Size of pixel patch removed when source is ID (def=3)', $
            "Name of filter function (def='min'==lower envelope)", $
            'Include map corners, selecting points w/mouse (def=yes)', $
            'Supress plotting (not valid if interactive)', $
            'Create a final source subtracted map (must specify cut)', $
            'Create a final background subtracted map', $
            'HELP', 'DONE']
inparm=UMENU(parmenu, title=0, init=inparm)
CASE inparm OF
  1: BEGIN
   READ,'Enter order for polynomial surface fit: ',fit
   pix_parm(0)=1
   GOTO, get_parms
   END
  2: BEGIN
   READ,'Enter size of spatial filter: ',fsize
   pix_parm(1)=1
   GOTO, get_parms
   END
  3: BEGIN
   READ,'Enter number of standard deviations used to ID source: ',cut 
   pix_parm(2)=1
   GOTO, get_parms
   END
  4: BEGIN
   READ,'Enter size of pixel patch removed: ',srcwidth
   pix_parm(3)=1
   GOTO, get_parms
   END
  5: BEGIN
   filter=''
   READ,'Enter name of filter function: ',filter
   pix_parm(4)=1
   GOTO, get_parms
   END
  6: BEGIN
   PRINT,'Excluding corners.'
   pix_parm(5)=1
   GOTO, get_parms
   END
  7: BEGIN
   PRINT,' No plotting.'
   pix_parm(6)=1
   GOTO, get_parms
   END
  8: BEGIN
   PRINT,'Creating source subtracted map.'
   pix_parm(7)=1
   GOTO, get_parms
   END
  9:BEGIN
   PRINT,'Creating background subtracted map.'
   pix_parm(8)=1
   GOTO, get_parms
   END
  10: BEGIN
     uimage_help, menu_title
     GOTO, get_parms
    END
  11: PRINT,'Done defining parameters, proceeding...'
ENDCASE
IF (pix_parm(7)) THEN BEGIN
 IF NOT(pix_parm(2)) THEN BEGIN
  PRINT,'Asked for source subtracted map w/o specifying cutoff value.'
  GOTO, get_parms
 ENDIF
ENDIF
backfit:
IF ((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
  MESSAGE,'The /noplot option will be set when not on an X-window.',/CONT
  pix_parm(6)=1
ENDIF
;
; Build call to MKBGRMOD
IF (using_zoomed EQ 1) THEN $
  j=EXECUTE('data =' + name + '.data(start_x:stop_x,start_y:stop_y)') $
  ELSE j=EXECUTE('data = ' + name +'.data')
exstr= 'bgr=MKBGRMOD(data'
IF (pix_parm(7)) THEN exstr=exstr+',dsrcmap'
IF (pix_parm(0)) THEN exstr=exstr+',fit=fit'
IF (pix_parm(1)) THEN exstr=exstr+',fsize=fsize'
IF (pix_parm(2)) THEN exstr=exstr+',cut=cut'
IF (pix_parm(3)) THEN exstr=exstr+',srcwidth=srcwidth'
IF (pix_parm(4)) THEN exstr=exstr+',filter='+filter
IF (pix_parm(5)) THEN exstr=exstr+',/nocorn'
IF (pix_parm(6)) THEN exstr=exstr+',/noplot'
exstr=exstr+',badval=badval)'
;
;  Call MKBGRMOD
;  ------------
  j = EXECUTE(exstr)
  err=SIZE(bgr)
  IF (err(0) EQ 0) THEN IF (bgr EQ -1) THEN GOTO,done
;
;  Store the result in a UIMAGE data-object.
;  -----------------------------------------
IF (using_zoomed EQ 0) THEN BEGIN
  j = EXECUTE('title = ' + name + '.title')
  titleout = 'Background Model for' + title
  IF(STRLEN(titleout) GT 40) THEN $
    titleout = append_number('Background Model')
  good = WHERE(bgr NE 0.)
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(bgr(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = 0.
    scale_max = 0.
  ENDELSE
  nameout = setup_image(data = bgr, title = titleout, badpixval = badval, $
    scale_min = scale_min, hidden=1, scale_max = scale_max, parent = parent)
;Data not displayed since MKBGRMOD automatically displays it.
;
; If asked for source subtracted map display object.
 IF (pix_parm(7)) THEN BEGIN
  titlesub = 'Source subtracted map for' + title
  IF(STRLEN(titlesub) GT 40) THEN $
    titlesub = append_number('Source subtracted map')
  good = WHERE(dsrcmap NE 0.)
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(dsrcmap(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = 0.
    scale_max = 0.
  ENDELSE
  nameout = setup_image(data = dsrcmap, title = titlesub, badpixval = badval, $
    scale_min = scale_min, scale_max = scale_max, parent = parent)
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, nameout
 ENDIF
 IF (pix_parm(8)) THEN BEGIN
  titleb = 'Background subtracted map for' + title
  bsub=data-bgr
  IF(STRLEN(titleb) GT 40) THEN $
    titleb = append_number('Background subtracted map')
  good = WHERE(bgr NE 0.)
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(bsub(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = 0.
    scale_max = 0.
  ENDELSE
  nameout = setup_image(data = bsub, title = titleb, badpixval = badval, $
    scale_min = scale_min, scale_max = scale_max, parent = parent)
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay, nameout
 ENDIF
ENDIF ELSE BEGIN
;
;create a zoomed object for background model
;
  w_available, zoomwin
  zbgr=bgr
  nblocks_x = 5
  nblocks_y = 5
  sz = SIZE(zoomed)
  ic = sz(3)
  zoomed = [zoomed, {zm_struct}]
  zoomed(ic).zoomflag=1
  zoomed(ic).hidden=-1
  zoomed(ic).window = zoomwin
  zoomed(ic).win_orig = -(win_orig)
  titleout = 'Background Model for ' + ztitle
  IF(STRLEN(titleout) GT 40) THEN $
    titleout = append_number('Background Model')
  zoomed(ic).title = titleout
  zoomed(ic).nblocks_x = nblocks_x
  zoomed(ic).nblocks_y = nblocks_y
  zname = 'ZOOMED(' + STRTRIM(STRING(ic), 2) + ')'
  nblocks_x_abs = nblocks_x*(2.^zoom_index)
  nblocks_y_abs = nblocks_y*(2.^zoom_index)
  make_window, zoomwin, nblocks_x_abs, nblocks_y_abs, titleout
  zoomed(ic).start_x = start_x
  zoomed(ic).start_y = start_y
  zoomed(ic).stop_x = stop_x
  zoomed(ic).stop_y = stop_y
  zoomed(ic).specific_zoom = specific_zoom
  good = WHERE(bgr NE 0.)
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(bgr(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = 0.
    scale_max = 0.
  ENDELSE
  zoomed(ic).scale_min = scale_min
  zoomed(ic).scale_max = scale_max
  wset,zoomwin
  total_zoom=specific_zoom*2^zoom_index
  dim1=1 + stop_x - start_x
  dim2=1 + stop_y - start_y
  scimage= BYTSCL(bgr,min=scale_min,max=scale_max,top=255-c_scalemin)+$
       c_scalemin
  vimage=CONGRID(scimage,dim1*total_zoom,dim2*total_zoom)
  tv,vimage
  IF (pix_parm(7)) THEN BEGIN
;
;create a zoomed object for source subtracted model
;
  zoomwin=zoomwin+1
  zdsrcmap=dsrcmap
  nblocks_x = 5
  nblocks_y = 5
  ic = ic +1
  zoomed = [zoomed, {zm_struct}]
  zoomed(ic).zoomflag=2
  zoomed(ic).hidden=-1
  zoomed(ic).window = zoomwin
  zoomed(ic).win_orig = -(win_orig)
  titlesub = 'Source Subtracted Map for ' + ztitle
  IF(STRLEN(titlesub) GT 40) THEN $
    titlesub = append_number('Source Subtracted Map')
  zoomed(ic).title = titlesub
  zoomed(ic).nblocks_x = nblocks_x
  zoomed(ic).nblocks_y = nblocks_y
  zname = 'ZOOMED(' + STRTRIM(STRING(ic), 2) + ')'
  nblocks_x_abs = nblocks_x*(2.^zoom_index)
  nblocks_y_abs = nblocks_y*(2.^zoom_index)
  make_window, zoomwin, nblocks_x_abs, nblocks_y_abs, titlesub
  zoomed(ic).start_x = start_x
  zoomed(ic).start_y = start_y
  zoomed(ic).stop_x = stop_x
  zoomed(ic).stop_y = stop_y
  zoomed(ic).specific_zoom = specific_zoom
  good = WHERE(bgr NE 0.)
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(dsrcmap(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = 0.
    scale_max = 0.
  ENDELSE
  zoomed(ic).scale_min = scale_min
  zoomed(ic).scale_max = scale_max
  wset,zoomwin
  total_zoom=specific_zoom*2^zoom_index
  dim1=1 + stop_x - start_x
  dim2=1 + stop_y - start_y
  scimage= BYTSCL(dsrcmap,min=scale_min,max=scale_max,top=255-c_scalemin)+$
       c_scalemin
  vimage=CONGRID(scimage,dim1*total_zoom,dim2*total_zoom)
  tv,vimage
  ENDIF
  IF (pix_parm(8)) THEN BEGIN
;
;create a zoomed object for background subtracted model
;
  zoomwin=zoomwin+1
  zbsub=data-bgr
  nblocks_x = 5
  nblocks_y = 5
  ic = ic +1
  zoomed = [zoomed, {zm_struct}]
  zoomed(ic).zoomflag=3
  zoomed(ic).hidden=-1
  zoomed(ic).window = zoomwin
  zoomed(ic).win_orig = -(win_orig)
  titleb = 'Background Subtracted Map for ' + ztitle
  IF(STRLEN(titleb) GT 40) THEN $
    titleb = append_number('Background Subtracted Map')
  zoomed(ic).title = titleb
  zoomed(ic).nblocks_x = nblocks_x
  zoomed(ic).nblocks_y = nblocks_y
  zname = 'ZOOMED(' + STRTRIM(STRING(ic), 2) + ')'
  nblocks_x_abs = nblocks_x*(2.^zoom_index)
  nblocks_y_abs = nblocks_y*(2.^zoom_index)
  make_window, zoomwin, nblocks_x_abs, nblocks_y_abs, titleb
  zoomed(ic).start_x = start_x
  zoomed(ic).start_y = start_y
  zoomed(ic).stop_x = stop_x
  zoomed(ic).stop_y = stop_y
  zoomed(ic).specific_zoom = specific_zoom
  good = WHERE(bgr NE 0.)
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(zbsub(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = 0.
    scale_max = 0.
  ENDELSE
  zoomed(ic).scale_min = scale_min
  zoomed(ic).scale_max = scale_max
  wset,zoomwin
  total_zoom=specific_zoom*2^zoom_index
  dim1=1 + stop_x - start_x
  dim2=1 + stop_y - start_y
  scimage= BYTSCL(zbsub,min=scale_min,max=scale_max,top=255-c_scalemin)+$
       c_scalemin
  vimage=CONGRID(scimage,dim1*total_zoom,dim2*total_zoom)
  tv,vimage
  ENDIF
ENDELSE
;
;  If journaling is enabled, then write info to the journal file.
;  --------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title
    PRINTF, luj, '  the following statement was executed:'
    PRINTF, luj, '  ' + STRMID(exstr, 0, 76)
    exstr = STRMID(exstr, 76, STRLEN(exstr))
    WHILE(STRLEN(exstr) GE 1) DO BEGIN
      PRINTF, luj, '    ' + STRMID(exstr, 0, 74)
      exstr = STRMID(exstr, 74, STRLEN(exstr))
    ENDWHILE
    PRINTF, luj, '  Background model --> ' + titleout
    IF (pix_parm(7)) THEN PRINTF,luj,'  Source subtracted map--> ' +titlesub
    IF (pix_parm(8)) THEN PRINTF,luj,'  Background subtracted map--> ' +titleb
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
  ENDIF
done:
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


