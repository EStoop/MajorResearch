PRO mkface, name
;+
;  MKFACE - a UIMAGE-specific routine.  This is UIMAGE's driver routine
;  for the GRABFACE procedure. This routine extracts a single FACE from
;  a 2D or 3D skymap. If there are weights (both 2D and 3D) associated
;  with the data then these are also extracted. 
;#
;  Written by J Newmark
;  SPR 10912  10 May 93 Creation.
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON info_3D,object3d,data3D_0,freq3D_0,data3D_1,freq3D_1,data3D_2,$
    freq3D_2,wt3d_0,wt3d_1,wt3d_2
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = 'Extract Face from Skymap'
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Put up a menu of potential operands.
;  ------------------------------------
show_menu:
  name = select_object(/map6,/map7,/map8,/map9,/object3d,/help,$
         /exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no appropriate objects available.', /CONT
    RETURN
  ENDIF
;
;  An operand has been selected - make sure it is valid.
;  -----------------------------------------------------
have_name:
  first_letter=STRMID(name,0,1)
  j = EXECUTE('faceno = ' +name + '.faceno')
  IF (faceno NE -1) THEN BEGIN
    PRINT,'Selected object is already a face'
    Return
  ENDIF
  READ,'Which face would you like extracted (0-5) ?',facenum
  parent = name
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('badpixval = ' + name + '.badpixval')
  j=EXECUTE('linkweight= ' +name+ '.linkweight')
  wsel=0
;        
; Create a new 2-D weight if have linked weights.
;
IF( linkweight NE -1) THEN BEGIN
  wsel=1
  namew=get_name(linkweight)
  wparent=namew
  j = EXECUTE('titlew = ' + namew + '.title')
  PRINT, 'Weights will be taken from ' + bold(titlew)
  j = EXECUTE('badvalw = ' + namew + '.badpixval')
  j = EXECUTE('weights = ' + namew + '.data')
  outwt=GRABFACE(weights,facenum)
  good = WHERE(outwt NE badvalw) 
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(outwt(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = 0.
    scale_max = 0.
  ENDELSE
  ntitlew='Face ' + STRMID(STRTRIM(facenum,1),0,1) + ' of ' + titlew
  wname = setup_image(data = outwt, title = ntitlew, faceno=facenum, $
    scale_min = scale_min, scale_max = scale_max, parent = wparent, $
    hidden=1)
  j=EXECUTE('linkweight=' + wname + '.window')
ENDIF
;
; Handle Case where selected object is a 3-dim one.
;
 IF (first_letter EQ 'O') THEN BEGIN
   str_index3d=STRMID(name,9,1)
   index3d=FIX(str_index3d)
   j=EXECUTE('freq = freq3d_' +str_index3d)
   j=EXECUTE('data = data3d_' + str_index3d)
   j=EXECUTE('havewt= wt3d_' + str_index3d)
   sz=SIZE(havewt)
   outwt=0
   IF (sz(0) EQ 3) THEN BEGIN
    j=EXECUTE('weight=wt3d_'+str_index3d) 
    outwt=GRABFACE(weight,facenum)
   ENDIF
;
; Extract face and set up UIMAGE object
   outface=GRABFACE(data,facenum)
   ntitle='Face ' + STRMID(STRTRIM(facenum,1),0,1) + ' of ' + title
   name=setup_object3d(data=outface,badpixval=badpixval,frequency=freq,$
     title=ntitle,units=object3d(index3d).units,$
     frequnits=object3d(index3d).frequnits,linkweight=linkweight,$
     instrume=object3d(index3d).instrume, weight3d=outwt,$
     orient=object3d(index3d).orient, faceno=facenum,$
     coordinate_system=object3d(index3d).coordinate_system)
   IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay,name
;
; Handle case where selected object is 2-dim.
;
 ENDIF ELSE BEGIN
  j=EXECUTE('data = ' +name + '.data')
; Extract face and set up UIMAGE object
  outface=GRABFACE(data,facenum)
  good = WHERE(outface NE badpixval)
  IF(good(0) NE -1) THEN BEGIN
    scale_min = MIN(outface(good), MAX = scale_max)
  ENDIF ELSE BEGIN
    scale_min = 0.
    scale_max = 0.
  ENDELSE
  ntitle='Face ' + STRMID(STRTRIM(facenum,1),0,1) + ' of ' + title
  name = setup_image(data = outface, title = ntitle, faceno=facenum, $
    scale_min = scale_min, scale_max = scale_max, linkweight=linkweight,$
    parent = parent)
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay,name
ENDELSE
;
;  If journaling is enabled, then write info to the journal file.
;  --------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    exstr='outface=GRABFACE(data,'+STRTRIM(facenum,1)+')'
    PRINTF, luj, menu_title
    PRINTF, luj, '  operand:  ' + title
    PRINTF, luj, '  the following statement was executed:'
    PRINTF, luj, '  ' + STRMID(exstr, 0, 76)
    exstr = STRMID(exstr, 76, STRLEN(exstr))
    WHILE(STRLEN(exstr) GE 1) DO BEGIN
      PRINTF, luj, '    ' + STRMID(exstr, 0, 74)
      exstr = STRMID(exstr, 74, STRLEN(exstr))
    ENDWHILE
    IF(wsel EQ 1) THEN BEGIN
      PRINTF, luj, '    weights came from:  ' + titlew
    ENDIF
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


