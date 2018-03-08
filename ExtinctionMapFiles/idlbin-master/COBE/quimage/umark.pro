PRO umark, name
;+
; NAME:	
;	UMARK
; PURPOSE:
;       For UIMAGE users, marks positions on chosen data object.
; CATEGORY:
;	Data access.
; CALLING SEQUENCE:
;       umark
; INPUTS:
;       None.
; OUTPUTS:
;       None.
;#
; COMMON BLOCKS:
;	uimage_data
; RESTRICTIONS:
;       For use with UIMAGE only.
; PROCEDURE:
;       Have the user identify an object via SELECT_OBJECT.  Then 
;       choose options for MARKMAP and call it.
; SUBROUTINES CALLED:
;       SELECT_OBJECT, MARKMAP, CONGRID
; MODIFICATION HISTORY:
;       SPR 10908 Creation:  J Newmark 05-May-93
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;----------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
     block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  esc = STRING(27b)
  cuu = esc + '[A'
  sentname=1
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
      MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
      RETURN
  ENDIF
  IF(NOT KEYWORD_SET(name)) THEN sentname=0
   menu_title = 'Markmap'
   IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Put up a menu so the user can select an operand.
;  ------------------------------------------------
obj_menu:
   name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
           /face8,/face9,/proj_map,/proj2_map,/help,/exit,$
           title=menu_title)
   IF(name EQ 'EXIT') THEN RETURN
   IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, obj_menu
   ENDIF
   IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no objects available.', /CONT
    RETURN
   ENDIF
have_name:
  j=EXECUTE('data = ' + name + '.data')
  sz=SIZE(data)
  dim1=sz(1)
  dim2=sz(2)
  zoom=2.^zoom_index
; 
; If the display is zoomed must resample data to match zoom.
;
  IF (zoom_index NE 0) THEN data=CONGRID(data,dim1*zoom,dim2*zoom)
  j=EXECUTE('window = ' + name + '.window')
  j=EXECUTE('facenum= ' +name+'.faceno')
  first_letter=STRMID(name,0,1)
  etype=['File','Keyboard','Cursor']
  menu=[' Choose coordinate and label entry mode: ','File','Keyboard',$
        'Cursor']
  eswitch=one_line_menu(menu,init=3)
  entry=etype(eswitch-1)
  color=0
  psym=1
  charsize=1
  charthick=1
  labelpos='T'
  j=EXECUTE('llx = ' + name + '.pos_x')
  j=EXECUTE('lly = ' + name + '.pos_y')
  get_parms:
  parmenu=['Define options if not using defaults',$
           'Color to plot labels','IDL plot symbol (def=plus)',$
           'Character size(def=1)','Character thickness(def=1)',$
           'Label position relative to plot symbol(def=top)',$
           'HELP','DONE']
  inparm=UMENU(parmenu,title=0,init=inparm)
  CASE inparm OF
   1:Begin
     READ,'Enter new color (0-255): ',color
     GOTO, get_parms
    END
   2:Begin
     READ,'Enter new plot symbol (1=+,2=*,3=.,4=diam,5=tri,6=squa,7=X)',psym
     GOTO, get_parms
    END
   3:Begin
     READ,'Enter new character size: ',charsize
     GOTO, get_parms
    END
   4:Begin
     READ,'Enter new character thickness: ',charthick
     GOTO, get_parms
    END
   5:Begin
     labchoice=['T','B','L','R']
     menu=['Choose label position relative to plot symbol','Top','Bottom',$
      'Left','Right']
     possw=one_line_menu(menu,init=1)
     labelpos=labchoice(possw-1)
     GOTO, get_parms
    END
   6:Begin
     uimage_help, menu_title
     GOTO, get_parms
    END
   7: PRINT,'Done defining parameters.'
  ENDCASE
  IF (first_letter EQ 'P') THEN BEGIN
    j=EXECUTE('proj= ' +name + '.projection')
    j=EXECUTE('pcoords= ' +name+ '.coordinate_system')
    IF (pcoords EQ 'ECLIPTIC') THEN pcoords='E'
    IF (pcoords EQ 'EQUATORIAL') THEN pcoords='Q'
    IF (pcoords EQ 'GALACTIC') THEN pcoords='G'
    IF (proj EQ 'AITOFF') THEN proj='A'
    IF (proj EQ 'MOLLWEIDE') THEN proj='M'
    IF (proj EQ 'GLOBAL SINUSOIDAL') THEN proj='S'
    markmap,data,win=window,entry=entry,color=color,psym=psym,llx=llx,lly=lly,$
     charsize=charsize,charthick=charthick,labelpos=labelpos,proj=proj,$
     pcoords=pcoords,facenum=facenum
  ENDIF ELSE BEGIN
    markmap,data,win=window,entry=entry,color=color,psym=psym,llx=llx,lly=lly,$
     charsize=charsize,charthick=charthick,labelpos=labelpos,facenum=facenum
  ENDELSE
IF (sentname EQ 1) THEN RETURN
IF (uimage_version EQ 2) THEN GOTO, obj_menu
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


