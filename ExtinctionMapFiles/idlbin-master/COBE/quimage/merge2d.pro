PRO merge2d
;+
;  MERGE2D - a UIMAGE-specific routine. This routine combines multiple
;            2D images into 1 3D object. The user has the choice of
;            using any weights associated with data to form a set of
;            weights for the 3d object. If only some of the 2D objects
;            have weights then the others are assigned a weight of 1.
;#
;  SPR 11103  Jun 28 93 Creation: J Newmark. 
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON info_3d,object3d,data3d_0,freq3d_0,data3d_1,freq3d_1,data3d_2,$
    freq3d_2,wt3d_0,wt3d_1,wt3d_2
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
;
;  Have the user select the number of operands.
;  --------------------------------------------
  READ,'Enter the number of 2-D objects to merge into a 3-D object.',num_images
  wmenu=['Use weights associated with 2-D objects?','Yes','No']
  wchoice=one_line_menu(wmenu,init=2)
  badval_in = FLTARR(num_images)
  names = STRARR(num_images)
  links = INTARR(num_images)
  weight = STRARR(num_images)
  freq = FLTARR(num_images)
  noweight=INTARR(num_images)
;
;  Have the user select the first operand.
;  ---------------------------------------
identify_1:
  menu_title='Select 2D object to combine'
  names(0) = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
    /face8,/face9,/help,/EXIT,title=menu_title)
  IF(names(0) EQ 'HELP') THEN BEGIN
    UIMAGE_HELP, menu_title
    GOTO, identify_1
  ENDIF
  IF(names(0) EQ 'EXIT') THEN RETURN
  parent = names(0)
  faceno=-1
  first_letter=STRMID(names(0),0,1)
  IF (first_letter EQ 'F') THEN j=EXECUTE('faceno='+names(0)+'.faceno')
  j = EXECUTE('badval_in(0) = ' + names(0) + '.badpixval')
  j = EXECUTE('temptitle = ' + names(0) + '.title')
  j = EXECUTE('freq(0) = ' + names(0) + '.frequency')
  j = EXECUTE('units = ' + names(0) + '.units')
  j = EXECUTE('orient = ' + names(0) + '.orient')
  j = EXECUTE('instrume = ' + names(0) + '.instrume')
  j = EXECUTE('coordinate_system = ' + names(0) + '.coordinate_system')
  PRINT, bold('X1') + ':  ' + temptitle
  IF (wchoice EQ 1) THEN BEGIN
    j=EXECUTE('links(0)=' + names(0) + '.linkweight')
    IF (links(0) EQ -1) THEN BEGIN
       weight(0)='none'
       PRINT,'There is no weight object associated with the chosen 2-D image.'
       PRINT,'The 3-D weights for this image will be set to 1.'
       noweight(0)=1
    ENDIF ELSE BEGIN
       weight(0)=get_name(links(0))
    ENDELSE
  ENDIF
  types = STRMID(parent, 0, STRPOS(parent, '('))
;
;  Have the user select all other operands such that their data arrays
;  match the size of the data array for the first operand.
;  -------------------------------------------------------------------
  IF(num_images GT 1) THEN BEGIN
    FOR i = 2, num_images DO BEGIN
identify_n:
      str_i = STRMID(STRCOMPRESS(STRING(i)), 1, 1)
      titlen = 'Select ' + str_i +' object to combine'
      j = EXECUTE('names(i-1)=select_object(/' + types + ',/help,/EXIT,' + $
        'title=titlen)') 
      IF(names(i - 1) EQ 'HELP') THEN BEGIN
        UIMAGE_HELP, menu_title
        GOTO, identify_n
      ENDIF
      IF(names(i - 1) EQ 'EXIT') THEN RETURN
      parentn = names(i - 1)
      j = EXECUTE('badval_in(i - 1) = ' + names(i - 1) + '.badpixval')
      j = EXECUTE('temptitle = ' + parentn + '.title')
      j = EXECUTE('freq(i-1) = ' + names(i-1) + '.frequency')
      PRINT, bold('X' + str_i) + ':  ' + temptitle
      IF (wchoice EQ 1) THEN BEGIN
         j=EXECUTE('links(i-1)=' + names(i-1) + '.linkweight')
         IF (links(i-1) EQ -1) THEN BEGIN
           weight(i-1)='none'
           PRINT,'There is no weight object associated with the chosen '+$
             '2-D image.'
           PRINT,'The 3-D weights for this image will be set to 1.'
           noweight(i-1)=1
         ENDIF ELSE BEGIN
           weight(i-1)=get_name(links(i-1))
         ENDELSE
      ENDIF
    ENDFOR
  ENDIF
;
j=EXECUTE('sz = SIZE('+names(0)+'.data)')
ndata=FLTARR(sz(1),sz(2),num_images)
nweight=FLTARR(sz(1),sz(2),num_images)
FOR i=0,num_images-1 DO BEGIN
 j=EXECUTE('ndata(*,*,i) =' + names(i) + '.data')
 IF (wchoice EQ 1) THEN BEGIN
   IF (noweight(i) EQ 1) THEN nweight(*,*,i)=ndata(*,*,i)*0.+1 ELSE $
     j=EXECUTE('nweight(*,*,i) =' + weight(i) + '.data')
 ENDIF
ENDFOR
newtitle=''
get_title:
READ,'Enter a title for the 3D image',newtitle
IF (newtitle EQ '') THEN BEGIN
  PRINT,'You must enter a character string for the title.'
  GOTO,get_title
ENDIF
badpixval=MIN(badval_in,max=max)
bad=WHERE(ndata EQ badpixval,count)
IF (count NE 0) AND (badpixval NE max) THEN BEGIN
  badpixval=MIN(ndata)-1
  ndata(bad)=badpixval
ENDIF
IF (wchoice EQ 1) THEN $
name3d=setup_object3d(data=ndata,frequency=freq,weight3d=nweight,$
 badpixval=badpixval,faceno=faceno,$
 instrume=instrume,orient=orient,units=units,$
 coordinate_system=coordinate_system,title=newtitle) ELSE $
name3d=setup_object3d(data=ndata,frequency=freq,badpixval=badpixval,$
 faceno=faceno,instrume=instrume,orient=orient,units=units,$
 coordinate_system=coordinate_system,title=newtitle)

IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN xdisplay,name3d 
;
;  If journaling is enabled, then send out info to the journal file.
;  -----------------------------------------------------------------
  IF(journal_on) THEN BEGIN
    exstr = ' ndata(*,*,i) = names(i).data'
    PRINTF, luj, '  the following statement was executed:'
    PRINTF, luj, '  ' + STRMID(exstr, 0, 76)
    exstr = STRMID(exstr, 76, STRLEN(exstr))
    WHILE(STRLEN(exstr) GE 1) DO BEGIN
      PRINTF, luj, '    ' + STRMID(exstr, 0, 74)
      exstr = STRMID(exstr, 74, STRLEN(exstr))
    ENDWHILE
    PRINTF, luj, '  where:'
    FOR i = 0, num_images - 1 DO BEGIN
     PRINTF, luj, '       using the data objects below'
     j = EXECUTE('tit = ' + names(i) + '.title')
     PRINTF, luj, '    ' + names(i) + ' --> ' + tit
     PRINTF, luj, '      and using the weight objects below'
     j = EXECUTE('titl = ' + weight(i) + '.title')
     PRINTF, luj, '    weight(i) --> ' + titl
    ENDFOR
    PRINTF, luj, '    result --> ' + newtitle
    PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
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


