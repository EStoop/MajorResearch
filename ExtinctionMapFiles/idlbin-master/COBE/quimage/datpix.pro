PRO datpix,name
;+
;  DATPIX - a UIMAGE-specific routine.  This routine allows an
;  user to input points from the keyboard on a selected image, and to have the
;  coordinates and data values for the selected points be printed out
;  in a table. This routine is called from IDENTIFY_PIXELS.
;#
;  SPR 11104 Written by J. Newmark
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON journal,journal_on,luj
  COMMON history,uimage_version
  COMMON zback,zbgr,zdsrcmap,zbsub
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = 'Keyboard pixel information'
  IF(N_PARAMS(0) GE 1) THEN GOTO, have_name
show_menu:
;
;  Put up a menu in which the user can select the image he wants to
;  work with.
;  ----------------------------------------------------------------
  name = select_object(/map6,/map7,/map8,/map9,/face6,/face7,$
           /face8,/face9,/proj_map,/proj2_map,/zoomed,/help,$
           /exit,title=menu_title)
  IF(name EQ 'EXIT') THEN RETURN
  IF(name EQ 'HELP') THEN BEGIN
    uimage_help, menu_title
    GOTO, show_menu
  ENDIF
  IF(name EQ 'NO_OBJECTS') THEN BEGIN
    MESSAGE, 'There are no objects available.', /CONT
    RETURN
  ENDIF
have_name:
  first_letter = STRMID(name, 0, 1)
  j = EXECUTE('title = ' + name + '.title')
  j = EXECUTE('window = ' + name + '.window')
  IF(first_letter EQ 'F') THEN j = EXECUTE('faceno = ' + name + '.faceno') $
                          ELSE faceno = -1
  ic = STRPOS(name, '(')
;
;  Handle complexity for different types of data-objects.
;  ------------------------------------------------------
  stype = STRUPCASE(STRMID(name, 0, ic))
  CASE stype OF
    'ZOOMED' : BEGIN
       j = EXECUTE('win_orig = ' + name + '.win_orig')
       bflag=0
       IF (win_orig LT 0) THEN BEGIN
         win_orig=ABS(win_orig)
         bflag=1
       ENDIF
       j=EXECUTE('zoomflag = ' + name + '.zoomflag')
       IF (zoomflag NE 0) THEN bflag=1
       name_orig = get_name(win_orig)
       IF(name_orig EQ 'UNDEFINED') THEN BEGIN
         MESSAGE, $
         'The user has deleted the data associated with a zoom window.', /CONT
         GOTO, show_menu
       ENDIF
       IF(STRMID(name_orig, 0, 1) EQ 'F') THEN $
         j = EXECUTE('faceno = ' + name_orig + '.faceno')
       j = EXECUTE('data = ' + name_orig + '.data')
       IF (bflag EQ 1) THEN BEGIN
         j = EXECUTE('start_x = ' + name + '.start_x')
         j = EXECUTE('start_y = ' + name + '.start_y')
         j = EXECUTE('stop_x = ' + name + '.stop_x')
         j = EXECUTE('stop_y = ' + name + '.stop_y')
         IF (zoomflag EQ 1) THEN data(start_x:stop_x,start_y:stop_y)=zbgr
         IF (zoomflag EQ 2) THEN data(start_x:stop_x,start_y:stop_y)=zdsrcmap
         IF (zoomflag EQ 3) THEN data(start_x:stop_x,start_y:stop_y)=zbsub
       ENDIF
       ic = STRPOS(name_orig, '(')
       dtype = STRUPCASE(STRMID(name_orig, 0, ic))
       IF((dtype EQ 'PROJ_MAP') OR (dtype EQ 'PROJ2_MAP')) THEN BEGIN
         PRINT,'WARNING...information from zoomed images of reprojections'+$
            ' is not photometrically accurate.'
         j = EXECUTE('win_orig_orig = ' + name_orig + '.win_orig')
         name_orig_orig = get_name(win_orig_orig)
         IF(STRMID(name_orig_orig, 0, 1) EQ 'F') THEN $
           j = EXECUTE('faceno = ' + name_orig_orig + '.faceno')
         j = EXECUTE('data = ' + name_orig_orig + '.data')
         IF (bflag EQ 1) THEN BEGIN
          IF (zoomflag EQ 1) THEN data(start_x:stop_x,start_y:stop_y)=zbgr
          IF (zoomflag EQ 2) THEN data(start_x:stop_x,start_y:stop_y)=zdsrcmap
          IF (zoomflag EQ 3) THEN data(start_x:stop_x,start_y:stop_y)=zbsub
         ENDIF
         ic = STRPOS(name_orig_orig, '(')
         dtype = STRUPCASE(STRMID(name_orig_orig, 0, ic))
       ENDIF 
     END
    'PROJ_MAP': BEGIN
       j = EXECUTE('win_orig = ' + name + '.win_orig')
       name_orig = get_name(win_orig)
       IF(name_orig EQ 'UNDEFINED') THEN BEGIN
         MESSAGE, 'The user has deleted the data associated with a ' + $
           'reprojected image.', /CONT
         GOTO, show_menu
       ENDIF
       j=EXECUTE('projection = ' + name + '.projection')
       IF(STRMID(name_orig, 0, 1) EQ 'F') THEN $
        j = EXECUTE('faceno = ' + name_orig + '.faceno')
       j = EXECUTE('data = ' + name_orig + '.data')
       CASE projection OF
         'AITOFF'            : proj = 'A'
         'GLOBAL SINUSOIDAL' : proj = 'S'
         'MOLLWEIDE'         : proj = 'M'
         ELSE                :
       ENDCASE
       ic = STRPOS(name_orig, '(')
       dtype = STRUPCASE(STRMID(name_orig, 0, ic))
     END
    'PROJ2_MAP': BEGIN
       j = EXECUTE('win_orig = ' + name + '.win_orig')
       name_orig = get_name(win_orig)
       IF(name_orig EQ 'UNDEFINED') THEN BEGIN
         MESSAGE, 'The user has deleted the data associated with a ' + $
           'reprojected image.', /CONT
         GOTO, show_menu
       ENDIF
       j=EXECUTE('projection = ' + name + '.projection')
       IF(STRMID(name_orig, 0, 1) EQ 'F') THEN $
        j = EXECUTE('faceno = ' + name_orig + '.faceno')
       j = EXECUTE('data = ' + name_orig + '.data')
       CASE projection OF
         'AITOFF'            : proj = 'A'
         'GLOBAL SINUSOIDAL' : proj = 'S'
         'MOLLWEIDE'         : proj = 'M'
         ELSE                :
       ENDCASE
       ic = STRPOS(name_orig, '(')
       dtype = STRUPCASE(STRMID(name_orig, 0, ic))
     END
     ELSE : BEGIN
       j = EXECUTE('data = ' + name + '.data')
       ic = STRPOS(name, '(')
       dtype = STRUPCASE(STRMID(name, 0, ic))
     END
  ENDCASE
  sz = SIZE(data)
  dim1_orig = sz(1)
  dim2_orig = sz(2)
;
;  Store the resolution number in RES, for the call to XY2PIX.
;  -----------------------------------------------------------
  CASE STRUPCASE(dtype) OF
    'MAP6'  : res = 6
    'MAP7'  : res = 7
    'MAP8'  : res = 8
    'MAP9'  : res = 9
    'FACE6' : res = 6
    'FACE7' : res = 7
    'FACE8' : res = 8
    'FACE9' : res = 9
     ELSE   : BEGIN
                MESSAGE, 'Unsupported resolution.', /CONT
                GOTO, show_menu
              END
  ENDCASE
 show_menu1:
  menu_title1 = 'Select input format'
  menu = [menu_title1, 'Pixel', 'Decimal hr/dec', 'Lon-lat', $
    ' ', 'Return to previous menu']
  sel = umenu(menu, title = 0)
  CASE sel OF
    1: infmt = 'P'
    2: infmt = 'H'
    3: infmt = 'L'
    5: GOTO, show_menu
    ELSE: GOTO, show_menu1
  ENDCASE
pix=lonarr(1)
show_menu2:
IF (infmt NE 'P') THEN BEGIN
  menu_title2 = 'Select input coordinate system'
  menu = [menu_title2, 'Ecliptic', 'Equatorial', 'Galactic', $
    ' ', 'Return to previous menu']
  sel = umenu(menu, title = 0)
  CASE sel OF
    1: inco = 'E'
    2: inco = 'Q'
    3: inco = 'G'
    5: GOTO, show_menu
    ELSE: GOTO, show_menu2
  ENDCASE
  inval=fltarr(2)
  str=''
read_lonlat:
  READ,'Enter input coordinates (dhr,dec or lon,lat)',str
  str = STRTRIM(str, 2)
  IF(str EQ '') THEN RETURN
  ic = STRPOS(str, ' ')
  IF(ic EQ -1) THEN ic = STRPOS(str, ',')
  IF(ic EQ -1) THEN BEGIN
    PRINT, 'Invalid delimiter (should be " " or ",").  Please enter again.'
    GOTO, read_lonlat
  ENDIF
  str1=STRMID(str,0,ic)
  str2=STRMID(str,ic+1,STRLEN(str))
  IF(validnum(str1) + validnum(str1) NE 2) THEN BEGIN
    PRINT, 'Invalid numerical value(s).  Please enter again.'
    GOTO, read_lonlat
  ENDIF
  inval(0) = FLOAT(str1)
  inval(1) = FLOAT(str2)
  ocoord=STRTRIM('R',2)+STRTRIM(STRING(res),2)
  pix = coorconv(inval, infmt = infmt, outfmt = 'P', $
     inco = inco, outco = ocoord)
ENDIF ELSE BEGIN
  inco=STRTRIM('R',2)+STRTRIM(STRING(res),2)
  READ,'Enter pixel number: ',pix
  inval=pix
ENDELSE
IF(STRMID(dtype,0,4) EQ 'FACE') THEN BEGIN
  onface = FIX(pix(0) / (4. ^ (res - 1)))
  IF(onface NE faceno) THEN BEGIN
    PRINT, 'Selected point does not lie on the appropriate face.' + $
      '  Please enter again.'
    GOTO, show_menu1
  ENDIF
 ENDIF 
 show_menu4:
  menu_title4 = 'Select output format'
  menu = [menu_title4, 'Pixel', 'Decimal hr/dec', 'Lon-lat', $
    ' ', 'Return to previous menu']
  sel = umenu(menu, title = 0)
  CASE sel OF
    1: outfmt = 'P'
    2: outfmt = 'H'
    3: outfmt = 'L'
    5: GOTO, show_menu
    ELSE: GOTO, show_menu4
  ENDCASE
IF (outfmt NE 'P') THEN BEGIN
show_menu3:
  menu_title3 = 'Select output coordinate system'
  menu = [menu_title3, 'Ecliptic', 'Equatorial', 'Galactic', $
    ' ', 'Return to previous menu']
  sel = umenu(menu, title = 0)
  CASE sel OF
    1: outco = 'E'
    2: outco = 'Q'
    3: outco = 'G'
    5: GOTO, show_menu
    ELSE: GOTO, show_menu3
  ENDCASE
ENDIF ELSE outco=STRTRIM('R',2)+STRTRIM(STRING(res),2)
val=pix2dat(pixel=pix,raster=data)
ll=coorconv(inval,infmt=infmt,inco=inco,outfmt=outfmt,outco=outco)
PRINT,'    '
PRINT,'   pixel #',pix
PRINT,'   data value =',val
IF (outfmt NE 'P') THEN $
 PRINT,format='("   output coordinates = ",2(A6,2x))',STRTRIM(STRING(ll),2) 
IF (journal_on) THEN BEGIN
  PRINTF,luj, '----------------------------------------' + $
             '--------------------------------------'
  PRINTF,luj,'   ' + menu_title
  PRINTF,luj,'   operand: ' + title
  PRINTF,luj,'   input format=' + infmt
  IF (infmt NE 'P') THEN PRINTF,luj,'   input coordinate system =' + inco
  PRINTF,luj,'   output format=' + outfmt
  IF (outfmt NE 'P') THEN PRINTF,luj,'   output coordinate system = '+ outfmt
  PRINTF,luj,'   pixel #',pix
  PRINTF,luj,'   data value =',val
  IF (outfmt NE 'P') THEN $
  PRINTF,luj,format=$
     '("   output coordinates = ",2(A6,2x))',STRTRIM(STRING(ll),2) 
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


