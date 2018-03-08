PRO gpixel,name
;+
;  GPIXEL - a UIMAGE-specific routine.  This routine displays information
;           on a selected pixel in a UIMAGE graph. 
;#
;  Written by J Newmark
;  SPR 11027  Jun 08 93 Creation.
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON journal,journal_on,luj
  COMMON xscreen,magnify,block_usage,scrsiz
  COMMON history,uimage_version
  IF(NOT defined(uimage_version)) THEN uimage_version = 2
  IF(NOT defined(journal_on)) THEN journal_on = 0
  menu_title = '1-D Pixel Information'
  IF ((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
     MESSAGE,'This routine can only be run on an X-windowing terminal.',/CONT
     RETURN
  ENDIF
  IF(N_PARAMS(0) EQ 1) THEN GOTO, have_name
;
;  Put up a menu of potential operands.
;  ------------------------------------
show_menu:
  name = select_object(/graph,/help,/exit,title=menu_title)
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
parent = name
j = EXECUTE('multi = ' + name + '.multi')
IF (multi NE 0) THEN BEGIN
  PRINT,'This routine should be run on the original graphs rather than ' + $
    'multiple graphs to avoid confusion. Please select the original graphs.'
  GOTO, show_menu
ENDIF
j = EXECUTE('title = ' + name + '.title')
j = EXECUTE('badpixval = ' + name + '.badpixval')
j = EXECUTE('linkweight= ' + name + '.linkweight')
j = EXECUTE('scatter = ' + name + '.scatter')
j = EXECUTE('num = ' +name+ '.num')
j = EXECUTE('window = ' + name + '.window')
IF (scatter EQ 1) THEN BEGIN
  j=EXECUTE('win_orig1 = ' + name + '.win_orig1')
  j=EXECUTE('win_orig2 = ' + name + '.win_orig2')
  name1=get_name(win_orig1)
  name2=get_name(win_orig2)
  j=EXECUTE('freq = ' + name1 +'.data')
  j=EXECUTE('spec = ' + name2 +'.data')
ENDIF ELSE BEGIN
 j=EXECUTE('data = ' + name + '.data')
 freq=data(0:num-1,0)
 spec=data(0:num-1,1)
ENDELSE
;
wset,window
wshow,window
plot_graph,name
PRINT,'Click left mouse button (mb1) on desired point or rightmost (mb3)' + $
      ' to exit'
select:
CURSOR,screenx,screeny,/data
CASE !ERR OF
 4: BEGIN
      IF(uimage_version EQ 2) THEN GOTO, show_menu
      RETURN
    END
 2: GOTO, select
 1: BEGIN
    IF (scatter EQ 1) THEN BEGIN
     dummy=MIN(ABS(freq-screenx),xelem)
     dummy=MIN(ABS(spec-screeny),yelem)
     PRINT,' '
     PRINT,'The index of the  X array element is: ',xelem
     PRINT,'The corresponding X-value is: ',freq(xelem)
     PRINT,'The index of the  Y array element is: ',yelem
     PRINT,'The corresponding Y-value is: ',spec(yelem)
    ENDIF ELSE BEGIN
     dummy=MIN(ABS(freq-screenx),elem)
     PRINT,' '
     PRINT,'The X-value chosen is: ',screenx
     PRINT,'The corresponding Y-value is: ',spec(elem)
     PRINT,'The index of the corresponding array element is: ',elem
    ENDELSE
;
;  If journaling is enabled, then write info to the journal file.
;  --------------------------------------------------------------
    IF(journal_on) THEN BEGIN
      PRINTF, luj, menu_title
      PRINTF, luj, '  operand:  ' + title
      IF (scatter EQ 1) THEN BEGIN
        PRINTF, luj, ' '
        PRINTF, luj, 'The index of the  X array element is: ',xelem
        PRINTF, luj, 'The corresponding X-value is: ',freq(xelem)
        PRINTF, luj, 'The index of the  Y array element is: ',yelem
        PRINTF, luj, 'The corresponding Y-value is: ',spec(yelem)
      ENDIF ELSE BEGIN
        PRINTF, luj, ' '
        PRINTF, luj, 'The X-value chosen is: ',screenx
        PRINTF, luj, 'The corresponding Y-value is: ',spec(elem)
        PRINTF, luj, 'The index of the corresponding array is: ',elem
      ENDELSE
      PRINTF, luj, '----------------------------------------' + $
                 '--------------------------------------'
    ENDIF
    END
  ENDCASE
  GOTO, select
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


