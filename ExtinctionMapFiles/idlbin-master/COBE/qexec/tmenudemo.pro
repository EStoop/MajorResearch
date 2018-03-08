PRO tmenudemo
;
;+
; NAME
;	TMENUDEMO
; PURPOSE:
;       This provides a demonstration of the TMENU function.
; CATEGORY:
;	User interface, Menu.
; CALLING SEQUENCE:
;       IDL>  tmenudemo
; INPUTS:
;       None.
; OUTPUTS:
;       None.
; COMMON BLOCKS:
;	None.
; RESTRICTIONS:
;       Only works on Tektronix-compatible terminals.
; MODIFICATION HISTORY:
;       Creation:  John Ewing, ARC, July 1992.
;       Changed some of the wording.  John Ewing, October 23, 1992.
;-
  esc = STRING(27b)
  sel = 2
screen1:
  PRINT,esc+'[2J'
  str = 'Demonstration of the TMENU routine'
  PRINT,esc+'[1;3f'+esc+'#3'+str
  PRINT,esc+'[2;3f'+esc+'#4'+str
  PRINT,esc+'[7f'+ $
    'TMENU can put up a simple vertical menu.'
  PRINT,esc+'[11f'+'  Please note that a cntl-W will'
  PRINT,esc+'[12f'+'    refresh the entire screen.'
  PRINT,esc+'[17f'+'Use the arrow & return keys to make'
  PRINT,esc+'[18f'+'any selection in order to continue.'
  sel = tmenu(['a menu','option1','option2','option3',' ','help','exit'],$
    xpos=50,ypos=10,title=0,init=sel,/exit_on_refresh)
  IF(sel LT 0) THEN BEGIN
    sel = - (sel + 1)
    GOTO, screen1
  ENDIF
  sel = 4
screen2:
  PRINT,esc+'[2J'
  str = 'Demonstration of the TMENU routine'
  PRINT,esc+'[1;3f'+esc+'#3'+str
  PRINT,esc+'[2;3f'+esc+'#4'+str
  PRINT,esc+'[7f'+ $
    'TMENU can also produce scrollable menus as shown below.'
  PRINT,esc+'[10f'+"'...'s mean that there's more to the menu."
  PRINT,esc+'[11f'+'Please note what happens when you try to'
  PRINT,esc+'[12f'+"select a '...'"
  PRINT,esc+'[14f'+'There is a YSIZE keyword.  For this example'
  PRINT,esc+'[15f'+'I set YSIZE equal to 8.  I could have made'
  PRINT,esc+'[16f'+'the menu longer or shorter in the vertical'
  PRINT,esc+'[17f'+'dimension.'
  PRINT,esc+'[21f'+'Use the arrow & return keys to '+ $
    'make any selection in order to continue.'
  menu = ['food-item','apple','banana','carrot','donut', $
    'eclaire','fritos','grape','ham','ice cream','jello','kitkat',$
    'lasagne','mars bar','nut','orange']
  sel = tmenu(menu,xpos=50,ypos=9,title=0,init=sel,ysize=8,$
    /exit_on_refresh)
  IF(sel LT 0) THEN BEGIN
    sel = - (sel + 1)
    GOTO, screen2
  ENDIF
  sel = 4
screen3:
  PRINT,esc+'[2J'
  str = 'Demonstration of the TMENU routine'
  PRINT,esc+'[1;3f'+esc+'#3'+str
  PRINT,esc+'[2;3f'+esc+'#4'+str
  PRINT,esc+'[5f'+ $
    "Here's another scollable menu, but it's bigger than the last one."
  PRINT,esc+'[10f'+"There's no limit as to how many options"
  PRINT,esc+'[11f'+'can go into a scrollable menu.'
  PRINT,esc+'[14f'+'Please note that you can hit "Next Screen"'
  PRINT,esc+'[15f'+'or "Prev Screen" to scroll through menus'
  PRINT,esc+'[16f'+'quickly.'
  PRINT,esc+'[22f'+'Use the arrow & return keys to '+ $
    'make any selection in order to continue.'
  menu = STRARR(100)
  menu(0) = 'Index   Frequency'
  FOR i=1,99 DO menu(i) = STRING(i,'(i3)')+STRING(i*.12435e-8,'(E15.5)')
  sel = tmenu(menu,xpos=50,ypos=7,title=0,init=sel,ysize=12,$
    /exit_on_refresh)
  IF(sel LT 0) THEN BEGIN
    sel = - (sel + 1)
    GOTO, screen3
  ENDIF
  sel1 = 2
  sel2 = 0
  use = 1
screen4:
  PRINT,esc+'[2J'
  str = 'Demonstration of the TMENU routine'
  PRINT,esc+'[1;3f'+esc+'#3'+str
  PRINT,esc+'[2;3f'+esc+'#4'+str
  PRINT,esc+'[5f'+ $
    'Some additional changes to TMENU were made so it can help to put ' + $
      'forth multiple'
  PRINT,esc+'[6f'+ $
    "menus at once.
  PRINT,esc+'[19f'+ $
   "Press the TAB key to jump between menus."
  PRINT,esc+'[22f'+'Use the arrow & return keys to '+ $
    'make any selection in order to exit.'
  menu1 = ['Images','DMR 53A','DMR 53B','DIRBE BAND 1','DIRBE BAND 2',$
    'FIRAS 32','ABS("DMR 53A")','SQR("DMR 53A")']
  menu2 = ['Graphs','Histogram of "DMR 53A"','Cross section [1]',$
    'Cross section [2]','Spectum [4395]']
  dum2 = tmenu(menu2,xpos=42,ypos=8,title=0,init=sel2,$
    /exit_on_refresh,/just_reg)
  dum1 = tmenu(menu1,xpos=12,ypos=8,title=0,init=sel1,ysize=6,$
    /exit_on_refresh,/just_reg)
  sel = -1
  WHILE(sel LT 0) DO BEGIN
    IF(use EQ 0) THEN BEGIN
      sel = tmenu(menu2,xpos=42,ypos=8,title=0,init=sel2,$
        /exit_on_refresh,/tababort)
      IF(sel LT -100) THEN BEGIN
        sel2 = sel + 1000
        use = 1 - use
      ENDIF ELSE BEGIN
        IF(sel LT 0) THEN BEGIN
          sel2 = - (sel + 1)
          GOTO, screen4
        ENDIF
      ENDELSE
    ENDIF ELSE BEGIN
      sel = tmenu(menu1,xpos=12,ypos=8,title=0,init=sel1,ysize=6,$
        /exit_on_refresh,/tababort)
      IF(sel LT -100) THEN BEGIN
        sel1 = sel + 1000
        use = 1 - use
      ENDIF ELSE BEGIN
        IF(sel LT 0) THEN BEGIN
          sel1 = - (sel + 1)
          GOTO, screen4
        ENDIF
      ENDELSE
    ENDELSE
  ENDWHILE
  IF(sel LT 0) THEN BEGIN
    sel = - (sel + 1)
    GOTO, screen4
  ENDIF
  PRINT,esc+'[2J'+esc+'[0f'
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


