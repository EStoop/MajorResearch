FUNCTION qmenu,strings,title=title,init=init, retain=retain, $
	xpos=xpos, ypos=ypos, demowait=demowait
;+NAME/ONE LINE DESCRIPTION:
;    QMENU puts up a one-column scrolling menu on an X-window terminal.
;
;DESCRIPTION:
;    QMENU puts up a one-column scrolling menu on an X-window terminal.
;    This routine uses a list widget for the menu.
;    It replaces BIGWMENU which did not scroll.
;
;CALLING SEQUENCE:
;    sel_index = qmenu( strings, title=[...], init=[...] )
;       
;     
;ARGUMENTS (I=input, O=output, []=optional)
;    strings    I   1-D arr  str        A 1-D array of option strings
;                                       (and possibly a title too).
;    title     [I]  keyword  int        Index of the element in STRINGS
;                                       which will be the title.
;    init      [I]  keyword  int        Index of the element in STRINGS
;                                       which will be the default choice.
;                                       of the menu window.
;    sel_index  O   scalar   int        Index of the selected option.
;
;WARNINGS:
;    1. A value of -1 is returned in the event of any error-conditions.
;
;EXAMPLE:
;    sel_index = qmenu(['title','option 1','option 2'],title=0)
;#
;COMMON BLOCKS:  none.
;
;LIBRARY CALLS:  none.
;
;PROCEDURE:
;    If this routine is invoked from a non-X-window terminal, then put
;    out a message and exit.  Check if a valid argument was supplied
;    for STRINGS, if not then put out a message and exit.  Determine
;    what will be the title of the menu (if any) and what will be the
;    options.  
;
;    Calculate an appropriate X-size for the menu window.  This does
;    not currently work due to the List Widget's insistance on 
;    setting the width based on the items in the list, w/out concern
;    for the title's width.
;
;    Put up the menu window, then monitor and respond to mouse activity.  
;    Exit when a mouse button is pressed and the option is not blank.
;
;REVISION HISTORY:
;	Created 23-Dec-1992 Peter Kryszak:
;		parameter checking from Ewing's BIGQMENU
;		use of IDL's scrolling list widget from Turpie's WDGM
;
; Prgmr    SPR    Date       Reason
; -------- ------ ---------- -------------------------------------------
; Ewing    10619  Feb 25 93  Supply a left margin.
; Turpie   unknwn Mar 19 93  Added large nonproportional font,
;                            fixed inconsistencies in margins.
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;.TITLE
;Routine QMENU
;-
  IF ((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE,'This routine can only be run on an X-window terminal.',/CONT
    RETURN,-1
  ENDIF
;
;  Make sure an array of strings was supplied via STRINGS.
;  -------------------------------------------------------
  IF(N_PARAMS(0) NE 1) THEN BEGIN
    MESSAGE,'One parameter (an array of strings) must be supplied.',/CONT
    RETURN,-1
  ENDIF
  sz = SIZE(strings)
  IF((sz(0) NE 1) OR (sz(2) NE 7)) THEN BEGIN
    MESSAGE,'The supplied parameter is not an array of strings.',/CONT
    RETURN,-1
  ENDIF
  num_str = sz(1)
;
;  Set up title string (TITLE_STR) and array of options (OPTIONS).
;  ---------------------------------------------------------------
  sz = SIZE(title)
  IF(sz(0)+sz(1) EQ 0) THEN BEGIN
    title_str = ''
    options = strings
    xsize_title = 0
    title = -1
    num_options = num_str
  ENDIF ELSE BEGIN
    IF((title LT 0) OR (title GE num_str)) THEN BEGIN
      MESSAGE,'An invalid value for TITLE was supplied.',/CONT
      RETURN,-1
    ENDIF
    IF(num_str LE 1) THEN BEGIN
      MESSAGE,'Too few strings were supplied.',/CONT
      RETURN,-1
    ENDIF

    title_str = strings(title)
    num_options = num_str - 1

    If (Title gt 0) then begin
      options = [ strings( 0:(Title-1) ), strings( (Title+1):* ) ]
      endif $
    else begin
      options = [ strings( 1:* ) ]
      endelse

    ENDELSE

  options = STRTRIM(options, 2)
  FOR i = 0, N_ELEMENTS(options) - 1 DO options(i) = '  ' + options(i)

  xsize_title =   fix( STRLEN(title_str) )
  xsize_options = max( STRLEN(options) )
  xsize = max( [xsize_title, xsize_options] )

  extra = fix(xsize-strlen(options(0)))
  padding = '  '
  if extra gt 0 then padding = string(replicate(32B,extra)) + padding
  options(0) = options(0) + padding


; default position for the menu is near the center
Device, Get_Screen_Size = Dim
x0 = Fix( Dim(0)*0.025 )
y0 = Fix( Dim(1)*0.61 )

; but check for explicit positioning info
if keyword_set(xpos) then x0 = xpos
if keyword_set(ypos) then y0 = ypos


menu = Widget_Base( Title = title_str, xoffset=x0, yoffset=y0 )


;   display no more than twenty elements of the menu at a time
;   if there are less then only display those
ysize = min( [20, n_elements(options)] ) ; in character lines not pixels
l = Widget_List( menu, Value=options, YSize = ysize )

item = ''

Widget_Control, menu, /REALIZE                                     

While (item eq '') do begin
                                                                          
             IF(KEYWORD_SET(demowait)) THEN BEGIN                               
                 WAIT, demowait                                                 
                 IF(NOT KEYWORD_SET(retain)) THEN Widget_Control, menu, /DESTROY
                 RETURN, init                                                   
                 ENDIF                                                          
                                                                                
             event = Widget_Event( menu )                                       
                                                                                
             index = event.index                                                
             item  = StrTrim( options(index), 2 )                               
                                                                                
             EndWhile



IF(NOT KEYWORD_SET(retain)) THEN Widget_Control, menu, /DESTROY

IF(title NE -1) THEN IF(index GE title) THEN index=index+1

RETURN,index

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


