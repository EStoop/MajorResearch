FUNCTION bigwmenu,strings,title=title,init=init,xpos=xpos,ypos=ypos, $
  retain=retain,demowait=demowait
;+NAME/ONE LINE DESCRIPTION:
;    BIGWMENU puts up a one-column menu on an X-window terminal.
;
;DESCRIPTION:
;    This routine is analogous to the intrinsic WMENU (they both put
;    up a one-column menu on an X-window terminal), only BIGWMENU
;    makes a menu-window that has bigger characters.  (The algorithm
;    came from the STAR package).
;
;CALLING SEQUENCE:
;    sel_index = bigwmenu( strings, title=[...], init=[...], $
;       xpos=[...], ypos=[...], retain=[...])
;
;ARGUMENTS (I=input, O=output, []=optional)
;    strings    I   1-D arr  str        A 1-D array of option strings
;                                       (and possibly a title too).
;    title     [I]  keyword  int        Index of the element in STRINGS
;                                       which will be the title.
;    init      [I]  keyword  int        Index of the element in STRINGS
;                                       which will be the default choice.
;    xpos      [I]  keyword  int        X position (in screen coords)
;                                       of the menu window.
;    ypos      [I]  keyword  int        Y position (in screen coords)
;                                       of the menu window.
;    retain    [I]  keyword  int        If set, then the menu window
;                                       will remain on the screen after
;                                       the user makes his selection.
;    sel_index  O   scalar   int        Index of the selected option.
;
;WARNINGS:
;    1. A value of -1 is returned in the event of any error-conditions.
;    2. If a color table is not set beforehand, then this routine will
;          establish the RED TEMPERATURE color table.
;
;EXAMPLE:
;    sel_index = bigwmenu(['title','option 1','option 2'],title=0)
;#
;COMMON BLOCKS:  colors.
;
;LIBRARY CALLS:  none.
;
;PROCEDURE:
;    If this routine is invoked from a non-X-window terminal, then put
;    out a message and exit.  Check if a valid argument was supplied
;    for STRINGS, if not then put out a message and exit.  Determine
;    what will be the title of the menu (if any) and what will be the
;    options.  If more options were suppled than how many can fit on
;    the screen, then put out a message and exit.  Calculate an appro-
;    priate X-size for the menu window.  If a color table is not defined,
;    then define one (use RED TEMPERATURE by default).  Put up the menu
;    window, then monitor and respond to mouse activity.  Exit when a
;    mouse button is pressed.
;
;REVISION HISTORY:
;    Originally written by the authors of the STAR package, to whom the
;    credit belongs.
;    Put in current form (small changes made) J Ewing, ARC, March 1992.
;    Changed to not draw boxes for blanks. J Ewing  Sept. 17, 1992.
; SPR 10264  Dec 03,92  Special color values for UIMAGE.  J Ewing
;
;.TITLE
;Routine BIGWMENU
;-
  COMMON colors,r_orig,g_orig,b_orig,r_curr,g_curr,b_curr
  IF(!D.NAME NE 'X') THEN BEGIN
    MESSAGE,'This routine can only be run on an X-window terminal.',/CONT
    RETURN,-1
  ENDIF
  locolor = 120
  hicolor = 200
  letcolor = 255
;
;  Redefine those values if run from UIMAGE.
;  -----------------------------------------
  COMMON history,uimage_version
  sz = SIZE(uimage_version)
  IF(sz(0)+sz(1) NE 0) THEN IF(uimage_version GT 0) THEN BEGIN
    locolor = 7
    hicolor = 8
    letcolor = 9
  ENDIF

  old_window = !D.WINDOW
  IF(!D.NAME NE 'vms') THEN DEVICE, PSEUDO_COLOR = 8
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
  height = 30
  csize = 1.6
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
    xsize_title = 74+STRLEN(title_str)*9
    num_options = num_str - 1
    options = STRARR(num_options)
    ii = 0
    FOR i=0,num_str-1 DO IF(i NE title) THEN BEGIN
      options(ii) = strings(i)
      ii = ii+1
    ENDIF
  ENDELSE
  options = STRTRIM(options, 2)
;
;  Check if too many strings were supplied.
;  ----------------------------------------
  DEVICE,GET_SCREEN_SIZE=scrsiz
  IF((16+height*num_options) GT scrsiz(1)) THEN BEGIN
    MESSAGE,'Too many strings were supplied.',/CONT
    RETURN,-1
  ENDIF
;
;  Calculate the X-size of the menu window.
;  ----------------------------------------
  WINDOW,/free,xsiz=500,/pixmap
  maxlen = 0
  FOR i=0,num_options-1 DO BEGIN
    XYOUTS,0,0,'!3'+options(num_options-1-i),size=csize,/device,width=width
    width = width * 500
    IF(width GT maxlen) THEN maxlen=width
  ENDFOR
  WDELETE
  xsize_options = (maxlen + (8 * 2)) > 80
  xsize = MAX([xsize_title, xsize_options])
;
;  Define a color table if one isn't already defined.
;  --------------------------------------------------
  so = SIZE(r_orig)
  sc = SIZE(r_curr)
  IF((so(0)+so(1) EQ 0) AND (sc(0)+sc(1) EQ 0)) THEN BEGIN
    quiet_orig = !quiet
    !quiet = 1
    LOADCT,3
    !quiet = quiet_orig
  ENDIF
;
;  Put up the menu window.
;  -----------------------
  sz = SIZE(xpos)
  IF(sz(0)+sz(1) EQ 0) THEN xpos = (scrsiz(0)-xsize)/2
  sz = SIZE(ypos)
  IF(sz(0)+sz(1) EQ 0) THEN ypos = (scrsiz(1)-(16+height*num_options))/2
  WINDOW,/free,xpos=xpos,ypos=ypos,xsize=xsize,ysize=num_options*height,$
    title=title_str,retain=2
  FOR i=0,num_options-1 DO XYOUTS,8,i*height+10,'!3'+options(num_options-1-i),$
    size=csize,/device,color=letcolor
  xx=[2,2,2,xsize-2,xsize-2,xsize-2,xsize-2,2]
  sz = SIZE(init)
  IF(sz(0)+sz(1) NE 0) THEN BEGIN
    sel = init
    IF((title EQ -1) OR (init LE title)) THEN BEGIN
      yinit = height*(num_options-init-1) + height/2
    ENDIF ELSE BEGIN
      yinit = height*(num_options-init) + height/2
      sel = init - 1
    ENDELSE
    TVCRS,xsize/2,yinit
  ENDIF ELSE sel = 0
  osel = 0
  WHILE(options(osel) EQ '') DO osel = osel + 1
  ybox = INTARR(8,num_options)
  FOR i=0,num_options-1 DO BEGIN
    y1 = (num_options - i) * height - (height - 2)
    y2 = (num_options - i) * height - 2
    ybox(0,i) = [y1,y2,y2,y2,y2,y1,y1,y1]
  ENDFOR
  IF(KEYWORD_SET(demowait)) THEN BEGIN
    PLOTS,xx,ybox(*,sel),/DEVICE,THICK=2.5,COLOR=hicolor
    FOR i=0,num_options-1 DO IF((i NE sel) AND (options(i) NE '')) THEN $
      PLOTS,xx,ybox(*,i),/DEVICE,THICK=2.5,COLOR=locolor
    WAIT, demowait
    IF(NOT KEYWORD_SET(retain)) THEN WDELETE
    IF(old_window NE -1) THEN WSET,old_window
    RETURN, init
  ENDIF
;
;  Monitor mouse movements.
;  ------------------------
  first = 1
  !ERR=0
  WHILE(!ERR EQ 0) DO BEGIN
    CURSOR,x,y,0,/DEVICE
    IF((y NE -1) OR (first EQ 1)) THEN BEGIN
      sel = num_options - FIX(y / height) - 1
      IF(options(sel) EQ '') THEN sel = osel
      PLOTS,xx,ybox(*,sel),/DEVICE,THICK=2.5,COLOR=hicolor
      FOR i=0,num_options-1 DO IF((i NE sel) AND (options(i) NE '')) THEN $
        PLOTS,xx,ybox(*,i),/DEVICE,THICK=2.5,COLOR=locolor
    ENDIF
    first = 0
    osel = sel
  ENDWHILE
  IF(title NE -1) THEN IF(sel GE title) THEN sel=sel+1
  IF(NOT KEYWORD_SET(retain)) THEN WDELETE
  IF(old_window NE -1) THEN WSET,old_window
  RETURN,sel
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


