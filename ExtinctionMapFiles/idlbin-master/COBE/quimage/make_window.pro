PRO make_window, window, nblocks_x, nblocks_y, title, error_status
;+
;  MAKE_WINDOW - a UIMAGE-specific routine.  This routine creates a window
;  at a size and position in conformance with the UIMAGE X-window screen
;  management scheme.
;#
;  Written by John Ewing.
;  SPR 10294  Dec 03 92  Changed BLOCKSIZE_X.  J Ewing
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;--------------------------------------------------------------------------
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  IF ((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    error_status = -1
    RETURN
  ENDIF
  IF(N_PARAMS(0) LT 4) THEN BEGIN
    MESSAGE, 'All four arguments must be supplied.', /CONT
    error_status = -1
    RETURN
  ENDIF
  error_status = 0
;
;  If the screen size was not previously determined, then determine it now.
;  ------------------------------------------------------------------------
  sz = SIZE(screensize_x)
  IF(sz(0) + sz(1) EQ 0) THEN BEGIN
    DEVICE, RETAIN = 2
    DEVICE, GET_SCREEN_SIZE = ss
    screensize_x = ss(0) + 2
    screensize_y = ss(1) - 2
  ENDIF
;
;  Specify the size of a UIMAGE screen "block".
;  --------------------------------------------
  blocksize_x = 58
  blocksize_y = 47
  maxblocks_x = screensize_x / blocksize_x
  maxblocks_y = screensize_y / blocksize_y
  IF(maxblocks_y * blocksize_y EQ screensize_y) THEN $
    maxblocks_y = maxblocks_y + 1
;
;  Handle the situation in which a requested window would be too big
;  to fit on the screen.
;  -----------------------------------------------------------------
  IF((nblocks_x GT maxblocks_x) OR (nblocks_y GT maxblocks_y)) THEN BEGIN
    error_status = -1
    IF(zoom_index GT 0) THEN BEGIN
      PRINT, 'An image is too big to be displayed with the current ' + $
        'window magnification.'
      zf = STRING(2^zoom_index, '(i1)')
      PRINT, '(The window magnification is currently ' + bold(zf) + ').'
      menu = ['Set the window magnification to 1?', 'Yes', 'No']
      sel = one_line_menu(menu)
      IF(sel EQ 2) THEN BEGIN
        PRINT, 'Image is not displayed.'
        RETURN
      ENDIF
      zoom_index = 0
      xdisplay_all
      RETURN
    ENDIF
    MESSAGE, /CONT, $
      'Image is too big to be displayed with the current zoom factor.'
    RETURN
  ENDIF
  sz = SIZE(block_usage)
  IF(sz(0) EQ 0) THEN block_usage = INTARR(maxblocks_y, maxblocks_x)
  xspace = blocksize_x * nblocks_x
  yspace = blocksize_y * nblocks_y
  windowsize_x = xspace - 32
  windowsize_y = yspace - 47
  sz = SIZE(block_usage)
  IF(sz(0) EQ 2) THEN nlevels = 1 ELSE nlevels = sz(3)
;
;  The following two statements create an area named WORK which has a value
;  of 0 for those elements which correspond to blocks on the screen which
;  are vacant.
;  ------------------------------------------------------------------------
  work = INTARR(maxblocks_y, maxblocks_x)
  FOR i = 0, nlevels - 1 DO work = work + block_usage(*, *, i)
  level = 0
  first = 1
find_x1y1:
;
;  If an all-zero subsection of WORK is identified, of the appropriate size,
;  (i.e., an available space of the right size on the screen has been
;  identified) then branch to label HAVE_X1Y1.
;  -------------------------------------------------------------------------
  FOR x1 = 0, maxblocks_x - 1 DO BEGIN
    x2 = x1 + nblocks_x - 1
    IF(x2 LE maxblocks_x - 1) THEN BEGIN
      FOR y1 = 0, maxblocks_y - 1 DO BEGIN
        y2 = y1 + nblocks_y - 1
        IF(y2 LE maxblocks_y - 1) THEN $
          IF(max(work(y1:y2, x1:x2)) EQ 0) THEN GOTO, have_x1y1
      ENDFOR
    ENDIF
  ENDFOR
;
;  Set WORK to indicate which blocks are being used at a given level of
;  overlay (i.e., for the third dimension of BLOCK_USAGE being constant).
;  First time through, use the highest level, then successively go through
;  lower levels.
;  -----------------------------------------------------------------------
  IF(first EQ 1) THEN BEGIN
    level = nlevels - 1
    first = 0
  ENDIF ELSE BEGIN
    IF(level GT 0) THEN BEGIN
      level = level - 1
    ENDIF ELSE BEGIN
;
;  Add another layer to BLOCK_USAGE.  The requested space will be given
;  from within this new layer (overlapping lower layers).
;  --------------------------------------------------------------------
      level = nlevels
      temp = block_usage
      block_usage = INTARR(maxblocks_y, maxblocks_x, nlevels + 1)
      FOR i = 0, nlevels - 1 DO block_usage(0, 0, i) = temp(*, *, i)
      nlevels = nlevels + 1
    ENDELSE
  ENDELSE
  work = block_usage(*, *, level)
  GOTO, find_x1y1
have_x1y1:
;
;  At this point, the values of X1 and Y1 indicate where on the screen
;  the new screen-object should be placed.  Translate this information
;  into units of screen coordinates and RETURN those coordinates.
;  -------------------------------------------------------------------
  xp = screensize_x - x1 * blocksize_x - xspace
  yp = screensize_y - y1 * blocksize_y - yspace
;
; funny business on DECstation
  IF (!version.os eq 'ultrix') and (windowsize_x lt 100) then windowsize_x=100
  window, window, xsiz = windowsize_x, ysiz = windowsize_y, xpos = xp, $
    ypos = yp, title = title
  block_usage(y1:(y1 + nblocks_y - 1), x1:(x1 + nblocks_x - 1), level) = $
    window + 1
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


