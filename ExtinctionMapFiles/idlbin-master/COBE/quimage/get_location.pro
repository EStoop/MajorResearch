FUNCTION get_location, name, nblocks_x = nblocks_x, nblocks_y = nblocks_y
;+
;  GET_LOCATION - a UIMAGE-specific routine.  This routine is only
;  actually invoked when UIMAGE runs at what is currently referenced
;  as version #4.  This function returns screen coordinates at which
;  a screen-object should be placed.  It calculates those screen
;  coordinates by the same algorithm as is used by MAKE_WINDOW (which
;  is not invoked under UIMAGE version #4).  There is one required
;  argument to this routine, NAME, which is a string which identifies
;  a UIMAGE data-object (such as "MAP6(0)").
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
;--------------------------------------------------------------------------
  COMMON xscreen,magnify,block_usage,scrsiz
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  IF(NOT defined(magnify)) THEN magnify = 2
;
;  If BLOCK_USAGE is undefined, then create it.  Determine its dimensions.
;  -----------------------------------------------------------------------
  sz = SIZE(block_usage)
  IF(sz(0) + sz(1) EQ 0) THEN BEGIN
    DEVICE, GET_SCREEN_SIZE = scrsiz
    maxblocks_x = (scrsiz(0) - 16) / 32
    maxblocks_y = (scrsiz(1) - 48) / 32
    block_usage = LONARR(maxblocks_y, maxblocks_x)
  ENDIF ELSE BEGIN
    maxblocks_y = sz(1)
    maxblocks_x = sz(2)
  ENDELSE
  IF(sz(0) LE 2) THEN nlevels = 1 ELSE nlevels = sz(3)
;
;  If the size for the new screen-object was not specified (via the NBLOCK*
;  parameters, then determine the appropriate size now (in units of blocks).
;  -------------------------------------------------------------------------
  IF(NOT KEYWORD_SET(nblocks_y)) THEN BEGIN
    IF(STRMID(name, 0, 5) EQ 'GRAPH') THEN BEGIN
      j = EXECUTE('nblocks_x = ' + name + '.nblocks_x * magnify')
      j = EXECUTE('nblocks_y = ' + name + '.nblocks_y * magnify')
    ENDIF ELSE BEGIN
      j = EXECUTE('sz = SIZE(' + name + '.data)')
      nblocks_x = sz(1)/32 * magnify
      nblocks_y = sz(2)/32 * magnify
    ENDELSE
    nblocks_x = nblocks_x + 1
    nblocks_y = nblocks_y + 3
  ENDIF
;
;  The following two statements create an area named WORK which has a value
;  of 0 for those elements which correspond to blocks on the screen which
;  are vacant.
;  ------------------------------------------------------------------------
  work = LONARR(maxblocks_y, maxblocks_x)
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
          IF(TOTAL(work(y1:y2, x1:x2)) EQ 0) THEN GOTO, have_x1y1
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
      block_usage = LONARR(maxblocks_y, maxblocks_x, nlevels + 1)
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
  xoffset = scrsiz(0) - (x1 * 32) - (nblocks_x * 32)
  yoffset = 56 + (y1 * 32)
  block_usage(y1:(y1 + nblocks_y - 1), x1:(x1 + nblocks_x - 1), level) = -1
  RETURN, [xoffset, yoffset]
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


