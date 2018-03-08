PRO show_zoomed, name
;+
;  SHOW_ZOOMED - a UIMAGE-specific routine.  This routine will display
;  a zoomed section of an image on the currently set X-window.
;#
;  Written by John Ewing.
;  SPR 10383  Jan 07 93  Add more comments.  J Ewing
;  SPR 10456  Jan 14 93  Bound prologue comments by ";+" and ";-".  J Ewing
;  SPR 11169  Jul 21 93  Work with zoomed output for bckgrnd fits. J Newmark
;  SPR 11226  Aug 18 93  Add large reprojections. J. Newmark
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;--------------------------------------------------------------------------
  COMMON uimage_data,map6,map7,map8,map9,face6,face7,face8,face9,graph,zoomed
  COMMON uimage_data2,proj_map,proj2_map
  COMMON xwindow,screensize_x,screensize_y,blocksize_x,blocksize_y,$
    block_usage,zoom_index
  COMMON color_values,c_badpix,c_draw,c_scalemin
  COMMON zback,zbgr,zdsrcmap,zbsub
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine can only be run on an X-window terminal.', /CONT
    RETURN
  ENDIF
;
;  Extract appropriate info from the ZOOMED structure, and from the structure
;  that is associated with the original unzoomed image.
;  --------------------------------------------------------------------------
  j = EXECUTE('zoomwin = ' + name + '.window')
  WSET, zoomwin
  ERASE
  j = EXECUTE('ztitle = ' + name + '.title')
  j = EXECUTE('nblocks_x = ' + name + '.nblocks_x')
  j = EXECUTE('nblocks_y = ' + name + '.nblocks_y')
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
    MESSAGE,'The user has deleted the data associated with a zoom window.',/CONT
    RETURN
  ENDIF
  j = EXECUTE('start_x = ' + name + '.start_x')
  j = EXECUTE('start_y = ' + name + '.start_y')
  j = EXECUTE('stop_x = ' + name + '.stop_x')
  j = EXECUTE('stop_y = ' + name + '.stop_y')
  j = EXECUTE('specific_zoom = ' + name + '.specific_zoom')
  j = EXECUTE('scale_min = ' + name + '.scale_min')
  j = EXECUTE('scale_max = ' + name + '.scale_max')
  j = EXECUTE('bad = WHERE('+name_orig+'.data EQ '+name_orig+'.badpixval)')
  j = EXECUTE('data = ' + name_orig + '.data')
  sz=SIZE(data)
  dim1_orig = sz(1)
  dim2_orig = sz(2)
  total_zoom = specific_zoom * 2^zoom_index
  xp = 0
  yp = 0
  IF(start_x LT 0) THEN BEGIN
    xp = (-1) * start_x * total_zoom
    start_x = 0
  ENDIF
  IF(start_y LT 0) THEN BEGIN
    yp = (-1) * start_y * total_zoom
    start_y = 0
  ENDIF
  stop_x = stop_x < (dim1_orig - 1)
  stop_y = stop_y < (dim2_orig - 1)
  IF (bflag EQ 1) THEN BEGIN
    IF (zoomflag EQ 1) THEN data(start_x:stop_x,start_y:stop_y)=zbgr
    IF (zoomflag EQ 2) THEN data(start_x:stop_x,start_y:stop_y)=zdsrcmap
    IF (zoomflag EQ 3) THEN data(start_x:stop_x,start_y:stop_y)=zbsub
  ENDIF
;
;  Set IMAGE to be a byte-scaled version of the original unzoomed data.
;  --------------------------------------------------------------------
  IF(scale_min LT scale_max) THEN BEGIN
    j = EXECUTE('vimage=BYTSCL(data,' + $
      'min=scale_min,max=scale_max,top=255-c_scalemin)+c_scalemin')
  ENDIF ELSE BEGIN
    j = EXECUTE('vimage=BYTE(data*0.)+255')
  ENDELSE
  IF(bad(0) NE -1) THEN vimage(bad) = c_badpix
;
;  Set IMAGE to be only that part which was zoomed in on.
;  ------------------------------------------------------
  vimage = vimage(start_x:stop_x, start_y:stop_y)
;
;  Use CONGRID to get IMAGE to be the right size, and then display it.
;  -------------------------------------------------------------------
  dim1 = 1 + stop_x - start_x
  dim2 = 1 + stop_y - start_y
  vimage = CONGRID(vimage, dim1 * total_zoom, dim2 * total_zoom)
  TV, vimage, xp, yp
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


