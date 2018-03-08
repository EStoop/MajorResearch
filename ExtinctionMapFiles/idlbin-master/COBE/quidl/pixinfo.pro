        FUNCTION pixinfo, data, outcoord=outcoord, win=win, $
        proj=proj, pcoords=pcoords, face=face, wpos=wpos, $
        pixels=pixels, xout=xout, yout=yout, datval=datval
;
;+NAME/ONE-LINE DESCRIPTION:
;  PIXINFO returns lon/lat, x/y, pixel number, and data value 
;
;DESCRIPTION: 
;  PIXINFO allows the user to choose different points on an image
;  with the mouse and find out the points' x/y coordinates,
;  longitude/latitude in any coordinate system, pixel number, and
;  data value.  Longitudes and latitude are returned in an (n,2)
;  array.  X/Y coordinates, pixel number, and data value may also
;  be returned by using the "xout", "yout", "pixels", and "datval"
;  keywords.  Users have the option of removing the marks 
;  (refreshing the image).
;
;CALLING SEQUENCE:
;  info = PIXINFO (pixinfo, data, outcoord=outcoord, win=win, $
;         proj=proj, pcoords=pcoords, face=face, wpos=wpos, $
;         pixels=pixels, xout=xout, yout=yout, datval=datval)
;	
;
; ARGUMENTS: (I=input, O=output, []=optional)
;  data	      I     2-D or 3-D data array
; [outcoords] I     Coordinates that the longitudes and latitudes
;                   will be given in.
; [win]       I     Window number where the display lives. 
; [proj]      I     Projection type if data array is a reprojection.
;                   Use proj='a'  for Aitoff
;                       proj='s'  for Global Sinusoidal
;                       proj='m'  for Mollweide
; [pcoords]   I     Coordinates of the projection.  If data array
;                   is a skycube or face, pcoords is assumed to be
;                   ecliptic.  For reprojections, use
;                       pcoords='g'  for galactic coords
;                       pcoords='e'  for ecliptic
;                       pcoords='q'  for equatorial
; [face]      I     If data array is a face, you may specify the 
;                   face number in the call.
; [wpos]      I     2-element array indicating Window POSitioning
;                   of the image.  Use form wpos=[x,y] where x is the
;                   X coordinate of the lower left corner of the
;                   image and y is the Y coordinate of the lower left 
;                   corner of the image.  i.e., if the image has a 32 
;                   screen pixel border on each side, wpos=[32,32]
;                   Default=[0,0]
; [pixels]    O     Pixel list of the pixels selected.
; [xout]      O     List of X coordinates of selected points.
; [yout]      O     List of Y coordinates of selected points.
; [datval]    O     List of data values of selected points.
;			
;	
; EXAMPLES:
;   To return the galactic longitude/latitude of some points on a 
;   skycube (data array=CUBE) displayed in IDL window 1, type
;
;        UIDL>  lons_lats=PIXINFO (cube, outcoord='g', win=1)
;
;   Click left mouse button to mark the points and right mouse button
;   to quit.  The galactic coordinates are held in LONS_LATS.
;	
;
;#
; SUBROUTINES CALLED:
;   u_menu, coorconv, proj2uv, xy2pix, validnum
;
;
; COMMON BLOCKS:  None
;
; LIBRARY CALLS:  None
;
; WARNINGS:
;   Will only run on certain reprojections: 
;       Aitoff:  Ecl (Sm,Lg),  Gal (Sm,Lg),  Equ (Sm) 
;       Mollwd:  Ecl (Sm,Lg),  Gal (Sm,Lg),  Equ (Sm) 
;       Gl.Sin:  Ecl (Sm)   ,  Gal (Sm)   ,  Equ (Sm)
;	
;   Skycubes and faces are assumed to be in ecliptic coordinates
;
; PROGRAMMING NOTES:
;
; MODIFICATION HISTORY:
;   Written by Celine Groden, USRA    April 1993  (SPR 10499)
;   Modified documentation CMG, USRA  May 1993 
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;
;.TITLE
; Routine PIXINFO
;
;-
;
; Establish error condition to return to calling routine in case of
; trouble.
;
ON_ERROR, 2
;
;  Miscellaneous errors.
;
IF (!D.NAME NE 'X' AND !D.NAME NE 'WIN') THEN MESSAGE, $
'Must be on an X Windows terminal.'
projsz=SIZE(data)
IF (projsz(0) EQ 0) THEN MESSAGE, 'No input array was supplied'
IF (projsz(0) NE 2 AND projsz(0) NE 3) THEN MESSAGE, $
   'Input array must be 2-D or 3-D'
IF (projsz(1) NE projsz (2) AND $
    projsz(1) NE 2*projsz(2) AND $
    3*projsz(1) NE 4*projsz(2)) THEN MESSAGE, $
    'Input array must be a skycube, face, or reprojection.'
;
; Get window number and activate proper window.  Set WPOS default if
; keyword was not specified.
; 
IF (N_ELEMENTS(win) EQ 0) THEN BEGIN
   win=''
   READ,'Enter IDL window number: ', win
   IF (VALIDNUM(win) EQ 0) THEN MESSAGE, 'Invalid window number.'
   win=FIX(win)
ENDIF
DEVICE, WINDOW_STATE=ws
IF (ws(win) EQ 0) THEN MESSAGE, 'Window is not available.'
WSET, win
IF (N_ELEMENTS(wpos) EQ 0) THEN wpos=[0,0]
;
;  Save original image in variable in case the user wants to remove
;  the marks at the end of the program.
;
vimage = TVRD (0, 0, !D.X_VSIZE, !D.Y_VSIZE)
;
;  Get reprojection type and coordinate system if not specified.
;
IF (projsz(1) EQ 2*projsz(2)) THEN BEGIN
   IF (N_ELEMENTS(proj) EQ 0) THEN BEGIN
      proj_types=['A', 'S', 'M']
      menu=['Enter Projection Type', 'Aitoff', 'Global Sinusoidal', $
	    'Mollweide'] 
      answer=UMENU(menu,title=0,init=1)
      proj=proj_types(answer-1)
   ENDIF ELSE proj=STRUPCASE(proj)
   IF (N_ELEMENTS(pcoords) EQ 0) THEN BEGIN
      pcoords_types=['E', 'G', 'Q']
      menu=['Enter Projection Coordinate System', 'Ecliptic', 'Galactic', $
	    'Equatorial']
      answer=UMENU(menu,title=0,init=1)
      pcoords=pcoords_types(answer-1)
   ENDIF ELSE pcoords=STRUPCASE(pcoords)
ENDIF
;
; For skycubes faces, determine resolution and set coord system to
; ecliptic.  
;
IF (3*projsz(1) EQ 4*projsz(2)) THEN BEGIN
   proj='SKYCUBE'  &  pcoords='E'
   res=FIX(2*ALOG((projsz(1)/4))/ALOG(4) + 1)
   res_cc='R'+STRTRIM(STRING(res),1)
ENDIF
IF (projsz(1) EQ projsz(2)) THEN BEGIN
   proj='FACE'  &  pcoords='E'
   res=FIX(2*ALOG((projsz(1)))/ALOG(4) + 1)
   res_cc='R'+STRTRIM(STRING(res),1)
   IF (N_ELEMENTS(face) EQ 0) THEN BEGIN
      face=''
      READ,'Enter face number: ', face
      IF (VALIDNUM(face) EQ 0) THEN MESSAGE, 'Invalid face number.'
      face=FIX(face)
      IF (face LT 0 OR face GT 5) THEN MESSAGE, 'Invalid face number.'
   ENDIF
ENDIF
;
; Determine the coordinate system of the output if not specified in
; the call.
;
IF (N_ELEMENTS(outcoord) EQ 0) THEN BEGIN
   outcoord_types=['E', 'G', 'Q']
   menu=['Enter Coordinate System of Output', 'Ecliptic', 'Galactic', $
	 'Equatorial']
   answer=UMENU(menu,title=0,init=1)
   outcoord=outcoord_types(answer-1)
ENDIF ELSE outcoord=STRUPCASE(outcoord)
;
;
; Mask arrays:
; Set up masking arrays for determining if points are actually on the
; image or are on the background.  Set array to zero where there is no 
; data.  For reprojections, use Joel's masks.
;
mk=''
mask=INTARR(!D.X_VSIZE, !D.Y_VSIZE)
IF (proj EQ 'SKYCUBE' OR proj EQ 'FACE') THEN BEGIN
   side_len = FIX(SQRT(4^(FLOAT(res)-1)))
   IF (proj EQ 'SKYCUBE') THEN BEGIN
      mask(0+wpos(0):4*side_len-1+wpos(0), 0+wpos(1):3*side_len-1+wpos(1))=1
      mask(0+wpos(0):3*side_len-1+wpos(0), 0+wpos(0):side_len-1+wpos(1))=0
      mask(0+wpos(0):3*side_len-1+wpos(0), $
      2*side_len+wpos(1):3*side_len-1+wpos(1))=0
   ENDIF ELSE mask(0+wpos(0):side_len-1+wpos(0),0+wpos(1):side_len-1+wpos(1))=1
ENDIF ELSE BEGIN
   res_cc=''
   READ,'Resolution of original skycube: ', res_cc
   IF (VALIDNUM(res_cc) EQ 0) THEN MESSAGE, 'Invalid resolution.'
   res_cc='R'+STRTRIM(res_cc,1)
   IF (projsz(1) EQ 1024) THEN BEGIN
      repjsz='big_' 
      mk='msk'
      xdim=1024  &  ydim=512
   ENDIF ELSE BEGIN
      repjsz=''
      mk='mask'
      xdim=512  &  ydim=256
   ENDELSE
   IF (pcoords EQ 'G') THEN cd='gal'
   IF (pcoords EQ 'Q') THEN cd='equ'
   IF (pcoords EQ 'E') THEN cd='ecl'
   IF (proj EQ 'A') THEN prj='ait'
   IF (proj EQ 'S') THEN prj='gs'
   IF (proj EQ 'M') THEN prj='mwd'
   dirspec=getenv('CGIS_DATA')
   filename=dirspec+repjsz+prj+'_'+cd+'.lut'
   big_ait_gal_msk=0&big_mwd_gal_msk=0&big_ait_ecl_msk=0&big_mwd_ecl_msk=0
   big_gs_gal_msk=0 & big_gs_ecl_msk=0
   ait_gal_mask=0 & ait_ecl_mask=0 & ait_equ_mask=0
   mwd_gal_mask=0 & mwd_ecl_mask=0 & mwd_equ_mask=0
   gs_gal_mask=0  & gs_ecl_mask=0  & gs_equ_mask=0
   ait_gal_i=0    & ait_ecl_i=0    & ait_equ_i=0
   ait_gal_j=0    & ait_ecl_j=0    & ait_equ_j=0
   mwd_gal_i=0    & mwd_ecl_i=0    & mwd_equ_i=0
   mwd_gal_j=0    & mwd_ecl_j=0    & mwd_equ_j=0
   gs_gal_i=0     & gs_ecl_i=0     & gs_equ_i=0
   gs_gal_j=0     & gs_ecl_j=0     & gs_equ_j=0
   RESTORE, filename
   i=EXECUTE("mask(wpos(0):xdim+wpos(0)-1, wpos(1):ydim+wpos(1)-1)="+repjsz+prj+'_'+cd+'_'+mk)
   one=WHERE(mask NE 255)
   zero=WHERE(mask EQ 255)
   mask(one)=1
   mask(zero)=0
   IF (TOTAL(wpos) NE 0) THEN BEGIN
      mask(0:!D.X_VSIZE-1, 0:wpos(1)-1)=0
      mask(0:wpos(0)-1, 0:!D.Y_VSIZE-1)=0
      mask(!D.X_VSIZE-wpos(0):!D.X_VSIZE-1, wpos(1):!D.Y_VSIZE-1)=0
      mask(wpos(0):!D.X_VSIZE-1, !D.Y_VSIZE-wpos(1):!D.Y_VSIZE-1)=0
   ENDIF
ENDELSE
;
; Get the points, display the information, and save the info in output
; arrays.
;
!ERR=0 & i=0
lonarr=0 & latarr=0
PRINT, ' '
PRINT, 'Mark points with <MB1>. Press <MB3> to exit.'
PRINT, ' '
WHILE (!ERR NE 4) DO BEGIN
   CURSOR, x, y, /DOWN, /DEVICE
   PLOTS, x, y, PSYM=1, /DEVICE
   IF (mask(x,y) EQ 0) THEN PRINT, 'Point is not on map.' $
   ELSE BEGIN
      x=x-wpos(0)
      y=y-wpos(1)
      IF (proj EQ 'SKYCUBE') THEN BEGIN
	 pixel = XY2PIX (x,y,res=res)
	 ll = COORCONV (pixel,infmt='P',inco=res_cc,outfmt='L',outco=outcoord)
      ENDIF
      IF (proj EQ 'FACE') THEN BEGIN
	 pixel = XY2PIX (x,y,res=res,face=face)
	 ll = COORCONV (pixel,infmt='P',inco=res_cc,outfmt='L',outco=outcoord)
      ENDIF
      IF (projsz(1) EQ 2*projsz(2)) THEN BEGIN
	 xy=FLTARR(1,2)
	 xy(0,0)=x  &  xy(0,1)=y
	 uvec = PROJ2UV (xy, proj, projsz)
	 pixel=COORCONV(uvec,infmt='U',inco=pcoords,outfmt='P',outco=res_cc)
	 ll = COORCONV (uvec,infmt='U',inco=pcoords,outfmt='L',outco=outcoord)
      ENDIF
      lon=STRTRIM(ll(0),1)
      lat=STRTRIM(ll(1),1)
      If (i EQ 0) THEN BEGIN
         lonarr = lon
         latarr = lat
	 xout   = x
	 yout   = y
         pixels = pixel
	 datval = data (x, y)
      ENDIF ELSE BEGIN
         lonarr = [lonarr, lon]
         latarr = [latarr, lat]
	 xout   = [xout, x]
	 yout   = [yout, y]
         pixels = [pixels, pixel]
	 datval = [datval, data (x, y)]
      ENDELSE
      PRINT, 'x= ',x,'          y= ',y, FORMAT='(a3,i3,a13,i3)'
      PRINT, 'lon= ',lon,'     lat= ',lat, FORMAT='(a5,f6.2,a10,f6.2)'
      PRINT,'pixel=',pixel,'   Data Value=',datval(i),FORMAT='(a6,i7,a14,e)'
      PRINT, ' '
      i=i+1
   ENDELSE
ENDWHILE
!ERR=0
;
; If the user wants to refresh the image, do so.
;
YN=['Y', 'N']
menu=['Remove Points from image?', 'Yes', 'No']
answer=UMENU(menu,title=0,init=1)
IF (YN(answer-1) EQ 'Y') THEN TV, vimage
;
; Put lon/lat coordinates into one array to be returned.
;
lonlat=FLTARR(N_ELEMENTS(lonarr),2)
lonlat(*,0)=lonarr
lonlat(*,1)=latarr
;
;
RETURN, lonlat
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


