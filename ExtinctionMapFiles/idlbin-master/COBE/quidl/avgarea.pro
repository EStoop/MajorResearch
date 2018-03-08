  	FUNCTION AVGAREA, data, weights=weights, $
	circle=circle, annulus=annulus,	facenum=facenum, $
	in_coords=in_coords, wpos=wpos, badval=badval, win=win, $
	pixels=pixels, regmin=regmin, regmax=regmax, standev=standev,$
	noplot=noplot
;
;+NAME/ONE-LINE DESCRIPTION:
;	AVGAREA finds the average over a circle or annulus
;
; DESCRIPTION: 
;	AVGAREA calculates a weighted or unweighted average over a given
;	circle or annulus.  A circle can be defined by an input
;	array containing the longitude, latitude, and radius of the
;	circle or X Window users can mark the circle center with a
;	mouse and then enter the radius from the keyboard.  An annulus
;	can be defined in the two ways similar to the ways a circle is
;	defined except two radii, the inner and the outer, are entered. 
;	The average over the specified area is reported, and the min,
;	max, standard deviation, and pixels in the region can also be
;	returned.
;
; CALLING SEQUENCE:
;	average = avgarea (data, circle=circle, $
;		  annulus=annulus, [weights=weights], [facenum=facenum],$
;		  [in_coords=in_coords], [wpos=wpos], [badval=badval], $
;		  [win=win], [pixels=pixels], [regmin=regmin], $
;		  [regmax=regmax], [standev=standev], [noplot=noplot])
;
;	Use only ONE of the following in call:  circle, annulus
;	
;
; ARGUMENTS: (I=input, O=output, []=optional)
;	data	   I	2-D or 3-D data array
;	circle     I	May be used one of two ways.  If keyword is
;			entered simply as "/circle", the center of the
;			circle will be marked with the mouse and the radius
;			will be entered from the keyboard.  The radius 
;			must be specified in degrees away from the circle
;			center, which acts as a pole around which the 
;			circle is drawn.  To enter the coordinates of the
;			center and the radius from the command line, use
;			an array in this form:
;			circle=[center_lon, center_lat, radius(in deg)].
;	annulus    I 	Specified in a manner similar to specifying a
;			circle, except that two radii are entered (an
;			inner then an outer) instead of a single radius.
;			Note: inner radius must be less than outer radius.
;			To enter the annulus info from the command line,
;			use the following form:
;			annulus=[cntr_lon,cntr_lat,in_radius,out_radius].
;	average	   O	Averaged data of the specified region.  For 2-D
;			input, one value is returned.  For 3-D input, 
;			an average over each slice is calculated and 
;			returned as an element in the average array.
;      [weights]	Array containing weights.  If no weight array is
;			specified, the weights are all set equal to 1.
;			Weight arrays are only 2-D, even if the array
;			is 3-D.
;      [in_coords] I	Coordinates that the polygon vertices or circle
;			center are given in.
;      [facenum]   I	If data array is a face, you may specify the 
;			facenum in the call.
;      [badval]    I	Bad pixel value.  These pixels will not be used
;			in the averaging.  Default value is 0.
;      [win]	   I	Window number where the display lives. 
;			N/A for non-X Window users.
;      [wpos]	   I	2-element array indicating Window POSitioning
;			of the image.  Use form wpos=[x,y] where x is the
;			X coordinate of the lower left corner of the
;			image and y is the Y coordinate of the lower left 
;			corner of the image.  i.e., if the image has a 32 
;			screen pixel border on each side, wpos=[32,32]
;			Default=[0,0]
;      [pixels]    O	Pixel list of the pixels in the polygon/circle/
;			annulus. 
;      [regmin]    O	Region (polygon,circle,annulus) minimum.
;      [regmax]    O	Region (polygon,circle,annulus) maximum.
;      [standev]   O	Standard deviation over region.
;      [noplot]	   I	Use the /noplot keyword if you do not want the
;			polygon or circle plotted on the image.  Non X
;			window users must use this option.
;			
;	
; EXAMPLES:
;
; To average a circle on a skycube called DATA using a weight array 
; called WTS and marking the center with the mouse, the call would 
; look like this:
;
;	UIDL>   circ_avg = avgarea (data, /circle, weights=wts)
;
; To do the same as above except specify the coordinates of the center
; and the radius size in the call, you would type:
;
;	UIDL>  circ_avg = avgarea (data, circle=[0,0,5], $
;	UIDL>  in_coords='g', weights=wts)
;
;
;
;#
; SUBROUTINES CALLED:
;	one_line_menu, coorconv, pix2xy, uv2proj, pix2fij, proj2uv, 
;	xy2pix, latcirc, validnum, pxincirc
;
;
; COMMON BLOCKS:  None
;
; LIBRARY CALLS:  None
;
; WARNINGS:
;   (1) If the program crashes, chances are good that your originally
;	2-D data array will now be a 3-D array. (see Programming Note
;	number 1.)
;
; PROGRAMMING NOTES:
;   (1) The main body of the code is written to operate on a 3-D object,
;	so if a 2-D array is input into DATA, AVGAREA converts it to a
;	degenerate 3-D array and then changes it back to a 2-D array 
;	before the program ends.  For this reason, all IDL MESSAGES are
;	preceded by an "If ThreeD Then REFORM" statement.  Also, in the
;	unlikely event of a program crash (cough,cough) the data array
;	will most likely be in a 3-D array.
;
;   (2) A "Mask" array is set up which is either the same size as the 
;	input array or the same size of the window where the image is.
;	The mask array is given	a value of 1 where the data is; else, 
;	mask is set to 0.  This is a quick way to check	if circle/annulus 
;	centers chosen with the mouse are actually on the image.
;
;
; MODIFICATION HISTORY:
;	Written by Celine Groden, USRA    April 1993
;
;.TITLE
; Routine AVGAREA
;
;-
;
; Establish error condition to return to calling routine in case of
; trouble.
;
ON_ERROR, 2
;
;
; Check that a proper input data array is being used.
;
projsz=SIZE(data)
IF (projsz(0) EQ 0) THEN MESSAGE, 'No input array was supplied'
IF (projsz(0) NE 2 AND projsz(0) NE 3) THEN MESSAGE, $
   'Input array must be 2-D or 3-D'
IF (projsz(1) NE projsz (2) AND $
    3*projsz(1) NE 4*projsz(2)) THEN MESSAGE, $
    'Input array must be a skycube or a face.'
IF (projsz(0) EQ 2) THEN BEGIN
   ThreeD=0
   slices=1
ENDIF
IF (projsz(0) EQ 3) THEN BEGIN
   ThreeD=1 
   slices=projsz(3)
ENDIF
;
;
; Set up a bogus weight array if one was not specified.  If an array
; was specified, check that it's ok.
;
IF (KEYWORD_SET(weights) EQ 0) THEN BEGIN
   weights=FLTARR(projsz(1),projsz(2),1)
   weights(*)=1.0
ENDIF
IF (KEYWORD_SET(weights)) THEN BEGIN
   wtsz=SIZE(weights)
   IF (wtsz(0) NE 2) THEN MESSAGE, 'Weight array must be 2-D.'
   IF (wtsz(1) NE projsz(1) OR wtsz(2) NE projsz(2)) THEN MESSAGE, $
   'Weight array must have same x,y dimensions as the data array.'
ENDIF
;
;
; Miscellaneous messages
;
IF (NOT KEYWORD_SET(circle) AND $
    NOT KEYWORD_SET(annulus)) THEN MESSAGE, $
   'Must specify keyword circle OR keyword annulus.'
IF (KEYWORD_SET(wpos)) THEN BEGIN
   IF (N_ELEMENTS(wpos) NE 2) THEN MESSAGE, 'Use WPOS=[x,y]'
ENDIF
IF (!D.NAME NE 'X' AND !D.NAME NE 'WIN') THEN BEGIN
   IF (NOT KEYWORD_SET(noplot)) THEN MESSAGE, $
   'Non-X Window users must specify noplot keyword in call.'
ENDIF
;
;
;
; Create size arrays for the polygon, circle, or annulus.  Terminate the 
; routine if a non-X Window user tries to use cursor input for polygon 
; vertices or circ/ann center.  Terminate program if (1) a 3-D data object was 
; used as input AND (2) an image is not present and the user wants to mark the
; region with the cursor (i.e., vertices or circle/annulus center and
; radii were not given on the comand line).  For polygons, check that an
; even number of coords were entered and at least 3 vertices are given.
; For annuli, check that the inner radius is smaller than the outer.
;
;
IF (KEYWORD_SET(circle )) THEN BEGIN
   circsz=SIZE(circle)
   IF (N_ELEMENTS(circle) EQ 1) THEN BEGIN
      IF (KEYWORD_SET(noplot)) THEN BEGIN
        MESSAGE,/CONT,'If NOPLOT keyword, enter circle info from command line.'
	MESSAGE,'e.g., circle = [center_lon, center_lat, radius]'
      ENDIF
      IF (!D.NAME NE 'X' AND !D.NAME NE 'WIN') THEN BEGIN
         MESSAGE,/CONT,'Non-X Window users must enter center and radius from'
         MESSAGE,/CONT,'command line.  Use the keyword CIRCLE in call, ' 
	 MESSAGE,'e.g., circle = [center_lon, center_lat, radius]'
      ENDIF
      IF (ThreeD AND NOT KEYWORD_SET(win)) THEN BEGIN
	 PRINT,' '
	 PRINT,'The input data array is a 3-D object.'
	 PRINT,'In order to mark the circle center with the cursor,'
	 PRINT,'a 2-D frequency slice of the object must be available and'
	 PRINT,'displayed in an IDL window.'
	 menu=['Is an image displayed?','Yes','No']
	 answer=one_line_menu(menu)
	 IF (answer EQ 1) THEN BEGIN
	    win=''
	    READ,'Enter IDL window number: ', win
	    IF (VALIDNUM(win) EQ 0) THEN MESSAGE, 'Invalid window number.'
	    win=FIX(win)
	    DEVICE, WINDOW_STATE=ws
	    IF (ws(win) EQ 0) THEN MESSAGE, 'Window is not available.'
	    WSET, win
	 ENDIF
	 IF (answer EQ 2) THEN BEGIN
	    MESSAGE,/CONT,'An image must be displayed in order to mark points'
	    MESSAGE,/CONT,'with the cursor.  Use command line input for circle'
            MESSAGE,'center and radius if no image is available.'
	 ENDIF
      ENDIF
   ENDIF
ENDIF
;
IF (KEYWORD_SET(annulus )) THEN BEGIN
   ansz=SIZE(annulus)
   IF (N_ELEMENTS(annulus) EQ 1) THEN BEGIN
      IF (KEYWORD_SET(noplot)) THEN BEGIN
        MESSAGE,/CONT,'If NOPLOT keyword, enter annulus info from command line'
	MESSAGE,'eg.,annulus=[center_lon,center_lat,inner_radius,outer_radius]'
      ENDIF
      IF (!D.NAME NE 'X' AND !D.NAME NE 'WIN') THEN BEGIN
         MESSAGE,/CON,'Non-X Window users must enter center and 2 radii from' 
         MESSAGE,/CONT,'command line.  Use the keyword ANNULUS in call, e.g., '
	 MESSAGE,'annulus=[center_lon, center_lat, inner_radius, outer_radius]'
      ENDIF
      IF (ThreeD AND NOT KEYWORD_SET(win)) THEN BEGIN
	 PRINT,' '
	 PRINT,'The input data array is a 3-D object.'
	 PRINT,'In order to mark the annulus center with the cursor,'
	 PRINT,'a 2-D frequency slice of the object must be available and'
	 PRINT,'displayed in an IDL window.'
	 menu=['Is an image displayed?','Yes','No']
	 answer=one_line_menu(menu)
	 IF (answer EQ 1) THEN BEGIN
	    win=''
	    READ,'Enter IDL window number: ', win
	    IF (VALIDNUM(win) EQ 0) THEN MESSAGE, 'Invalid window number.'
	    win=FIX(win)
	    DEVICE, WINDOW_STATE=ws
	    IF (ws(win) EQ 0) THEN MESSAGE, 'Window is not available.'
	    WSET, win
	 ENDIF
	 IF (answer EQ 2) THEN BEGIN 
	    MESSAGE,/CONT,'An image must be displayed in order to mark points '
	    MESSAGE,/CONT,'with the cursor. Use command line input for annulus'
	    MESSSAGE, 'center and two radii if no image is available.'
         ENDIF
      ENDIF
   ENDIF
   IF (N_ELEMENTS(annulus) GT 1) THEN BEGIN
      IF (N_ELEMENTS(annulus) NE 4) THEN MESSAGE, $
      'Invalid annulus vector. Use annulus=[center_lon, center_lat, inner_radius, outer_radius].'
      IF (annulus(2) EQ 0 OR annulus(3) EQ 0) THEN MESSAGE, $
      'Radius cannot be equal to zero.'
      IF (annulus(2) GE annulus(3)) THEN MESSAGE, $
      'Annulus inner radius must be less than its outer radius.'
   ENDIF
ENDIF
;
;
; Activate the proper window.  If the data array is a 3-D object and the 
; polygon/circle/annulus information was specified at the command line,
; there is no need for a window set (it was set before).  In this case, 
; the program skips the WSET line and continues onward. 
;
IF (!D.NAME EQ 'X' OR !D.NAME EQ 'WIN') THEN BEGIN
   IF (NOT KEYWORD_SET(noplot)) THEN BEGIN
    IF (N_ELEMENTS(win) EQ 0) THEN BEGIN
      win=''
      READ,'Enter IDL window number: ', win
      IF (VALIDNUM(win) EQ 0) THEN MESSAGE, 'Invalid window number.'
      win=FIX(win)
    ENDIF
    DEVICE, WINDOW_STATE=ws
    IF (ws(win) EQ 0) THEN MESSAGE, 'Window is not available.'
    WSET, win
   ENDIF
ENDIF
;
;
; Set defaults for some keywords not specified in the call.
;
IF (KEYWORD_SET (badval) EQ 0) THEN badval=0.
IF (KEYWORD_SET (wpos) EQ 0) THEN wpos=[0,0]
;
;
; If the data array is a 2-D object, REFORM it to a 3-D object for the
; execution of the program.
;
IF (not ThreeD) THEN data = REFORM(data,projsz(1),projsz(2),1)
;
;
; Skycubes and faces are assumed to be in an ecliptic coordinate system.
;
IF (3*projsz(1) EQ 4*projsz(2)) THEN BEGIN
   pcoords='E'
   proj='SKYCUBE'
   facenum=-1
ENDIF
IF (projsz(1) EQ projsz(2)) THEN BEGIN
   pcoords='E'
   proj='FACE'
   IF (N_ELEMENTS(facenum) EQ 0) THEN  BEGIN
	facenum=''
	READ, 'Enter Face Number: ', facenum
	IF (VALIDNUM(facenum) EQ 0) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Invalid face number.'
	ENDIF
	facenum=FIX(facenum)
   ENDIF
   IF (facenum lt 0 or facenum gt 5) THEN BEGIN
	IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))	
	MESSAGE, 'Face number must be between 0 and 5.'
   ENDIF
ENDIF
;
;
; 
; Find resolution of skycubes and faces.  Create RES_CC from the RES
; for later CoorConv input. 
;
IF (proj EQ 'SKYCUBE') THEN BEGIN
   res=FIX(2*ALOG((projsz(1)/4))/ALOG(4) + 1)
   res_cc='R'+STRTRIM(STRING(res),1)
ENDIF
IF (proj EQ 'FACE') THEN BEGIN
   res=FIX(2*ALOG((projsz(1)))/ALOG(4) + 1)
   res_cc='R'+STRTRIM(STRING(res),1)
ENDIF
;
;
;
; Set up masking arrays for determining if points are actually on the
; image or are on the background.  Set array to zero where there is no 
; data.
;
IF (!D.NAME EQ 'WIN' OR !D.NAME EQ 'X') THEN BEGIN
   side_len = FIX(SQRT(4^(FLOAT(res)-1)))
   mask=INTARR(!D.X_VSIZE, !D.Y_VSIZE)
   IF (proj EQ 'SKYCUBE') THEN BEGIN
      mask(0+wpos(0):4*side_len-1+wpos(0), 0+wpos(1):3*side_len-1+wpos(1))=1
      mask(0+wpos(0):3*side_len-1+wpos(0), 0+wpos(0):side_len-1+wpos(1))=0
      mask(0+wpos(0):3*side_len-1+wpos(0), $
           2*side_len+wpos(1):3*side_len-1+wpos(1))=0
   ENDIF ELSE mask(0+wpos(0):side_len-1+wpos(0),0+wpos(1):side_len-1+wpos(1))=1
ENDIF
;
;
;
; CIRCLES   CIRCLES   CIRCLES   CIRCLES   CIRCLES   CIRCLES   CIRCLES
; ===================================================================
;
;                     CIRCLES FROM MOUSE INPUT
;
; This segment makes the circle from mouse input for center with the radius
; entered from the keyboard.  The outline of the circle is plotted if 
; desired.
; 
;
IF (KEYWORD_SET(circle)) THEN BEGIN
   IF (N_ELEMENTS(circle) EQ 1) THEN BEGIN
	PRINT, ' '
	PRINT, 'Click on the circle center first, '
	PRINT, 'then enter the radius from keyboard.'
	PRINT, ' '
	CURSOR, x, y, /DEVICE, /DOWN
	IF (NOT KEYWORD_SET(noplot)) THEN PLOTS, x, y, /DEVICE, PSYM=3
	x_center=x-wpos(0)   &   y_center=y-wpos(1)
	IF (mask(x,y) EQ 0) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Chosen center point is not on the image.'
	ENDIF
	radius=''
	READ, 'Enter radius in degrees: ', radius
	IF (VALIDNUM(radius) EQ 0) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Radius must be a number.'
	ENDIF
	radius=FLOAT(radius)
        IF (radius LE 0. OR radius GT 90.) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Radius of circle must be between zero and 90.'
	ENDIF
;
; Create the circle outline using LATCIRC.
;
	IF (proj EQ 'SKYCUBE') THEN BEGIN
		center_pix= XY2PIX(x_center, y_center, res=res)
		pix_outline=LATCIRC(pole=center_pix, angle=radius, res=res)
		PIX2XY,pix_outline,x_outline,y_outline,res=res
	ENDIF ELSE BEGIN
		center_pix= XY2PIX(x_center, y_center, res=res, face=facenum)
		pix_outline=LATCIRC(pole=center_pix,angle=radius,res=res)
		check=PIX2FIJ(pix_outline, res)
		use=WHERE(check(*,0) EQ facenum, count)
		IF (count EQ 0) THEN BEGIN
		   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
		   MESSAGE, 'No part of the circle falls on face.'
		ENDIF
		pix_outline=pix_outline(use)
	   	PIX2XY,pix_outline,x_outline,y_outline,res=res,/face
	ENDELSE
   ENDIF
;
;
;
;                 CIRCLES FROM COMMAND LINE INPUT
;
; Make a circle from coordinatess which were given on the command line.
; First, make sure the circle array is in the proper form, then create
; an array of longitudes and latitudes for the circle using LATCIRC. 
;
   IF (N_ELEMENTS(circle) GT 1) THEN BEGIN
	IF (N_ELEMENTS(circle) NE 3) THEN BEGIN
	   MESSAGE,/CONT,'Circle input array must be in the following form:'
	   MESSAGE,'Circle = [center_lon, center_lat, radius]'
	ENDIF
 	IF (circle(2) LE 0 OR circle(2) GT 90.) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Radius of circle must be between zero and 90.'
	ENDIF
	IF (NOT KEYWORD_SET (in_coords)) THEN BEGIN
	   in_coords_types=['E', 'G', 'Q']
	   menu=['Enter Coordinate System of Circle Center', $
		 'Ecliptic', 'Galactic', 'Equatorial']
	   answer=one_line_menu(menu)
	   in_coords=in_coords_types(answer-1)
	ENDIF ELSE in_coords=STRUPCASE(in_coords)
	ll_cntr=[FLOAT(circle(0)),FLOAT(circle(1))]
	radius= circle(2)
;
; Check that input coordinates for circle center are valid.
;
	IF (ll_cntr(0) GT 360. OR ll_cntr(0) LT -360.) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Decimal longitude must be in range (-360, 360).'
	ENDIF
	IF (ll_cntr(1) GT 90. OR ll_cntr(1) LT -90.) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Decimal latitude must be in range (-90, 90).'
	ENDIF
;
; Convert the circle outline to XY coordinates.
;
	z=WHERE(ll_cntr EQ 0.)
	IF (z(0) NE -1) THEN ll_cntr(z)=1.e-6
	center_pix=COORCONV(ll_cntr,inco=in_coords,infmt='L',$
	outco=res_cc,outfmt='P')
	ll_outline=LATCIRC(pole=ll_cntr,angle=radius,res=res)
	z=WHERE(ll_outline EQ 0.)
	IF (z(0) NE -1) THEN ll_outline(z)=1.e-6
	pix_outline=COORCONV(ll_outline,infmt='L',outfmt='P', $
	inco=in_coords,outco=res_cc)
	IF (proj EQ 'SKYCUBE') THEN $
	PIX2XY,pix_outline,x_outline,y_outline,res=res $
	ELSE BEGIN
	   check=PIX2FIJ(pix_outline,res)
	   use=WHERE(check(*,0) EQ facenum, count)
	   IF (count EQ 0) THEN BEGIN
	      IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	      MESSAGE, 'No part of the circle lies on the face.'
	   ENDIF
	   pix_outline=pix_outline(use)
	   PIX2XY,pix_outline,x_outline,y_outline,res=res,/face
	ENDELSE
   ENDIF
;
;
; Plot circle if keyword NOPLOT has not been set.
;
   IF (NOT KEYWORD_SET(noplot)) THEN BEGIN
      fij_input=PIX2FIJ(pix_outline,res)
      vec_arc=[0,0,0]
      re_proj=0
      color=255
      DRAWCIRC, fij_input, N_ELEMENTS(x_outline), x_outline, y_outline, $
      vec_arc, re_proj, wpos, facenum, color
   ENDIF
;
;
;
; Find pixels which are inside the circle using PXINCIRC.  Of these,
; use only the pixels that contain values not equal to the set BADVAL.
; Find the average using only the inside, good pixels.
;
   pixels=PXINCIRC(cen=center_pix,ang=radius,res=res)
   IF (N_ELEMENTS(pixels) EQ 0) THEN BEGIN
      IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
      MESSAGE, 'No pixels inside the circle.'
   ENDIF
   PRINT, 'Number of pixels is circle is ',N_ELEMENTS(pixels)
   IF (proj NE 'FACE') THEN PIX2XY, pixels, xin, yin, res=res $
   ELSE PIX2XY, pixels, xin, yin, res=res, /face
   in_indices=FLTARR(N_ELEMENTS(xin)*2)
   in_indices=xin+yin*projsz(1)
   regmin=FLTARR(slices)
   regmax=FLTARR(slices)
   average=FLTARR(slices)
   standev=FLTARR(slices)
   FOR i=0, slices-1 DO BEGIN 
      temp=data(*,*,i) 
      good=WHERE(temp(in_indices) NE badval,ngood)
      IF(ngood GT 0) THEN in_indices=in_indices(good) ELSE BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'All pixels in circle have the bad value.'
      ENDELSE
      use=FLTARR(N_ELEMENTS(good),slices)
      use(*,i)=temp(in_indices)
      wts=weights(in_indices)
      regmin(i)=min(use(*,i))
      regmax(i)=max(use(*,i))
      average(i)=TOTAL(use(*,i)*wts) / TOTAL(wts)
      N=N_ELEMENTS(good)
      IF ( N GT 1) THEN BEGIN
	 diffsq = (use(*,i)-average(i))^2
	 standev(i)=SQRT((N*TOTAL(wts*diffsq))/((N-1)*TOTAL(wts)))
      ENDIF ELSE PRINT, 'Only 1 point.  No standard deviation calculated.'
      IF (slices EQ 1) THEN BEGIN
	 PRINT, 'Number of Good pixels is ', ngood
	 PRINT, 'Minimum over circle is ', regmin(i)
	 PRINT, 'Maximum over circle is ', regmax(i)
	 PRINT, 'Average over circle is ', average(i)
	 IF (N GT 1) THEN PRINT,'Standard Deviation over circle is ',standev(i)
      ENDIF ELSE BEGIN
	 index = STRTRIM(STRING(i+1),1)
	 PRINT, 'Number of good pixels in circle on frequency level '$
		+index+' is', N_ELEMENTS(good)
	 PRINT,'Minimum over circle on frequency level '+index+' is',regmin(i)
	 PRINT,'Maximum over circle on frequency level '+index+' is',regmax(i)
	 PRINT,'Average over circle on frequency level '+index+' is',average(i)
	 IF (N GT 1) THEN PRINT, 'Standard Deviation over circle on frequency level '+index+' is', standev(i)
	 PRINT, ' '
      ENDELSE
   ENDFOR
   ;
   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
   RETURN, average
ENDIF
;
;
;
;
; ANNULI   ANNULI   ANNULI   ANNULI   ANNULI   ANNULI   ANNULI
; ============================================================
;
;                     ANNULUS FROM MOUSE INPUT
;
; This segment makes the annulus from mouse input for center with the radii
; entered from the keyboard. 
;
IF (KEYWORD_SET(annulus)) THEN BEGIN
   IF (N_ELEMENTS(annulus) EQ 1) THEN BEGIN
	PRINT, ' '
	PRINT, 'Click on the annulus center first, ' 
	PRINT, 'then enter the radii from keyboard.'
	PRINT, ' '
	CURSOR, x, y, /DEVICE, /DOWN
	IF (NOT KEYWORD_SET(noplot)) THEN PLOTS, x, y, /DEVICE, PSYM=3
	x_center=x-wpos(0)   &   y_center=y-wpos(1)
	IF (mask(x,y) EQ 0) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Chosen center point is not on the image.'
	ENDIF
	in_rad=''	 &   out_rad=''
	READ, 'Enter inner radius in degrees: ', in_rad
	READ, 'Enter outer radius in degrees: ', out_rad
	IF (VALIDNUM(in_rad) EQ 0 OR VALIDNUM(out_rad) EQ 0) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Both radii must be numbers.'
	ENDIF
	in_rad=FLOAT(in_rad)
	out_rad=FLOAT(out_rad)
	IF (in_rad GE out_rad) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Inner radius must be less than outer radius.'
	ENDIF
	IF (in_rad LE 0 OR in_rad GT 90. OR $
            out_rad LE 0 OR out_rad GT 90.) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Radius must be between zero and 90.'
	ENDIF
;
; Create the annulus outline using LATCIRC.
;
	IF (proj EQ 'SKYCUBE' OR proj EQ 'FACE') THEN BEGIN
	   IF (proj EQ 'SKYCUBE') THEN BEGIN
		center_pix= XY2PIX(x_center, y_center, res=res)
		inpix_outln=LATCIRC(pole=center_pix, angle=in_rad, res=res)
		outpix_outln=LATCIRC(pole=center_pix, angle=out_rad, res=res)
		PIX2XY,inpix_outln,inx_outln,iny_outln,res=res
		PIX2XY,outpix_outln,outx_outln,outy_outln,res=res
	   ENDIF ELSE BEGIN
		center_pix= XY2PIX(x_center, y_center, res=res, face=facenum)
		inpix_outln=LATCIRC(pole=center_pix,angle=in_rad,res=res)
		outpix_outln=LATCIRC(pole=center_pix,angle=out_rad,res=res)
		incheck=PIX2FIJ(inpix_outln, res)
		inuse=WHERE(incheck(*,0) EQ facenum, incount)
		IF (incount EQ 0) THEN BEGIN
		   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
		   MESSAGE, 'No part of the inner circle falls on face.'
		ENDIF
		inpix_outln=inpix_outln(inuse)
		outcheck=PIX2FIJ(outpix_outln, res)
		outuse=WHERE(outcheck(*,0) EQ facenum, outcount)
		IF (outcount EQ 0) THEN BEGIN
		   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
		   MESSAGE, 'No part of the outer circle falls on face.'
		ENDIF
		outpix_outln=outpix_outln(outuse)
	   	PIX2XY,inpix_outln,inx_outln,iny_outln,res=res,/face
	   	PIX2XY,outpix_outln,outx_outln,outy_outln,res=res,/face
	   ENDELSE
	ENDIF
   ENDIF
;
;
;
;                 ANNULI FROM COMMAND LINE INPUT
;
; Make an annulus from coordinates which were given on the command line.
; First, make sure the annulus array is in the proper form, then create
; an array of longitudes and latitudes for the circles using LATCIRC. 
;
   IF (N_ELEMENTS(annulus) GT 1) THEN BEGIN
	IF (N_ELEMENTS(annulus) NE 4) THEN BEGIN
	   MESSAGE,/CONT,'Annulus input array must be in the following form:'
	   MESSAGE,'Annulus=[center_lon,center_lat,inner_radius,outer_radius]'
	ENDIF
	IF (NOT KEYWORD_SET (in_coords)) THEN BEGIN
	   in_coords_types=['E', 'G', 'Q']
	   menu=['Coordinate System of Annulus Center', $
		 'Ecliptic', 'Galactic', 'Equatorial']
	   answer=one_line_menu(menu)
	   in_coords=in_coords_types(answer-1)
	ENDIF ELSE in_coords=STRUPCASE(in_coords)
	ll_cntr=[FLOAT(annulus(0)),FLOAT(annulus(1))]
	in_rad= annulus(2)
	out_rad= annulus(3)
	IF (in_rad LE 0 OR in_rad GT 90. OR $
            out_rad LE 0 OR out_rad GT 90.) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Radius must be between zero and 90.'
	ENDIF
;
; Check that input coordinates for circle center are valid.
;
	IF (ll_cntr(0) GT 360. OR ll_cntr(0) LT -360.) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Decimal longitude must be in range (-360, 360).'
	ENDIF
	IF (ll_cntr(1) GT 90. OR ll_cntr(1) LT -90.) THEN BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'Decimal latitude must be in range (-90, 90).'
	ENDIF
;
; Convert the annulus outline arrays to XY coordinates.
;
	IF (proj EQ 'SKYCUBE' OR proj EQ 'FACE') THEN BEGIN
	   z=WHERE(ll_cntr EQ 0.)
	   IF (z(0) NE -1) THEN ll_cntr(z)=1.e-6
	   center_pix=COORCONV(ll_cntr,infmt='l',inco=in_coords,$
	   outco=res_cc,outfmt='p')
	   inll_outln=LATCIRC(pole=ll_cntr,angle=in_rad,res=res)
	   outll_outln=LATCIRC(pole=ll_cntr,angle=out_rad,res=res)
	   z=WHERE(inll_outln EQ 0.)
	   IF (z(0) NE -1) THEN inll_outln(z)=1.e-6
	   inpix_outln=COORCONV(inll_outln,infmt='L',outfmt='P', $
	   inco=in_coords,outco=res_cc)
	   z=WHERE(outll_outln EQ 0.)
	   IF (z(0) NE -1) THEN outll_outln(z)=1.e-6
	   outpix_outln=COORCONV(outll_outln,infmt='L',outfmt='P', $
	   inco=in_coords,outco=res_cc)
	   IF (proj EQ 'SKYCUBE') THEN BEGIN
		PIX2XY,inpix_outln,inx_outln,iny_outln,res=res
		PIX2XY,outpix_outln,outx_outln,outy_outln,res=res
	   ENDIF ELSE BEGIN
		incheck=PIX2FIJ(inpix_outln,res)
		outcheck=PIX2FIJ(outpix_outln,res)
		inuse=WHERE(incheck(*,0) EQ facenum, count)
		outuse=WHERE(outcheck(*,0) EQ facenum, count)
	        IF (inuse(0) EQ -1 OR outuse(0) EQ -1) THEN BEGIN
		   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
		   MESSAGE,'No pixels on face.'
	        ENDIF
		inpix_outln=inpix_outln(inuse)
		outpix_outln=outpix_outln(outuse)
		PIX2XY,inpix_outln,inx_outln,iny_outln,res=res, /face
		PIX2XY,outpix_outln,outx_outln,outy_outln,res=res, /face
	   ENDELSE
	ENDIF
   ENDIF
;
;
;  Plot annulus if keyword NOPLOT has not been set.
;
   IF (NOT KEYWORD_SET(noplot)) THEN BEGIN
      xplotin=inx_outln    &  yplotin=iny_outln
      xplotout=outx_outln  &  yplotout=outy_outln
      fij_in=PIX2FIJ(inpix_outln,res)
      fij_out=PIX2FIJ(outpix_outln,res)
      vec_arc=[0,0,0]
      re_proj=0
      color=255
      DRAWCIRC, fij_in, N_ELEMENTS(xplotin), xplotin, yplotin, $
      vec_arc, re_proj, wpos, facenum, color
      DRAWCIRC, fij_out, N_ELEMENTS(xplotout), xplotout, yplotout, $
      vec_arc, re_proj, wpos, facenum, color
   ENDIF
;
;
;
; Find indices of the data array which are inside both circles.  Of these,
; use only the pixels that contain values not equal to the set BADVAL.
; Eliminate the pixels inside the inner circle and you're left with the
; pixels in the annulus.  The zero_arr variable  creates a dummy array 
; that is the same size as the input data array. The elements contained 
; in the outer circle are set to -100, then the elements of the inner 
; circle are set to -101.  The elements of the dummy array that remain 
; at value -100 are the elements of the annulus. 
;
   in_pixels=PXINCIRC(cen=center_pix,ang=in_rad,res=res)
   out_pixels=PXINCIRC(cen=center_pix,ang=out_rad,res=res)
   IF (proj NE 'FACE') THEN BEGIN
      PIX2XY, in_pixels, xin, yin, res=res
      PIX2XY, out_pixels, xout, yout, res=res
   ENDIF ELSE BEGIN
      PIX2XY, in_pixels, xin, yin, res=res, /face
      PIX2XY, out_pixels, xout, yout, res=res, /face
   ENDELSE
   in_indices=xin+yin*projsz(1)
   out_indices=xout+yout*projsz(1)
   zero_arr=data*0 
   zero_arr(out_indices)= -100
   zero_arr(in_indices) = -101
   ring_indices=where(zero_arr EQ -100)
   IF (ring_indices(0) EQ -1) THEN BEGIN
      IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
      MESSAGE, 'No pixels found in annulus.'   
   ENDIF
   y = ring_indices / projsz(1)
   x = ring_indices-y*projsz(1)
   IF (proj NE 'FACE') THEN pixels = XY2PIX(x, y, res=res) $
   ELSE pixels = XY2PIX(x, y, res=res, face=facenum)
   PRINT, 'Number of pixels in annulus is ', N_ELEMENTS(pixels)
   regmin=FLTARR(slices)
   regmax=FLTARR(slices)
   average=FLTARR(slices)
   standev=FLTARR(slices)
   FOR i=0, slices-1 DO BEGIN
      temp=data(*,*,i) 
      good=WHERE(temp(ring_indices) NE badval,ngood)
      IF (ngood GT 0) THEN ring_indices=ring_indices(good) $
      ELSE BEGIN
	   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
	   MESSAGE, 'All pixels in annulus have badval value.'
      ENDELSE
      wts=weights(ring_indices)
      use=FLTARR(N_ELEMENTS(ring_indices),slices)
      use(*,i)=temp(ring_indices)
      regmin(i)=min(use(*,i))
      regmax(i)=max(use(*,i))
      average(i)=TOTAL(use(*,i)*wts) / TOTAL(wts)
      N=N_ELEMENTS(ring_indices)
      IF (N GT 1) THEN BEGIN
	 diffsq = (use(*,i)-average(i))^2
	 standev(i)=SQRT((N*TOTAL(wts*diffsq))/((N-1)*TOTAL(wts)))
      ENDIF ELSE PRINT, 'Only 1 point.  No standard deviation calculated.'
      IF (slices EQ 1) THEN BEGIN
	 PRINT, 'Number of Good pixels is ', ngood
	 PRINT, 'Minimum over annulus is ', regmin(i)
	 PRINT, 'Maximum over annulus is ', regmax(i)
	 PRINT, 'Average over annulus is ', average(i)
	 IF (N GT 1) THEN PRINT, 'Standard Deviation over annulus is ', $
	 standev(i)
      ENDIF ELSE BEGIN
	 index = STRTRIM(STRING(i+1),1)
	 PRINT, 'Number of good pixels in annulus on frequency level '$
		+index+' is', N_ELEMENTS(ring_indices)
	 PRINT,'Minimum over annulus on frequency level '+index+' is',regmin(i)
	 PRINT,'Maximum over annulus on frequency level '+index+' is',regmax(i)
   
	 PRINT,'Average over annulus on frequency level '+index+' is',average(i)
	 IF (N GT 1) THEN PRINT, 'Standard Deviation over annulus on frequency level ' +index+' is', standev(i)
	 PRINT, ' '
      ENDELSE
   ENDFOR
   ;
   ;
   IF (NOT ThreeD) THEN data=REFORM(data,projsz(1),projsz(2))
   RETURN, average
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


