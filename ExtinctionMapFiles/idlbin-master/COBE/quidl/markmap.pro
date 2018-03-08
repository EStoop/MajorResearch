   PRO markmap, image, proj=proj, pcoords=pcoords, file=file, $
   fcoords=fcoords, facenum=facenum, win=win, entry=entry, $
   color=color, psym=psym, charsize=charsize, $
   charthick=charthick, labelpos=labelpos, llx=llx, lly=lly
;
;+NAME/ONE-LINE DESCRIPTION:
;  MARKMAP places labels on images given coordinates in any system
;
; DESCRIPTION: 
;  MARKMAP places labels on an image given a file of coordinates
;  and corresponding labels.  The image may be a skycube, face, or
;  a reprojection in any coordinate system.  The coordinates given
;  in the file may also be in any coordinate system.  Coordinates
;  and labels may also be entered from the keyboard or positions
;  may be marked with the mouse.
;
; CALLING SEQUENCE:
;  markmap, image, [proj=proj], [pcoords=pcoords], [file=file], $
;  [fcoords=fcoords], [facenum=facenum], [win=win], [entry=entry], $
;  [color=color], [psym=psym], [charsize=charsize], $
;  [charthick=charthick], [labelpos=labelpos], [llx=llx], [lly=lly]
;
; ARGUMENTS: (I=input, O=output, []=optional)
;  image      I   2-D data array
;  [proj]     I   Projection type.  If data is a skycube or face,
;                 proj need not be specified.  Reprojections
;                 must be 'A', 'M', or 'S' for Aitoff, Mollweide,
;                 or Global Sinusoidal.  Need not be in capitals.
;  [pcoords]  I   Projection coordinate system.  Must be
;                 'E', 'G', or 'Q' for Ecliptic, Gal, or Equatorial.
;                 If the projection type is either 'skycube' or 
;                 'face', pcoords is set to 'E'.  Need not be in caps.
;  [file]     I   Name of the file holding the labels and coordinates.
;                 File needs to be of one of these two forms:
;                 HR MIN SEC DEG DEC_MIN DEC_SEC LABEL, i.e., 
;                       7 columns for coords in RA Dec Sexigesimal or
;                 LONG LAT LABEL, i.e.,
;                       3 columns for coords in long/lat format.
;                 Labels cannot have spaces in them if they are read
;                 in from a file;  labels entered from the keyboard
;                 may contain spaces.  Labels 
;                 Accepted values for label coordinates are:
;                       All latitudes: (-90, +90)
;                       Decimal longitudes:  (-360, +360)
;                       RA hr:  ( 0, 24)
;                 If the program finds coordinates that do not match
;                 one of these criterion, this line will not be 
;                 plotted;  the others, however, will be displayed.
;  [fcoords]  I   File coordinate system. Must be
;                 'E', 'G', or 'Q' for Ecliptic, Gal, or Equatorial.
;                 Equatorial coords can be specified as "QD" for
;                 Decimal coords or "QS" for Sexigesimal coords.
;                 Need not be in caps.  
;  [facenum]  I   If data array is a face, you may specify the 
;                 facenum in the call.
;  [win]      I   Window number where the display lives. 
;  [entry]    I   Specifies coordinate and label entry mode.
;                 'file' or 'keyboard' gets the input from a
;                 file or keyboard entry;  'cursor' allows you to
;                 type the label in and then click the spot on the
;                 the image where you want it to appear.  If the
;                 keyword FILE is set, you DO NOT have to specify
;                 that entry='file'.
;  [color]    I   Color to plot labels with. 0-255
;  [psym]     I   IDL plot symbol.  1=plus sign (default),
;                 2=asterisk, 3=period, 4=diamond, 5=triangle,
;                 6=square, 7=X
;  [charsize]  I  Size of the printed labels and symbols. Default=1
;  [charthick] I  Thickness of the printed labels.  Default=1
;  [labelpos]  I  The position of the printed label relative to 
;                 the plot symbol.  Possible choices are 'T' for top,
;                 'B' for bottom, 'L'/'R' for left and right. The 
;                  default is to place the label above the symbol (T).
;  [llx]       I  X coordinate of the lower left corner of the
;                 image.  i.e., if the image has a 32 screen
;                 pixel border on each side, llx=32. Def=0
;  [lly]       I  Y coordinate of the lower left corner of the
;                 image.  Def=0
;			
;	
; EXAMPLE:
;  If you have an Aitoff/Galactic reprojection in a data array called 
;  DATA and a file with ecliptic coords, and you would like the labels
;  placed to the right of the points, the call would look like
;
;  UIDL>  markmap, data, proj='a', pcoords='g', file='myfile' $
;  UIDL>  fcoords='e', labelpos='r'
;
;#
; SUBROUTINES CALLED:
;  umenu, readcol, ten, coorconv, pix2xy, uv2proj, pix2fij
;
; COMMON BLOCKS:  None
;
; LIBRARY CALLS:  None
;
; WARNINGS:  None
;
; PROGRAMMING NOTES:
;  There is a minor "kludge" in this program.  COORCONV doesn't like
;  an ecliptic longitude=0 when the projection is also in ecliptic
;  coordinates.  To get around this, MARKMAP looks for this situation
;  and will offset the latitude by 1.e-6 so that COORCONV will not
;  display an error message.
;
; MODIFICATION HISTORY:
;  Written by Celine Groden, USRA    Feb 1993    SPR # 10498
;  Modified by C. Groden, USRA       May 1993    SPR # 10916
;
;.TITLE
;Routine MARKMAP
;
;-
;
; Establish error condition to return to calling routine in case of
; trouble.
;
ON_ERROR, 2
;
; Check that a proper input data array is being used.
;
projsz=SIZE(image)
IF (projsz(0) NE 2) THEN MESSAGE, 'Input array must be 2-D'
;
IF (projsz(1) NE projsz(2) AND $
    projsz(1) NE 2*projsz(2) AND $
    3*projsz(1) NE 4*projsz(2)) THEN MESSAGE, $
    'Input array must be a skycube, face, or a reprojection'
;
;
; Set defaults for keywords not specified in the call.
;
IF (KEYWORD_SET(psym) EQ 0) THEN psym=1
IF (KEYWORD_SET(color) EQ 0) THEN color = 255
IF (KEYWORD_SET(charsize) EQ 0) THEN charsize=1
IF (KEYWORD_SET(charthick) EQ 0) THEN charthick=1
IF (KEYWORD_SET(labelpos) EQ 0) THEN labelpos='T'
IF (KEYWORD_SET (llx) EQ 0) THEN llx=0
IF (KEYWORD_SET (lly) EQ 0) THEN lly=0
IF (3*projsz(1) EQ 4*projsz(2)) THEN BEGIN
	pcoords='E'
	proj='SKYCUBE'
ENDIF
IF (projsz(1) EQ projsz(2)) THEN BEGIN
	pcoords='E'
	proj='FACE'
	IF (KEYWORD_SET(facenum) EQ 0) THEN  BEGIN
		facenum=' '
		READ, 'Enter Face Number: ', facenum
	ENDIF
	IF (facenum lt 0 or facenum gt 5) THEN MESSAGE, $
		'Face number must be between 0 and 5.'
ENDIF
;
;
skip=0  ; SKIP is turned on if the menu to differentiate between
	; equatorial sexigesimal and decimal is not needed (i.e.,
	; if equatorial coordinates are read in from the keyboard)
;
;
; Activate proper window.
;
IF (N_ELEMENTS(win) EQ 0) THEN BEGIN
    win=''
    READ,'Enter IDL window number: ', win
    IF (VALIDNUM(win) EQ 0) THEN MESSAGE, 'Invalid window number.'
ENDIF
win=FIX(win)
DEVICE, WINDOW_STATE=ws
IF (ws(win) EQ 0) THEN BEGIN
   MESSAGE, 'Window is not available.',/CONTINUE
   RETURN
ENDIF
WSET, win
;
;
; If the plotting color is changed, save the old !p.color so it can
; be restored at the end of the program.
;
old_color=!p.color
IF (KEYWORD_SET (color)) THEN !p.color=color
;
;
; Use keyboard entry to specify the name of the reprojection
; if it was not supplied in the call and image is not a skycube or face.
;
IF (KEYWORD_SET(proj) EQ 0) THEN BEGIN
	proj_types=['A', 'S', 'M']
	menu=UMENU(['Enter Projection Type', $
		'Aitoff', $
		'Global Sinusoidal', $
		'Mollweide'], $
		 title=0, init=1)
	proj=proj_types(menu-1)
ENDIF ELSE proj=STRUPCASE(proj)
;
;
; Determine the coord system of the projection.  If it was specified in
; call, change the string value to all caps.
;
IF (KEYWORD_SET(pcoords) EQ 0) THEN BEGIN
	pcoords_types=['E', 'G', 'Q']
	menu=UMENU(['Enter Projection Coordinate System', $
		'Ecliptic', $
		'Galactic', $
		'Equatorial'], $
		 title=0, init=1)
	pcoords=pcoords_types(menu-1)
ENDIF ELSE pcoords=STRUPCASE(pcoords)
;
;
; Enter coordinates and labels from the keyboard if a file is not 
; available.  Flag unacceptable coordinates (see 'file' in
; prologue).  Equatorial coordinates need to be entered in 
; sexigesimal format.
;
IF (KEYWORD_SET(file) EQ 0) THEN BEGIN
    file=''
    IF (KEYWORD_SET(entry)) THEN entry=STRUPCASE(entry)
    IF (KEYWORD_SET(entry) EQ 0) THEN BEGIN
	entry_mode=['CURSOR', 'KEYBOARD']
	menu=UMENU(['Select Entry Mode', $
		'Cursor', $
		'File or Keyboard'], $
		 title=0, init=1)
	   entry=entry_mode(menu-1)
   ENDIF
   IF (entry EQ 'KEYBOARD' OR entry EQ 'FILE') THEN BEGIN
	READ, 'Enter filename (enter <CR> for keyboard entry): ', file
	IF (file EQ '') THEN BEGIN
		READ, 'Number of entries: ', numline
		IF (numline le 0) THEN MESSAGE, $
			'Number of entries must be greater than zero.'
		numline=FIX(numline) & fcoords=' '
		READ, 'Coordinate system of entries (E,G, or Q): ', fcoords
		fcoords=STRUPCASE(fcoords)
		IF (fcoords NE 'E' AND fcoords NE 'G' AND fcoords NE 'Q') $
			THEN MESSAGE, 'Invalid coordinate system'
		tmpcoords=FLTARR(2,numline)
		label=STRARR(numline)
		IF (fcoords EQ 'E' OR fcoords EQ 'G') THEN BEGIN
			FOR i=0, numline-1 DO BEGIN
				number= STRTRIM(STRING(i+1),1)
				PRINT, ' '
				PRINT,'Separate numbers with a space or comma.'
				READ, 'Lon,Lat of Point '+number+': ',lon,lat
				IF (lon gt 360. OR lon lt -360.) THEN MESSAGE, $
				   'Longitude must be in degree range (-360, +360)'
				IF (lat gt 90. OR lat lt -90.) THEN MESSAGE, $
				   'Latitudes must be in degree range (-90, +90)'
				lab=''
				READ, 'Label for Point '+number+': ', lab
				tmpcoords(0,i) =lon
				tmpcoords(1,i) =lat
				label(i)=STRING(lab)
			ENDFOR
		ENDIF
		IF (fcoords EQ 'Q') THEN BEGIN
			PRINT, ' '
			PRINT, 'Enter points in sexigesimal format.'
			PRINT,'Separate numbers with a space or comma.'
			PRINT, 'Six numbers must be entered.'
			PRINT, 'RA in (hr, min, sec), Dec in (deg, min, sec)'
			FOR i=0, numline-1 DO BEGIN
				number= STRTRIM(STRING(i+1),1)
				READ, 'RA of Point '+number+': ',hr,mn,sec
				IF (hr gt 24 OR hr lt 0) THEN MESSAGE, $
				   'Right Ascension must be in hour range (0, 24)'
				READ, 'Dec of Point '+number+': ',deg,dmin,dsec
				lab=''
				READ, 'Label for Point '+number+': ', lab
				tmpcoords(0,i)= 15*TEN(hr,mn,sec)
				tmpcoords(1,i)= TEN(deg,dmin,dsec)
				IF (tmpcoords(1,i) gt 90. OR $
				    tmpcoords(1,i) lt -90.) THEN MESSAGE,$
				    'Declination must be in degree range (-90, +90)'
				label(i)=STRING(lab)
			ENDFOR
		skip=1
		ENDIF
	ENDIF
   ENDIF 
;
;
; Mouse position entry.  CASE segment defines the label positioning.
; Then a label is read in from the keyboard, followed by a mouse click
; for positioning.  Hit third mouse button to quit.
;
   IF (entry EQ 'CURSOR') THEN BEGIN
	label=''
	multfact= 5.*charsize
	IF (KEYWORD_SET(labelpos)) THEN BEGIN
	  CASE (STRUPCASE(labelpos)) OF
		'T':	BEGIN
			xmove= 0  &  ymove= multfact + 3
			alignment=0.5
			END
		'B':	BEGIN
			xmove= 0  &  ymove= - 2*multfact - 3
			alignment=0.5
			END
		'L':	BEGIN
			xmove = 0  &  ymove= - multfact/2
			alignment = 1.0
			END
		'R':	BEGIN
			xmove = multfact  &  ymove= - multfact/2
			alignment = 0.0
			END
	  ENDCASE
	ENDIF 
	PRINT, 'Click on points with <MB1>. Click the LAST POINT with <MB3>.'
	!ERR=0
	WHILE (!ERR NE 4) DO BEGIN
		READ, 'Enter Label: ', label
		IF (STRUPCASE(labelpos) EQ 'L') THEN label=label+' '
		PRINT, 'Click corresponding position on image with <MB1>.'
		PRINT, 'Click <MB3> if this is the LAST point.'
		PRINT, ' '
		CURSOR, x_coord, y_coord, /DEVICE, /DOWN
		PLOTS, x_coord, y_coord, /DEVICE, PSYM=psym, $
		THICK=charthick, SYMSIZE=charsize
		XYOUTS, x_coord+xmove, y_coord+ymove, label, /DEVICE, $
		ALIGNMENT=alignment, CHARSIZE=charsize, CHARTHICK=charthick
	ENDWHILE
	!err=0
	RETURN
   ENDIF
ENDIF ELSE file=file
;
;  
; If a file is being used, determine the coord system of the coords 
; in the file.  If it was specified in call, change the string value to 
; all caps.
;
IF (N_ELEMENTS(fcoords) EQ 0) THEN BEGIN
	fcoords_types=['E', 'G', 'QS', 'QD']
	menu=UMENU(['File Coordinates', $
			'Ecliptic', $
			'Galactic', $
			'Equatorial-Sexigesimal Format', $
			'Equatorial-Decimal Format'], $
			 title=0, init=1)
	fcoords=fcoords_types(menu-1)
ENDIF 
fcoords=STRUPCASE(fcoords)
;
;
; If the file coords are Equatorial, determine their format, either
; decimal or sexigesimal.
;
IF (fcoords EQ 'Q' AND skip NE 1) THEN BEGIN
	fcoords_types_2=['QS', 'QD']
	menu=UMENU(['File Coordinates', $
			'Equatorial-Sexigesimal Format', $
			'Equatorial-Decimal Format'], $
			 title=0, init=1)
	fcoords=fcoords_types_2(menu-1)
ENDIF ELSE fcoords=STRUPCASE(fcoords)
;
;
; Read in the coordinates from the specified file.  Coordinates that
; in equatorial sexigesimal need to be in six columns of hour, min, sec, 
; deg, degmin, degsec and a seventh column of string labels; all other 
; coordinate systems are assumed to be in two columns of longitude
; and latitude with a third column of string labels.  Checks are done
; to check that RA hr is between (0,24), decimal longitudes are
; between (-360,360), and latitudes are between (-90, 90).  If a bad
; line is detected, this line is excluded from the array of coordinates
; to be plotted.  If ALL lines are detected as bad, the program does
; not plot anything and the routine ends.
;
IF (file NE '' ) THEN BEGIN
	IF (fcoords EQ 'QS') THEN BEGIN
		READCOL, file, FORMAT='f,f,f,f,f,f,a', $
		hr, mn, sec, deg, dmin, dsec, label
		filesz=SIZE(hr)
		numline=filesz(1)
		radec=FLTARR(2,numline)
		FOR i=0,numline-1 DO BEGIN
			radec(0,i)= 15*TEN(hr(i),mn(i),sec(i))
			radec(1,i)= TEN(deg(i),dmin(i),dsec(i))
		ENDFOR
; Check RA.
		toobig=where(hr gt 24. OR hr lt 0., countbig)
		keep=where(hr le 24. AND hr ge 0., countkeep)
		IF (countbig EQ numline) THEN BEGIN
			PRINT, ' '
			PRINT, 'Right Ascension must be in hour range (0, 24)'
			PRINT, 'No valid coordinates found.'
			RETURN
		ENDIF
		IF (countbig GT 0) THEN BEGIN
			PRINT, ' '
			PRINT, 'Right Ascension must be in hour range (0, 24)'
			IF (countbig EQ 1) THEN BEGIN
			   PRINT, 'Markmap will not use the following line:'
			   PRINT, hr(toobig), mn(toobig), sec(toobig), $
			   deg(toobig), dmin(toobig), dsec(toobig), $
			   label(toobig), FORMAT='(6f8.2,a20)'
			ENDIF ELSE BEGIN
			   PRINT, 'Markmap will not use the following lines:'
			   FOR a=0, countbig-1 DO BEGIN
			      i=toobig(a)
			      PRINT, hr(i), mn(i), sec(i), deg(i), dmin(i), $
			      dsec(i), label(i), FORMAT='(6f8.2,a20)'
			   ENDFOR
			ENDELSE
			;  Count the number of accepted lines
			numline=numline-countbig
		ENDIF
		IF (keep(0) NE -1) THEN BEGIN
			tmpcoords=FLTARR(2,numline)
			tmpcoords(0,*)=radec(0,keep)
			tmpcoords(1,*)=radec(1,keep)
			label=label(keep)
		ENDIF ELSE BEGIN
			tmpcoords=FLTARR(2,numline)
			tmpcoords(0,*)=radec(0,*)
			tmpcoords(1,*)=radec(1,*)
			label=label(*)
		ENDELSE
	ENDIF ELSE BEGIN
;
; Read in three column files.
;
		READCOL, file, FORMAT='f,f,a', dlon, dlat, label
		filesz=SIZE(dlon)
		numline=filesz(1)
		lonlat=FLTARR(2,numline)
		lonlat(0,*)=dlon & lonlat(1,*)=dlat
		keep=where(dlon le 360. AND dlon ge -360.    $
		       AND dlat le  90. AND dlat ge  -90.)
; Check longitudes.
		lon_out=where(dlon gt 360. OR dlon lt -360., countout)
		lon_in=where(dlon le 360. OR dlon ge -360., countin)
		IF (countout EQ numline) THEN BEGIN
			PRINT, ' '
			IF (fcoords EQ 'E' or fcoords EQ 'G') THEN PRINT, $
			'Longitudes must be in range (-360, 360)'
			IF (fcoords EQ 'QD') THEN PRINT, $
			'RA must be in range (-360, 360)'
			PRINT, 'No valid coordinates found.'
			RETURN
		ENDIF
		IF (countout GT 0) THEN BEGIN
			PRINT, ' '
			IF (fcoords EQ 'E' or fcoords EQ 'G') THEN PRINT, $
			'Longitudes must be in range (-360, 360)'
			IF (fcoords EQ 'QD') THEN PRINT, $
			'Decimal RA must be in range (-360, 360)'
			IF (countout EQ 1) THEN BEGIN
			   PRINT, 'Markmap will not use the following line:'
			   PRINT, dlon(lon_out), dlat(lon_out), $
			   label(lon_out), FORMAT='(2f10.2,a20)'
			ENDIF ELSE BEGIN
			   PRINT, 'Markmap will not use the following lines:'
			   FOR a=0, countout-1 DO BEGIN
			      i=lon_out(a)
			      PRINT, dlon(i), dlat(i), label(i), $
			      FORMAT='(2f10.2,a20)'
			   ENDFOR
			ENDELSE
			;  Count the number of accepted lines
		ENDIF
; Check latitudes.
		lat_out=where(dlat gt 90. OR dlat lt -90., countout2)
		lat_in=where(dlat le 90. OR dlat ge -90., countin2)
		IF (countout2 EQ numline) THEN BEGIN
			PRINT, ' '
			IF (fcoords EQ 'E' or fcoords EQ 'G') THEN PRINT, $
			'Latitudes must be in range (-90, 90)'
			IF (fcoords EQ 'QD' OR fcoords EQ 'QS') THEN PRINT, $
			'Declinations must be in range (-90, 90)'
			PRINT, 'No valid coordinates found.'
			RETURN
		ENDIF
		IF (countout2 GT 0) THEN BEGIN
			PRINT, ' '
			IF (fcoords EQ 'E' or fcoords EQ 'G') THEN PRINT, $
			'Latitudes must be in range (-90, 90)'
			IF (fcoords EQ 'QD') THEN PRINT, $
			'Declinations must be in range (-90, 90)'
			IF (countout2 EQ 1) THEN BEGIN
			   PRINT, 'Markmap will not use the following line:'
			   PRINT, dlon(lat_out), dlat(lat_out), $
			   label(lat_out), FORMAT='(2f10.2,a20)'
			ENDIF ELSE BEGIN
			   PRINT, 'Markmap will not use the following lines:'
			   FOR a=0, countout2-1 DO BEGIN
			      i=lat_out(a)
			      PRINT, dlon(i), dlat(i), label(i), $
			      FORMAT='(2f10.2,a20)'
			   ENDFOR
			ENDELSE
			;  Count the number of accepted lines
		ENDIF
		numline=numline-countout-countout2
		tmpcoords=FLTARR(2,numline)
		IF (keep(0) NE -1) THEN BEGIN
			tmpcoords(0,*)=lonlat(0,keep)
			tmpcoords(1,*)=lonlat(1,keep)
			label=label(keep)
		ENDIF ELSE BEGIN
			tmpcoords(0,*)=lonlat(0,*)
			tmpcoords(1,*)=lonlat(1,*)
			label=label(*)
		ENDELSE
	ENDELSE
ENDIF
;
;
IF (fcoords EQ 'QS' OR fcoords EQ 'QD') THEN fcoords='Q'
;
;
;
; Check for the "Ecliptic longitude=0 feature" in COORCONV.  If 
; projection is in ecliptic coords and an input coordinate is ecliptic
; longitude=0, offset the longitude by a tiny bit so that COORCONV 
; doesn't give an error.
;
IF (pcoords EQ 'E' AND fcoords EQ 'E') THEN BEGIN
	lonzero=where(tmpcoords(0,*) EQ 0, countlon)
	IF (countlon NE 0) THEN tmpcoords(0,lonzero)=1.e-6
ENDIF
;
;
; Determine screen pixel adjustments for label positioning and place
; these shift values in Xmove and Ymove.  Top and Bottom center the
; text and vertically move them above or below the symbol.  Left/Right 
; move the text vertically by half the height of the symbol (so label
; and symbol appear to be on the same line) and right/left justify,
; respectively, the label.  Left positioning, which is right justified,
; requires a space to be added to each of the labels so a space is 
; left between the label and symbol.
;
xmove=fltarr(numline)
ymove=fltarr(numline)
multfact= 5.*charsize
;
CASE STRUPCASE(labelpos) OF
	'T':	BEGIN
		xmove(*)=  0.
		ymove(*)= + multfact+3
		alignment=0.5    ;  set XYOUTS keyword for centered text
		END
	'B':	BEGIN
		xmove(*)=  0.
		ymove(*)= -2*multfact-3.
		alignment=0.5    ;  set XYOUTS keyword for centered text
		END
	'L':	BEGIN
		ymove(*)= - multfact/2
		label=label + ' '   ; add spaces at the end of the label
				    ; array to create a space
		alignment= 1.0   ;  set XYOUTS keyword for rt justified
		END
	'R':	BEGIN
		ymove(*)= - multfact/2
		xmove(*)=   multfact
		alignment = 0.0  ;  set XYOUTS keyword for left justified
		END
ENDCASE
;
;
; If array is a skycube, determine plot coordinates in xy coords
; and plot labels.  Skycube resolution must be supplied if it was
; not set in the call.
;
IF (3*projsz(1) EQ 4*projsz(2)) THEN BEGIN
	res=FIX(2*alog((projsz(1)/4))/alog(4) + 1)
	outco='R'+STRTRIM(STRING(res),1)
	tmpcoords=TRANSPOSE(tmpcoords)
	pixels=COORCONV(tmpcoords,infmt='L',outfmt='P', $
	inco=fcoords, outco=outco)
	PIX2XY,pixels,xout,yout,res=res
	XRANGE=[0,projsz(1)]
	YRANGE=[0,projsz(2)]
	FOR j=0, numline-1 DO BEGIN
		PLOTS,xout(j)+llx,yout(j)+lly, $
		THICK=charthick, SYMSIZE=charsize, PSYM=psym,/DEVICE
		XYOUTS,xout(j)+llx+xmove(j),yout(j)+lly+ymove(j),$
		label(j), CHARSIZE=charsize, CHARTHICK=charthick,$
		ALIGNMENT=alignment, /DEVICE
	ENDFOR
ENDIF
;
;
; If array is a face, calculate resolution and determine plot 
; coordinates in xy coords and plot labels.
;
IF (projsz(1) EQ projsz(2)) THEN BEGIN
	res=FIX(2*alog(projsz(1))/alog(4) + 1)
	outco='R'+STRTRIM(STRING(res),1)
	tmpcoords=TRANSPOSE(tmpcoords)
	pixels=COORCONV(tmpcoords,infmt='L',outfmt='P', $
	inco=fcoords, outco=outco)
;
; Plot only the labels that fall in this face.
;
	facecheck=PIX2FIJ(pixels,res)
	use=WHERE(facecheck(*,0) EQ facenum, count_pts)
	IF (count_pts EQ 0) THEN MESSAGE, 'No labels fall on given face.'
	count_pts=STRTRIM(STRING(count_pts),1)
	PRINT, ''+count_pts+' label(s) fall on the given face.'
	label=label(use)	
	sz_leftovers=SIZE(use)
	numline=sz_leftovers(1)
	PIX2XY,pixels(use),xout,yout,res=res,/face
	XRANGE=[0,projsz(1)]
	YRANGE=[0,projsz(2)]
	FOR j=0, numline-1 DO BEGIN
		PLOTS,xout(j)+llx,yout(j)+lly, $
		thick=charthick, symsize=charsize, psym=psym,/DEVICE
		XYOUTS,xout(j)+llx+xmove(j),yout(j)+lly+ymove(j),$
		label(j), charsize=charsize, charthick=charthick,$
		ALIGNMENT=alignment, /DEVICE
	ENDFOR
ENDIF
;
;
; If image is a reprojection, determine the plot coordinates in 
; xy screen coords and plot labels.
;
IF (projsz(1) EQ 2*projsz(2)) THEN BEGIN
	tmpcoords=TRANSPOSE(tmpcoords)
	uvects=COORCONV(tmpcoords,infmt='L',outfmt='U', $
	inco=fcoords, outco=pcoords)
	xy=UV2PROJ(uvects, proj, projsz)
	XRANGE=[0,projsz(1)]
	YRANGE=[0,projsz(2)]
	FOR j=0, numline-1 DO BEGIN
		PLOTS, xy(j,0)+llx, xy(j,1)+lly, thick=charthick, $
		symsize=charsize, psym=psym,/DEVICE
		XYOUTS, xy(j,0)+llx+xmove(j), xy(j,1)+lly+ymove(j),$
		label(j), charsize=charsize, charthick=charthick,$
		ALIGNMENT=alignment, /DEVICE
	ENDFOR
ENDIF
;
;
IF (KEYWORD_SET(color)) THEN !p.color = old_color
;
RETURN
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


