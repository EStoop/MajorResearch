	FUNCTION LATCUT, lat0, lat1, data=data, res=res, $
			pixels=pixels, coords=coords, noabs=noabs
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     LATCUT selects data from a map in a specified latitude range.
;
;DESCRIPTION: 
;     LATCUT takes as input a pair of latitudes that it uses to
;     define a band (or 2 bands) on the sky.  If a data array is
;     supplied it returns a list of data from that array that
;     fall in the specified latitude range.  If, instead of a data 
;     array a resolution number is supplied, LATCUT will
;     return a list of pixel numbers in the specified latitude
;     range.  The function does not care what in order the latitudes
;     are specified; it sorts out which is the larger one.  The
;     default condition is that absolute value of the latitudes is
;     used so that data from both the northern and southern latitude
;     bands is returned; the /NOABS keyword turns this effect off so
;     that the user can specify a range, e.g., [-5, +15], or a region
;     around one pole.  
;    
;CALLING SEQUENCE:  
;     DATALIST = LATCUT (lat0, lat1, [data=data], [res=res], $
;			[/noabs], [coords=coords],[pixels=pixels]) 
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     DATALIST   O   arr       List of latitude-selected data (if
;                                 "data" keyword used, or pixel list
;                                 (if "res" keyword used).
;     LAT0       I   flt       Lower or upper latitude limit
;     LAT1       I    "        Other latitude boundary
;    [DATA]     [I}  array     Skymap, either in sixpack format or
;                                 unfolded T.
;    [RES}      [I}  int       Resolution of skycube; required only
;                                 if DATA array is not supplied.
;    [/NOABS]   [I]  keywd     If specified, only the specified
;                                 latitude range is used; otherwise
;                                 absolute value is taken so that both
;                                 northern and southern bands are
;                                 returned.
;    [COORDS]   [I]  string    Coordinate system in which to interpret
;                                 latitudes.  Only allowed values are:
;                                 'G':  galactic  (Default)
;                                 'E':  ecliptic
;                                 'Q':  equatorial
;    [PIXELS]   [O]  int arr   Optional output latitude-selected
;                                 pixel list.  Identical to function
;                                 output if DATA array not specified.
;
;WARNINGS:
;     1.  One and only one of DATA or RES must be specified.
;
;EXAMPLES:
;     For an array MAP from which we wish to extract all data within
;     |b|<15, we could say
;
;             output = LATCUT(0,15,data=map)             or
;             output = LATCUT(-15,15,data=map,/noabs)
;
;     We could also have added pixels=<variable name> to the call
;     to get the list of pixel numbers corresponding to the data. 
;     If we wish to know which pixels fall within 10 deg of the
;     southern ecliptic pole on a DIRBE-resolution skycube, we say
;
;             pixlist = LATCUT(-80,-90,/noabs,res=9,coords='e')
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Standard coordinate conversion calls.  Does assorted error
;    checking to make sure data is a skycube or sixpack, and does a
;    little footwork to sort out which latitude is the northern
;    boundary (depending on absolute value.)
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     Calls COORCONV and PIX2DAT from the UIDL library.
;  
;MODIFICATION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  18 Aug 1993
;.TITLE
; Function LATCUT
;-
	ON_ERROR, 2
;
;  First figure out which is the lower latitude and which is the
;  higher; this depends on whether we are taking the absolute value.
;
	latpair = [lat0, lat1]
	IF (KEYWORD_SET(noabs)) THEN BEGIN
	   latmin = MIN(latpair)
	   latmax = MAX(latpair)
	ENDIF ELSE BEGIN
	   latmin = MIN(ABS(latpair))
	   latmax = MAX(ABS(latpair))
	ENDELSE
;
;  Check that latitudes and coordinate system are physically 
;  reasonable.  Default to galactic coordinates if unspecified.
;
	IF (latmax EQ latmin OR ABS(latmax GT 90) OR $
		ABS(latmin) GT 90) THEN MESSAGE, $
		'You must supply sensible min and max latitudes.'
;
	IF (NOT KEYWORD_SET(coords)) THEN coords = '?'
	coordsys = STRUPCASE (coords)
	IF (coordsys NE 'Q' AND $
	    coordsys NE 'G' AND $
	    coordsys NE 'E') THEN BEGIN
		PRINT, "Galactic coordinates assumed."
		coordsys = 'G'
	ENDIF
;
;  If a data array is present then figure out the cube resolution by
;  looking at the size of a single face.  If no data array present,
;  get the resolution from the "res" keyword.
;
	IF (N_ELEMENTS(data) GT 0) THEN BEGIN
	   dims = SIZE(data)
	   IF (2*dims(1) EQ 3*dims(2)) THEN sixpak = 1 ELSE $
	      IF (3*dims(1) EQ 4*dims(2)) THEN sixpak = 0 ELSE $
	         MESSAGE,'Data must be a T or sixpack.'
	   facesize = dims(2)/(3-sixpak)
	   res = FIX(3.3*ALOG10(facesize) + 1.5)
	ENDIF ELSE IF (NOT KEYWORD_SET(res)) THEN MESSAGE, $
 	   'You must supply either a data array or a cube resolution.'
	IF (res NE 6 AND res NE 9) THEN MESSAGE, $
		'Skycube resolution must be 6 or 9.'
;
;  Generate an array of all pixel numbers, then convert them all
;  to longitudes and latitudes in the desired coordinate system.
;  Use the latitude list to filter out the ones in the desired
;  latitude range. If the /NOABS keyword is set then grab only
;  the band [latmin --> latmax]; otherwise, get the corresponding
;  band on the other side of the equator as well.
;
	allpix = LINDGEN(3*2l^(2*res-1))
	lonlat = coorconv(allpix, infmt='p', $
			  inco='r'+STRTRIM(STRING(res),1), $
			  outco=coordsys, outfmt='L')
	lats = lonlat(*,1)
	IF (NOT KEYWORD_SET(noabs)) THEN pixels = $
			WHERE (ABS(lats) GE latmin AND $
			       ABS(lats) LE latmax) $
	ELSE pixels = WHERE (lats GE latmin AND lats LE latmax)
;
;  If no data array is present then just return the pixel list
;  as the function value.  Otherwise, return the selected data
;  in a list. (The pixel list is also returned in the PIXELS 
;  keyword.)
;	
	IF (N_ELEMENTS(data) EQ 0) THEN RETURN, pixels
	RETURN, pix2dat(pix=pixels, ras=data)
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


