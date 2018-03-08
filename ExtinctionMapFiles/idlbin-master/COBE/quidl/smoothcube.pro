	function smoothcube, map, kernel, badval=badval
;+NAME/ONE LINE DESCRIPTION:
;    SMOOTHCUBE convolves an unfolded skycube with an arbitrary kernel.
;
;DESCRIPTION:
;    This smoothing routine operates on unfolded skycubes or single
;    faces, properly accounting in the former case for the way the 
;    edges of the six faces match up.  It will work with an arbitrary 
;    kernel.  For unfolded cubes,  it requires that 
;    the input array be a cube unfolding of the "right-T" type, i.e.
;                                  0
;                               4321
;                                  5
;    If the input array is square, it assumes it to be a single face.
;    A flag value can be specified so that pixels with bad or missing 
;    data are not included in the smoothing process.  The routine 
;    ignores such pixels and renormalizes the part of the kernel 
;    overlying "good" data.
;    When smoothing a single face, all pixels less than a kernel width
;    from the edge are replaced by the flag value.
;
;CALLING SEQUENCE:
;    outmap = smoothcube (inmap, kernel, [badval=...])
;
;ARGUMENTS (I=input, o=output, []=optional): 
;    inmap      I   2-D arr, any type   Input map
;    kernel     I   2-D arr, any type   convolution kernel
;    badval    [I]  keyword  flt        Flag value that identifies empty 
;                                         pixels.  If not present, take
;                                         value from the unused "face" 
;                                         array adjacent to face 5.  If
;                                         input is a single face, use 0.
;    outmap     O   2-D arr             Smoothed map, same type as input
;
;WARNINGS:
;    1. Input map must be an unfolded cube (3:4 aspect ratio) or a
;          square face.
;    2. Convolving kernel need not be square but must be odd 
;          in both dimensions in order to preserve position info in the map.
;    3. Mathematically meaningless results may be obtained if the
;          smoothing kernel is not rotationally symmetric!!
;    4. Finally, it is always safest to specify the badval if known.
;
;EXAMPLE:
;    Consider a map with e.g., -999 in the unfilled pixels, which one 
;    wishes to smooth with the following array:
;                                               1  2  1
;                                               2  5  2
;                                               1  2  1
;    The invocation would be:
;    IDL> kernel = [[1,2,1],[2,5,2],[1,2,1]]
;    IDL> outmap = smoothcube (inmap, kernel, badval=-999)
;
;               ...or, IF the map is well-behaved...
;
;    IDL> outmap = smoothcube (inmap, kernel)
;#
;COMMON BLOCKS: none
;
;PROCEDURE:
;    IF "badval" keyword not present, and the input data is a full cube,
;    SMOOTHCUBE identifies the unused pixel value by looking at face "6", 
;    where there is no data.  If input is a single face, the default 
;    flag value is taken to be zero.  The program uses this flag value 
;    to weight the convolution so as not to include unused pixels in
;    the smoothing process. In the full sky case, the edges are matched
;    up by filling in the borders of the map with the appropriate pieces
;    of connecting faces prior to smoothing.  For a single face, pixels
;    less than a kernel width from the boundary are simply set to the
;    flag value on output. 
;
;LIBRARY CALLS:  None.
;
;REVISION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  31 December 1991
;    RBI fixed incorrect handling of right-T maps.     17 Jan 92  SPR 9407
;    RBI restricted use to right-T only.               28 Feb 92
;    RBI enhanced to handle single faces.              17 Jun 92 SPR 9773
; SPR 9616
;    Dave Bazell added check for WHERE = -1             4 Mar 93 SPR 10642
;
;.TITLE
;Routine SMOOTHCUBE
;-
;
;  First do some checks:  the kernel must be of odd size, and the map 
;  must have a 4:3 (or 1:1) aspect ratio to be an unfolded cube (or cube 
;  face).
;
	ON_ERROR, 2
	kersize = SIZE(kernel)
	IF kersize(0) ne 2 THEN MESSAGE, $
		'Smoothing kernel must be 2-dimensional!'
	IF (kersize(1)/2. - FIX(kersize(1)/2.) lt .49) or $ 
	   (kersize(2)/2. - FIX(kersize(2)/2.) lt .49) THEN MESSAGE,$
		 'Smoothing kernel must be of odd size!'
;
	mapsize = SIZE(map)
	IF mapsize(0) ne 2 THEN MESSAGE, $
		  'Input map must be 2-dimensional.'
	IF (3*mapsize(1) ne 4*mapsize(2)) and $
	     (mapsize(1) ne   mapsize(2)) THEN MESSAGE, $
		  'Input map is not a 3:4 sky cube or square cube face!'
;
;  Branch here if we are only dealing with a single face
;
	IF mapsize(1) eq mapsize(2) THEN GOTO, singleface
;
;  We are dealing with a full sky map (unfolded cube).
;  Flip everything since we figured out the convolution for a left-T...
;
	fsize = mapsize(1)/4
	fs2 = 2 * fsize
	fs3 = 3 * fsize
	fs4 = 4 * fsize
	map = map(mapsize(1)-INDGEN(mapsize(1))-1,*)
;
;  Now get to work.  We will take one face at a time and form a 
;  "bigface" by attaching to the face's vertices the bordering pieces of
;  the appropriate neighboring faces.  Then we will do the convolution.
;
	IF KEYWORD_SET(badval) THEN blank = badval $
	ELSE blank = map(fsize,0)               ; a blank pixel for reference
	outmap = map
	borderx = FIX(kersize(1)/2)
	bordery = FIX(kersize(2)/2)
	bigface = FLTARR(fsize+2*borderx,fsize+2*bordery)
;
	FOR faceno = 0,5 DO BEGIN
	   bigface = 0*bigface + blank
	   CASE faceno OF
		0: BEGIN
			face   = map(0:fsize-1,fs2:fs3-1)
			top    = ROTATE (map(fs2:fs3-1, fs2-bordery:fs2-1),2)
			left   = ROTATE (map(fs3:fs4-1, fs2-bordery:fs2-1),3)
			right  = ROTATE (map(fsize:fs2-1,fs2-bordery:fs2-1),1)
			bottom = map(0:fsize-1, fs2-bordery:fs2-1)
		   END
		1: BEGIN
			face   = map(0:fsize-1,fsize:fs2-1)
			top    = map(0:fsize-1, fs2:fs2+bordery-1)
			left   = map(fs4-borderx:fs4-1, fsize:fs2-1)
			right  = map(fsize:fsize+borderx-1, fsize:fs2-1)
			bottom = map(0:fsize-1, fsize-bordery:fsize-1)
		   END
		2: BEGIN
			face   = map(fsize:fs2-1,fsize:fs2-1)
			top    = ROTATE (map(fsize-bordery:fsize-1, $
							fs2:fs3-1),3)
			left   = map(fsize-borderx:fsize-1, fsize:fs2-1)
			right  = map(fs2:fs2+borderx-1, fsize:fs2-1)
			bottom = ROTATE (map(fsize-bordery:fsize-1, $
							0:fsize-1),1)
		   END    
		3: BEGIN
			face   = map(fs2:fs3-1,fsize:fs2-1)
			top    = ROTATE (map(0:fsize-1, fs3-bordery:fs3-1),2)
			left   = map(fs2-borderx:fs2-1, fsize:fs2-1)
			right  = map(fs3:fs3+borderx-1, fsize:fs2-1)
			bottom = ROTATE (map(0:fsize-1, 0:bordery-1),2)
		   END
		4: BEGIN
			face   = map(fs3:fs4-1,fsize:fs2-1)
			top    = ROTATE (map(0:bordery-1, fs2:fs3-1),1)
			left   = map(fs3-borderx:fs3-1, fsize:fs2-1)
			right  = map(0:borderx-1, fsize:fs2-1)
			bottom = ROTATE (map(0:bordery-1,0:fsize-1),3)
		   END
		5: BEGIN
			face   = map(0:fsize-1,0:fsize-1)
			top    = map(0:fsize-1, fsize:fsize+bordery-1)
			left   = ROTATE (map(fs3:fs4-1, $
					fsize:fsize+borderx-1),1)
			right  = ROTATE (map(fsize:fs2-1, $
						fsize:fsize+borderx-1),3)
			bottom = ROTATE (map(fs2:fs3-1, $
					     fsize:fsize+bordery-1),2)
		   END
	   ENDCASE
;
;       Build the "bigface" and create a weights array corresponding to 
;       blank pixels.
;
	   bigface(borderx:borderx+fsize-1,bordery:bordery+fsize-1) = face
	   bigface(0:borderx-1,bordery:bordery+fsize-1) = left
	   bigface(borderx:borderx+fsize-1, $
				bordery+fsize:2*bordery+fsize-1) = top
	   bigface(borderx+fsize:2*borderx+fsize-1,  $
				bordery:bordery+fsize-1) = right
	   bigface(borderx:borderx+fsize-1,0:bordery-1) = bottom
;
	   weights = bigface
	   index = where(bigface eq blank, cnt)
	   if (cnt ne 0) then weights(index) = 0.
	   index = where(bigface ne blank, cnt)
	   if (cnt ne 0) then weights(index) = 1.
;
;       Do the convolution of the data and the weights, and normalize. 
;
	   index = where(bigface eq blank, cnt)
	   if (cnt ne 0) then bigface(index) = 0.
	   bigface = CONVOL (bigface, kernel)
	   weights = CONVOL (weights, kernel)
	   index = where(weights eq 0., cnt)
	   if (cnt ne 0) then weights(index) = 1.0e15
	   bigface = bigface/weights
	   index = where(weights gt 1.e14, cnt)
	   if (cnt ne 0) then bigface(index) = blank
;
;       Bigface now holds the smoothed face, plus the stuff at the 
;       edges.  Trim off the edges and stick it back into the map array.
;
	   newface = $
		bigface(borderx:borderx+fsize-1,bordery:bordery+fsize-1)
	   CASE faceno of 
		0: outmap(0:fsize-1,fs2:fs3-1) = newface
		1: outmap(0:fsize-1,fsize:fs2-1) = newface
		2: outmap(fsize:fs2-1,fsize:fs2-1) = newface
		3: outmap(fs2:fs3-1,fsize:fs2-1) = newface
		4: outmap(fs3:fs4-1,fsize:fs2-1) = newface
		5: outmap(0:fsize-1,0:fsize-1) = newface
	   ENDCASE
	ENDFOR
;
;  Now flip everything back to a right-T before returning.
;
	outmap = outmap(mapsize(1)-INDGEN(mapsize(1))-1,*)
	map = map(mapsize(1)-INDGEN(mapsize(1))-1,*)
	RETURN, outmap
;
;  Everything from this point on is for a single face only!
;
	singleface: $
;
;  As for a whole skycube we will form a "bigface" by attaching a 
;  border to the outside edges.  Then we will do the convolution.
;
	if KEYWORD_SET(badval) THEN blank = badval $
	ELSE blank = 0.
	fsize = mapsize(1)
	outmap = map
	borderx = FIX(kersize(1)/2)
	bordery = FIX(kersize(2)/2)
	bigface = FLTARR(fsize+2*borderx,fsize+2*bordery) + blank
;
;  Build the "bigface" and create a weights array corresponding to 
;  blank pixels.
;
	   bigface(borderx:borderx+fsize-1,bordery:bordery+fsize-1) = $
								 map
;
	   weights = bigface
	   index = where(bigface eq blank, cnt)
	   if (cnt ne 0) then weights(index) = 0.
	   index = where(bigface ne blank, cnt)
	   if (cnt ne 0) then weights(index) = 1.
;
;  Do the convolution of the data and the weights, and normalize.
;
	   index = where(bigface eq blank, cnt)
	   if (cnt ne 0) then bigface(index) = 0.
	   bigface = CONVOL (bigface, kernel)
	   weights = CONVOL (weights, kernel)
	   index = where(weights eq 0., cnt)
	   if (cnt ne 0) then weights(index) = 1.0e15
	   bigface = bigface/weights
	   index = where(weights gt 1.e14, cnt)
	   if (cnt ne 0) then bigface(index) = blank
;
;  Bigface now holds the smoothed face, plus the stuff at the 
;  edges.  Trim off the edges and stick it back into the map array.
;
	   outmap = $
		bigface(borderx:borderx+fsize-1,bordery:bordery+fsize-1)
;
;  Now zap the pixels that are less than a kernel width from the edge
;
	outmap(0:fsize-1,0:kersize(2)-2) = blank
	outmap(0:fsize-1,fsize-kersize(2)+1:fsize-1) = blank
	outmap(0:kersize(1)-2,0:fsize-1) = blank
	outmap(fsize-kersize(1)+1:fsize-1,0:fsize-1) = blank
	RETURN, outmap
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


