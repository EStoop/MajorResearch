	function downres, inmap, nsteps, badval=badval
;+NAME/ONE-LINE DESCRIPTION:
;    DOWNRES "coarsens" the resolution of a skycube or cube face.
;
;DESCRIPTION:
;    "Downsizes" a skycube or cube face NSTEPS in resolution, i.e.
;    turns an order-N cube into one of order N-NSTEPS.  Each unit of
;    nsteps is a factor of 2 in spatial resolution.  Both 2-D and 3-D
;    maps are handled; if 3-D, the third dimension is assumed to be
;    non-spatial so its resolution is unchanged. 
;
;CALLING SEQUENCE:
;    outmap = DOWNRES (inmap, nsteps, [badval=bad_pixel_value])
;
;ARGUMENTS (I=input, O=output, []=optional): 
;    outmap   O   2D/3D arr       Smoothed map, smaller than inmap by
;				    a factor of 2^nsteps in x and y
;    inmap    I   2D/3D arr       Input map: byte, real or integer
;    nsteps   I   int             Number of steps to decrease res'n
;    badval  [I]  flt   keyword   The value of pixels assumed to
;				     be empty or bad; defaults to 0.
;
;WARNINGS:
;    1. Pixels whose value is equal to badval are not included in the 
;          averaging process.  Thus, stripes whose width is smaller than
;          2^nsteps pixels will appear to have disappeared in the low
;          resolution output map.
;    2. Input map must have either a 3:4 aspect ratio in y:x (unfolded 
;          cube), a 2:3 ratio ("sixpacked" cube), 
;          or be square.  (In the latter case, we assume 
;          we are looking at a single cube face.)
;    3. Only a single bad value is accepted.  Maps with multiple
;          sentinel values (e.g. some DIRBE maps) will be confused
;          unless changed to a single value prior to being "DOWNRES"ed.
;
;EXAMPLE:
;    outmap = DOWNRES(inmap, 3, badval=-999) 
;
;    This would convert e.g., a 1024x768 DIRBE map (i.e. order 9) 
;    into a (128x96) FIRAS map (order 6).
;#
;COMMON BLOCKS:  none
;
;PROCEDURE:
;    If, for example, nsteps=1 then the resulting map will be half as 
;    big (in numbers of pixels) as the input map in x and y, with pixels
;    that are 4 times as large in area.  If a bad pixel is 
;    present in the original map, the resulting big pixel there will 
;    be the average of the remaining three good pixels.
;
;    Information in the "z" direction is unchanged; each slice at 
;    constant z is handled as a separate map and downsized 
;    separately.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS:  none
;
;MODIFICATION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  2 January 1992
;
; SPR 9616
;.TITLE
;Routine DOWNRES
;-
;
;  Do some error checking for aspect ratio and figure out 
;  whether we are looking at a whole unfolded cube or just a face.
;
	dims = size (inmap)
	ndims = dims(0)
	if (ndims lt 2) or (ndims ge 4) then begin
	   print, 'DOWNRES error: Can only handle 2-D or 3-D maps'
	   return, 0
	endif

	mapsize = size (inmap)
	if 3*mapsize(1) eq 4*mapsize(2) then face = mapsize(1)/4    $
	   else if mapsize(1) eq mapsize(2) then face = mapsize(1)  $
	   else if 2*mapsize(1) eq 3*mapsize(2) then face = mapsize(1)/3  $
	   else begin
	   	   print, $
		'DOWNRES error:  Map is neither a whole cube nor a face.'
		   return, 0
	endelse

	newsize = 2^nsteps
	if (nsteps lt 0) or (newsize gt face) then begin
	   print, 'DOWNRES error: Resolution step size makes no sense'
	   return, 0
	endif
;
;  We will always work with a 3-D map, with a degenerate z-dimension 
;  if necessary.
;
	case ndims of
	   2:  begin
		 zrange = 0
		 inmap = reform (inmap, mapsize(1), mapsize(2), 1)
	       end
	   3:  zrange = dims(3) - 1
	endcase
	outmap = fltarr (mapsize(1)/newsize, mapsize(2)/newsize, zrange+1)
	if (not keyword_set(badval)) then badval = 0.
;
;  Now the main loop.  Since the new pixel size is 2^nstep, the DO
;  loops hopscotch through each dimension with steps of that size.  At 
;  each point, the square of data is averaged to make the new big pixel.
;  The averaging is done by totalling the number of non-bad pixels in 
;  the region and dividing by that number.  If there are no good pixels, 
;  add a one so you don't divide zero by zero.
;
	for z=0,zrange do begin
	   xpos = -1
	   for x=0,mapsize(1)-newsize,newsize do begin
	      xpos = xpos + 1
	      ypos = -1
	      for y=0,mapsize(2)-newsize,newsize do begin
		 ypos = ypos + 1
		 newdat = inmap(x:x+newsize-1,y:y+newsize-1,z)
		 goodpix = newdat ne badval
		 newdat = newdat * goodpix			; Mask out bad guys (set them = 0)
		 area = total(goodpix) + (total(goodpix) eq 0.) ; Total good pixels; add one if they're all bad
		 outmap(xpos,ypos,z) = total(newdat)/area
	      endfor
	   endfor
	endfor
;
;  Get rid of degenerate 3rd dimension if necesssary, and return.
;
	if zrange eq 0 then begin 
	   outmap = reform (outmap, mapsize(1)/newsize, mapsize(2)/newsize)
	   inmap = reform (inmap, mapsize(1), mapsize(2))
	endif
	return, outmap
	end
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


