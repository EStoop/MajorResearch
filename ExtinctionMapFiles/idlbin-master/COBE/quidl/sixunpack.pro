	function sixunpack, box_in, badval=badval
;+NAME/ONE-LINE DESCRIPTION:
;    SIXUNPACK "uncompresses" the packed skycube created by SIXPACK.
;
;DESCRIPTION:
;    Data "uncompression" routine takes the space-saving form of 
;    the unfolded sky cube, i.e.
;
;                       445500
;                       445500
;                       332211
;                       332211
;
;    and unpacks it into the "standard" unfolded right-T:
;
;                           00
;                           00
;                     44332211
;                     44332211
;                           55
;                           55
;
;    It will work for both single "sheets" (2-D arrays) and spectral 
;    cubes (3-D).  The empty areas of the array will be filled with
;    zeros or, optionally, a specified bad pixel value.
;
;CALLING SEQUENCE:
;    t_out = SIXUNPACK (box_in [,badval=value] )
;
;ARGUMENTS: 
;    box_in	   2-D or 3-D input data array in packed unfolded cube
;    badval	   keyword used to specify bad pixel or fill value
;    t_out	   2-D or 3-D output array (same type as input) in
;                     right-T-shaped unfolded cube, if an error occurs
;		      then this is equal to box_in
;    !err	
;    !error	     0 if no errors, 1 if an error occurs
;
;
;WARNINGS:
;    1. Sixunpack checks to see that the input array has the proper x/y
;         aspect ratio for an packed cube, and quits if it does not.
;    2. One may run into virtual memory problems for large 
;         spectral cubes since the output array is twice as large as the 
;         original.  
;    3. It is possible to do the uncompression "in place" by using the 
;	  same variable as input and output in the calling sequence.  
;
;EXAMPLE:
;    tdata = SIXUNPACK(pdata, badval=-9999 )
;	
;    The data array will be reformatted from the packed cube into an
;    unfolded cube with the empty portions set to -9999. 
;
;#
;COMMON BLOCKS: none
;
;PROCEDURE:
;    Creates a degenerate 3rd dimension for two dimensional arrays,
;    and does all its work on 3-D arrays.  (The 3rd dimension is
;    removed prior to output.)  
;
;PERTINENT ALGORITHMS, LIBRARY CALLS:  none
;
;MODIFICATION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  23 July 1991
;    Revised by RBI to restrict output to a right-T.  28 Feb 92
;    RBI fixed bug regarding variable "nullface"  9 Apr 92
;    Pete Kryszak-Servin, added BADVAL keyword for fill and corrected
;			  documentation to allow "in-place" unpacking.
;			  SER 9838 
;    Pete Kryszak-Servin, error handling thru !err and message, returns
;			  original array in those cases also!  SPR 9841
;
; SPR 9616
;.TITLE
;Routine SIXUNPACK
;-
;
;  guilty of errors until successful 
   !err = 1
   !error = 1
;
;  first check dimensionality and create degenerate 3rd dimension if 
;  needed
;
	ndims = size(box_in)
	if ndims(0) eq 2 then begin
	   box_in = reform (box_in, ndims(1), ndims(2), 1)
	   depth = 1
	endif else begin
	   depth = ndims(3)
	endelse
	xsize = ndims(1)
	ysize = ndims(2)
;
;  make sure that this is in fact a compressed unfolded cube
;
	if 2*xsize ne 3*ysize then begin
	   Message, /continue, $
		"This is not a compressed sky cube!  SIXUNPACK dies horribly!"
	   return, box_in
	endif
;
;  create the output array with appropriate data type
;
	fsize = xsize/3
	case ndims(ndims(0)+1) of
		1:  t_out = bytarr(4*fsize, 3*fsize, depth)
		2:  t_out = intarr(4*fsize, 3*fsize, depth)
		3:  t_out = lonarr(4*fsize, 3*fsize, depth)
		4:  t_out = fltarr(4*fsize, 3*fsize, depth)
		5:  t_out = dblarr(4*fsize, 3*fsize, depth)
		6:  t_out = complexarr(4*fsize, 3*fsize, depth)
		7:  t_out = strarr(4*fsize, 3*fsize, depth)
		else:  begin
		          Message,/continue, $
				'Flaky data type!  SIXUNPACK croaks!'
			  return, box_in
		       end
	endcase
;
;  initialize t_out to the bad pixel value if given, otherwise
;  let the array stand as initialized by IDL
if Keyword_Set( badval ) then t_out(*,*,*) = badval

;
;  fill in the upper left 3 x 2 corner of the T.  Will have to zero out
;  box faces 4 and 5 later, so create a null face.
;
	t_out(fsize:4*fsize-1, fsize:3*fsize-1, 0:depth-1) = box_in
	nullface = t_out(0:2*fsize-1, 0:fsize-1,*)
;
;  fill in the 4 and 5 faces, then blank them out of the T
;
	t_out(3*fsize:4*fsize-1, 0:fsize-1, 0:depth-1) = $
	   box_in(fsize:2*fsize-1, fsize:2*fsize-1, *)   ;face 5
	t_out(0:fsize-1, fsize:2*fsize-1, 0:depth-1) = $
	   box_in(0:fsize-1, fsize:2*fsize-1, *) ;face 4
	t_out(fsize:3*fsize-1, 2*fsize:3*fsize-1, 0:depth-1) = $
	   nullface
;
;  eliminate the degenerate dimension if this was a 2-D array
;
	if depth eq 1 then begin
	   box_in = reform (box_in, ndims(1), ndims(2))
	   t_out = reform (t_out, 4*fsize, 3*fsize)
	endif
	!err = 0
	!error = 0
	return, t_out
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


