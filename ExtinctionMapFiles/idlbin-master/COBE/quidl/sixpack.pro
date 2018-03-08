	function sixpack, t_in, t_orient
;+NAME/ONE-LINE DESCRIPTION:
;    SIXPACK packs an unfolded skycube into 3x2 format (no wasted space).
;
;DESCRIPTION:
;    Data compression routine takes a standard unfolded skycube, i.e.
;
;       00                                    00
;       00                                    00
;       11223344       or               44332211
;       11223344                        44332211
;       55                                    55
;       55                                    55
;
;       "Left T"                       "Right T" 
;
;    and compresses it into the space-saving form:
;
;                        445500
;                        445500
;                        332211
;                        332211
;
;    It will work for both single "sheets" (2-D arrays) and spectral 
;    cubes (3-D).  Note that the output corresponds to a "right-T"
;    regardless of the orientation of the input.
;
;CALLING SEQUENCE:
;    box_out = SIXPACK (t_in, [t_orient])
;
;ARGUMENTS: 
;    t_in	     2-D or 3-D data array in "standard" unfolded cube
;    t_orient        optional character string 'L' or 'R' to describe 
;		         orientation of input T.  Default is 'R'.
;    box_out	     2-D or 3-D array (same as input) compressed 
;                        unfolded cube, if an error occurs then
; 			 this is equal to t_in
;    !err	
;    !error	     0 if no errors, 1 if an error occurs
;
;WARNINGS:
;    1. Sixpack checks to see that the input array has the proper x/y
;         aspect ratio for an unfolded cube, and quits if it does not.
;    2. One may run into virtual memory problems for large 
;         spectral cubes since the output array is half as large as the 
;         original.  
;    3. It is possible to do the compression "in place" by using the 
;	  same variable as input and output in the calling sequence.  
;
;EXAMPLE:
;    outmap = sixpack(incube,'R')
;    outmap = sixpack(incube)
;#
;COMMON BLOCKS: none
;
;PROCEDURE:
;    Creates a degenerate 3rd dimension for two dimensional arrays,
;    and does all its work on 3-D arrays.  (The 3rd dimension is
;    removed prior to output.)  For historical reasons, the routine
;    flips a right-T to the left before compressing, then "unflips" 
;    it for output. 
;
;PERTINENT ALGORITHMS, LIBRARY CALLS:  none
;
;MODIFICATION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  23 July 1991
;    Revised by RBI 28 Feb 92 to use a right-T packing and change 
;        default input to right-T
;    Revised by Pete Kryszak-Servin, GSC, 27 Jul 1992 to document
;	ability to do reformatting 'in-place'.  SPR 9841.
;    Revised by Pete Kryszak-Servin, GSC, 7 Aug 1992 to use !err and
;	message in the case of errors, returns original in case of error
;	SPR 9838.
;
; SPR 9616
;.TITLE
;Routine SIXPACK
;-
;
; guilty of errors until successful
!err = 1
!error = 1
;
;
;  first check dimensionality and create degenerate 3rd dimension if 
;  needed
;
	ndims = size(t_in)
	if ndims(0) eq 2 then begin
	   t_in = reform (t_in, ndims(1), ndims(2), 1)
	   depth = 1
	endif else begin
	   depth = ndims(3)
	endelse
	xsize = ndims(1)
	ysize = ndims(2)
;
;  make sure that this is in fact a standard unfolded cube
;
	if 3*xsize ne 4*ysize then begin
	   Message, /continue, "This is not a standard unfolded cube!"  
	   Message, /continue, "SIXPACK fails messily, all over the terminal!"
	   return, t_in
	endif
;
;  now flip cube to right-T t_orientation if necessary.  'R' is default.
;
	if n_elements(t_orient) eq 0 then t_orient = 'R'
	if strupcase(t_orient) eq 'L' then  	    $
	   t_in(0:xsize-1, 0:ysize-1, 0:depth-1) = $
	      t_in(xsize-1-indgen(xsize), *, *)
;
;  create the output array
;
	fsize = xsize/4
	box_out = t_in(fsize:4*fsize-1, fsize:3*fsize-1, *)  ;faces 0,1,2,3
	box_out(fsize:2*fsize-1,fsize:2*fsize-1,0:depth-1) = $
		t_in(3*fsize:4*fsize-1,0:fsize-1,*)		  ;face 5
	box_out(0:fsize-1,fsize:2*fsize-1,0:depth-1) = $
		t_in(0:fsize-1,fsize:2*fsize-1,*)  		;face 4
;
;  restore the input map to its original t_orientation and depth 
;
	if strupcase(t_orient) eq 'L' then begin
	   t_in(0:xsize-1, 0:ysize-1, 0:depth-1) = $
	      t_in(xsize-1-indgen(xsize), *, *)
	endif
	if depth eq 1 then begin
	   t_in = reform (t_in, ndims(1), ndims(2))
	   box_out = reform (box_out, 3*fsize, 2*fsize)
	endif
	!err = 0
	!error = 0
	return, box_out
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


