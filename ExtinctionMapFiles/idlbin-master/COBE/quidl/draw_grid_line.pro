PRO draw_grid_line, cube_input, vec_arc, proj, sz_proj, w_pos, color, $
		    zoom=zoom
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    draw_grid_line - Draws coordinate grid lines on skycubes and projections.
;
;DESCRIPTION:
;    This IDL procedure draws the coordinate grid lines (meridians and 
;    parallels) on either unfolded skycubes or reprojected images.  It 
;    is called by JPRO and is not a user routine.
;
;CALLING SEQUENCE:
;    draw_grid_line, cube_input, vec_arc, proj, sz_proj, w_pos, color, $
;                    zoom=zoom
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    cube_input  I     int arr             face number input,
;					   skycube column of grid line,
;					   skycube row of grid line
;					   (quad cube projection)
;    vec_arc     I     flt arr             unit vector array of grid line
;					   (non-cube projection)
;    proj        I     string              projection type
;    sz_proj     I     long arr            projection size array
;    w_pos       I     int arr (2 elem)    position of projection in window
;    color       I     int                 color of grid line
;    zoom        I     [flt]               zoom factor
;
;WARNINGS:
;    None
;
;EXAMPLE: 
;
;    None.  Not user routine.
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    If sky cube overlay then break the grid line into separate
;    face units (to eliminate drawing within "empty space" of T.)
;    If projection overlay, use UV2PROJ to convert from
;    spatial unit vector to projection (x,y) coordinate.  In both
;    cases, useS the IDL PLOTS command to overlay grid.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Calls procedure UV2PROJ
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Apr 92
;
;    SPR 12144  Pass projection type and size through command line
;               Remove common block
;               Pass quad-cube input via single array CUBE_INPUT
;
;    16-Mar-95  J.M. Gales
;
;.TITLE
;Routine DRAW_GRID_LINE
;-

IF (NOT KEYWORD_SET(zoom)) THEN zoom = 1

IF (N_ELEMENTS(cube_input) NE 0) THEN BEGIN

	n_pix_arc = N_ELEMENTS(cube_input(*,0))
	n_brk = 0
	prev_face = cube_input(0,0)
	face_break = intarr(16)
	FOR i=1,n_pix_arc-1 DO BEGIN

		IF (cube_input(i,0) NE prev_face) THEN BEGIN

			n_brk = n_brk + 1
			face_break(n_brk) = i
			prev_face = cube_input(i,0)

		ENDIF

	ENDFOR

	face_break(n_brk+1)= n_pix_arc-1

	FOR i=0,n_brk DO BEGIN

	ival = cube_input(face_break(i):face_break(i+1)-1,1)
	jval = cube_input(face_break(i):face_break(i+1)-1,2)
        good=where((ival ge 0) and (jval ge 0), ngood)
        if (ngood gt 0) then $
	PLOTS,(ival(good)+w_pos(0))*zoom,(jval(good)+w_pos(1))*zoom,$
           COLOR=color,/DEVICE,LINESTYLE=0,THICK=1

	ENDFOR

ENDIF ELSE BEGIN

	n_pix_arc = N_ELEMENTS(vec_arc(*,0))
	n_brk = 0
	gc = uv2proj(vec_arc,proj,sz_proj)

	prev_seg_i = gc(0,0)
	prev_seg_j = gc(0,1)
	seg_break = intarr(16)
	FOR i=1,n_pix_arc-1 DO BEGIN

		IF ((ABS(gc(i,0)-prev_seg_i) GT 20) OR $
		    (ABS(gc(i,1)-prev_seg_j) GT 20)) THEN BEGIN

			n_brk = n_brk + 1
			seg_break(n_brk) = i

		ENDIF

		prev_seg_i = gc(i,0)
		prev_seg_j = gc(i,1)

	ENDFOR

	seg_break(n_brk+1)= n_pix_arc-1

	FOR i=0,n_brk DO BEGIN

		IF (seg_break(i) LT seg_break(i+1)) THEN BEGIN
		        ival=gc(seg_break(i):seg_break(i+1)-(i NE n_brk),0)
        		jval=gc(seg_break(i):seg_break(i+1)-(i NE n_brk),1)

	        good=where((ival ge 0) and (jval ge 0), ngood)
        	if (ngood gt 0) then $
		PLOTS,(ival(good) + w_pos(0)) * zoom, $
		      (jval(good) + w_pos(1)) * zoom, $
		      COLOR=color,/DEVICE,LINESTYLE=0,THICK=1

		ENDIF

	ENDFOR

ENDELSE

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


