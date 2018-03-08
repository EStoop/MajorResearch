	FUNCTION gapfill, skycube, badval=badval
;+NAME/ONE-LINE DESCRIPTION:
;   GAPFILL interpolates across gaps in an unfolded skycube.
;
;DESCRIPTION:
;   GAPFILL uses Delaunay triangulation to interpolate across the gaps
;   in an unfolded skycube.  It will work in principle on any array;
;   the only significance of the cube is that the shape of the good data
;   area is known so that the array points outside of the T are ignored.
;   This routine is NOT intended for rigorous work, since the 
;   interpolation is only linear and no statistical weighting is 
;   employed; it is primarily for creating nicely filled-in maps for
;   presentations and other glitz.
;   
;CALLING SEQUENCE:
;   outmap = GAPFILL (inmap, [BADVAL=flag_value])
;
;ARGUMENTS (I=input, O=output, []=optional):
;   OUTMAP    O    2-D array    Output interpolated map
;   INMAP     I    2-D array    Input map with gaps    
;  [BADVAL]  [I]   scalar       Flag value, i.e. the value in the gaps.
;                               Default = 0
;
;WARNINGS:
;   Don't try to do anything serious with the output.
;
;EXAMPLE:
;   If a DIRBE map "dmap" has gaps containing a sentinel value -16000 
;   then the call   
;                     Nice_Map = GAPFILL (dmap, badval=-16000)
;
;   will produce a completely filled-in map in the variable Nice_Map.
;#
;COMMON BLOCKS:  None
;
;PROCEDURE:
;   Uses the WHERE function to identify the points with valid data and
;   triangulates on those using IDL's TRIANGULATE routine.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS:
;   UPKWHERE:  This UIDL routine is called to unpack the WHERE vector
;              into x and y coordinates.
;
;MODIFICATION HISTORY:
;   Written by Rich Isaacman, General Sciences Corp.    26 Apr 93 
;-
	ON_ERROR,2
	IF NOT KEYWORD_SET(badval) THEN badval = 0        ; Set default
	dims = SIZE (skycube)
	IF (4*dims(2) NE 3*dims(1)) THEN MESSAGE, $
		'Input map is not an unfolded cube.'
	fsize = dims(1)/4
	outmap = skycube*0 + badval                       ; Make output map
;
;  Identify the valid data locations, then unpack them into x-y pairs and
;  use the pairs to create the Delaunay triangles
;
	goodpts = WHERE(skycube NE badval)
        IF (N_ELEMENTS(goodpts) LT 3) THEN MESSAGE, $
             'There is too little good data in this map!'
	UPKWHERE,skycube,goodpts,xok,yok
	TRIANGULATE,xok,yok,triangles
;
;  Create an output grid with the same spacing as the input map, then
;  do the interpolation and fill it.
;
	spacing = [1,1]
	newgrid = TRIGRID(xok,yok,skycube(goodpts),triangles,spacing)
	outmap(*,fsize:2*fsize-1) = newgrid(*,fsize:2*fsize-1)
	outmap(3*fsize:4*fsize-1,*) = newgrid(3*fsize:4*fsize-1,*)
	RETURN,outmap
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


