pro deproj, proj_arr, qcube, res, proj=proj, coord=coord, $
    min=min,max=max,win=win,noshow=noshow
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    DEPROJ is a COBE specific reprojection facility. This does not 
;    preserve photometric accuracy.
;
;DESCRIPTION:
;    This is an IDL procedure that creates COBE quadrilateralized
;    skycubes from various projections. The actual display of the 
;    projection can be suppressed if desired in which case the rasterized 
;    array can be returned as output.  The output skycube may be in either 
;    unfolded or sixpacked format. The output resolutions must be 6 or 9.  
;    If any of the following parameters are not specified on
;    the command line then the user is prompted for them:  projection 
;    type (Aitoff, Global Sinusoidal, Mollweide), coordinate
;    system (Ecliptic, Galactic, Equatorial) and low and high values to 
;    be used in the display scaling (the defaults are the input min and max). 
;
;    This routine is meant to be used for display purposes only or for
;    large scale structures. Due to the lack of a 1 to 1 mapping in
;    this reprojection small scale structures and photometry can be
;    distorted!
;
;    Only right oriented, ecliptic coordinate images are built.
;
;CALLING SEQUENCE:
;    pro deproj,proj_arr,qcube,res,[coord=coord] [,proj=proj] [,/sixpack]
;        [,min=min] [,max=max] [,/noshow]
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    proj_arr    I     flt arr            "Projected" input map
;    qcube       O     flt arr             Rasterized output array
;    res         I     integer             Resolution= 6 or 9
;    [proj]      I     string              Projection Type:
;                                          ('A','S','M')
;    [coord]     I     string              Coordinate System:
;                                          ('E','G','Q')       
;    [sixpack]   I     qualifier           Output array in sizpack format
;    [min]       I     flt                 Image Scale Minimum
;                                          Default: Input Min
;    [max]       I     flt                 Image Scale Maximum
;                                          Default: Input Max
;    [win]       I     integer             display window
;    [noshow]    I     qualifier           no display proj switch
;
;WARNINGS:
;   1) THIS ROUTINE DOES NOT PRESERVE PHOTOMETRIC ACCURACY!!!!!!
;    NO PARTIAL PIXEL SAMPLING IS PERFORMED!!!!!!!!!!!
;    In order to test accuracy the computed skycube can be used
;    as an input to the procedure REPROJ. This will re-create
;    the original projection and can be compared to the original
;    to note the errors.
;
;EXAMPLE: 
;
; To reproject an Aitoff projection in galactic coordinates (stored
; in the IDL array 'inproj') into a right handed quadrilateralized
; skycube (out_raster) of resolution 9, using the default min and 
; max and displaying the image in window 1 of X-terminal:
;
; deproj, inproj,out_raster,9,proj='a',coord='g',win=1
;
; ____
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Compute the unit vectors for a skycube of designated resolution
;      and proper projected coordinate system
;    Compute x and y projection (screen) coordinates for these
;      unit vectors and projection type
;    Obtain data values at these coordinates
;    Produce rasterized skycube
;    Display image if enabled
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: coorconv,uv2proj,pix2xy
;
;MODIFICATION HISTORY
;    Written by J. Newmark, Applied Research Corp.  Oct. 1994 
;    In response to SPR's 10407, 10506. Note, no direct link
;     into UIMAGE, 2 step process by design.
;
;.TITLE
;Routine DEPROJ
;-

IF (keyword_set(noshow) EQ 0) THEN noshow = 0
IF (NOT KEYWORD_SET(sixpack)) THEN sixpack = 0
IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN noshow = 1
CASE res OF
   9 :  begin
         pixnum=393216.
         inco='r9'
        end
   6 :  begin
          pixnum=6144.
          inco='r6'
        end
 ELSE:	begin
	str = 'Unknown or improper resolution, (6 or 9):'
	str = str + '"' + string(res) + '"'
	message,str,/cont
	return
	end
ENDCASE
;obtain unit vectors for all the points on a quad-cube in
;    chosen projection:
pix=findgen(pixnum)
pix=long(pix)
coor_table = ['e','g','q']
IF (keyword_set(coord) EQ 0) THEN BEGIN
    i_coor = umenu(['Enter Projection Coordinate System', $
		  'Ecliptic', $
                  'Galactic', $
	          'Equatorial'], title=0,init=1)
    coord = coor_table(i_coor-1)
ENDIF
uv=coorconv(pix,infmt='p',outfmt='u',inco=inco,outco=coord)

;compute xy positions of unit vectors in reprojected coords
;NOTE: there may not be the same number of unit vectors as pixels
;in the projection
sz=size(proj_arr)
proj_table = ['A','S','M']
IF (keyword_set(proj) EQ 0) THEN BEGIN
	i_proj = umenu(['Enter Projection Type', $
                	'Aitoff', $
	                'Global Sinusoidal', $
        	        'Mollweide'],title=0,init=1) 
	proj = proj_table(i_proj-1)
ENDIF 
xy=uv2proj(uv,proj,sz)

;grab each data value at x,y from the projection and assign it
;     to a pixel number. THIS IS THE STEP THAT DOES NOT DO
;     PROPER SAMPLING!!
quad_array=proj_arr(xy(*,0),xy(*,1))

;create new quad-cube
IF (sixpack ne 0) then $
  pix2xy,pix,res=res,data=quad_array,raster=qcube,/sixpack $
ELSE pix2xy,pix,res=res,data=quad_array,raster=qcube

;clean up temporary array
quad_array=0

IF (noshow NE 1) THEN BEGIN
        if (keyword_set(win) eq 0) then win=0
	IF (N_ELEMENTS(min) EQ 0) THEN $
		min = MIN(proj_arr) ELSE min = FLOAT(min)
	IF (N_ELEMENTS(max) EQ 0) THEN $
		max = MAX(proj_arr) ELSE max = FLOAT(max)
        szq=SIZE(qcube)
	bim = BYTSCL(qcube,min = min,max = max)
	WINDOW,win,xsize=szq(1)+64, ysize=szq(2)+64,retain=2
	w_pos = [32,32]
	TV, bim, w_pos(0),w_pos(1)
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


