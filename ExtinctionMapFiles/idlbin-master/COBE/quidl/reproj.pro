pro reproj,inmap,raster,noshow=noshow,proj=proj,coord=coord, $
	   psize=psize,face=face_num,win=win,merd=merd,para=para, $
           gcoord=gcoord,min=min,max=max
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    REPROJ is a general reprojection facility.
;
;DESCRIPTION:
;    This is an IDL procedure that creates reprojections of the
;    native skycubes.  The actual display of the projection can be
;    suppressed if desired in which case the rasterized array can be
;    returned as output.  The skycube may be in either unfolded or 
;    sixpacked format.  Note that a single face may also be used for 
;    input.   If any of the following parameters are not specified on
;    the command line then the user is prompted for them:  projection 
;    type (Quad Cube, Aitoff, Global Sinusoidal, Mollweide), coordinate
;    system (Ecliptic, Galactic, Equatorial), projection size (Small:
;    512 x 256, Large:  1024 x 512), and low and high values to be used
;    in the display scaling (the defaults are the input min  and max). 
;    Meridians and parallels can be overlaid on either the
;    reprojections or the skycubes.  These must be specified on the
;    command line.  The meridians and parallels are NOT stored as part
;    of the rasterized image.
;
;CALLING SEQUENCE:
;    pro reproj,inmap,[raster],[proj=proj],[coord=coord], $
;	   [psize=psize],[face=face_num],[win=win],[merd=merd], $
;          [para=para],[gcoord=gcoord],[min=min],[max=max],[/noshow]
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    inmap       I/O   flt arr             Sky Cube/Face (Res 6-9)
;    [raster]    O     flt arr             Rasterized output array
;    [proj]      I     string              Projection Type:
;                                          ('Q','A','S','M')
;    [coord]     I     string              Coordinate System:
;                                          ('E','G','Q')
;    [psize]     I     string              Projection Size ('S','L')
;    [face]      I     int                 face number (face inpt only)
;    [win]       I     int                 Window # (default:0)
;    [merd]      I     flt arr             Meridian (longitude) array
;                                          Enter "-1" to suppress
;    [para]      I     flt arr             Parallel (latitude) array
;                                          Enter "-1" to suppress
;    [gcoord]    I     string              Grid coordinate system
;                                          ('E','G','Q','N')
;                                          'N' suppresses overlays
;    [min]       I     flt                 Image Scale Minimum
;                                          Default: Input Min
;    [max]       I     flt                 Image Scale Maximum
;                                          Default: Input Max
;    [noshow]    I     qualifier           no display proj switch
;
;WARNINGS:
;    If displaying an image, an X-windows terminal must be used.
;
;EXAMPLE: 
;
; To reproject a skycube (stored in the IDL array called 'incube') into
; a small, Aitoff projection in galactic coordinates, store the
; rasterized image in the IDL array 'ras_1' and display it in window 2
; with the default meridians (0,60,120,180) and parallels (0,30,60)
; and a "forced" minimum of 0:
;
; reproj, incube,out_raster,proj='a',coord='g',psize='s',win=2,min=0
;
; It is not necessary to designate whether the skycube is in unfolded
; or sixpacked format.
; ____
;
; To reproject a skycube (stored in the IDL array called 'incube')
; and then store the rasterized image in the IDL array 'ras_1' without 
; displaying it:
;
; reproj, incube,ras_1,/noshow
;
; It is not necessary to designate the meridian and parallel 
; information in this case, indeed it will be ignored if provided.
; ____
;
; To reproject a skycube (stored in the IDL array called 'incube')
; and display the rasterized image (in the default window) with the
; meridians (40,50,60) and no parallels, without storing the raster
; image:
;
; reproj, incube,merd=[40,50,60],para=-1,/nofile
;
; In this case the user will be prompted for projection type,
; coordinate system, projection size, and minimum and maximum 
; for image display scaling.
; ____
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Essentially a routing routine.  Querys for info and uses 
;    if-then-else and case statements to make the proper call
;    to JPRO which actually generates and displays the reprojection.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: jpro
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Apr 92
;
;    SPR 9770	1) noshow set to 1 for non-X windows terminals
;		2) grid coordinate system prompt bypassed for
;		   proj = -1 & merd = -1
;
;    17 JUN 1992
;
; SER 9957    Sixpack input capabilities
;             Allow sixpacked cube input in addition to unfolded
;             skycubes and single faces
; 9-SEP-92    J.M. Gales
;
; SPR 10043   '/no_c' command line switch to allow overide of
;             cproj call
; 24-SEP-92   J.M. Gales
;
; SPR 10212   Eliminate '/no_c' command line switch
; 13-NOV-92   J.M. Gales
;
; SPR 11151   1) Set image scaling min/max to input min/max if not supplied
;                on the command line.
;             2) Suppress coordinate grid overlays by "gcoord='n'"
; 13-JUL-93   J.M. Gales
;
;.TITLE
;Routine REPROJ
;-
;
if (N_ELEMENTS(face_num) eq 0) then face_num = -1
if (keyword_set(win) eq 0) then win=0
if (keyword_set(noshow) eq 0) then noshow = 0
if (keyword_set(merd) eq 0) then merd = [0,60,120,180]
if (keyword_set(para) eq 0) then para = [0,30,60]

proj_table = ['q','a','s','m']
size_table = ['s','l']
coor_table = ['e','g','q']
orient_table= ['l','r']


dims = SIZE(inmap)
IF (dims(0) NE 2) THEN BEGIN
	MESSAGE,'Input Array Must be Two-Dimensional',/CONT
	GOTO, exit
ENDIF

orient = 'R'

IF (!d.name NE 'X' and !d.name ne 'WIN') THEN noshow = 1

IF (keyword_set(proj) EQ 0) THEN BEGIN
	i_proj = umenu(['Enter Projection Type', $
		        'Quad Cube/Face', $
                	'Aitoff', $
	                'Global Sinusoidal', $
        	        'Mollweide'],title=0,init=1) 

	proj = proj_table(i_proj-1)
ENDIF

IF ((noshow EQ 1) AND (STRUPCASE(STRMID(proj,0,1)) EQ 'Q')) THEN GOTO, exit

IF (keyword_set(psize) EQ 0) THEN BEGIN
	psize = ' '
	IF (STRUPCASE(proj) NE 'Q') THEN BEGIN

	i_size = umenu(['Enter Projection Size', $
		        'Small (512 x 256)', $
                	'Large (1024 x 512)'],title=0,init=1) 

	psize = size_table(i_size-1)

	ENDIF

ENDIF

IF ((keyword_set(coord) EQ 0) AND $
    (STRUPCASE(STRMID(proj,0,1)) NE 'Q')) THEN BEGIN

	i_coor = umenu(['Enter Projection Coordinate System', $
		        'Ecliptic', $
                	'Galactic', $
	                'Equatorial'], title=0,init=1)

	coord = coor_table(i_coor-1)

ENDIF



IF (keyword_set(gcoord) NE 0) THEN BEGIN
	IF (STRUPCASE(STRMID(gcoord,0,1)) EQ 'N') THEN BEGIN
		merd = -1
		para = -1
	ENDIF
ENDIF


sz_m = SIZE(merd)
sz_p = SIZE(para)

IF ((keyword_set(gcoord) EQ 0) AND (noshow EQ 0) AND $
    (sz_m(0)+sz_p(0) GT 0)) THEN BEGIN

	i_coor = umenu(['Enter Grid Coordinate System', $
		        'Ecliptic', $
                	'Galactic', $
	                'Equatorial', 'None'], title=0,init=1)

	IF (i_coor LE 3) THEN BEGIN
		gcoord = coor_table(i_coor-1)
	ENDIF ELSE BEGIN 
		gcoord = ' '
		merd = -1
		para = -1
	ENDELSE

ENDIF ELSE IF (keyword_set(gcoord) EQ 0) THEN gcoord = ' '



IF (noshow NE 1) THEN BEGIN

	IF (N_ELEMENTS(min) EQ 0) THEN $
		min = MIN(inmap) ELSE min = FLOAT(min)
	IF (N_ELEMENTS(max) EQ 0) THEN $
		max = MAX(inmap) ELSE max = FLOAT(max)

ENDIF	; noshow


IF (STRUPCASE(proj) NE 'Q') THEN BEGIN

case face_num of -1: begin

jpro,inmap,proj=proj,coord=coord,psize=psize,min=min,max=max,win=win, $
     merd=merd,para=para,gcoord=gcoord,noshow=noshow,raster

		end

else: begin

jpro,inmap,proj=proj,coord=coord,psize=psize,min=min,max=max, $
     win=win,merd=merd,para=para,gcoord=gcoord,noshow=noshow, $
     face=face_num,raster

	end

endcase

endif else begin ; quad cube display

case face_num of -1: begin

jpro,inmap,coord=coord,psize=psize,min=min,max=max,win=win, $
     merd=merd,para=para,gcoord=gcoord,noshow=noshow,raster

		end

else: begin

jpro,inmap,coord=coord,psize=psize,min=min,max=max, $
     win=win,merd=merd,para=para,gcoord=gcoord,noshow=noshow, $
     face=face_num,raster

		end
endcase
endelse

exit:

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


