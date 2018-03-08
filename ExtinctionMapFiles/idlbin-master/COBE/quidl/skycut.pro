PRO skycut, input, nat_cube=nat_cube, proj=proj, coord=coord, $
            win=win, w_pos=w_pos, pwin=pwin, noplot=noplot, $
            sclmin=sclmin, sclmax=sclmax, face=face_num, $
            color=color,n_pnts=n_pnts, cs_plot=cs_plot, $
            cs_titl=cs_titl
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    SKYCUT extracts sky data along a cross-sectional cut.
;
;DESCRIPTION:
;    This procedure extracts and displays sky data along a cross
;    sectional cut.  Because this program can be run on both X- and 
;    non-X windows terminals, and the cut specified in a number of
;    ways on various input formats, the user should consult the 
;    EXAMPLES section below for a fuller description of this facility.
;    The first time user should probably read the entire section.        
;
;
;CALLING SEQUENCE:
;    PRO skycut, input,[nat_cube=nat_cube],[proj=proj],[coord=coord],
;                [win=win],[w_pos=w_pos],[pwin=pwin],[noplot=noplot],
;                [sclmin=sclmin],[sclmax=sclmax],[face=face_num],
;                [color=color],[n_pnts=n_pnts],[cs_plot=cs_plot],
;                [cs_titl=cs_titl]
;
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    input       I     flt arr             skycube/face/projection input
;    [nat_cube]  I     flt arr             Native skycube/face input
;                                          (unfolded or sixpack)
;    [proj]      I     string              Projection Type:
;                                          ('A' - Aitoff
;                                           'S' - Global Sinusoidal
;                                           'M' - Mollweide)
;    [coord]     I     string              Coordinate System:
;                                          ('E' - Ecliptic
;                                           'G' - Galactic
;                                           'Q' - Equatorial)
;    [win]       I     int                 Window # of image
;    [w_pos]     I     int arr             Image offset within window
;    [pwin]      I     int                 Plot window #
;    [noplot]    I     int                 noplot flag for non-X
;                                          windows terminals
;    [sclmin]    I     flt                 Image Scale Minimum
;    [sclmax]    I     flt                 Image Scale Maximum
;    [face]      I     int                 face number (face inpt only)
;    [color]     I     int                 skycut line color
;    [n_pnts]    O     int arr             # points in each plot
;    [cs_plot]   O     flt arr             storage array for plots
;    [cs-titl]   O     str arr             description array for plots
;
;WARNINGS:
;    If displaying an image, an X-windows terminal must be used.
;    It is recommended when accessing DIRBE data to use either a
;    single face or a small (512 by 256) projection.  Use of an
;    entire DIRBE sky cube or large projection will slow the
;    routine.
;
;EXAMPLES: 
;
; To extract and store a skycut from a skycube called 'cube' using a 
; non-X window terminal:
;
; skycut,cube,n_pnts=npnts,cs_plot=plts,cs_titl=titl
;
; The user will be prompted for the longitude and latitude of the end
; points of the desired cut.  SKYCUT then generates a plot with the 
; distance in degrees along the (great circle) cut on the x-axis and 
; the values of the skymap along the cut on the y-axis.  This plot is
; displayed in the graphics window of the terminal.  The user is then
; given the following options: (1) to specify another cross sectional
; cut, (2) to plot the "complement" of the cut, that is the arc 
; connecting the endpoints with length greater than 180 degress, 
; (3) to store the sky cut, or (4) to exit the program.  In the third 
; case, the relevant data is stored in the array specified by the 
; 'cs_plot' keyword which has dimensions: (720,4,16).  This data array
; contains, for each cut, the x and y "vectors" displayed in the plot
; (the (*,0,*) and (*,1,*) elements, respectively), plus the longitude
; and latitude values of the arc (the (*,2,*) and (*,3,*) elements).
; Upon exiting the program, SKYCUT truncates the output arrays to 
; minimize the number of unfilled elements.  The actual number
; of points for a particular cut is stored in an integer vector of 
; dimension 16 specified by the 'n_plts' keyword.  In addition the 
; description of the cut, containing the coordinates of the endpoints is 
; stored in the string array specified by the 'cs_titl' keyword. All
; three storage keywords must be provided on the command line otherwise
; SKYCUT will not allow plots to be saved for later access.
;
; ____
;
; To extract from a single face, 'in_face', on a non-X windows terminal,
; the user should supply the face number on the command line:
;
; skycut,in_face,face=4
;
; If the face number is not given on the command line, the user will
; be prompted for it.  Note that in this case, the plots cannot be
; stored.
;
; ____
;
; To extract from a reprojected image, for example, a galactic
; Aitoff projection called, 'gal_ait', the user must supply the 
; projection type and the coordinate system:
;
; skycut,gal_ait,proj='a',coord='g'
;
; If either of these are not provided on the command line, the user
; will be prompted for them.  The facility has no way of checking 
; whether these values actually describe the projection provided so
; care should be taken.  In addition, because the reprojection
; procedure introduces small spatial distortions into the data, the
; resulting cuts will not be as accurate as those from a skycube,
; therefore this type of invocation should be avoided if complete data
; "fidelity" is necessary.
;
;
; The three invocations above will also work on an X-windows terminal
; as long as the user gives the window for the skycut plot in the
; 'pwin' keyword.  If not provided, the user will be prompted for it.
; In the following examples, it is assumed that the user is running
; SKYCUT from an X-windows terminal.
;
; ____
;
; To suppress the display of the cross-sectional plot on a non-X 
; windows terminal, simply supply the /noplot qualifier on the
; command line.
;
; skycut,cube,/noplot
;
; ____
;
; If the user provides a displayed image (unfolded skycube, single face
; reprojection), they may use the cursor to specify the cross section.
; The user must provide two additional keywords on the command line:
; (1) 'win' which gives the window number of the displayed image and
; (2) 'w_pos' which gives the offset within the window of the image.
; This is a 2 element integer vector.  If either of the keywords is
; missing, the user will be prompted for them.
;
; To extract from a skycut from a skycube called, 'cube':
;
; skycut,cube,win=0,w_pos=[32,32]
;
; The image of 'cube' should already exist with window 0.  This can be
; accomplished with the commands: 'WINDOW,0' and 'TVSCL,cube,32,32'.
; 
; To specify a cross section, the user should 'point and click' within
; the image window using the left button of the mouse.  The arc will
; be displayed on the window overlayed on the image.  The user can also
; specify the cut using the longitude and latitude of the arc endpoints.
;
; For single face and reprojection images enter, respectively:
;
; skycut,in_face,face=4,win=1,w_pos=[0,0]
;
; skycut,gal_ait,proj='a',coord='g',win=2,w_pos=[16,32]
;
; If the user supplies the storage keywords, the cuts can be saved
; for latter use.
;
; ____
;
; If the user has access to both a skycube (or single face) and a
; reprojected image of it, SKYCUT allows the user to specify the
; cross section using the reprojected image for better orientation,
; while extracting the data from the cube (or face) thereby retaining
; complete accuracy.
;
; The invocation is as follows:
;
; skycut,gal_ait,proj='a',coord='g',win=2,w_pos=[16,32],nat_cube=cube
;
; The 'nat_cube' keyword contains the name to the skycube from which
; the galactic Aitoff image, 'gal_ait', was created.  The 'win' and
; 'w_pos' keywords refer to the window containing the reprojection.
;
; ____
;
;               Rules for specifying longitude and latitude
;               -------------------------------------------
;
;  The longitude and latitude for a given point are specified by
;  two numbers separated by spaces or commas.
;
;  To specify a cross section with endpoints 30 long -25 lat and
;  260 long 23.5 lat:
;
;  First point : 30 -25
;  Second point: 260  23.5
;
;
;  To specify a cross section along the 42nd parallel from 20 long
;  to 180 long:
;
;  First point : 20 42.0
;  Second point: 180
;
;
;  To specify a cross section along the 42nd parallel from 270 long
;  to 10 long:
;
;  First point : 270 42.0
;  Second point: 10
;
;  Latitude cross-sections go in the direction of increasing longitude.
;
;
;  To specify a cross section along the 90th meridian from pole to pole:
;
;  First point : 90
;  Second point: <CR>
;
;
;#
;COMMON BLOCKS:
;    image_parms,arc
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Get input parameters.
;    Get entry mode (cursor of lon/lat).
;    Specify arc, get data along arc (getarc).
;    Display arc along image (if X-windows terminal)
;    Display cross-sectional plot.
;    Draw arc on input.
;    Get next command.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: getinpar, getemode, getarc, drawcirc,
;                        getnext
;
;MODIFICATION HISTORY
;
;    Written by J.M. Gales, Applied Research Corp.   May 92
;
; SPR 9805 7/2/92 JMG
; Eliminate KEYWORD_SET check of storage arrays.
; Eliminate 'nostore' flag.
;
; SER 9953     Sixpack Input Capabilities
;              Sixpacked skycube can now be used as input
; 9-SEP-92     J.M. Gales
;
; SER 9954     Incorporation of GETNEXT into code
; 9-SEP-92     J.M. Gales
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;
;.TITLE
;Routine SKYCUT
;-
;
COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
		    cube_side, proj0, coord0, sz_proj

COMMON arc, x_plot,y_plot,win_title,fij_input,n_pix_arc,i_arc,j_arc, $
            vec_arc,end_pixs


reply = ' '
mask = 0
noinpt = 0
n_plot = 0

If (N_ELEMENTS(sclmin) EQ 0) THEN sclmin = min(input)
If (N_ELEMENTS(sclmax) EQ 0) THEN sclmax = max(input)

IF (((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) OR (N_ELEMENTS(win) EQ 0)) $
     THEN noshow = 1 ELSE noshow = 0

IF (((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) AND (N_ELEMENTS(pwin) EQ 0)) $
     THEN READ,'Please enter Plot window: ',pwin

IF ((noshow EQ 0) AND (N_ELEMENTS(w_pos) EQ 0)) THEN BEGIN
	w_pos = intarr(2)
	READ,'Please enter image offsets within window. (Ex: 32 32): ',w_pos
ENDIF

IF (KEYWORD_SET(color) EQ 0) THEN color = 255

IF (KEYWORD_SET(noplot) EQ 0) THEN noplot = 0

n_pnts = INTARR(16)
cs_plot = FLTARR(768,4,16)
cs_titl = STRARR(16)

sz = SIZE(input)
IF (sz(1) EQ 2*sz(2)) THEN BEGIN

	IF (KEYWORD_SET(proj) EQ 0) THEN BEGIN
		entry = umenu(['Projection Type','Aitoff', $
        	               'Global Sinusoidal', 'Mollweide'], $
                	        title=0,init=1)

		proj = STRMID('ASM',entry-1,1)
	ENDIF

	IF (KEYWORD_SET(coord) EQ 0) THEN BEGIN
		entry = umenu(['Projection coordinate System', $
			       'Ecliptic', 'Galactic', 'Equatorial'], $
                	        title=0,init=1)

		coord = STRMID('EGQ',entry-1,1)
	ENDIF

ENDIF

IF (KEYWORD_SET(coord)) THEN coord0 = STRUPCASE(STRMID(coord,0,1))
IF (KEYWORD_SET(proj)) THEN proj0 = STRUPCASE(STRMID(proj,0,1))

ret_stat = getinpar(input,nat_cube,res,re_proj,psize,mask,face_num)
IF (ret_stat NE 0) THEN GOTO, exit


e_mode:

arc_type = 'S'
getemode, entry_mode,noshow,coord0,infmt,re_proj


loop:

IF (noshow EQ 0) THEN WSET, win

ret_stat = getarc(entry_mode,input,nat_cube,mask,res,face_num,w_pos, $
                  re_proj,arc_type,lon_flag,lat_flag,noinpt,color)
IF (ret_stat NE 0) THEN GOTO, exit


IF (noshow EQ 0) THEN BEGIN
	ERASE
	TVSCL, sclmax<input>sclmin,   w_pos(0),w_pos(1) 
	drawcirc, fij_input, n_pix_arc, i_arc, j_arc, vec_arc, $
		  re_proj, w_pos, face_num, color
ENDIF


IF ((!D.NAME EQ 'X') or (!D.NAME EQ 'WIN')) THEN BEGIN
	WINDOW,pwin,XSIZE=320,YSIZE=256,TITLE=win_title
	PLOT,x_plot,y_plot,TITLE=win_title
ENDIF ELSE IF (noplot EQ 0) THEN BEGIN
	PLOT,x_plot,y_plot,POSITION=[0.25,0.25,0.75,0.75],TITLE=win_title
	PRINT,'Press <CR> to continue'
	READ,reply
	PRINT,STRING(27b) + '2'		; return to text screen (GRAPH-ON)
ENDIF


str2 = 'This operation not available for single faces and latitude cuts.'
str3 = 'Cursor entry mode not allowed for this invocation or terminal.'

wish:	wish = umenu(['Do you wish to:', $
		      'Specify new cross-section?', $
		      'View "complementary" cross-section?', $
		      'Change Entry Mode?', $
		      'Save cross-sectional plot?', $
		      'Exit procedure?', $
		      'HELP'],title=0,init=1)

	CASE wish OF

	1 :	BEGIN
		arc_type = 'S'
		GOTO, loop
		END

	2 :	BEGIN

		IF ((lat_flag EQ 0) AND (lon_flag EQ 0) AND $
		(face_num EQ -1)) THEN BEGIN
			IF (arc_type EQ 'S') THEN arc_type = 'L' $
                        ELSE arc_type = 'S'
			noinpt = 1
			GOTO, loop
		ENDIF ELSE BEGIN
			PRINT,str2
			GOTO, wish
		ENDELSE

		END

	3 :	BEGIN

		IF (noshow EQ 1) THEN BEGIN 
			PRINT,str3
			GOTO, wish
		ENDIF ELSE GOTO, e_mode

		END

	4 :	BEGIN

		n_pnts(n_plot) = n_pix_arc
		PRINT,STRCOMPRESS('Plot'+STRING(n_plot+1)+' Stored')
		ll = coorconv(vec_arc,infmt='U',outfmt='L')
		FOR i=0,n_pix_arc-1 DO BEGIN
			cs_plot(i,0,n_plot) = x_plot(i)
			cs_plot(i,1,n_plot) = y_plot(i)
			cs_plot(i,2,n_plot) = ll(i,0)
			cs_plot(i,3,n_plot) = ll(i,1)
		ENDFOR
		cs_titl(n_plot) = win_title
		n_plot = n_plot + 1
		GOTO, wish
		END

	5 :	GOTO, exit

	6 :	BEGIN
		PRINT,'No help available at this time'
		GOTO, wish
		END

	ENDCASE



exit:

IF (n_plot GT 0) THEN BEGIN
	max_n = MAX(n_pnts)
	n_pnts = n_pnts(0:n_plot-1)
	cs_plot = cs_plot(0:max_n,*,0:n_plot-1)
	cs_titl = cs_titl(0:n_plot-1)
ENDIF


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


