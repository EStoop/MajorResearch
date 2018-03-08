PRO batchcut, input, coord=coord,proj=proj,ll_in=ll_in, $
              arc_flag=arc_flag,lat_flag=lat_flag,face=face_num, $
              n_pnts=n_pnts,cs_plot=cs_plot,cs_titl=cs_titl
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    BATCHCUT extracts skycuts defined by a lon/lat array.
;
;DESCRIPTION:
;    This procedure extracts and stores skycuts defined by a user-
;    supplied array containing the longitude and latitude of the
;    endpoints.  Both great-circle and constant latitude cuts can
;    be specified as well as short (< 180 degrees) and long (> 180)
;    arcs.  The input array can be either a native skycube (unfolded
;    or sixpacked format), a single face, or a reprojection.
;
;
;CALLING SEQUENCE:
;    PRO batchcut, input,coord=coord,[proj=proj],ll_in=ll_in,
;                  [lat_flag=lat_flag],[arc_flag=arc_flag],
;                  [face=face_num],[n_pnts=n_pnts],[cs_plot=cs_plot],
;                  [cs_titl=cs_titl]
;
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    input       I     flt arr             skycube/face/projection input
;    coord       I     string              Coordinate System:
;                                          ('E' - Ecliptic
;                                           'G' - Galactic
;                                           'Q' - Equatorial)
;    [proj]      I     string              Projection Type:
;                                          ('A' - Aitoff
;                                           'S' - Global Sinusoidal
;                                           'M' - Mollweide)
;    ll_in       I     flt arr (2,2*n)     Lon/lat input array
;    [lat_flag]  I     int arr (n elem)    Latitude cut flag
;    [arc_flag]  I     int arr (n elem)    Arc type flag
;    [face]      I     int                 face number (face inpt only)
;    [n_pnts]    O     int arr             # points in each plot
;    [cs_plot]   O     flt arr             storage array for plots
;    [cs-titl]   O     str arr             description array for plots
;
;WARNINGS:
;    Care be taken not to specify points diametrically opposed as the
;    great circle is then not uniquely defined.
;
;EXAMPLES: 
;
; The skycut specification array, "ll_in", is a 2 by 2*n array, where
; n is the number of cuts.  The first column contains the longitudes
; in degrees (0-360), the second contains the latitudes (-90 - 90).
;
; To specify three cuts along the longitude lines -20, 0, 20 from -25
; to 25 latitude, enter:
;
; IDL> ll_in = FLTARR(2,6)
; IDL> READ, ll_in
; : -20 -25 -20 25  0 -25 0 25  20 -25 20 25 <CR>
;
; To extract and store these cuts from a input skycube, 'sky_input',
; (in either unfolded or sixpacked format), enter at the IDL prompt:
;
; batchcut,sky_input,ll_in=ll_in,coord='g', $
;          n_pnts=npnts,cs_plot=plts,cs_titl=titl
;
;
; The 'coord' keyword tells the procedure that these endpoints are in
; galactic coordinates.  The retrieved data is stored in the array
; specified by the 'cs_plot' keyword.  This data array contains,
; for each cut, the x and y "vectors" displayed in the plot
; (the (*,0,*) and (*,1,*) elements, respectively), plus the longitude
; and latitude values of the arc (the (*,2,*) and (*,3,*) elements).
; Upon exiting the program, SKYCUT truncates the output arrays to 
; minimize the number of unfilled elements.  The actual number
; of points for a particular cut is stored in an integer vector 
; specified by the 'n_plts' keyword.  In addition, the description of
; the cut, containing the coordinates of the endpoints is stored in
; the string array specified by the 'cs_titl' keyword. All three
; storage keywords must be provided on the command line otherwise
; SKYCUT will not allow plots to be saved for later access.
;
; If the input array is a single face then the face number must be
; specified through the 'face' keyword.  If the input array is a
; reprojection, then the 'proj' keyword must be specified and the
; 'coord' keyword must correspond to the coordinate system of the
; projection.  Cuts from reprojections are designated by the string
; 'PR' in the title, all others with a 'SC'.
;
;
; "Long" Arcs
; -----------
; By default great circle skycuts use the smaller of the two great
; circles defined by the endpoints.  To extract the larger of the two,
; the user passes an 'n' element integer array via the 'arc_flag' 
; keyword with a zero specifying the default "short" arc and a '1',
; the "long" ( > 180 deg) arc.
;
; For example, if the user desires a cut along the 0 longitude
; excluding points within 25 degrees of the galactic plane, as well
; as the two galactic cuts at -20 and 20 longitude, they would enter:
; 
; batchcut,sky_input,ll_in=ll_in,coord='g',arc_flag=[0,1,0], $
;          n_pnts=npnts,cs_plot=plts,cs_titl=titl
;
;
; Latitude Cuts
; -------------
; To specify a latitude cut, the user supplies an 'n' element integer
; array passed through the 'lat_flag' keyword where '0' specifies the
; default great circle arc and '1' a latitude cut.  To specify a cut
; along the ecliptic latitude 15 from 90 longitude to 270 longitude
; in addition to a regular great circle cut (specifed as the third and
; fourth row of 'll'), enter:
; 
; IDL> ll = [[90.,15],[270,15],[10,-45],[340,36]]
;
; batchcut,sky_input,ll_in=ll,coord='e',lat_flag=[1,0], $
;          n_pnts=npnts,cs_plot=plts,cs_titl=titl
;
;
; The latitude arc goes from 90 longitude to 270.  To specify the
; arc from 270 longitude to 90 simply exchange the numbers in 'll':
;
; IDL> ll = [[270.,15],[90,15],[10,-45],[340,36]]
;
;#
;COMMON BLOCKS:
;    image_parms
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Get input parameters.
;    Get pixel numbers and unit vectors along arc
;    Get data along arc
;    Store data
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: getinpar,grtcirc,pix2fij,coorconv
;            
;
;MODIFICATION HISTORY
;
;    Written by J.M. Gales, Applied Research Corp.   Sep 92
;
;.TITLE
;Routine BATCHCUT
;-
;
COMMON image_parms, input_l,input_h, cur_2_face, offx, offy, $
		    cube_side, proj0, coord0, sz_proj

d2r = !pi / 180

re_proj = 0

IF (N_ELEMENTS(input) EQ 0) THEN BEGIN
	MESSAGE,'Input array not specified or does not exist',/CONT
	GOTO, bail
ENDIF

n = N_ELEMENTS(ll_in) / 4
		; get number of cuts specified

IF (n EQ 0) THEN BEGIN
	MESSAGE,'LL_IN array not specified or does not exist',/CONT
	GOTO, bail
ENDIF


IF (N_ELEMENTS(lat_flag) EQ 0) THEN lat_flag = INTARR(n)
		; set lat_flag to 0 vector if not present

IF (N_ELEMENTS(arc_flag) EQ 0) THEN arc_flag = INTARR(n)
		; set arc_flag to 0 vector if not present


arc = ['S','L']		; arc type array


n_pnts = INTARR(n)
cs_plot = FLTARR(768,4,n)
cs_titl = STRARR(n)
		; allocate storage arrays


sz = SIZE(input)
IF (sz(1) EQ 2*sz(2)) THEN BEGIN
	re_proj = 2
	IF (KEYWORD_SET(proj) EQ 0) THEN BEGIN
		MESSAGE,'Projection Type Must be Supplied',/CONT
		GOTO, exit
	ENDIF
ENDIF
		; if reprojection cut then re_proj value
		; make sure proj keyword is supplied



IF (KEYWORD_SET(coord) EQ 0) THEN BEGIN
	MESSAGE,'Coordinate System Must be Supplied',/CONT
	GOTO, exit
ENDIF
		; make sure coord keyword is supplied



ll_in0 = TRANSPOSE(FLOAT(ll_in))
		; "rotate" ll_in0 to proper orientation for coorconv



; Adjust "problem" points
; -----------------------

i = where(ll_in0 EQ 0.0)
IF (i(0) NE -1) THEN ll_in0(i) = ll_in0(i) + 1.e-4

i = where(ll_in0(*,1) EQ 90.0)
IF (i(0) NE -1) THEN ll_in0(i,1) = ll_in0(i,1) - 1.e-4

i = where(ll_in0(*,1) EQ -90.0)
IF (i(0) NE -1) THEN ll_in0(i,1) = ll_in0(i,1) + 1.e-4




ret_stat = getinpar(input,nat_cube,res,re_proj,psize,mask,face_num)
		; get input parameters

out_res = STRCOMPRESS('R' + STRING(res),REMOVE_ALL=1)
		; build out_res string

end_pix = coorconv(ll_in0,infmt='L',inco=coord,outfmt='P',outco=out_res)
		; convert endpoints for lon/lat to pixel number




; Get pixels and unit vectors along arc
; -------------------------------------

FOR n_plot=0,n-1 DO BEGIN	; loop through cuts

	arc_deg = 0.5

	IF (lat_flag(n_plot) EQ 0 AND $
	    ll_in0(2*n_plot,0) NE ll_in0(2*n_plot+1,0)) THEN BEGIN

		arc_type = arc(arc_flag(n_plot))

		ep = end_pix(2*n_plot:2*n_plot+1)

		pixel_arc = grtcirc(ep,vec_arc,input_type='P', $
			    res=res,arc_inc=arc_deg,arc_type=arc_type, $
			    arc_actl)
				; get pixel #'s along great circle


		vec_arc = coorconv(vec_arc,infmt='U',inco='E',outco=coord)
				; convert path in unit vectors 
				; from ecliptic (skycube system)
				; to user specified system

	ENDIF


	IF (lat_flag(n_plot) EQ 0 AND $
	    (ll_in0(2*n_plot,0) EQ ll_in0(2*n_plot+1,0))) THEN BEGIN

		lon_rad = d2r * FLOAT(ll_in0(2*n_plot,0))
		IF (lon_rad EQ 0.0) THEN lon_rad = d2r*1.e-6

		arc_len = ABS(ll_in0(2*n_plot,1)-ll_in0(2*n_plot+1,1))
		lat_min = MIN([ll_in0(2*n_plot,1),ll_in0(2*n_plot+1,1)])
		n_arc = arc_len / 0.5
		arc_actl = arc_len / n_arc
		lat_rad = d2r*(arc_actl*indgen(n_arc+1)+lat_min)

		vec_arc = fltarr(n_arc+1,3)
		vec_arc(*,0) = cos(lon_rad)*cos(lat_rad)
		vec_arc(*,1) = sin(lon_rad)*cos(lat_rad)
		vec_arc(*,2) = sin(lat_rad)

		pixel_arc = coorconv(vec_arc,infmt='U',inco=coord, $
				     outfmt='P',outco=out_res)

	ENDIF



	IF (lat_flag(n_plot) NE 0) THEN BEGIN	; latitude cut

		long_1 = ll_in0(2*n_plot,0)
		long_2 = ll_in0(2*n_plot+1,0)
		IF (long_2 LT long_1) THEN long_2 = long_2 + 360
		lat_rad = d2r*ll_in0(2*n_plot,1)
		lat_len = (long_2-long_1) * cos(lat_rad)
		n_arc = FIX(ABS(lat_len) / arc_deg) + 1
		del_lat = (long_2-long_1) / n_arc
		vec_arc = FLTARR(n_arc+1,3)
		arc_actl = del_lat

		FOR i=0,n_arc DO BEGIN
		  vec_arc(i,0) = cos(d2r*(long_1+(i*del_lat)))*cos(lat_rad)
		  vec_arc(i,1) = sin(d2r*(long_1+(i*del_lat)))*cos(lat_rad)
		  vec_arc(i,2) = sin(lat_rad)
		ENDFOR

		pixel_arc = coorconv(vec_arc,infmt='U',inco=coord, $
				     outfmt='P',outco=out_res)


	ENDIF



; Extract data along arc
; ----------------------

n_pix_arc = n_elements(pixel_arc)

IF (re_proj EQ 2) THEN BEGIN	; reprojection cut

	proj_xy = uv2proj(vec_arc,proj,sz_proj)
	i_arc = proj_xy(*,0)
	j_arc = proj_xy(*,1)
	win_title = 'PR  '
			; determine projection coordinates along arc

ENDIF ELSE BEGIN		; cube cut

	fij_input = pix2fij(pixel_arc,res)

	i_arc = offx(fij_input(*,0)) * cube_side + fij_input(*,1)
	i_arc = input_l - (i_arc+1) 
	j_arc = offy(fij_input(*,0)) * cube_side + fij_input(*,2)

	win_title = 'SC  '
			; determine skycube coordinates along arc

ENDELSE


; Build window title
; ------------------
IF (lat_flag(n_plot) NE 0) THEN win_title = win_title + ' (lat) '

ll_pos = coorconv([vec_arc(0,*),vec_arc(n_pix_arc-1,*)], $
		  infmt='U',outfmt='L')


w_str = string(ll_pos,format='(f6.1)')
win_title = win_title + $
  STRCOMPRESS('(' +w_str(0)+ ',' +w_str(2)+ ')',REMOVE_ALL=1) + '  to  ' + $
  STRCOMPRESS('(' +w_str(1)+ ',' +w_str(3)+ ')',REMOVE_ALL=1) + ' ' + coord

d_mask = intarr(n_pix_arc)
IF (face_num EQ -1) THEN d_mask(*) = 1 ELSE d_mask = (fij_input(*,0) $
					    EQ face_num)
				; mask out points not in active face

x_plot = FINDGEN(n_pix_arc)*arc_actl
IF (re_proj EQ 1) THEN y_plot = nat_cube(i_arc,j_arc)*d_mask $
		  ELSE y_plot = input(i_arc,j_arc)*d_mask
				; extract data along path

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

ENDFOR

exit:

n_pnts = n_pnts(0:n-1)
cs_plot = cs_plot(0:MAX(n_pnts)-1,*,0:n-1)
cs_titl = cs_titl(0:n-1)
		; truncate storage arrays

bail:

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


