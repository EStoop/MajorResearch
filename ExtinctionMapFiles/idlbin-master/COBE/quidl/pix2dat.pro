function pix2dat, pixel=pixel,x_in=x_in,y_in=y_in,raster=raster
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    pix2dat extracts data for a given set of raster coordinates.
;
;DESCRIPTION:
;    This IDL function creates a data array given either a list of
;    pixels or a set of x and y raster coordinates and a raster image
;    (sky cube or face).  The skycube can be in either unfolded or six
;    pack formate.  This routine is the "complement" to PIX2XY. 
;    The program assumes a right oriented, ecliptic coordinate input 
;    raster.
;
;CALLING SEQUENCE:
;    data =  pix2dat(pixel=pixel,raster=raster) or
;    data =  pix2dat(x_in=x_in,y_in=y_in,raster=raster)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    pixel       I     int/long arr        Input pixel list
;    x_in        I     int/long arr        x coord list
;    y_in        I     int/long arr        y coord list
;    raster      I     int/flt/dbl arr     input raster array
;    data        O     int/flt/dbl arr     output data array
;                                          
;    Either the pixel keyword or the x_in AND y_in keywords must
;    be specified.
;
;    The x and y coordinates for a single face are those within the
;    face, not those for the face positioned within a unfolded sky
;    cube.  For example, for a FIRAS face the coordinates range from
;    0 to 31.
;
;WARNINGS:
;
;EXAMPLE: 
;
; To extract the data from a raster for all pixels with galactic 
; latitude less than 30 degrees:
;
; pixel = indgen(6144)
;
; ll = coorconv(pixel,infmt='p',inco='f',outfmt='l',outco='g')
;
; pix30 = pixel(where(abs(ll(*,1)) lt 30.0))
;
; dat30 = pix2dat(pixel=pix30,raster=sky_cube)
;
; Alternatively, we can use x and y coordinate arrays generated with
; pix2xy:
;
; pix2xy,pix30,x30,y30,res=6
;
; dat30 = pix2dat(x_in=x30,y_in=y30,raster=sky_cube)
;
;
;
; To extract from a face:
;
; face4 = sky_cube(0:31,32:63) ; extract face 4 from sky cube
;
; fij = pix2fij(pix30,6)	; get face and coordinates
;
; pix4 = pix30(where(fij(*,0) eq 4)) ; extract face 4 pixels
;
; dat4 = pix2dat(pixel=pix4,raster=face4)
;
;
; To extract along a galactic meridian of a DIRBE map:
;
; ll = fltarr(361,2)	; allocate long-lat array
;
; ll(*,0) = 0.		; chose central meridian 
;
; ll(*,1) = 90 - findgen(361)*0.5 ; generate points along meridian at
;                                   half degree intervals.
;
; pix = coorconv(ll,infmt='l',inco='g',outfmt='p',outco='b')
;			; get pixels along meridian
;
; dat = pix2dat(pixel=pix,raster=band_01)
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Uses PIX2FIJ to "decode" pixel list into face,i,j vector.
;    Fills data array from i,j vector and raster array
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: pix2fij, pix2xy
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   May 92
;
; SER 9857  2-AUG-1992  JMG
; Input skycube files can now be in sixpacked format.
;
; SER 9983    Call PIX2XY to generate column & row arrays (x_in,y_in)
;             (via CRASTER c module)
; 14-SEP-92   JMG;
;
; SPR 10346   Remove obsolete 'no_c' keyword in call to PIX2XY
; 15-DEC-92   JMG
;
;.TITLE
;Routine PIX2DAT
;-


; Get size of input raster
; ------------------------
sz = SIZE(raster)
input_l = sz(1)
input_h = sz(2)

IF (input_l EQ input_h) THEN cube_side = input_l
IF (3*input_l EQ 4*input_h) THEN cube_side = input_l/4
IF (2*input_l EQ 3*input_h) THEN cube_side = input_l/3


; Determine resolution of quad cube
; ---------------------------------
FOR bit=0,15 DO BEGIN
	IF ((cube_side XOR 2^bit) EQ 0) THEN BEGIN
		res = bit+1
		GOTO, lbl_1
	ENDIF
ENDFOR
MESSAGE,'Improper Image Size',/cont
goto, exit

lbl_1:


; Determine number of pixels / Get column & row numbers if pixel entry
; --------------------------------------------------------------------
IF (N_ELEMENTS(pixel) NE 0) THEN BEGIN
	num_pix = N_ELEMENTS(pixel)

	IF (input_l EQ input_h) THEN $		; face
	pix2xy,pixel,x_in,y_in,res=res,/face

	IF (3*input_l EQ 4*input_h) THEN $	; skycube
	pix2xy,pixel,x_in,y_in,res=res

	IF (2*input_l EQ 3*input_h) THEN $	; sixpack
	pix2xy,pixel,x_in,y_in,res=res,/sixpack
	
ENDIF ELSE BEGIN

	num_pix = N_ELEMENTS(x_in)
	IF (N_ELEMENTS(x_in) NE N_ELEMENTS(y_in)) THEN BEGIN
		str = 'Column and Row arrays have incompatible sizes'
		MESSAGE,str,/CONT
	ENDIF

ENDELSE


; Build data array
; ----------------
ras_sz = SIZE(raster)
IF (ras_sz(0) EQ 2) THEN num_ras = 1 ELSE num_ras = ras_sz(3)

arr_type = ['BYT','INT','LON','FLT','DBL','COMPLEX']
arr = arr_type(ras_sz(ras_sz(0)+1)-1)
stat = EXECUTE('data = ' + arr + 'arr(num_pix,num_ras)')


; Load data array
; ---------------
IF (num_ras EQ 1) THEN data = raster(x_in,y_in) ELSE BEGIN

	FOR i=0,num_ras-1 DO BEGIN
		ras1 = raster(*,*,i)
		data(*,i) = ras1(x_in,y_in)
	ENDFOR

ENDELSE

exit:

RETURN, data
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


