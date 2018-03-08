function buildmap,file=in_file,fields=fields,face=face,multobs=multobs, $
                  x_out=x_out,y_out=y_out,pixel=pixel,data=data, $
                  real=real,imag=imag,sixpack=sixpack,sent_val=sent_val
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    BUILDMAP provides user-friendly access to COBE archive maps.
;
;DESCRIPTION:
;    This IDL function creates an array of rasterized images (sky
;    cubes or faces) using the READ_SKYMAP and PIX2XY facilities.
;    If not supplied on the command line, the user is prompted for
;    the skymap data file name (including the archive directory), the 
;    data fields to be extracted (up to 7 if the AV1 or AV2 switches
;    are specified, up to 8 otherwise), and the face number (-1 if
;    entire cube).  The output consists of the raster image array,
;    and if desired, the x and y raster coordinate arrays and the 
;    pixel and data arrays (from READ_SKYMAP).
;
;    For DIRBE skymaps, this facility can handle multiple observations
;    per pixel in a number of ways.  They are (1) a simple average,
;    (2) an average for a given approach vector value, (3) the first 
;    observation for the pixel, or (4) the last observation for the
;    pixel.
;
;
;CALLING SEQUENCE:
;    outmap = buildmap([file=file],[fields=fields],[face=face], $
;                      [multobs=multobs],[x_out=x_out],[y_out=y_out], $
;                      [pixel=pixel],[data=data],/real,/imag,/sixpack)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    [outmap]    O     int/flt/dbl/cpx arr Rasterized output array
;    [file]      I     str                 Skymap file
;    [fields]    I     str                 Skymap data fields
;    [face]      I     int                 Skymap face
;    [multobs]   I     str                 Multuple observations
;                                          per pixel qualifier
;    [x_out]     O     int arr             x-coordinate array
;    [y_out]     O     int arr             y-coordinate array
;    [pixel]     O     int arr             returned pixel vector
;    [data]      O     int/flt/dbl/cpx arr Returned data list
;    [real]      I     qualifier           takes real part of complex
;                                          data
;    [imag]      I     qualifier           takes imag part of complex
;                                          data
;    [sixpack]   I     qualifier           sixpack output switch
;    [sent_val]  I     float               Sentinal value upper bound
;
;
;    Multiple Observations per Pixel (MULOPTS) Values
;    ------------------------------------------------
;
;    'avg'    - Average multiple values in pixel
;    'av1'    - Average multiple values in pixel with approach vector 1
;    'av2'    - Average multiple values in pixel with approach vector 2
;    'first'  - Use first observation in pixel
;    'last'   - Use last observation in pixel (DEFAULT)
;
;
;WARNINGS:
;    Data fields must be all of the same data type.  BUILDMAP uses
;    the first data field specified to determine the data type.
;
;    BUILDMAP uses a lot of memory.  The user should free up as
;    much space as possible before running.
;
;    If any non-recoverable error occurs in BUILDMAP, the outmap
;    return variable will be set to 0.
;
;    Because of the special nature of the ADT times (I*8), the user
;    should not use BUILDMAP to produce an averaged value.
;
;EXAMPLE: 
;
; To build an array of DIRBE skymaps for the photometry channels 1, 2
; and 7, from a weekly skymap file for approach vector 1 data:
;
; outarr = buildmap(file='dirbe_edit:bpw_wf.week_avg_913611607', $
;                   fields='photometry(1:2),photometry(7)', $
;                   face=-1,multobs='av1',pixel=pixout)
;
; The output array, 'outarr' will have dimensions (1024,768,3).  The
; pixel list is stored in the vector, 'pixout'.
;
; Note that the FORTRAN convention is used (subscripts starting with 1)
; when specifying array fields.
;
;
; To build a high frequency FIRAS (real) spectral skymap in sixpack
; format:
;
; spcmap = buildmap(file='firsky:fcs_ccmsp_rh.ed_8934302_9026408', $
;                   fields='chan.spec(5:171)',face=-1,/real/sixpack)
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Gets Skymap file name, fields, and face number if not supplied.
;    Gets resolution using GET_RESOLUTION.
;    Reads data using READ_SKYMAP.
;    Combines separate data arrays into single array.
;    If DIRBE then combine multiple observations.
;    Rasterizes data using PIX2XY.
;
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: get_resolution, read_skymap, pix2xy
;                        read_flds, pixavg (c routine)
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Apr 92
;
; SPR 9780  Temporary fix for IDL TRANSPOSE function
; 19-JUN-92 JMG
;
; SER 9782  Addition of /real and /imag qualifiers
; 22-JUN-92 JMG
;
; SER 9787  Use auxilary fields variable when calling READ_SKYMAP.
;           This will eliminate changes to input fields variable
;           when appending approach vector field for DIRBE reads.
; 24-JUN-92 JMG
;
; SER 9833  Allow user to designate sixpack output using sixpack
;           qualifier.
; 24-JUL-92 JMG
;
; SER 9981    Implement PIXAVG (C module pixel averaging)
; 9-SEP-1992  JMG
;
; SPR 10827   Remove DIRBE sentinal values from averaging.
; 16-APR-1993 JMG
;
; SPR 10962   Add SENT_VAL keyword to pass to pixavg C module
; 27-MAY-1993 JMG
;
; SPR 10914  4-May-1994 Improper error handling for MULTOBS input. J. Newark
;
;.TITLE
;Routine BUILDMAP
;-
;
arr_type = ['BYT','INT','LON','FLT','DBL','COMPLEX']

IF (N_ELEMENTS(in_file) EQ 0) THEN in_file = ''
IF (NOT KEYWORD_SET(sixpack)) THEN sixpack = 0
IF (NOT KEYWORD_SET(sent_val)) THEN sent_val = -1.e-35

get_sky:

read_sky = 0
IF (in_file EQ '') THEN BEGIN
	PRINT, ' '
	read_sky = 1
	in_dir = ''
	READ, 'Skymap Directory or Archive <Default>: ',in_dir

	IF ((in_dir NE '') AND $
            (STRMID(in_dir,STRLEN(in_dir)-1,1) NE ']') AND $
            (STRMID(in_dir,STRLEN(in_dir)-1,1) NE ':')) THEN $
	in_dir = in_dir + ':'

	PRINT,' '
	READ, 'File Name (<CR> to exit BUILDMAP): ',in_file
	IF (in_file EQ '') THEN BEGIN
		outmap = 0
		GOTO, exit
	ENDIF
	PRINT,' '
	in_file = in_dir + in_file
ENDIF

PRINT, 'Getting Resolution of Skymap'
stat = EXECUTE('res = get_resolution(in_file)')

IF (stat EQ 0) THEN BEGIN
	MESSAGE,'Skymap file not found',/CONT
	IF (read_sky EQ 0) THEN BEGIN
		outmap = 0
		GOTO, exit
	ENDIF ELSE BEGIN
		in_file = ''
		GOTO,get_sky
	ENDELSE
ENDIF

cube_side = 2l ^ (res-1)
PRINT, ' '

IF (N_ELEMENTS(face) EQ 0) THEN BEGIN
	face = umenu(['Cube/Face Selection','Entire Cube', $
                      'Face 0 - North Ecliptic Pole', $
                      'Face 1 - South Galactic Pole', $
                      'Face 2 - Galactic Anti-Center', $
                      'Face 3 - North Galactic Pole', $
                      'Face 4 - Galactic Center', $
                      'Face 5 - South Ecliptic Pole'], $
		      title=0,init=1)

	face = face - 2
ENDIF


IF (face EQ -1) THEN BEGIN
	first_pixel = 0
	num_pix = 6 * (cube_side ^ 2)
ENDIF ELSE BEGIN
	first_pixel = face * (cube_side ^2)
	num_pix  = cube_side ^ 2
ENDELSE


IF (KEYWORD_SET(multobs) EQ 0) THEN multobs = 'LAS'

multobs = STRUPCASE(STRMID(multobs,0,3))

CASE multobs OF
	'AV1' : max_flds = 7
	'AV2' : max_flds = 7
	'AVG' : max_flds = 8
	'FIR' : max_flds = 8
	'LAS' : max_flds = 8
	ELSE : BEGIN
		MESSAGE,'Unknown Value of MULTOBS Qualifier,/CONT
		outmap = 0
		GOTO, exit
		END
ENDCASE



IF (N_ELEMENTS(fields) EQ 0) THEN BEGIN

	read_flds, max_flds,n_flds,fld_rd

ENDIF ELSE BEGIN

	n_flds = 1
	FOR i=0,STRLEN(fields)-1 DO $
		IF (STRMID(fields,i,1) EQ ',') THEN n_flds = n_flds + 1
	fld_rd = fields

ENDELSE

get_flds: 


IF (n_flds EQ 0) THEN BEGIN
	outmap = 0
	GOTO, exit
ENDIF



IF (multobs EQ 'AV1' OR multobs EQ 'AV2') THEN $
    fld_rd = fld_rd + ',approach_vector'



data_str = ''
dat_str = STRARR(8)
FOR i=0,n_flds-1 DO BEGIN
	dat_str(i) = STRCOMPRESS('dat' + STRING(i),/remove_all)
	data_str = data_str + ',' + dat_str(i)
ENDFOR
IF (multobs EQ 'AV1' OR multobs EQ 'AV2') THEN data_str = data_str + ',av'

PRINT,'Reading Data From Skymap File'
stat = EXECUTE('status = read_skymap(in_file,fld_rd,pixel' + data_str + $
	       ',startpixel=first_pixel,count=num_pix)')

IF (stat EQ 0) THEN BEGIN
	PRINT, ' '
	read_flds, max_flds,n_flds,fld_rd
	GOTO, get_flds
ENDIF


tot_row = 0
n_row = intarr(8)
FOR i=0,n_flds-1 DO BEGIN

	IF (KEYWORD_SET(real) NE 0) THEN $
	stat = EXECUTE(dat_str(i) + '= FLOAT(' + dat_str(i) + ')')

	IF (KEYWORD_SET(imag) NE 0) THEN $
	stat = EXECUTE(dat_str(i) + '= IMAGINARY(' + dat_str(i) + ')')

	stat = EXECUTE('sz = SIZE(' + dat_str(i) + ')')
	CASE sz(0) OF
		1: n_row(i) = 1

		2: BEGIN

			n_row(i) = sz(1)
;			stat = EXECUTE(dat_str(i) + '=TRANSPOSE(' + $
;			               dat_str(i) + ')')

; Begin IDL TRANSPOSE bug fix

		CASE sz(sz(0)+1) OF
			1: temp = bytarr(sz(2),sz(1))
			2: temp = intarr(sz(2),sz(1))
			3: temp = lonarr(sz(2),sz(1))
			4: temp = fltarr(sz(2),sz(1))
			5: temp = dblarr(sz(2),sz(1))
			6: temp = complexarr(sz(2),sz(1))
		ENDCASE

			FOR j=0,sz(1)-1 DO BEGIN
			str = 'temp(0,j)=REFORM('+dat_str(i)+'(j,*))'
			stat = EXECUTE(str)
			ENDFOR
			stat = EXECUTE(dat_str(i)+'=temp')
			temp = 0b

		   END
; End IDL TRANSPOSE bug fix

	ENDCASE

	tot_row = tot_row + n_row(i)

ENDFOR

data_info = SIZE(dat0)
arr = arr_type(data_info(data_info(0)+1)-1)
stat = EXECUTE('data = ' + arr + 'arr(N_ELEMENTS(pixel),tot_row)')

start_row = 0
FOR i=0,n_flds-1 DO BEGIN
	str = 'data(*,' + STRING(start_row) + ':' + $
	       STRING(start_row+n_row(i)-1) + ') = ' + dat_str(i)
	stat = EXECUTE(STRCOMPRESS(str,/remove_all))
	stat = EXECUTE(dat_str(i) + ' = 0b')
	start_row = start_row + n_row(i)
ENDFOR


PRINT,' '

IF (multobs NE 'LAS') THEN BEGIN

	CASE multobs OF

	'AV1' : BEGIN
		PRINT, 'Extracting Approach Vector 1 Data'
		pixel = pixel(WHERE(av EQ 1))
		data = data(WHERE(av EQ 1),*)
		PRINT,STRING(N_ELEMENTS(pixel))+' pixels extracted'
		END

	'AV2' : BEGIN
		PRINT, 'Extracting Approach Vector 2 Data'
		pixel = pixel(WHERE(av EQ 2))
		data = data(WHERE(av EQ 2),*)
		PRINT,STRING(N_ELEMENTS(pixel))+' pixels extracted'
		END

	else :

	ENDCASE

	IF ((multobs EQ 'AVG') OR (multobs EQ 'AV1') OR (multobs EQ 'AV2')) $
	THEN BEGIN

		PRINT, ' '
		PRINT,'Averaging Multiple Observations'

		pixavg,pixel,data,FLOAT(sent_val)
			; C pixel averaging routine
			; Condensed pixel and data arrays
			; passed back to their corresponding
			; input arrays.

		PRINT,STRING(N_ELEMENTS(pixel))+' pixels remaining'

	ENDIF

	IF (multobs EQ 'FIR') THEN BEGIN
		flip = N_ELEMENTS(pixel) - 1 - LONG(FINDGEN(N_ELEMENTS(pixel)))
		data = data(flip,*)
		pixel = pixel(flip)
	ENDIF

ENDIF

PRINT,'Rasterizing Data'

IF (face EQ -1) THEN $
	pix2xy,pixel,x_out,y_out,res=res,data=data,raster=outmap, $
	       sixpack=sixpack $
ELSE $
	pix2xy,pixel,x_out,y_out,res=res,data=data,raster=outmap, $
	       sixpack=sixpack,/face


PRINT, ' '
PRINT, 'Build Completed'
PRINT, ' '

exit:

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


