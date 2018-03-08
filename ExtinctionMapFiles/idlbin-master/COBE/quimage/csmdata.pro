Pro csmdata, 		$
	dsname, 		$
	dataset_field, 		$
	subscript,		$
	facenum, 		$
	index_levels,		$
	data
	
;+
; NAME:
;	CSMDATA
;
; PURPOSE:
;       This procedure reads the specified data set and 
;	puts it into an unfolded cube format in sky looking
;	orientation (right-hand T).  It uses Read_Sky_Map.
;	
; CATEGORY:
;	User interface, Menu.
;
; CALLING SEQUENCE:
;	csmdata, dsname, dataset_field, subscript, $
;			facenum, index_levels, data
;	
; INPUTS:
;	dsname	the name of the data set
;	dataset_field	the field name 
;	subscript	the subscript of the array
;	facenum		the number of the face or 7 (All)
;	index_levels	indicates the pixelization or pixel size
;
; OUTPUTS:
;       data 		an array of data read from the file 
;
; SIDE EFFECTS:
;	None.
;
; COMMON BLOCKS:
;	None.
;
; RESTRICTIONS:
;	None.
;
; SUBROUTINES/FUNCTIONS CALLED:
;	Read_Skymap
;
; MODIFICATION HISTORY:
;       Creation:  Pete Kryszak, GSC, February 1992.
;	14-Jan-1993	10177	Kryszak	debug code removed
;       17-Mar-1993     10626   Turpie  change second call to Pix2xy to use
;                                       /Face instead of Face=FaceNum
;
;-
;
;!quiet = 0
message,  "this is csmdata...", /informational
message,  "the name of the data set "+dsname, /informational
message,  "the field name "+dataset_field, /informational
message,  "the subscript of the array "+subscript, /informational

;
; if there are subscripts specified then they must be included
; ------------------------------------------------------------
field_name = dataset_field
if subscript ne '' then begin
	field_name = dataset_field + '(' + strtrim(subscript,2) + ')'
        endif


;
; if only a map face, then find start pixel and count
; ---------------------------------------------------
pixcnt  = LONG( 2L^( 2 * (index_levels(0)-1) ) )
if facenum eq 7 then begin
  strtpix = 0          ; start at beginning
  pixcnt = pixcnt * 6  ; get six faces of data
  endif $
else begin
  strtpix = pixcnt * facenum ; start at first pixel in face
  endelse				  ; and get one face worth
;
;
; get the data via Read Skymap
;-----------------------------
print, 'reading sky map data'

message,  "this is csmdata again...", /informational
message,  "the name of the data set "+dsname, /informational
message,  "the field name "+dataset_field, /informational
message,  "the subscript of the array "+subscript, /informational

status = $
read_skymap(dsname,field_name,pixel,array, $
		startpixel=strtpix,count=pixcnt)

if status ne 1 then begin
	message, 'Read_Skymap did not complete successfully.', /continue 
	!err = -1
	!error = -1
	return 
	endif

;
; see if this needs to be a three dimensional object
; -------------------------------------------------------
dim3 = SIZE(array)
thrddim = 0 ; defaults to 2-D, third dimension is zero
IF dim3(0) EQ 3 THEN thrddim = dim3(3)	; size of third dimension
IF dim3(0) GT 3 THEN BEGIN
	MESSAGE, 'Too many dimensions returned by Read_Skymap.', /continue 
	return  
	ENDIF
;
; find the actual number of pixels returned
; there may be future problems here due to
; multiple pixels / record, but I don't care now.
; -----------------------------------------------
count = N_ELEMENTS( pixel )


dimns = 2^(index_levels(0)-1)	;32,64, etc pixels on a side

;print, 'index_levels = ', index_levels(0)
;print, 'dimns = ', dimns
;print, 'thrddim = ', thrddim

IF facenum EQ 7 THEN BEGIN	; full map in unfolded T format
  faceflag = 'U' ;for pix_to_xy, unfolded cube
  dimnsx = 4*dimns ; four faces wide
  dimnsy = 3*dimns ; three faces high
  IF thrddim ne 0 then  data = FLTARR( dimnsx, dimnsy, thrddim ) $
  ELSE data = FLTARR( dimnsx, dimnsy )
  ENDIF $
ELSE BEGIN; not all faces, but only one, so square
  faceflag = 'F' ;for pix_to_xy, face not unfolded cube
  IF thrddim ne 0 then  data = FLTARR( dimns, dimns, thrddim ) $
  ELSE data = FLTARR( dimns, dimns )
  ENDELSE

;
; reformat the data from pixel list to xy raster 
; ----------------------------------------------
print, 'reformatting pixel list to unfolded cube'

if faceflag eq 'U' then begin
   pix2xy, pixel, res=index_levels(0),data=array,raster=data
   endif $
  else begin 
   pix2xy, pixel, res=index_levels(0),data=array,raster=data,/face
   endelse

return

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


