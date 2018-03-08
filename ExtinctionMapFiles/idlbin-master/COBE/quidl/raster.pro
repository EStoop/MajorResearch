pro raster,pixel,res,face,sixpack,data,raster,x_out,y_out
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    RASTER generates rasterized images
;
;DESCRIPTION:
;    This IDL procedure is called by PIX2XY to create a raster image
;
;CALLING SEQUENCE:
;    raster,pixel,res,face,sixpack,data,raster,x_out,y_out
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    pixel       I     int/long arr        Input pixel list
;    res         I     int                 cube resolution
;    face        I     int                 Face # (-1 if cube)
;    sixpack     I     int                 Sixpack flag (1-six,0-not)
;    data        I     int/flt/dbl arr     input data array 
;                                          (=[-1] if no data)
;    raster      O     int/flt/dbl arr     rasterized output array
;    x_out       O     int arr             x-coordinate array
;    y_out       O     int arr             y-coordinate array
;
;WARNINGS:
;
;EXAMPLE: 
;
; Not user procedure
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Uses PIX2FIJ to "decode" pixel list into face,i,j vector.
;    Creates x_out and y_out arrays from this info and offsets for
;    right_T.
;    Fills raster array if data is provided.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: pix2fij
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Oct 92
;    Delivered 13-NOV-1992  SPR 10212
;-
;
	pix_sz = SIZE(pixel)
	data_info = SIZE(data)

	IF (sixpack EQ 1) THEN BEGIN	; sixpack

		i0 = 3
		j0 = 2
		offx = [0,0,1,2,2,1]
		offy = [1,0,0,0,1,1]
			; Note: These are packed Left_t offsets

	ENDIF ELSE BEGIN

		IF (face EQ -1) THEN BEGIN	; cube

			i0 = 4	; T len (# face sides)
			j0 = 3	; T hgt (# face sides)
			offx = [0,0,1,2,3,0]
			offy = [2,1,1,1,1,0]
				; Note: These are left_t offsets

		ENDIF ELSE BEGIN	; face

			i0 = 1
			j0 = 1
			offx = [0,0,0,0,0,0]
			offy = [0,0,0,0,0,0]

		ENDELSE

	ENDELSE

	fij = pix2fij(pixel,res)
		; get face,column,row info for pixels

	cube_side = 2^(res-1)

	len = i0 * cube_side

	x_out = offx(fij(*,0)) * cube_side + fij(*,1)
	x_out = len - (x_out+1)
			; flip to right orientation
	y_out = offy(fij(*,0)) * cube_side + fij(*,2)

	IF (n_elements(data) NE 1) THEN BEGIN

		thrd = data_info(data_info(0)+2) / pix_sz(1)
			; get third dimension size

		CASE data_info(data_info(0)+1) OF
			1: raster = bytarr(i0*cube_side,j0*cube_side,thrd)
			2: raster = intarr(i0*cube_side,j0*cube_side,thrd)
			3: raster = lonarr(i0*cube_side,j0*cube_side,thrd)
			4: raster = fltarr(i0*cube_side,j0*cube_side,thrd)
			5: raster = dblarr(i0*cube_side,j0*cube_side,thrd)
			6: raster = complexarr(i0*cube_side,j0*cube_side,thrd)
		ENDCASE

		IF (thrd EQ 1) THEN raster(x_out,y_out) = data $

		ELSE BEGIN

			FOR k=0,thrd-1 DO BEGIN

				temp_arr = raster(*,*,k)
				temp_arr(x_out,y_out) = data(*,k)
				raster(*,*,k) = temp_arr

			ENDFOR

		ENDELSE

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


