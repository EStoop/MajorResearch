;+
; NAME: Pixface
;
; PURPOSE: To extract a face from a pixel and data list (ie: before
;          rasterization.
;
;  CATEGORY: Utility
;
; CALLING SEQUENCE: Pixface, pixel, data, moredata=moredata, faceno=faceno,
;                            res=res
;
; INPUTS:
;          pixel	a pixel vector
;	   data		a data vector, data can be 1 or 2 dimensional
;			the routine will figure out which is the pixel
;			dimension.  (Note: there HAS to be a one-to-one
;			correspondence between the pixel array and one
;			dimension of the data vector.)
;	   moredata	an additional data array which EXACTLY matches up
;                       with the data array passed in earlier (this was 
;                       added to allow subsetting of the weights array at
;                       the same time).  If this array does not match up
;			EXACTLY then the results of this routine will be
;                       undefined.
;	   faceno       the face number to be extracted
;	   res          the resolution of the resulting skycube
;
; OUTPUTS:
;	   pixel	the reduced pixel vector
;	   data		the reduced data vector
;	   !err		will be set to 0 for successful completion of
;			this routine, and -1 if an error occured.
;
; WARNINGS:  This routine destroys the original data and pixel vectors
;	     so if you need to use them again, be sure to save them
;	     before calling this routine!
;#
; COMMON BLOCKS: None
;
; INCLUDE FILES: None
;
; RESTRICTIONS: None
;
; SUBROUTINES CALLED:
;
; MODIFICATION HISTORY: Created Dalroy Ward, GSC    May 20, 1993
;---------------------------------------------------------------------
;
	Pro pixface, pixel, data, moredata=moredata,faceno=faceno, res=res

	!err = 0
;	check to be sure that all of the required parameters were passed
;	in
	if n_elements(res) le 0 then begin
	   message, 'You must supply the resolution of the data',/informational
	   !err = -1
	endif

	if n_elements(faceno) le 0 then begin
	   message, 'You must supply the face number to extract',/informational
	   !err = -1
	endif

	if faceno lt 0 or faceno gt 5 then begin
	   message, 'You have supplied an invalid face number', /informational
	   !err = -1
	endif

	if n_elements(pixel) le 0 then begin
	   message, 'You must supply the pixel data!',/informational
	   !err = -1
	endif

	if n_elements(data) le 0 then begin
	   message, 'You must supply the data along with the pixel list',$
		    /informational
	   !err = -1
	endif

;	now check to see if any errors occured, if so then return
	if !err ne 0 then return

	; begin to do some real work here!

	; calculate the number of points on a side for this resolution
	res_calc = res - 1
	pface = (long(2)^res_calc)^2

	; calculate the starting and ending pixel numbers for the face
	end_face = (pface * (faceno+1)) - 1
	start_face = end_face - pface + 1

	pixel_size = size(pixel)
	if pixel_size(0) gt 1 then begin
	   message, 'Pixel array must be 1 dimensional',/informational
	   !err = -1
	   return
	endif

	pixel_index = where(pixel ge start_face and pixel le end_face)

;	check to see if we found any data
	if pixel_index(0) eq -1 then begin
	   message, /continue, 'No data in selected face'
	   !err = -1
	   return
	endif

;	now, resize the pixel array
	pixel = pixel(pixel_index)

	data_size = size(data)

	if data_size(0) eq 1 then begin
	   data = data(pixel_index)
	endif

	if data_size(0) eq 2 then begin
	   if data_size(1) eq pixel_size(1) then data = data(pixel_index,*) $
           else data = data(*,pixel_index)
	endif

	if n_elements(moredata) gt 0 then begin
	   more_size = size(moredata)
	   if more_size(0) eq 1 then begin
	      moredata = moredata(pixel_index)
	   endif

	   if more_size(0) eq 2 then begin
	      if more_size(1) eq pixel_size(1) then moredata = $
					       moredata(pixel_index,*) $
              else moredata = moredata(*,pixel_index)
	   endif
	endif

	!err = 0
	return
	end
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


