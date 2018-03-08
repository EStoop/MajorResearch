Function lstflds, data_type, data_set_name, flist, fcount
;+
; NAME:
;	lstflds 
; PURPOSE:
;       This function returns a list of field names that
;	are contained in the data set.  It will hide the
;	fact that some files may be CSM, CISS, or FITS.
;	CSM files will require the Data Server to get this
;	information.
;	
; CATEGORY:
;	UIMAGE, files, utility, data base access
;
; CALLING SEQUENCE:
;	status = lstflds( type, name, flist, fcount )
;
; INPUTS:
;       type 	a string containing the data set type 
;	name	a string naming the data set
;
; OUTPUTS:
;       status 	-1 if there is a problem
;		 0 if there are none found
;		 number of fields normally (fcount)
;	flist	a string array with field names
;	fcount 	the number of elements in flist
;	
; COMMON BLOCKS:
;	None.
; RESTRICTIONS:
;	None.
; MODIFICATION HISTORY:
;       Creation:  Pete Kryszak, GSC, February 1992.
;-
;
;
fcount = 0
;                                      
flist = [ 'SIGNAL  - Signal, Primary Field', $
	  'ECLON   - Ecliptic Longitude', $
	  'ECLAT   - Ecliptic Latitude', $ 
	  'SERROR  - Signal Error (+/-)', $
	  'PIXEL   - Quad-Cube Pixel Number' $
	]

fcount = n_elements( flist )
;
;print, 'this is lstflds.pro'
;print, 'flist  = ', flist
;print, 'fcount = ', fcount
;
if fcount eq 0 then begin
	print, 'no fields found.'
	return, 0
	endif
;
; some debugging
;print, fcount
;for i = 0, fcount-1 do begin
;	print, flist(i)
;	endfor 
;
return, fcount
End
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


