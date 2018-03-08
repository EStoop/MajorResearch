;+
; NAME:
;       FITSHDR
; PURPOSE:
;       This Procedure will display the primary fits header in a fits file
; CALLING SEQUENCE:
;       FITSHDR, dsname
; INPUTS:
;       dsname           the name of the file to display
; OUTPUTS:
;         None, the file is displayed using uscroll
;
; COMMON BLOCKS:
;  None.
; RESTRICTIONS:
;  None.
; MODIFICATION HISTORY:

	Pro fitshdr, dsname

	files = findfile(dsname,count=fcount)
	if fcount eq 0 then begin
	   print, 'FITSHDR: cannot find file: ',dsname
	   goto, abort
	endif

	result = is_fits(dsname,extension)
	if result ne 1 then begin
	   print, 'FITSHDR: The file you have specified is not a FITS format file'
	   goto, abort
	endif

;	set data set name according to selection by user
	use_file = 'N'
	header = headfits(dsname)
	review_header, header, 'FITS Header = '+dsname, use_file, $
                       question='no q'

abort:
	Return

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


