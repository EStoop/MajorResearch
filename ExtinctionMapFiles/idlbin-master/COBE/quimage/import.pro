Pro import
;
; MODIFICATIONS
;	name	date		spr	purpose
;	ward	24-feb-1993	10819	modify call to selinput to add
;					fits_exten variable to support
;					reading FITS extensions
;
!err = 1
!error = 1
old_quiet=!quiet
;
SelInput, inform, intype, dsname, fits_exten, field, subscr, weights, fnum
;print, '!err = ', !err, ' ', !err_string
if !ERR Eq 99 then begin
	!ERR = 1
	!ERROR = 1
	return
	EndIf
if !err ne 1 then begin
        !quiet=0
	message, 'Input data selection invalid.', /informational
        !quiet=old_quiet
	return
	endif

; this is import, so output_format is always UIMAGE
output_format = 'UIMAGE'
outname = ''
;
;
dconvert, inform, intype, dsname, fits_exten, fnum, field, subscr, $
          weights, output_format, outname
;
;
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


