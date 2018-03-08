Pro Translat
;
;
!err = 1
!error = 1
message, '  ', /noprint, /continue
;
SelInput, inform, dsname, field, subscr, fnum
;print, '!err = ', !err, ' ', !err_string
if !err le 0 then begin
	message, '%TRANSLAT, Input data selection invalid.', /informational
	return
	endif
;
;
SelOutpt, outform, outname 
;print, '!err = ', !err, ' ', !err_string
if !err le 0 then begin
	message, '%TRANSLAT, Output data file selection invalid.', $ 
								/informational
	return
	endif
;
;print, 'input data set format  = ', inform
;print, 'input data set name    = ', dsname
;print, 'input data set face    = ', BYTE(fnum), '  (0-6, 7=All)'
;print, 'input data set field   = ', field
;print, 'input data subscript   = ', subscr
;print, 'output data set format = ', outform
;print, 'output data set name   = ', outname
;
;
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


