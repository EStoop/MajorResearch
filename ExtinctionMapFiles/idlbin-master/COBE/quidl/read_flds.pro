pro read_flds,max_flds,n_flds,fields

; This procedure is called by BUILDMAP

	fields = ''
	fld = ''

	FOR i=0,max_flds-1 DO BEGIN

	fn = STRCOMPRESS('Field Name' + STRING(i+1))

read_fld:

	IF (i EQ 0) THEN $
	READ, fn + ' (<CR> to exit BUILDMAP): ', fld ELSE $
	READ, fn + ' (<CR> to proceed)      : ', fld

		IF (fld EQ '') THEN BEGIN
			fields = STRMID(fields,1,STRLEN(fields)-1)
			n_flds = i
			PRINT,' '
			GOTO, lbl_1
		ENDIF ELSE fields = fields +  ',' + fld

	ENDFOR

lbl_1:

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


