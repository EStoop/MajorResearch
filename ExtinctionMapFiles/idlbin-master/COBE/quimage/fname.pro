	function fname, filepath, path=path

;+
; fname takes a string as input, determines if it's only a filename or
; both a filename and path and returns the filename along.  Optionally
; it will also return the path in the optional arugement path
; note: the path is always returned with the trailing seperator already
; attached (ie: ] for vms, / for unix, \ for windows)
;-
	if !cgis_os eq 'vms' then begin
	   pos = strpos(filepath,']')
	   if pos gt 0 then begin
	      path = strmid(filepath,0,pos+1)
	      file = strmid(filepath,pos+1,strlen(filepath))
	   endif else begin
	      file = filepath
	      path = ''
	   endelse
	endif

	if !cgis_os eq 'unix' then begin
	   pos = -1
           for i = 0,strlen(filepath) do $
	       if strmid(filepath,i,1) eq '/' then pos = i
	   if pos gt 0 then begin
	      path = strmid(filepath,0,pos+1)
	      file = strmid(filepath,pos+1,strlen(filepath))
	   endif else begin
	      file = filepath
	      path = ''
	   endelse
	endif

	if !cgis_os eq 'windows' then begin
	   pos = -1
           for i = 0,strlen(filepath) do $
	       if strmid(filepath,i,1) eq '\' then pos = i
	   if pos gt 0 then begin
	      path = strmid(filepath,0,pos+1)
	      file = strmid(filepath,pos+1,strlen(filepath))
	   endif else begin
	      file = filepath
	      path = ''
	   endelse
	endif

	return, file
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


