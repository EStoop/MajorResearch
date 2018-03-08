	Pro fieldlist,bin_header,data_set_type,fields,dims,fieldcnt
;+
;       FIELDLIST- This primitive routine was part of SELINPUT, but
;       was separated for future use in other parts of the code. This
;       routine gets the FITS fields from the header for user to choose.
;       A new field not contained in the header has been added, namely
;       SIGNAL(OR IFG)+WEIGHT.
;-
;       SPR 10829 created 15-Apr-1993 D. Ward, J. Newmark
;  SPR 11759 17-May-1994   Ingest FIRAS PDS's. J. Newmark
; SPR 11905 07-Sep-1994   Ingest DIRBE PDS's. J. Newmark
;  SPR 11981 08-Nov-1994  Read in skymap_info files.  J. Newmark 
; SPR 12140 16-Mar-1995   Ingest FIRAS line maps. J. Newmark
;
	tfields = fxpar(bin_header,'TTYPE*',count=fieldcnt)
	fields = strarr(fieldcnt+1)

	if strtrim(data_set_type) eq 'DMR' then begin
	    p = where('SIGNAL  ' Eq tfields)
	    fields = [tfields(0:p(0)), 'SIGNAL+WEIGHTS', $
		     tfields(p(0)+1:fieldcnt-1)]
	    dims = -1
	    fieldcnt = fieldcnt + 1
	endif

	if strtrim(data_set_type) eq 'FIRAS' then begin
	    tdims =  fxpar(bin_header,'TFORM*')
	    dims  = strtrim(tdims)
            vals = ['REAL_SPE','COADDED_','SPECTRUM','TEMP']
            found = lonarr(n_elements(vals))
	    for i = 0, n_elements(vals)-1 do found(i) = $
                  where(vals(i) Eq strtrim(tfields))
            vm = where(found ne -1)
            IF (vm(0) eq -1) THEN begin
                fields=tfields
                goto, abort
            endif else p = found(vm(0))
            if p(0) lt fieldcnt - 1 then begin
               fields = [tfields(0:p(0)),vals(vm(0))+'+WEIGHT',tfields(p(0)+1:fieldcnt-1)] 
               dims = [tdims(0:p(0)),tdims(p(0)),tdims(p(0)+1:fieldcnt-1)]
            endif else begin 
               fields = [tfields(0:p(0)),vals(vm(0))+'+WEIGHT']
               dims = [tdims(0:p(0)),tdims(p(0))]
	    endelse
            fieldcnt = fieldcnt + 1
	endif
	if data_set_type eq 'HLIN-MAP' then begin
	    tdims =  strarr(fieldcnt)
            tdims(3:4) = '(10)'
	    p = where('LINE_FLU' Eq tfields)
	    fields = [tfields(0:p(0)), 'LINE_FLU+WEIGHT', $
		     tfields(p(0)+1:fieldcnt-1)]
	    dims = [tdims(0:p(0)),tdims(p(0)),tdims(p(0)+1:fieldcnt-1)]
	    fieldcnt = fieldcnt + 1
	endif
	if data_set_type eq 'LLIN-MAP' then begin
	    tdims =  strarr(fieldcnt)
            tdims(3:4) = '(8)'
	    p = where('LINE_FLU' Eq tfields)
	    fields = [tfields(0:p(0)), 'LINE_FLU+WEIGHT', $
		     tfields(p(0)+1:fieldcnt-1)]
	    dims = [tdims(0:p(0)),tdims(p(0)),tdims(p(0)+1:fieldcnt-1)]
	    fieldcnt = fieldcnt + 1
	endif

	if strtrim(data_set_type) eq 'DIRBE' then begin
	    tdims =  fxpar(bin_header,'TFORM*')
	    fields = tfields
	    dims  = strtrim(tdims)
	endif

	If data_set_type eq 'PXINFO' then begin
	    fields = tfields
	    dims  = -1
	Endif
	If data_set_type eq '        ' then begin
	    fields = tfields
	    dims  = tdims
	Endif
abort:
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


