;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;	LSTDSNAM returns a list of files of specified format.
;
;DESCRIPTION:
;       This function returns a list of file names that
;	match the given specification.  It will hide the
;	fact that some files may have DECnet node names
;	included (this function is waiting on this 
;	capability from other utilities).  
;
; CALLING SEQUENCE:
;	status = lstdsnam( input_format, spec, flist, extension, fcount )
;
; ARGUMENTS: (I=input, O = output, [] = optional ):
;	input_format	I str		String containing the input
;                                       format 'CISS', 'FITS' or 'CSM'
;
;       spec		I str		A string containing the node, 
;	                                disk, directory, and file
;	                                specification 
;       status 		O int		-1 if there is a problem
;	                                0 if there are no files found
;	                                count of files normally (fcount)
;	flist		O str [arr]	a string array with file names
;	extension	O str [arr]	a string array with file extensions
;	fcount		O int           the number of elements in
; 	                                flist
;
; EXAMPLES:
;	Suppose you want to find the FITS files in the directory 
;	cgis$data, then spec = 'cgis$data:*.*' and input_format =
;	'FITS'.  Invoke as
;
;	    status = lstdsnam( input_format, spec, flist, extension, fcount )
;
;	Both status and fcount return the number of files found, valid
;	files are listed in flist.
;
; WARNINGS:
;	Doesn't work across network
;#
; COMMON BLOCKS:
;	None.
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;	Calls IS_SS and IS_FITS to determine valid save sets and fits
;	files
;
; MODIFICATION HISTORY:
;       Creation:  Pete Kryszak, GSC, February 1992.
;	Dave Bazell	4 Feb 1993	Added calles to is_ss and
; 	                                is_fits to check validity of
;                                       files.
;                                       sprs 10463, 10477
;	Dalroy Ward    22 Feb 1993	Added extension to variables
;					returned from this routine,
;					supports FITS extensions
;					spr 10819
;-
;
;
Function lstdsnam, input_format, spec, flist, extension, fcount
;
fcount = 0

files = findfile( spec, count=fcount )
flist = strarr(1)
extension = strarr(1)

; Check files only if fcount > 0, otherwise just return
if (fcount gt 0 ) then begin

    ; Copy CISS files into returned file list
    if (strupcase(input_format) eq 'CISS' ) then begin
        for i = 0,fcount-1 do begin
            if ( is_ss(files(i)) ) then begin
                flist = [flist,files(i)]
            endif
        endfor
    endif

    ; Copy FITS files into returned file list
    if ( strupcase(input_format) eq 'FITS' ) then begin
        for i = 0, fcount-1 do begin
            if ( is_fits(files(i), exten) ) then begin
                flist = [flist,files(i)]
	        extension = [extension,exten]
            endif
        endfor
    endif

    ; For CSM format just return all the appropriate files found
    if (strupcase(input_format) eq 'CSM') then begin
	flist = [flist,files]
    endif

    ; Check that some files were found.  Remember that flist was
    ; dimensioned to be 1 at the start so flist(0) is empty.  Same
    ; holds for extension, but only look at extension if we were
    ; looking at fits files

    fcount = n_elements(flist)
    if (fcount gt 1) then begin
        flist = flist(1:fcount - 1)
        if (strupcase(input_format) eq 'FITS') then $
            extension = extension(1:fcount-1)
    endif else begin
        fcount = 0
    endelse

endif

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


