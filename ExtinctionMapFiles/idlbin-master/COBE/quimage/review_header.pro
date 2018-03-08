	pro review_header, fits_header, title, use_file, question=question
;+
; NAME:
;	REVIEW_HEADER
; PURPOSE:
;	   This procedure allows easy review of a FITS header (passed in)
;	   if the user would like to see it.
; CATEGORY:
;	UIMAGE
;
; CALLING SEQUENCE:
;	review_header (fits_header, title, use_file)
;
; INPUTS: 
;	fits_header: the fits header as returned by readfits.pro or a 
;	             similiar routine
;	title      : the title to use on the display, should the user
;		     choose to review the header
;	use_file   : if passed in as Y then the user will be asked if
;		     they would like to use the file, if passed in as
;		     anything else then the user will not be asked.
;		     If not set, use_file defaults to Y.
;	question   : optional, if set this is the question that the user
;		     will be asked in place of "review the fits header?"
;		     If this parameter is set to 'no q' then no question
;		     will be asked and the file will be automatically 
;		     displayed.
;
; OUTPUTS:
;	use_file   : returned as Y if the user has selected the file,
;		     returned as N if not to use the file.  Note: if
;		     use_file is passed into the routine as anything
;		     but Y you should not check it's value upon return
;		     as it's value is undefined.
;	
; COMMON BLOCKS:
;	None.
; RESTRICTIONS:
;	None.
; MODIFICATION HISTORY:
;       Creation:  Dalroy Ward, GSC, March 1993   SPR 10819
;-
	if (n_elements(use_file) Eq 0) then use_file = 'Y'

        if n_elements(question) eq 1 then $
           if question eq 'no q' then ans = 1

	if (n_elements(question) Eq 0) Then begin
	    question = 'Would you like to review the FITS header?'
	endif

	; ask the user if they want to see the file, if not turned off.
	if question ne 'no q' then $
	   ans = one_line_menu([question,'Yes','No'], init=1)

	if ans eq 1 then begin
	    uscroll, fits_header, title=title
	    if use_file Eq 'Y' then begin
	        ans = one_line_menu( ['Use this FITS file?','Yes', 'No' ], $
		      init=1 )
	        if ans eq 2 then use_file = 'N'
	    Endif
	Endif

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


