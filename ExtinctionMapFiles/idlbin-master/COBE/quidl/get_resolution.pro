;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    GET_RESOLUTION returns the Skymap's resolution.
;
;DESCRIPTION:
;    GET_RESOLUTION is an IDL function which uses the Cobetrieve Data 
;    Server to determine a Skymap dataset's resolution.
;
;CALLING SEQUENCE:
;    resolution = GET_RESOLUTION ('dataset')
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    DATASET           I      string    The name of the Skymap dataset.
;    RESOLUTION        O      long      The resolution number.
;
;WARNINGS:
;    None.
;
;EXAMPLE:
;    The following call returns the resolution of a FIRCOADD:[FIRAS.
;    SKYMAP_2]FCS_CCMSP_LL.ED_8934308_9012014 Skymap dataset.
;
;    resolution = GET_RESOLUTION ('fircoadd:[firas.skymap_2]
;    fcs_ccmsp_ll.ed_8934308_9012014')
;
;    The following block of statements may be imbedded into a user's
;    IDL procedure.
;
;    index = 0
;    index = GET_RESOLUTION ('user_disk:[temp]test.skymap')
;    if (index ne 0) then begin
;       .
;       .
;       .
;    endif
;
;#
;COMMON BLOCKS:
;    None.
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    The user should set the resolution variable to 0 prior to invoking
;    this function if the variable is going to be used to check the status
;    of this function call in an IDL program.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS,ETC.:
;    None.
;
;MODIFICATION HISTORY:
;    Written by Song Yom, January 1992 (SER# 9616)
;
;.TITLE
;Routine GET_RESOLUTION
;-
FUNCTION GET_RESOLUTION
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


