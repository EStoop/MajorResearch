;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    READ_RMS returns data for multiple fields from an RMS dataset.
;
;DESCRIPTION:
;    READ_RMS is an IDL function which uses the Cobetrieve Data 
;    Server to extract data from an RMS archive dataset.  User must
;    define a Record Definition File by calling SET_RDF procedure
;    prior to invoking this function.
;
;CALLING SEQUENCE:
;    status = READ_RMS ('dataset', 'field[,...]', data[,...] )
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    DATASET     I   string     The RMS dataset specification.
;    FIELD[,...] I   string     Field in the dataset to be retrieved.
;                               Multiple fields must be separated by commas.
;    DATA[,...]  O   arr[,...]  IDL array into which READ_RMS stores the data 
;                               for the specified field.  Multiple data
;                               variables must be separated by commas, and 
;                               there must be enough data variables for 
;                               the number of fields specified.
;    STATUS      O   int        Integer condition value.  A value of 1 is 
;                               returned for success.
;
;WARNINGS:
;    READ_RMS retrieves up to 8 fields from an RMS archive dataset.
;
;EXAMPLE:
;    The following call retrieves NAME and VERSION fields data from
;    the CSDR$ANCIL_REFERENCE:XAC_SHARP_H2.HDR RMS archive dataset.
;
;    status = READ_RMS ('csdr$ancil_reference:xac_sharp_h2.hdr',
;                       'name,version', d1, d2)
;
;    The following block of statements may be imbedded into a user's
;    IDL procedure:
;
;    status = 0
;    status = READ_RMS ('csdr$this:and.that','one,two,three',d1,d2,d3)
;    if (status eq 1) then begin
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
;    The user should set the status variable to 0 prior to invoking
;    this function if the variable is going to be used to check the status
;    of this function call in an IDL program.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS,ETC.:
;    None.
;
;MODIFICATION HISTORY:
;    Written by Leon Herreid,
;               Dave Cottingham,
;               Song Yom, March 1992 (SER# 9616)
;
;.TITLE
;Routine READ_RMS
;-
FUNCTION READ_RMS
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


