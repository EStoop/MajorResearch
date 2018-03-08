;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  READ_RDF returns field descriptions from RDF files.
;
;DESCRIPTION:
;    An IDL function which uses the Cobetrieve Data Server to return an 
;    alphabetical list of field names from Record Definition Files.
;
;CALLING SEQUENCE:
;    names = read_rdf ('dataset')
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    DATASET     I   str         The name of the dataset.
;    NAMES       O   str array   The listing of the field names.
;
;EXAMPLE:
;    The following call lists the array of field names for the dataset
;    fcs_ccmsp:
;
;        arr = read_rdf('fcs_cmmsp')
;
;    The following block of statements may be imbedded into a user's
;    IDL procedure:
;
;        arr = read_rdf ('fcs_cmmsp')
;        print, arr
;#
;.TITLE
; Routine READ_RDF
;-
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


