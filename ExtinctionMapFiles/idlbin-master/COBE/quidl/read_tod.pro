;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    READ_TOD returns data for multiple fields from a Time-ordered dataset.
;
;DESCRIPTION:
;    READ_TOD is an IDL function which uses the Cobetrieve Data 
;    Server to extract data from a Time-ordered archive dataset.
;    User may override the default Record Definition File
;    specified in the DAFS database by calling the SET_RDF procedure.
;    Refer to the "CSDR Data Base Manual" for complete descriptions
;    of fields contained in COBE Time-ordered archive datasets.
;
;CALLING SEQUENCE:
;    status = READ_TOD ('dataset','field[,...]','start','stop',data[,...],
;                [maxrec=nrec] )
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    DATASET     I   string     The Time-ordered dataset.  If the dataset 
;                               parameter includes an archive logical name 
;                               specification, the first five characters of 
;                               the logical name must be 'CSDR$'.  Refer to 
;                               the "COBETRIEVE Programmer's Manual" for 
;                               details on archive dataset naming conven-
;                               tions.  If you omit the archive logical 
;                               name, the default archive location will be 
;                               determined from the DAFS database.
;    FIELD[,...] I   string     Field in the dataset to be retrieved.
;                               Multiple fields must be separated by commas.
;    START       I   string     The start time of the time-ordered dataset
;                               in either the GMT format (yydddhhmmssmmm) 
;                               or the ADT format (dd-mmm-yyyy hh:mm:ss.cc).
;                               If parts of the time fields are omitted,
;                               missing values default to zero.
;    STOP        I   string     The stop time of the time-ordered dataset 
;                               in the GMT format (yydddhhmmssmmm) or the
;                               ADT format (dd-mmm-yyyy hh:mm:ss.cc).  
;                               If parts of the time fields are omitted,
;                               missing values default to zero.
;    MAXREC=n    I   keyword    Maximum number of records to be fetched. 
;                               (Useful in cases where the default algorithm
;                               would have expected more records in the 
;                               requested time interval.)
;    DATA[,...]  O   arr[,...]  IDL array in which READ_TOD stores the data 
;                               for the specified field.  Multiple data
;                               variables must be separated by commas, and 
;                               there must be enough data variables for 
;                               the number of fields specified.
;    STATUS      O   int        Integer condition value.  A value of 1 is 
;                               returned for success.
;
;WARNINGS:
;    READ_TOD retrieves up to 8 fields from a Time-ordered archive 
;    dataset.
;
;EXAMPLE:
;    The following calls retrieve the CT_HEAD.GMT and CT_HEAD.TIME data
;    fields from the CSDR$FIRAS_EDIT:FDQ_SDF_RH time-ordered archive 
;    dataset, the second call specifying maximum records to retrieve.
;
;    status = READ_TOD ('csdr$firas_edit:fdq_sdf_rh',
;             'ct_head.gmt,ct_head.time', '901090903', '90109100', d1, d2)
;
;    status = READ_TOD ('csdr$firas_edit:fdq_sdf_rh', maxrec=200,
;             'ct_head.gmt,ct_head.time', '901090903', '90109100', d1, d2)
;
;    The following block of statements may be imbedded into a user's
;    IDL procedure.
;
;    status = 0
;    status = READ_TOD ('csdr$this:andwhat', 'who,what,when', '1-jan-1990',
;                       '2-jan-1990:20', this, that, huh)
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
;               Shirley Masiee, December 1991 (SER# 9616)
;
;.TITLE
;Routine READ_TOD
;-
FUNCTION READ_TOD
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


