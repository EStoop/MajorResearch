;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    READ_MAP_FACE returns data for a face from a Skymap dataset.
;
;DESCRIPTION:
;    READ_MAP_FACE is an IDL function which uses the Cobetrieve Data
;    Server to extract data from a Skymap archive datasets.  It reads 
;    a single field from a Skymap stored in quad-tree binary format,
;    and rasterizes those pixel observations which lie on a specific
;    map face into an output array.  Refer to the "CSDR Data Base
;    Manual" for complete descriptions of fields contained in the
;    COBE Skymap archive datasets.
;
;CALLING SEQUENCE:
;    data = READ_MAP_FACE ('dataset','rdf.field',face)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    DATASET           I      string    Skymap dataset to be accessed.  If
;                                       the dataset parameter does not
;                                       include an archive specification,
;                                       the default archive location will
;                                       be determined from the DAFS database.
;    RDF.FIELD         I      string    Record Defintion File and the field
;                                       of the specified Skymap dataset to
;                                       be accessed.
;    FACE              I      int       The particular face of the Skymap
;                                       dataset to be accessed.
;    DATA              O      arr       An IDL array. The X-Axis of the
;                                       array has the 'sky_looking' orientation.
;
;WARNINGS:
;    READ_MAP_FACE only retrieves a single field from a Skymap archive
;    dataset.  To access multiple fields, refer to the documentations
;    for the READ_SKYMAP function.  Also, VAX Absolute Date and Time
;    (ADT) and any string type (i.e. GMT) fields are not supported by
;    this routine.
;
;EXAMPLE:
;    The following call returns face 0 of CHAN.SPEC field data for the
;    CACFR:[FIRAS.SKYMAP]FCS_CCMSP_LL.TD_8934319_9012105 Skymap dataset.
;
;    data = READ_MAP_FACE ('cacfr:[firas.skymap]fcs_ccmsp_ll.td_
;    8934319_9012105','fcs_ccmsp.chan.spec',0)
;
;#
;COMMON BLOCKS:
;    None.
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    None.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS,ETC.:
;    None.
;
;MODIFICATION HISTORY:
;    Written by Leon Herreid,
;               Dave Cottingham,
;               Shirley Masiee, April 1991 (SER# 9616)
;
;.TITLE
;Routine READ_MAP_FACE
;-
FUNCTION READ_MAP_FACE
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


