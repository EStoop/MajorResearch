
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    SET_RDF defines an alternate Record Definition File.
;
;DESCRIPTION:
;    SET_RDF is an IDL procedure which allows the user to define an
;    alternate Record Definition File different from one in the DAFS
;    database.  Subsequent calls to READ_TOD and READ_SKYMAP IDL
;    functions are affected.  Users should note that this function
;    defines a either logical name (VMS) or a environment variable
;    (UNIX) which is translated during runtime of READ_TOD,
;    READ_SKYMAP, READ_RMS, and GET_RESOLUTION IDL functions.
;
;CALLING SEQUENCE:
;    SET_RDF,rdf
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    RDF               I      string    The Record Definition File
;                                       specification.
;
;WARNINGS:
;
;EXAMPLE:
;    The following call defines a Record Definition File to be used in
;    subsequent READ_TOD, READ_SKYMAP, READ_RMS, and GET_RESOLUTION IDL
;    function invocations.
;
;    SET_RDF,'user_disk:[test]temp.rdf'
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
;    Written by Song Yom, December 1991 (SER# 9616)
;
;.TITLE
;Routine SET_RDF
;-
PRO SET_RDF,RDF
;
;	If the operating system is VMS, then define a logical name.
;	Otherwise, add an environment variable.
;
IF !cgis_os EQ "vms"    THEN SETLOG, 'UFS_RDF', RDF
IF !cgis_os EQ "unix" THEN SETENV, 'UFS_RDF=' + RDF
IF !cgis_os EQ "windows" THEN SETENV, 'UFS_RDF=' + RDF
;
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


