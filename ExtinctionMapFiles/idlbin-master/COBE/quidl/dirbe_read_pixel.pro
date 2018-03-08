;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;  DIRBE_READ_PIXEL gets data for a single pixel within a time range
;
;DESCRIPTION:
;   Retrieves all observations of a user-specified
;   skymap quantity, within a specified time period, for a
;   single pixel.
;
;CALLING SEQUENCE:
;  data = DIRBE_READ_PIXEL(start,stop,fieldname,pixno[,filename])
;
;ARGUMENTS (I = input, O = output, [] = optional):
;
;  START     I  str     the time at which DIRBE_READ_PIXEL starts 
;                       grabbing data from the skymap(s).  'Start' 
;                       may be specified by the user either as a 
;                       'ZULU' time (e.g., 89365010000) or as a
;                       'VAX' time (e.g.,31-dec-1989:01:00:00).
;                       The user must enclose the start time within
;                       single quotes.
;  STOP      I  str     the time at which DIRBE_READ_PIXEL stops 
;                       grabbing data from the skymap(s). 
;                       Acceptable formats are the same as for
;                       'start'. The user must enclose the stop time 
;                       within single quotes.
;  FIELDNAME I  str     The skymap field of interest.  The field name
;                       must include both the RDL name, as well as the
;                       name of the field within the RDL.  For example,
;                       'bci_cirssm.photometry(1)' is a valid 'fieldname'
;                       but just 'photometry(1)' is not.  Only one 
;                       element of an array may be retrieved in a call.
;                       The user must enclose the fieldname within 
;                       single quotes.
;  PIXNO     I  int     The pixel number of interest.  Only data for 
;                       this pixel is retrieved.
;  FILENAME [I] str     An optional parameter.  Use ONLY if you wish to
;                       get data from a non-archive (RMS) file.  
;                       'Filename' is the name of the RMS file.  The 
;                       user must enclose the filename within single 
;                       quotes.
;  DATA      O  flt arr The output IDL variable into which the requested
;                       observations are stored.  'Data' can have any
;                       variable name.
;
;WARNINGS:
;  1. USE AT YOUR OWN RISK!  Use of the routine READ_SKYMAP is preferred
;     in most cases, since DIRBE_READ_PIXEL does not release
;     virtual memory correctly and can cause your IDL session to bomb
;     with repeated use within a session.  However, DIRBE_READ_PIXEL
;     possesses the unique capability of reading more than one
;     input skymap, which READ_SKYMAP cannot do.
;  2. In theory, DIRBE_READ_PIXEL can handle any arbitrary skymap.
;     For RDLs not in the default CSDR$RDF location, the user must
;     redefine CSDR$RDF to a location which contains the correct RDF.
;  3. If a user specifies a large enough time period to cause
;     DIRBE_READ_PIXEL to access more than one skymap, the user
;     usually will find that the time ranges of the observations
;     returned are not limited to the input time interval.  Because
;     of the way the skymap access works, DIRBE_READ_PIXEL returns
;     ALL observations from the accessed skymap files.
;  4. Byte, I*2,I*4,R*4 and ADT(time) fields are acceptable to
;     DIRBE_READ_PIXEL.  Character strings and R*8 are not.
;     ******NOTE!***** Currently, a bug exists which terminates
;     your session if you try to retrieve an ADT TIME field.
;  5. DIRBE_READ_PIXEL handles ADT (time) fields a little differently 
;     from other field types.  READ_PIXEL returns elapsed time (in 
;     seconds) from the first observation read.  The time of the first 
;     observation is stored in a separate IDL character string variable 
;     with the name of ZTIME.
;
;EXAMPLE:
;   The user wants to read the MUX process 1 photometry data for 
;   pixel 204353 from the DIRBE daily skymaps (BCI_CIRSSM) located in
;   CSDR$DIRBE_ARCHIVE.  All daily skymaps which contain data between
;   23-nov-1989 00:01:23 and 24-nov-1989 00:00:00 are accessed.
;   The output is stored in IDL variable PROC1.
;
;   proc1 = DIRBE_READ_PIXEL('23-nov-1989:00:01:23',$
;           '24-nov-1989:00:00:00','bci_cirssm.photometry(1)',204353)
;
;#
;COMMON BLOCKS:
;    There are no IDL common blocks.
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    This routine is written in FORTRAN and linkimaged into IDL.
;    The FORTRAN code makes use of the FIELD RETRIEVER and 
;    CSA capabilities in the CSDR environment.
;
;
;PERTINENT ALGORITHMS, LIBRARY CALLS,ETC.:
;    None.
;
;MODIFICATION HISTORY:
;    Written JP Weiland  August 1989 for IDL V1
;    Modified Vidya Sagar(ARC) Nov '90 for IDL V2
;    Documentation updated by Urmila Prasad (ARC) Aug' 91
;
;.title
;Routine DIRBE_READ_PIXEL
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


