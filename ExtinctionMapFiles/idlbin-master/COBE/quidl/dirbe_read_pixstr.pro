;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;   DIRBE_READ_PIXSTR gets time strings through a central pixel
;
;DESCRIPTION:
;   DIRBE_READ_PIXSTR retrieves time-strings of a user-specified             
;   skymap quantity, within a specified time period and radius,        
;   which passed through a specified 'central' pixel.                  
;
;
;CALLING SEQUENCE:
;   data = DIRBE_READ_PIXSTR(start,stop,fieldname,pixno,$
;                            radius[,filename])
;
;ARGUMENTS (I = input, O = output, [] = optional):
;  START     I  str     the time at which DIRBE_READ_PIXSTR starts 
;                       grabbing data from the skymap(s).  'Start' 
;                       may be specified by the user either as a 
;                       CSDR GMT (e.g., 89365010000) or as a
;                       'VAX' time (e.g.,31-dec-1989:01:00:00).
;                       The user must enclose the start time within
;                       single quotes.
;  STOP      I  str     the time at which DIRBE_READ_PIXSTR stops 
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
;                       CURRENTLY, ONLY BCI_CIRSSM-formatted SKYMAPS 
;                       ARE ACCEPTABLE!    
;  PIXNO     I  int     The pixel number of the 'bullseye' pixel.  All 
;                       timestrings returned by DIRBE_READ_PIXSTR will 
;                       pass through this pixel.
;  RADIUS    I  int     Radius of a roughly circular region on the sky 
;                       centered on the 'bullseye' pixel.  Only 
;                       observations within this radius are 
;                       used to construct the time-strings.  Units
;                       of the radius are specified in DIRBE pixels.  
;  FILENAME [I] str     An optional parameter.  Use ONLY if you wish to
;                       get data from a non-archive (RMS) file.  
;                       'Filename' is the name of the RMS file.  The 
;                       user must enclose the filename within single 
;                       quotes.
;  DATA      O  flt arr The output IDL variable into which the requested
;                       observations are stored.  'Data' can have any
;                       variable name.  The output data is stored as a        
;                       2-dimensional array: the first dimension is the       
;                       time-string number, and the second dimension 
;                       contains all the observations for each time 
;                       string.  Since not all time-strings are the same
;                       length, the second dimension is the length of 
;                       the longest time-string -- shorter time-strings 
;                       are padded out to 0.  Note: because the time
;                       string lengths are variable, the location of
;                       the 'bullseye' pixel is not the same from
;                       string to string.
;
;WARNINGS:
;  1. USE AT YOUR OWN RISK! DIRBE_READ_PIXSTR does not release
;     virtual memory correctly and can cause your IDL session to bomb
;     with repeated use within a session. 
;  2. DIRBE_READ_PIXSTR can only handle BCI_CIRSSM-formatted skymaps.  
;     Furthermore, it is limited to skymaps at DIRBE resolution ONLY.
;  3. Byte, I*2,I*4,R*4 and ADT(time) fields are acceptable to   
;     READ_PIXSTR.  Character strings and R*8 are not.
;     ******NOTE!***** Currently, a bug exists which terminates
;     your session if you try to retrieve an ADT TIME field.
;  4. DIRBE_READ_PIXSTR handles ADT (time) fields a little differently 
;     from other field types.  DIRBE_READ_PIXSTR returns elapsed time 
;     (in seconds) from the first observation read for each time-string. 
;     The time of the first observation of each time-string is stored 
;     in a separate IDL character string array with the name of ZTIME.
;  5. The data retrieved are strictly limited by the START and STOP
;     times requested.  (This is unlike DIRBE_READ_PIXEL).
;
;EXAMPLE:
;
;   The user wants to retrieve all PHOTOMETRY(1) timestrings which 
;   pass through DIRBE pixel 49152 during the time interval
;   9007209400377 to 900739500314.  The user has limited the radius
;   of the search region to 6 DIRBE pixels.
;
;    data = READ_PIXSTR('9007209400377', '9007309500314', 
;                       'BCI_CIRSSM.PHOTOMETRY(1)', 49152, 6)               
;
;    The timestrings are returned into IDL variable data.
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
;   Written   JP Weiland August 1989 for IDL V1
;   Modified  Vidya Sagar   (ARC)    Nov '90  for IDL V2
;   Documentation Update   Urmila Prasad  (ARC)  Aug 1991         
;
;.title
;Routine DIRBE_READ_PIXSTR
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


