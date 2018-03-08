;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    READ_SKYMAP returns data for multiple fields from a Skymap dataset.
;
;DESCRIPTION:
;    READ_SKYMAP is an IDL function which uses the Cobetrieve Data 
;    Server to extract data from a Skymap archive dataset, and return 
;    the data in pixel order.  This routine automatically returns
;    the pixel number.  User may override the default Record Definition 
;    File specified in the DAFS database by calling the SET_RDF procedure.
;    Refer to the "CSDR Data Base Manual" for complete descriptions of 
;    fields contained in COBE Skymap archive datasets.
;
;CALLING SEQUENCE:
;    status = READ_SKYMAP ('dataset','field[,...]',pixel,data[,...],
;               [startpixel=n],[count=n])
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    DATASET      I  string     The Skymap dataset specification.  If the 
;                               dataset parameter does not include an archive 
;                               specification, the default archive location 
;                               will be determined from the DAFS database.
;    FIELD[,...]  I  string     Field in the dataset to be retrieved.
;                               Multiple fields must be separated by commas.
;    STARTPIXEL=n I  keyword    Reference start pixel number.  The keyword para-
;                               meter 'n' must be an integer.  If this keyword
;                               is omitted, 'startpixel' defaults to zero.
;    COUNT=n      I  keyword    Total number of pixels.  The keyword parameter 
;                               'n' must be an integer.  If this keyword is 
;                               omitted, the total number of pixels defaults 
;                               to the number of pixels in the whole map.
;    FACE=n	  I  keyword	Number of face to be returned.  Overrides 
;				startpixel and count.
;    /RESOL	  I  keyword	Requests READ_SKYMAP to return the skymap's
;				resolution.  This requires an IDL variable to
;				receive the resolution.
;    PIXEL        O  arr        IDL variable into which READ_SKYMAP stores the 
;                               pixel numbers corresponding to the data returned
;    DATA[,...]   O  arr[,...]  IDL array(s) into which READ_SKYMAP stores the 
;                               data for the specified field.  Multiple data
;                               variables must be separated by commas, and 
;                               there must be enough data variables for the 
;                               number of fields specified.
;    [RES]	  O   int	IDL variable which READ_SKYMAP will set to the
;				skymap's resolution.  NOTE THAT THIS WILL BE 
;				THE FIRST VARIABLE AFTER THOSE USED FOR DATA. 
;    STATUS       O   int       Integer condition value.  A value of 1 is 
;                               returned for success.
;
;WARNINGS:
;    READ_SKYMAP retrieves up to 8 fields plus the pixel number field
;    from a Skymap archive dataset.
;
;EXAMPLE:
;    The following call retrieves PIXEL, TEMPERATURE, and TIME fields
;    from the CSDR$DMR_REFERENCE:DSK_SKY31A.MISSION Skymap dataset.
;
;    status = READ_SKYMAP ('csdr$dmr_reference:dsk_sky31a.mission', $
;	'temperature,time',pixel,data1,data2,startpixel=0,count=1024)
;
;    This call requests the return of the "photometry" field from Face
;    0 of a skymap, and that the skymap resolution be returned in IRES.
;
;    istat = READ_SKYMAP 
;	('dataset='dirbe_pds1:[dirbe_edit]bci_sad.skymp_cal_913620211', $
;	 'photometry', pix, data, ires, FACE=0, /RESOL)
;
;    The following block of statements may be imbedded into a user's
;    IDL procedure.
;
;    status = 0
;    status = READ_SKYMAP ('who:that.be', 'this,and,that', 
;                           pix, d1, d2, d3, startpixel=0, count=65536)
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
;               Shirley Masiee, December 1991 (SPR# 9616)
;
;.TITLE
;Routine READ_SKYMAP
;-
FUNCTION READ_SKYMAP
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


