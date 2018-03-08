 Function WEXName, FileName
;+
;   Description
;   -----------
;   WEXName searches for a WEX file whose name is specified in the input
;   parameter FileName and returns the fully qualified name and path.  The
;   function looks in the directories whose names are contained in the global
;   array WEXPaths (found in wexpath.inc).  Upon the first instance of the
;   file, the directory name and file name are concatenated and returned.
;   A null string is returned if the file is not found.
;
;#
;   Called by
;   ---------
;   LoadWEX
;
;   Routines Called
;   ---------------
;   InitWEX
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  November 1992
;
;   11/20/92  KRT  Added call to InitWEX
;   12/28/92  KRT  Modified to use arrays of paths and FindFile.
;   01/22/93  KRT  Improved documentation per SPR 10480.
;   SPR 11127  18 Oct 93  IDL for Windows Compatability. J. Newmark
;
;.Title
;   Function WEXName
;-
;
 On_Error, 2  ; Return to caller if an error occurs
;
@wexpath.inc  ; Access to WEX file names, paths, and extensions
;
;Begin,
    InitWEX
;
    i        = 0
    IF (!cgis_os EQ 'windows') THEN filename = filename + '*'
    LastPath = N_Elements( WEXPaths )-1
    FileList = [ "" ]
;
    While ((i le LastPath) and (FileList(0) eq "")) do begin
       Key      = StrTrim( WEXPaths(i), 2 ) + StrTrim( FileName, 2 )
       FileList = FindFile( Key )
       i        = i + 1
    EndWhile
;
    Return, FileList(0)
 End
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


