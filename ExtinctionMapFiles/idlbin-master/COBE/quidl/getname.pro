 Function GetName, FileName, Paths, IDLPath=IDLPath
;+
; NAME
;       GetName
;
; PURPOSE:
;       GetName searches for a file whose name is specified in the input
;       parameter FileName and returns the fully qualified name and path.
;       The function looks in the directories whose names are contained
;       in the input string array (or string) Paths.  Upon the first in-
;       stance of the file, the directory name and file name are conca-
;       tonated and returned.  A null string is returned if the file is
;       not found.  The keyword IDLPath causes GetName to search direct-
;       ories in the default directory and the IDL path (!path) first.
;
; CATEGORY:
;       I/O Utility
;
; CALLING SEQUENCE:
;
;       Name = GetName( FileName [, Paths ] [,/IDLPath ] )
;
;
; INPUTS:
;       FileName = a simple file name (no directory spec.).
;                  If this is the only parameter used, the
;                  default directory is searched.
;
;       Paths    = a directory name or an array of directory
;                  names.  The full file name specification
;                  of the first occurance of the file is re-
;                  turned.  The default directory can be
;                  specified using a blank element.
;
;       IDLPath  = will cause the program to search through
;                  the default directory first, then any
;                  directories in the IDL path (!path), and
;                  finally any directories in Paths (if spec-
;                  ified).
;
; OUTPUTS:
;       A string containing the fully qualified file or a null
;       if the file was not found.
;
; EXAMPLES:
;
;       Find and return the name of a file on !path:
;
;             UIDL>print, getname( "planck.pro", /idlpath )
;             CGIS$ROOT:[IDL]PLANCK.PRO;4
;
;       Identify a special path and find a file and its name:
;
;             UIDL>paths = [ "CGIS$UITREE:", "CSDR$BULLETIN:[IDL]" ]
;             UIDL>print, getname( "news.asc", paths )
;             CSDR$BULLETIN:[IDL]NEWS.PRO;4
;
;       If GetName can not find the file, it returns a null string:
;
;             UIDL>print, getname( "fudged.data", paths, /idlpath )
;
;             UIDL>
;
;#
; ROUTINES CALLED:
;       findfile (standard IDL)
;
; COMMON BLOCKS:
;	None.
;
; RESTRICTIONS:
;       None.
;
; MODIFICATION HISTORY:
;
;   Created by K. Turpie,  General Sciences Corporation,  January 1992
;
;.Title
;   Function GetName
;-
;
 On_Error, 2  ; Return to caller if an error occurs
;
;Begin,
;
    If  (N_Params()             eq  0)  then  $
       Message, "No parameters.  A file name must be specified at a mininum."
;
    If ((N_Elements( FileName ) gt  1)  or    $
        (N_Elements( FileName ) eq  0)) then  $
       Message, "No file name.  A file name must be specified at a mininum."
;
    If  (StrTrim( FileName )    eq "")  then  $
       Message, "Null file name."
;
    Using_Paths   = ((N_Params() gt 1) and (N_Elements( Paths ) gt 0))
    Using_IDLPath = Keyword_Set( IDLPath )
;
    If (Using_IDLPath) then begin
;
       PathList = StrParse( !Path, Delim=',' )
       indx = Where( StrPos( PathList, '@' ) lt 0 )
;
       If (indx(0) ne -1) then PathList = [ "", PathList( indx ) ] $
                          else PathList = [ "" ]
    EndIf
;
    If (Using_Paths) then                                                   $
       If (N_Elements( PathList ) gt 0) then PathList = [ PathList, Paths ] $
                                        else PathList = Paths               $
    Else                                                                    $
       If (N_Elements( PathList ) le 0) then PathList = [ "" ]
;
    LastPath = N_Elements( PathList )-1 > 0
;
    i        = 0
    FileList = [ "" ]
;
    While ((i le LastPath) and (FileList(0) eq "")) do begin
       Key      = StrTrim( PathList(i), 2 ) + StrTrim( FileName, 2 )
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


