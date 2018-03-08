 PRO Extract, InFileName, OutFileName, ErrStat
;+
;   Description
;   -----------
;   Extract copies the file specified by the name InFileName to the
;   file with the name given in OutFileName.
;
;#
;   Called by
;   ---------
;   UHelp
;
;   Routines Called
;   ---------------
;   LoadWEX
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  May 1992
;
;   11/16/92  KRT  Documentation added per SPR 10216.
;
;.Title
;   Procedure Extract
;-
 On_Error,2     ; Return to caller if an error occurs
;
;Begin,
;
    NameEnd = StrPos( InFileName, '.' )
    If (NameEnd eq -1) then NameEnd = StrLen( InFileName )
;
    LoadWEX, InFileName,  Text,      Heading,  $
             Subtopics,   SubFile,   ErrStat
;
    If (ErrStat eq 0) then begin
;
       Get_LUN, LUN
       OpenW, LUN, OutFileName
;
       StatVar = Size( Text )
       NLines  = StatVar(1)-1
;
       For i = 0, NLines do PrintF, LUN, Text(i)
;
       Close, LUN
       Free_LUN, LUN
;
       Print, " Help information for "+StrMid( InFileName, 0, NameEnd )+ $
              " was written to file "+OutFileName+"."
       Print, " "

    EndIf Else begin
;
       Message, " No help was found for routine "+ $
                StrMid( InFileName, 0, NameEnd )+ ".", /NoName, /NoPrefix
;
    EndElse
;
    Return
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


