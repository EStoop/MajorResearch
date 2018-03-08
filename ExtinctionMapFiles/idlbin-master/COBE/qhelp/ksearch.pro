 PRO KSearch, KeyWord, Text, Sentence, FirstWord
;+
;   Description
;   -----------
;   KSearch searches the string array Text for the first occurance of
;   the string in found in the argument KeyWord.  KSearch returns the
;   element of Text containing the string through the argument Sentence
;   and the first word of Sentence through the argument FirstWord.
;
;   [This routine was originally KeySearch created by K.Turpie in June
;    of 1992]
;#
;   Called by
;   ---------
;   KeyHelp
;
;   Routines Called
;   ---------------
;   <none>
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  November 1992
;
;.Title
;   Procedure KSearch
;-
 On_Error,2   ; Return to caller if an error occurs
;
;Begin,
;
    Sentence  = [ '' ]
    FirstWord = [ '' ]
;
    StartUp   = 1
    VarStat   = Size( Text )
;
    LastItem  = VarStat(1)-1
    LowerKey  = StrLowCase( KeyWord )
;
    For i = 1, LastItem do begin
;
       WorkStr  = StrTrim( Text(i), 2 )
       LowerStr = StrLowCase( WorkStr )
;
       If (StrPos( LowerStr, LowerKey ) gt -1) then begin
;
          First = StrMid( WorkStr, 0, StrPos( WorkStr, " " ) )
;
          If (StartUp) then begin
;
             Sentence  = [ WorkStr ]
             FirstWord = [ First   ]
             StartUp   = 0
;
          EndIf Else begin
;
             Sentence  = [ Sentence,  WorkStr ]
             FirstWord = [ FirstWord, First   ]
;
          EndElse
       EndIf
    EndFor
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


