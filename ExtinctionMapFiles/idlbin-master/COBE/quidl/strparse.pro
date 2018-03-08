 Function StrParse, Strings, Delim=Delim
;+
;   Description
;   -----------
;   StrParse returns an array of non-blank entities which are extracted
;   from the input string Strings.
;
;#
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  February 1992
;
;   Date          SPR
;   05-Jan-1993   10419	  Kryszak-Servin  allows any delimiter
;
;.Title
;   Function StrParse
;-
;Begin,
;
    If not(Keyword_Set(Delim)) then delim = ' '
;
    WorkString = StrTrim( Strings, 2 )
    WorkLen    = StrLen( WorkString )
;
    If (WorkLen gt 0) then begin
;
       DelimPos = StrPos( WorkString, delim )
;
       If (DelimPos gt 0) then begin
          FirstPart   = StrMid( WorkString, 0, DelimPos )
          Remaining   = StrMid( WorkString, DelimPos+1, WorkLen )
          StringArray = [ FirstPart, StrParse( Remaining, delim=delim ) ]
       EndIf Else begin
          StringArray = WorkString
       EndElse
;
    EndIf
;
    Return, StringArray
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


