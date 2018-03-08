 PRO SpcShell
;
 On_Error, 2 ; Return to caller if an error occurs
;
@wexdisp.inc
;
;Begin,
;   moving through tree using a user specified "path" :
;
    Print, EraseScrn+CursrHome
    Print, Format = '($," Enter specific topic ")'
;
    SpecificList = ''
    Read, SpecificList
    SpecificList = StrTrim( SpecificList, 2 )
;
    If (SpecificList ne '') then begin
                   TopicList = StrParse( SpecificList ) ; create an array of
                   SpecHelp, TopicList                  ; strings and process
    EndIf
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


