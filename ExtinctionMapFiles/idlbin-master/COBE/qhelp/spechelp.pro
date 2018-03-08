 PRO SpecHelp, Topics
;+
;   Description
;   -----------
;   SpecHelp takes a string of related topic files and loads them
;   by calling itself recursively.  Each time a parent topic is loaded,
;   the parent topic's related topic list search for the child topic.
;   If is found the recursion continues until the end of the string at
;   which DispHelp is called to display the help information.  If a
;   child topic is not found, SpecHelp indicates that such subtopic is
;   not found under the parent topic.
;
;   Note: It is assumed that the all parent topics are found in .Hlp
;         files, whereas the lowest level child is a .PRO file.  This
;         routine is used currently to read procedure source code prologs.
;#
;   Called by
;   ---------
;   SpcShell
;   SpecHelp (recursion)
;
;   Routines Called
;   ---------------
;   Action         LoadWEX       SpecHelp
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  February 1992
;
;   11/16/92  KRT  Documentation extended per SPR 10216 and integrated into
;                  the WEX system.
;
;   01/07/92  KRT  Changed DispHelp to the WEX routine Action.
;
;.Title
;   Procedure SpecHelp
;-
 On_Error,2   ; Return to caller if an error occurs
;
@wexdisp.inc  ; Common access to display control escape sequences
@wexpath.inc  ; Common access to file path and extensions
;
;Begin,
;
    TopicList = [ Topics ]
    VarStat   = Size( TopicList )
    NTopics   = VarStat(1)
;
    If (NTopics ge 2) then begin
;
       Parent  = TopicList(0)
       Child   = TopicList(1)
;
       LoadWEX, Parent+HlpExt, Text, Heading, Subtopics, SubFile, ErrStat
;   
       If (ErrStat ne 0) then begin
;         handling file error exception,
;
          Print, EraseScrn+CursrHome
;
          Print, " Information for "+Parent+" was not found."
          Print, " Please check the spelling."
          Print, " "
          Print, " Press return to continue."
          Hold = Get_Kbrd( 1 )
;
          Return
;
       EndIf
;
       Key      = StrLowCase( Child )
       KeyLen   = StrLen( Key )
       Selected = Where( Key eq StrLowCase( StrMid( SubFile, 0, KeyLen ) ) )
;
       VarStat  = Size( Selected )
       NSelectd = VarStat(1)
;
       Case 1 of
;
          (NSelectd eq 1) : $
             begin
                If (Selected(0) ge 0) then begin
;
;                  FirstTopic = Parent
                   Remaining  = TopicList( IndGen( NTopics-1 )+1 )
                   SpecHelp, Remaining
;
                EndIf Else begin
;
                   Print, " Subtopic " + StrUpCase( Child ) + $
                          " is not found under " + StrUpCase( Parent )
                   Print, " "
                   Print, " Press the return key to continue . . ."
;
                   Dummy = ''
                   Read, Dummy
;
                EndElse
             end
;
          (NSelectd gt 1) : $
             begin
;
                Print, " " + StrUpCase( Child ) + " is not a unique term."
                Print, " Please try again spelling out further."
                Print, " "
                Print, " Press the return key to continue . . ."
;
                Dummy = ''
                Read, Dummy
;
             end
       EndCase
;
    EndIf Else begin
;
       Action, TopicList(0)+ProExt
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


