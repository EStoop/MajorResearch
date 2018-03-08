 PRO  DispText,  Text,    Width,     Top,      Bottom
;+
;   Description
;   -----------
;   DispText display the text array Text on the screen.  Top and Bottom
;   are array indices which specify the range of Text that is to be
;   show to the user.
;
;#
;   Called by
;   ---------
;   TScroll
;
;   Routines Called
;   ---------------
;   MovToCol
;   MovToRow
;   InitWEX
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  February 1992
;
;   11/23/92  KRT  Documented and integrated into the WEX system (SPR 10216)
;
;.Title
;   Procedure DispText
;-
;
 On_Error,2   ; Return to caller if an error occurs
;
@wexdisp.inc  ; Access to screen control commands
@wexcntrl.inc ; Access to scroll window parts
;
;Begin,
;
    InitWEX
;
;   Abbreviate screen display command sequence variables :
;
    wb  = BoldText        ; Display bold characters
    wl  = Underline       ; Display underlined characters
    wr  = RevVideo        ; Display reverse video text
    wz  = NormVideo       ; Set all text attributes to normal
;
    cs  = GrphicsON       ; Graphics text on
    cz  = TextON          ; Normal text mode on
;
;   Abbreviate cursor movements command sequence variables :
;
    cuu = CursrUp         ; Move cursor up
    cud = CursrDwn        ; Move cursor down
;
    VarStat = Size( Text )
    NLines  = VarStat(1)
;
;   Set-up commands for left and top margins, operand for a horzontal
;   line command, and build blank lines for filling in the bottom of
;   the window (if necessary) :
;
    TM = MovtoRow( CurrRow )
    LM = MovtoCol( CurrCol )
;
    HL = 'q'
    For i = 2, Width+2 do HL = HL + 'q'
;
    Spaces = ' '
    For i = 2, Width do Spaces = Spaces + ' '
;
;   Center title string for display :
;
    TtlWdth = StrLen( WEXTitle )
;
    Padding = (Width - TtlWdth) / 2
    If (Padding lt 0) then Padding = 0
;
    TtlStr  = Spaces
    StrPut, TtlStr, WEXTitle, Padding
    TtlStr  = StrMid( TtlStr, 0, Width )
;
;   Draw top of window and window title :
;
    Print, TM+LM+cs+'l'+HL+'k'+cz, Format='(A)'
    Print, LM+cs+'x'+cz+' '+wr+TtlStr+wz+' '+cs+'x'+cz, Format='(A)'
;
;   Draw text with verticle lines on either side :
;
    For i = Top, Bottom do Print, LM+cs+'x'+cz+' '+Text(i)+' '+cs+'x'+cz, $
                                  Format='(A)'
;
;   Draw blank text lines to fill the remainder of the window :
;
    Remaindr = Bottom - Top - NLines + 1
    For i = 1, Remaindr do Print, LM+cs+'x'+cz+' '+Spaces+' '+cs+'x'+cz, $
                                  Format='(A)'
;
;   Center message for message bar :
;
    MsgWdth = StrLen( WEXMessage )
;
    Padding = (Width - MsgWdth) / 2
    If (Padding lt 0) then Padding = 0
;
    MsgStr  = Spaces
    StrPut, MsgStr, WEXMessage, Padding
    MsgStr  = StrMid( MsgStr, 0, Width )
;
;   Draw message bar and bottom of the window :
;
    Print, LM+cs+'x'+cz+' '+wr+MsgStr+wz+' '+cs+'x'+cz, Format='(A)'
    Print, LM+cs+'m'+HL+'j'+cz+cuu+cuu, TxtCursrOff, Format='(2A)'
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


