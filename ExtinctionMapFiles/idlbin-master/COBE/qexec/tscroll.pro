 PRO TScroll,   Text,    TITLE=Title,                 $
                         XPOS=XPos,   YPOS=YPos,      $
                         XSIZE=XSize, YSIZE=YSize,    $
                         INIT=Init,   QUIT=TQuit,     $
                         RETAIN=Retain,               $
                         SCROLL_HELP=Scroll_Help
;+
;   Description
;   -----------
;   TScroll creates a scrolling window for text at the user's terminal.
;   This is used primarily for VT200 or other devices without a GUI.
;
;   An example of the calling sequence is a follows:
;
;       TScroll, Text, Title=Title, Quit=TQuit, $
;                XPos=10, YPos=10, XSize=40, YSize=10
;
;   where,
;
;       INIT    = When set to an integer, TScroll will be display the
;                 text starting at the line specified by the value.
;                 INIT returns the value of the line at the top of
;                 screen when TScroll is exited.
;
;       RETAIN  = If set, the text display will remain after exiting
;                 TScroll.
;
;       TEXT    = A string array with each element containing the text
;
;       TITLE   = The title to appear on the top of the scrolling window
;
;       TQUIT   = Returns a signal to the calling routine that the user
;                 wishes to abort without clearing the screen.
;
;       XPOS    = The horizontal position of the left side of the window,
;                 as measured in characters for the left side of the
;                 screen.
;
;       YPOS    = The vertical position of the top of the window as
;                 measured in lines from the top of the screen.
;                  
;       XSIZE   =  The horizontal width of the window from side to side
;                  as measured in characters.
;
;       YSIZE   =  The vertical width of the window from top to bottom
;                  as measured in lines.
;
;       When the SCROLL_HELP keyword is set, TScroll will operate in help
;       mode.  This means that the ability to invoke another scrolling
;       help screen is turned off while the help is being viewed, prevent-
;       unnecessary recursion.
;
;#
;   Called by
;   ---------
;   DispHelp
;   DispWEX
;   TScroll (recursion)
;   UScroll
;
;   Routines Called
;   ---------------
;   DispText       InitWEX       LoadWEX       TErase       TScroll
;
;   History
;   -------
;   Created by K.Turpie,   General Sciences Corporation,  February 1992
;
;   09/14/92  KRT  New key actions added and default sizes changed.
;   10/26/92  KRT  Key actions change to no longer use the space bar.
;                  Message bell fixed.
;   11/22/92  KRT  Changed WEX globles Title and Message to WEXTitle and
;                  WEXMessage.
;                  Changed the keyword assignment Title=Heading to
;                  Title=Title.
;   12/02/92  KRT  Removed bell from TopMsg and BotMsg and send to terminal
;                  via Print command.  This eliminated the appearance of a
;                  reversed block character in the lower right corner of the
;                  screen every time the bell rang.
;   12/15/92  KRT  Reorganized code to use Bit7Command and Bit8Command
;                  booleans in checking keyboard entry.  This facilitated
;                  using the no-wait argument in Get_Kbrd and fixed a
;                  problem with responding to 8 bit commands.
;   12/21/92  KRT  Added /RETAIN and /INIT keywords.  Changed default width
;                  and height to fill a typical screen.
;   01/22/93  KRT  Put in error checking for file extraction per SPR 10478.
;                  Cleaned the output of all text commands before writing.
;
;.Title
;   Procedure TScroll
;-
 On_Error,2   ; Return to caller if an error occurs
;
@wexpath.inc  ; Access to help files names and file name parts
@wexdisp.inc  ; Access to screen display command sequences
@wexcntrl.inc ; Access to scrolling window parts
;
;Begin,
;
;   Initialize Variables :
;
    InitWEX
;
    Blank    = $
'                                                                              '
    ScrnWdth = 80              ; Width of screen
    ScrnHght = 24              ; Height of screen
;
    ExtraX   = 4               ; Account for spacing and sides of box
    ExtraY   = 4               ; Account for top, bottom, title bar, and
;                              ; message bar
;
    TxtMax   = ScrnWdth-ExtraX ; Maximum width of text
    DefWdth  = ScrnWdth-ExtraX ; Default width of scroll window
    DefHgt   = ScrnHght-ExtraY ; Default height of scroll window
;
    If (KeyWord_Set( Title )) then WEXTitle = Title $
                              else WEXTitle = " "
;
    GenMsg   = "Next Scrn=(Space)  Prev Scrn=(B)  Help=(H)  Exit=(Tab)"
    NormMsg  = GenMsg + "       "
    BottMsg  = GenMsg + "   END!"
    TopMsg   = GenMsg + "   TOP!"
;
;   Check the validity of what was supplied for the "TEXT" argument :
;
    S = Size(Text)
    If (S(0) ne 1) or (S(2) ne 7) then begin
       Message, 'An inappropriate argument was supplied.'
    EndIf
;
    NLines = S(1)
;
;   See if the vertical and horizontal dimensions of the window have
;   given (otherwise use defaults) :
;
    If (Keyword_Set( XSIZE )) then begin
       If (XSize gt ScrnWdth) then XSize = ScrnWdth
    EndIf Else begin
       XSize = DefWdth
    EndElse
;
    If (Keyword_Set( YSIZE )) then begin
       If (YSize gt ScrnHght) then YSize = ScrnHght
    EndIf Else begin
       YSize = DefHgt
       If (YSize gt (NLines+ExtraY)) then YSize = NLines+ExtraY
    EndElse
;
    If (XSize lt ExtraX+1) then XSize = ExtraX+1
    If (YSize lt ExtraY+1) then YSize = ExtraY+1
;
;   Use the value of "XPOS" (if supplied) to determine the amount of
;   initial spacing :
;
    If (Keyword_Set( XPOS )) then begin
       DiffX = ScrnWdth-XSize
       If (XPos gt DiffX) then XPos = DiffX
    EndIf Else begin
       XPos = 0
    EndElse
;
    CurrCol = XPos
;
;   Use the value of "YPOS" (if supplied) to determine the verticle place-
;   ment (YPOS=0 --> row=0, i.e., the upper row) :
;
    If (Keyword_Set( YPOS )) then begin
       DiffY = ScrnHght-YSize
       If (YPos gt DiffY) then YPos = DiffY
    EndIf Else begin
       YPos = 0
    EndElse
;
    CurrRow = YPos
;
    If (KeyWord_Set( INIT )) then Top = Init $   ;|
                             else Top = 0        ;|- text array pointers are
;                                                ;|  initialized here.
    Bottom = Top + YSize-ExtraY-1                ;|
;
    If (Top    lt        0) then Top    = 0
    If (Bottom gt Nlines-1) then Bottom = NLines-1
;
;
;   Fill up array "Q" with text, pad right with blanks to fill up window :
;
    Q = StrArr( NLines )
;
    MaxText = XSize - ExtraX
    For i = 0, NLines-1 do begin
;
       TransTab, Text(i),    TabFreeStr
       TransCom, TabFreeStr, SpecChrStr, Len, SpcLen
;
       Pad  = MaxText - Len
       Q(i) = StrMid( SpecChrStr, 0, MaxText+SpcLen ) $
            + StrMid( Blank,      0,            Pad )
    EndFor
;
;   Print the menu out on the terminal screen :
;
    Print, TxtCursrON
;
    WEXMessage = NormMsg
    DispText, Q, MaxText, Top, Bottom
;
;   Monitor and respond to key presses until a carriage return is sent :
;
    TQuit = 0
    Done  = 0
    While (not Done) do begin
;
       ArrowUp   = 0
       ArrowDown = 0
       PageUp    = 0
       PageDown  = 0
       GotoTop   = 0
       GotoEnd   = 0
       Help      = 0
       Extract   = 0
       Refresh   = 0
;
       Ch1 = StrUpCase( Get_Kbrd(1) )
;
       Bit8Command = (Ch1 eq CSI)
       If (Ch1 eq ESC) then Bit7Command = (Get_Kbrd(1) eq '[') $
                       else Bit7Command = 0
;
       If (Bit7Command or Bit8Command) then begin
;
          Ch2 = Get_Kbrd(1)
;
          Case 1 of
             (Ch2 eq 'A') : ArrowUp   = 1
             (Ch2 eq 'B') : ArrowDown = 1
             (Ch2 eq '5') : If (Get_Kbrd(1) eq '~') then PageUp   = 1
             (Ch2 eq '6') : If (Get_Kbrd(1) eq '~') then PageDown = 1
             Else:
          EndCase
;
       EndIf Else begin
;
          Case 1 of
             (Ch1 eq CntlB)      : PageUp    = 1
             (Ch1 eq 'B')        : PageUp    = 1
             (Ch1 eq ' ')        : PageDown  = 1
             (Ch1 eq  CR)        : PageDown  = 1
             (Ch1 eq 'T')        : GotoTop   = 1
             (Ch1 eq 'E')        : GotoEnd   = 1
             (Ch1 eq 'H')        : Help      = 1
             (Ch1 eq Tab)        : Done      = 1
             (Ch1 eq 'Q')        : TQuit     = 1
             (Ch1 eq CntlW) or $
             (Ch1 eq CntlR)      : Refresh   = 1
             (Ch1 eq CntlF)      : Extract   = 1
             Else:
          EndCase
;
       EndElse
;
       Case 1 of
          ArrowUp : Begin ; up one
             If (Top gt 0) then begin
                Top    = Top - 1
                Bottom = Top + YSize-ExtraY-1
                DispText, Q, MaxText, Top, Bottom
             EndIf Else begin
                WEXMessage = TopMsg
                DispText, Q, MaxText, Top, Bottom
                Print, Bell, Format="(A,$)"
                WEXMessage = NormMsg
             EndElse
          End
;
          ArrowDown : Begin ; down one
             If (Bottom lt NLines-1) then begin
                Bottom = Bottom + 1
                Top    = Bottom - (YSize-ExtraY-1)
                DispText, Q, MaxText, Top, Bottom
             EndIf Else begin
                WEXMessage = BottMsg
                DispText, Q, MaxText, Top, Bottom
                Print, Bell, Format="(A,$)"
                WEXMessage = NormMsg
             EndElse
          End
;
          PageUp : Begin ; Backward one page
             Bottom = Top
             Top    = Bottom - (YSize-ExtraY-1)
             If (Top lt 0) then begin
                Top        = 0
                Bottom     = (Top + (YSize-ExtraY-1)) < (NLines-1) ; JE 4/92
                Bottom     = Bottom
                WEXMessage = TopMsg
                DispText, Q, MaxText, Top, Bottom
                Print, Bell, Format="(A,$)"
                WEXMessage = NormMsg
             EndIf Else begin
                DispText, Q, MaxText, Top, Bottom
             EndElse
          End
;
          PageDown : Begin ; Forward one page
             Top    = Bottom
             Bottom = Top + (YSize-ExtraY-1)
             If (Bottom gt NLines-1) then begin
                Bottom     = NLines-1
                Top        = (Bottom - (YSize-ExtraY-1)) > 0       ; JE 4/92
                WEXMessage = BottMsg
                DispText, Q, MaxText, Top, Bottom
                Print, Bell, Format="(A,$)"
                WEXMessage = NormMsg
             EndIf Else begin
                DispText, Q, MaxText, Top, Bottom
             EndElse
          End
;
          GotoTop : Begin ; Go to top page
             Top        = 0
             Bottom     = YSize-ExtraY-1
             WEXMessage = TopMsg
             If (Bottom gt NLines-1) then Bottom = NLines-1
             DispText, Q, MaxText, Top, Bottom
             Print, Bell, Format="(A,$)"
             WEXMessage = NormMsg
          End
;
          GotoEnd : Begin ; Go to end page
             Bottom     = NLines-1
             Top        = Bottom - (YSize-ExtraY-1)
             WEXMessage = BottMsg
             If (Top lt 0) then Top = 0
             DispText, Q, MaxText, Top, Bottom
             Print, Bell, Format="(A,$)"
             WEXMessage = NormMsg
          End
;
          Refresh : Begin ; Repaint the display
             Print, TxtCursrON + CursrHome
             WEXMessage = NormMsg
             DispText, Q, MaxText, Top, Bottom
          End
;
          Extract : Begin ; Extract WEX file to user's directory
             Print, TxtCursrON + EraseScrn + CursrHome
;
             Print, Format = '($,"Enter output file name")'
             FileName = ""
             Read, FileName
;
             If (StrTrim( FileName, 2 ) ne "") then begin
;
                Get_LUN, LUN
                OpenW, LUN, FileName, Error=IOS
;
                If (IOS eq 0) then begin
;
                   For i = 0, NLines-1 do begin
                      TransTab, Text(i),    TabFreeStr
                      TransCom, TabFreeStr, PlainText, Len, SpcLen, /Plain
                      PrintF, LUN, PlainText
                   EndFor
;
                   Print, " "
                   Print, " Text was written to file "+FileName+"."
                   Print, " "
                   Print, Format = '($," Press any key to continue...")'
                   Answer = Get_Kbrd( 1 )
;
                EndIf Else begin
;
                   Print, " "
                   Print, " File " + FileName + " cannot be opened."
                   Print, " Sorry, but no text was written."
                   Print, " "
                   Print, Format = '($," Press any key to continue...")'
                   Answer = Get_Kbrd( 1 )
;
                EndElse
;
                Close, LUN
                Free_LUN, LUN
;
             EndIf
;
             Print, EraseScrn + TxtCursrON + CursrHome
             WEXMessage = NormMsg
             DispText, Q, MaxText, Top, Bottom
          End
;
          Help : Begin
;
             If (not KeyWord_Set( Scroll_Help )) then begin
;
                LoadWEX, ScrollHelp, HelpText, HelpTitle,  $
                         HelpTopics, HelpFiles, ErrStat
;
                TScroll, HelpText,                     $
                         Title=HelpTitle, Quit=HQuit,  $
                         XSize=XSize,     YSize=YSize, $
                         YPos=YPos,       /Scroll_Help
;
                If (KeyWord_Set( Title )) then WEXTitle = Title $
                                          else WEXTitle = " "
;
                WEXMessage = NormMsg
                DispText, Q, MaxText, Top, Bottom
             EndIf
          End
;
          TQuit : Begin
             Print, TxtCursrON  ; Enable text cursor
             Done = 1
          End
;
          Done : Begin ; Erase window and return
             If (not KeyWord_Set( Retain )) then TErase, XSize, YSize
             Print, TxtCursrON  ; Enable text cursor
          End
;
          Else :
       EndCase
;
    EndWhile
;
    Init = Top
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


