 PRO DispWEX, Name
;+
;   Description
;   -----------
;   DispWEX loads in a WEX format file with LoadWEX and displays any
;   text in it with UScroll.  If the file contains a menu, the menu
;   is presented to the user by using UMenu.  When the user select an
;   item from such a menu, Action is called to process the associated
;   command or file.
;
;   Calling Sequence:
;
;        DispWEX, Name
;
;   where Name is an input string variable containing the WEX com-
;   ponent name.
;
;#
;   Called by
;   ---------
;   Action
;
;   Routines Called
;   ---------------
;   Action        InitWEX      LoadWEX
;   UScroll       UMenu        WEXName
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  August 1992
;   01/07/93  KRT  Improved WEX component handling using index and state
;                  transition methods.
;   01/14/93  KRT  Change XSize, YSize, YPos, and XPos handling.  Default
;                  it used generally for XSize. (SPR 10449)
;                  Also changed I/O error message to "not found." instead
;                  of "not understood," giving full file name.
;                  Conditioned screening clearing to only occur if not
;                  in a GUI (inspired by SPR 10449).
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;
;.Title
;   Procedure DispWEX
;-
;
 On_Error,2   ; Upon an error, return to the calling routine.
;
@wexdisp.inc  ; Access to screen command and special characters
@wexpath.inc  ; Access to WEX file extensions, names, and areas
@wexcntrl.inc ; Access to the WEX stack and recursion control
;
;Begin,
    XEnv = (!D.Name eq 'X') OR (!D.Name eq 'WIN')
;
;   first put the current window onto the WEX stack:
;
    InitWEX
;
    Stat = Size( WEXStack )
    If (Stat(0) eq 0) then WEXStack = [ StrLowcase( Name ) ] $
                      else WEXStack = [ WEXStack, StrLowcase( Name ) ]
    NWEXs    = NWEXs + 1
;                                                             ; DEBUG
;    Print, "WEXStack in DispWEX: ("+Name+")"                 ; CODE
;    Print, "-----------------------------------"             ;
;    For ii = 0, NWexs-1 do Print, ii,".)", WEXStack(ii)      ;
;    Dum = ""                                                 ;
;    Read, Dum                                                ;
;
;   then load in the information from the WEX file,
    LoadWEX, Name, Text, Heading, Subtopics, SubFile, ErrStat
;
    If (ErrStat ne 0) then begin
;      handling file error exception,
;
       If (not XEnv) then Print, EraseScrn+CursrHome $
                     else Print, ' '
;
       BadName = WEXName( Name )
       If (StrTrim( BadName, 2 ) eq '') then BadName = Name
;
       Print, " "+BadName+" was not found."
       Print, " Please check the spelling of your specification."
       Print, " "
       Print, " Press return to continue."
       Hold = Get_Kbrd( 1 )
;
       NWEXs    = NWEXs - 1
       If (NWEXs gt 0) then WEXStack = WEXStack( IndGen( NWEXs ) ) $
                       else WEXStack = ""
;
       Return
    EndIf
;                                          THE PLAN
;   Initialize        Defaults for these variables will be eventually set up
;   Variables         in LoadWEX.  The user could specify values in WEX files.
;   ---------         These values would be interpreted by LoadWEX, which would
;                     load them in WEXStack in WEXCntrl.Inc.
;
    MaxHgt   = 24   ; Maximum height of a scrolling WEX in lines
;   ExtraX   =  4   ; Total size of right and left sides of WEX in characters
    ExtraY   =  4   ; Total size of top and bottom side of WEX in lines
;
    NLines   = N_Elements( Text )
    NItems   = N_Elements( Subtopics )
;
    VarStat  = Size( Heading )
    If (VarStat(1) ne 7) then Heading = ""
;
    Text_Comp = "T"   ; This is to be moved to InitWEX
    Menu_Comp = "M"   ; 
;
    If (NLines gt 0) then $
       If (N_Elements( WEXComps ) eq 0) then WEXComps = [ Text_Comp ] $
                                        else WEXComps = [ WEXComps, Text_Comp ]
    If (NItems gt 0) then $
       If (N_Elements( WEXComps ) eq 0) then WEXComps = [ Menu_Comp ] $
                                        else WEXComps = [ WEXComps, Menu_Comp ]
;
    LastComp    = N_Elements( WEXComps ) - 1
    Complex_WEX = (LastComp gt 0)
;
    XSize = IntArr( LastComp+1 )    ; Width of the WEX in characters
    YSize = IntArr( LastComp+1 )    ; Height of the WEX in lines
    XPos  = IntArr( LastComp+1 )+1  ; Col of the left side of the WEX
    YPos  = IntArr( LastComp+1 )+1  ; Row of the top of the WEX
;
    Menus   = Where( WEXComps eq Menu_Comp )
    Scrolls = Where( WEXComps eq Text_Comp )
;
    If (Complex_WEX and (NLines gt 0)) then begin
;      cascading menus on top of text...
       XPos( Menus ) = 5
       YPos( Menus ) = 5
    EndIf
;
    If (NLines gt 0) then $
       If (NLines gt MaxHgt-ExtraY) then YSize( Scrolls ) = MaxHgt $
                                    else YSize( Scrolls ) = NLines+ExtraY
;
;   (A similar algorithm as above should replaced by one in LoadWEX.
;    WEXComps would then be passed globally through WEXStack as an array.)
;   ---------
;
;   Realizing a WEX for this structure:
;
    If (not XEnv) then Print, EraseScrn+CursrHome, Format="(A,$)"
;
    indx     = 0
    EndLoop  = 0
;
    TQuit    = 0
    Top      = 0
;
    Sel      = 1
;
    While ((not EndLoop) and (not Quit) and (CollapseWEX eq "")) do begin
;      looping to receive and act upon selections from the user:
;
;      SCROLLING WINDOW
;      ----------------
       If ((WEXComps(indx) eq Text_Comp) and not TQuit) then begin
;         displaying text associate with this window structure:
;
          If (not XEnv) then $
             If (Complex_WEX) then Print, CursrHome, Format="(A,$)" $
                              else Print, EraseScrn+CursrHome, Format="(A,$)"
;
          UScroll, Text, Title=Heading,     Quit=TQuit,        $
                         YSize=YSize(indx), YPos=YPos(indx),   $
                         INIT=Top,                             $
                         RETAIN=Complex_WEX
;
          If (Complex_WEX) then indx = indx + 1 $
                           else EndLoop = 1
;
          If (TQuit) then begin
;            collapsing entire stack and return to main, since the user has
;            aborted in UScroll, 
             CollapseWEX = "quit"
             EndLoop = 1
          EndIf
       EndIf
;
;      MENU
;      ----
       If ((WEXComps(indx) eq Menu_Comp) and not TQuit) then begin
;
          If (NLines gt 0) then Items = [ "Subtopics", Subtopics ] $
                           else Items = [  Heading,    Subtopics ]
;
          If (not XEnv and not Complex_WEX) then  $
             Print, EraseScrn+CursrHome, Format="(A,$)"
;
          Sel = UMenu( Items, Title=0,         Init=Sel,         $
                              XPos=XPos(indx), YPos=YPos(indx) )
;
          If (not XEnv and not Complex_WEX) then  $
             Print, EraseScrn+CursrHome, Format="(A,$)"
;
          ActStr = StrLowcase( StrTrim( SubFile( Sel-1 ), 2 ) )
;
          If (Sel le 0) then begin
             If (not XEnv) then Print, EraseScrn+CursrHome, Format="(A,$)" $
                           else Print, ' '
             Message, " Error encountered in selection menu!"
          EndIf
;
          Case ActStr of
             "exit"   : Exit               ; return to the operating system!!!
             "quit"   : Quit    = 1        ; return to IDL!!
             "return" : EndLoop = 1        ; return to text or the last menu!
             "next"   : indx    = ((indx+1) mod LastComp)
;                                          ; jump to next WEX component
             "prev"   : indx    = ((LastComp+indx-1) mod LastComp)
;                                          ; jump to previous WEX component
              Else    : Action, ActStr     ; do action
          EndCase
;
          If (EndLoop and not XEnv) then  $
             Print, EraseScrn+CursrHome, Format="(A,$)"
;
;         The following statement stops the recursion collapse:
          If (CollapseWEX eq Name) then CollapseWEX = ""
       EndIf
    EndWhile
;
    NWEXs    = NWEXs - 1
    If (NWEXs gt 0) then WEXStack = WEXStack( IndGen( NWEXs ) ) $
                    else WEXStack = ""
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


