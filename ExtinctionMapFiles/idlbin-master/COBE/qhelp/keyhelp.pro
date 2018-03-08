 PRO KeyHelp, SearchKey
;+
;   Description
;   -----------
;   KeyHelp is the main driver for a key word search capability. It uses
;   the routine KSearch to perform a simple search of the file specified
;   by the global WEX variable Keyfiles.  The results are presented to the
;   user in the form of a list using UMenu.  If an item is selected from
;   the menu, then Action is used to display a file associated with the
;   key.
;
;   Note: This routine is specialized to search on a list of routines names.
;         The .PRO extention is appended automatically to the first name
;         found in the file record using the global WEX variable found in
;         WEXPath.Inc.
;#
;   Called by
;   ---------
;   UHelp
;
;   Routines Called
;   ---------------
;   Action         LoadWEX       KSearch     UMenu
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  June 1992
;
;   11/16/92  KRT  Integrated into the WEX system (SPR 10216).
;   12/21/92  KRT  Changed permanent menu items to lower case and QUIT to
;                  Exit UHelp
;   12/30/92  KRT  Reorganized to run as a stand-alone program from a
;                  parameter file.  Used Action instead of DispHelp.
;                  Modified search to look though multiple files.
;
;.Title
;   Procedure KeyHelp
;-
 On_Error,2   ; Return to caller if an error occurs
;
@wexdisp.inc  ; Access to display control command sequences
@wexpath.inc  ; Access to global WEX file name and path (e.g., KeyFile)
@wexcntrl.inc ; Access to global WEX recursion control
;
;Begin,
;
    LastFile = N_Elements( KeyFiles )-1
    If (LastFile lt 0) then Message, "The variable KEYFILES was not defined."
;
;   Load in the files to be searched :
;
    First = 1
    For i = 0, LastFile do begin
;
       LoadWEX, KeyFiles(i), ProDirectory, Heading, Subtopics, SubFile, ErrStat
;
       If (ErrStat ne 0) then begin
;         handling file error exception,
;
          Print, EraseScrn+CursrHome
          Print, " Error occurred routine directory file "+MainHelp+" ..."+  $
                 " check variable array KEYFILES, element (", i ,")."
          Print, " "
          Message, " Key word search terminated."
       EndIf
;
       If (First) then begin
          Text = [ ProDirectory ]
          First   = 0
       EndIf Else begin
          Text = [ Text, ProDirectory ]
       EndElse
    EndFor
;
    If (N_Params() le 0) then begin
;      getting the keyword from the user:
;
       Print, EraseScrn+CursrHome
       Print, Format = '($," Enter keyword ")'
;
       KeyWord = ''
       Read, KeyWord
       KeyWord = StrTrim( KeyWord, 2 )
;
    EndIf Else begin
;      taking the keyword from the parameter list:
;
       KeyWord = StrTrim( SearchKey, 2 )
;
    EndElse
;
;
    Quit     = 0
    Extra    = 3
;
    KeyItem  = "Try another key"
    BackItem = "Return to MAIN MENU"
    QuitItem = "Exit UHelp"
;
    Sel     = 1
    Back    = (KeyWord eq '') ; return if nothing was entered above.
    NewKey  = 0
;
    While ((not Quit) and (not Back)) do begin
;
       If (KeyWord ne '') then begin
;
          KSearch, KeyWord, Text, ProList, ProName
          KeyWord = ''
;
          If (ProList(0) ne '') then NItems = N_Elements( ProList )+Extra+1 $
                                else NItems = Extra
;
          If (NItems gt Extra) then begin
             Instruct = "KEY SEARCH: Select a Routine"
             Items    = [ Instruct, ProList, ' ', KeyItem, BackItem, QuitItem ]
          EndIf Else begin
             Instruct = "KEY SEARCH: No Routine Found"
             Items    = [ Instruct, KeyItem, BackItem, QuitItem ]
          EndElse
;
          SelKey  = NItems-2
          SelBack = NItems-1
          SelQuit = NItems
       EndIf
;
       Print, CursrHome + EraseScrn ; Clear the screen.
       Sel    = UMenu( Items, Title=0, Init=Sel )
;
       Back   = (Sel eq SelBack)
       Quit   = (Sel eq SelQuit)
       NewKey = (Sel eq SelKey)
;
       If (Sel le 0) then begin
          Print, CursrHome + EraseScrn
          Message, " Error encountered in Keyword Search Menu."
       EndIf
;
       Case 1 of
          Back   : ; return to calling procedure...
          Quit   : ; return with condition to end program...

          NewKey : begin
             Print, EraseScrn+CursrHome
             Print, Format = '($," Enter keyword ")'
;
             KeyWord = ''
             Read, KeyWord
             KeyWord = StrTrim( KeyWord, 2 )
          end
;
          Else   : begin
;            moving recursively through the help tree :
;
             Print, CursrHome + EraseScrn
             Name = ProName( Sel-1 )+ProExt
;
             Action, Name
          end
       EndCase
;
    EndWhile
;
    If (not Quit) then Print, EraseScrn+CursrHome
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


