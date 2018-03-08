 PRO UHelp, ProName, FILE=File, KEY=Key
;+
;   Description
;   -----------
;   UHelp is the main driver to the UHelp program.  When called from the
;   command line, the following options can be used.
;
;             IDL>UHELP
;                 this invokes the WEX system (Action) to provide menu
;                 supported help.
;
;             IDL>UHELP, "<procedure name>"
;                 this invokes WEX system (Action) directly to display the
;                 help for the procedure specified between the quotes.
;
;             IDL>UHELP, "<procedure name>", FILE="<output file name>"
;                 this loads the help information for the procedure
;                 specified between the quotes and prints it to the file
;                 specified by name, also between quotes.
;
;             IDL>UHELP, "<procedure name>", KEY="<keyword>"
;                 this does a keyword search on the UIDL procedure
;                 dictionary (found in oneliner.pro).
;#
;   Called by
;   ---------
;   Main (UIDL)
;
;   Routines Called
;   ---------------
;   Action       Extract      KeyHelp       SpecHelp
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  February 1992
;
;   11/16/92  KRT  Integrated into the WEX system and extended documentation
;                  per SPR 10216.
;   12/30/92  KRT  Removed dependency on the old routines DispHelp and MenuHelp,
;                  switching entirely to using the WEX routine Action.
;   01/14/93  KRT  Modified to clear screen when the user quits UHelp ONLY IF
;                  the user is not in a GUI. (SPR 10449)
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;
;.Title
;   Procedure UHelp
;-
 On_Error,2   ; Return to caller if an error occurs
;
@wexdisp.inc  ; Access to display control command sequences
@wexpath.inc  ; Access to global help file name and path
@wexcntrl.inc ; Access to WEX global control variables
;
;Begin,
    SaveQuiet  = !Quiet
    !Quiet     = 1
;
    InitWEX, /Again
;
;   Identify routine dictionary file,
;
    KeyFile    = "OneLiner.Pro"
;
    If (N_Elements( ProName )) then begin
;
       VarStat = Size( ProName )
;
       If (VarStat(0) ne 0) then begin
          !Quiet = SaveQuiet 
          Message, "ERROR: UHelp expects a scalar variable as an argument."
       EndIf
;
       If (VarStat(1) ne 7) then begin
          Print, "Please type UHelp, followed by a string with quotes."
          Print, " "
          Print, "Example:"
          Print, "--------"
          Print, " "
          Print, "IDL>UHelp, 'PLANCK'"
          Print, "... displays help for routine PLANCK"
          Print, " "
          Print, " "
          !Quiet = SaveQuiet 
          Return
       EndIf
;
       If (StrUpCase( ProName ) eq 'UHELP') then begin
          Print, "Syntax:"
          Print, "-------"
          Print, " "
          Print, "    UHELP, ['<routine name>' [,FILE='<output file name>'] ]"
          Print, " "
          Print, "Examples:"
          Print, "---------"
          Print, " "
          Print, "IDL>UHelp
          Print, "... puts you into help menu display."
          Print, " "
          Print, "IDL>UHelp, 'PLANCK'"
          Print, "... displays help for routine PLANCK"
          Print, " "
          Print, "IDL>UHelp, 'PLANCK', File='PLANCK.HLP"
          Print, "... sends help info to file PLANCK.HLP"
          Print, " "
          Print, " "
          !Quiet = SaveQuiet 
          Return
       EndIf
;
       Case 1 of
;
          (KeyWord_Set( FILE )) : begin
             Extract, ProName+ProExt, File, ErrStat
          end
;
          (KeyWord_Set( KEY )) : begin
             KeyHelp, ProName
          end
;
          Else : begin
             Print, EraseScrn+CursrHome
;
             SFName = ProName + ProExt
             Action, SFName
             If ((Quit) and (!D.Name ne 'X') AND (!D.NAME NE 'WIN')) $
                then Print, EraseScrn+CursrHome
          end
;
       EndCase
;
    EndIf Else begin
;
       If (KeyWord_Set( KEY )) then KeyHelp, Key $
                               else begin
                                  Action, "Main.Hlp"
                                  If ((Quit) and (!D.NAME NE 'X') AND $
                                    (!D.NAME NE 'WIN')) $
                                     then Print, EraseScrn+CursrHome
                               endelse
;
    EndElse
;
    !Quiet = SaveQuiet 
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


