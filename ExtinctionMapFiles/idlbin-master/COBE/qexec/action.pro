 PRO Action, Command
;+
;   Description
;   -----------
;
;   Calling Sequence:
;
;       Action, Command
;
;   where, Command is an input variable of type String.
;
;   Action process the input string found in the argument Command.  The
;   method of processing depends on the first character of the string.
;
;   If the first character is one of three special symbols ($,*, or #)
;   then the remained of the string is executed in a way associate with
;   the symbol.  The association goes as follows:
;
;   $ - spawns the string
;   * - executes the string with the IDL Execute routine
;   # - executes the string with the IDL Call_Procedure routine
;
;   If the first symbol is not one of these, the string is assumed to
;   be a WEX parameter file name and calls DispWEX to process it.
;
;#
;   Called by
;   ---------
;   CGIS
;   DispWEX
;
;   Routines Called
;   ---------------
;   DispWEX
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  August 1992
;   SPR 11003  Jun 02 93   Change !version.os to !cgis_os.  J Newmark
;   SPR 11127  Oct 18 93   IDL for Windows Compatability. J. Newmark
;
;.Title
;   Procedure Action
;-
 On_Error,2     ; Return to caller if an error occurs
;
@wexdisp.inc    ; Access to display control escape sequences
@wexcntrl.inc   ; Access to menu stack
;
;Begin,
;
    InitWEX ; Start up WEX system if not already started
;
;    Print, "WEXStack in Action: ("+Command+")"             ; DEBUG
;    Print, "-----------------------------------"           ; CODE
;    For ii = 0, NWexs-1 do Print, ii,".)", WEXStack(ii)    ;
;    Dum = ""                                               ;
;    Read, Dum                                              ;
;
    If (CollapseWEX eq "") then begin
;      by checking WEX stack,
;
       windex = where( StrLowcase( Command ) eq WEXStack )
;
       If (windex(0) lt 0) then begin
;         processing since the command doesn't appear on the WEX stack:
;
          SpwnMark = "$"
          ExecMark = "*"
          ProcMark = "#"
          ComStr   = StrTrim( Command, 2 )
          FirstCh  = StrMid( ComStr, 0, 1 )
;
          Case FirstCh of
;
             SpwnMark : Begin
                ComStr = StrMid( ComStr, 1, StrLen( ComStr ) )
                Case StrLowcase( !cgis_os ) of
                   "vms" :  Spawn, ComStr
                   "windows" :  Spawn, ComStr
                   "unix" :  Spawn, ComStr, /NoShell
                EndCase
             End
;
             ExecMark : Begin
                ComStr = StrMid( ComStr, 1, StrLen( ComStr ) )
                Status = Execute( ComStr )
             End
;
             ProcMark : Begin
                ComStr = StrMid( ComStr, 1, StrLen( ComStr ) )
                Call_Procedure, ComStr
             End
;
             Else : Begin
                ActStr = StrLowcase( ComStr )
;
                Case ActStr of
                   "RETURN" : Return                 ; return to the last menu!
                   "QUIT"   : Quit = 1               ; return to IDL!!
                   "EXIT"   : Exit                   ; return to the oper sys!!!
                    Else    : DispWEX, ActStr        ; display text!
                EndCase
;
             End
;
          EndCase
;
       EndIf Else begin
;         to collapse the recursion stack to where this command appeared before,
;
          CollapseWEX = Command
;
       EndElse
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


