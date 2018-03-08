 PRO InitWEX, AGAIN=Again
;+
;   Description
;   -----------
;   InitWEX initializes WEX system.  WEX stack and control variables are
;   defined and initialized.  File and directory names are defined for
;   global accessibility throughout the WEX system and its various compo-
;   nents.  Also, terminal control sequences are setup and passed through
;   a call to TextGrph.
;#
;   Called by
;   ---------
;   CGIS
;   UImage
;   UHelp
;   UScroll
;
;   Routines Called
;   ---------------
;   TextGrph
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  November 1992
;
;   12/21/92  KRT  Added keyword /AGAIN for forcing reinitialization of
;                  WEXStack and other control variables.
;
;   12/28/92  KRT  Changed the ____Path variables to arrays of paths and
;                  called them ____Paths.  WEXName uses each element of
;                  array WEXPaths to search for a WEX file.
;
;                  Changing the names facilitated moving the reinitialization
;                  block around the entire procedure (except Quit = 0).
;
;                  Used N_Elements to test for previous initialization instead
;                  of Size.
;
;   01/06/92  KRT  Change Oneliner.Pro to Oneliner.Hlp and put a / at the
;                  beginning of the default UNIX path specification.
;
;   01/03/93  DB   Change OSPaths to include doc subdirectories
;                  SPR 10630
;   04/07/93  DB   Change OSPaths to reference a relative path rather
;                  than an absolute path  SPR 10805
;
;   05/28/93  KRT  Added global variable defaults for mail and editor
;   SPR 11003  Jun 02 93   Change !version.os to !cgis_os.  J Newmark
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;
;.Title
;   Procedure InitWEX
;-
 On_Error,2    ; Return to caller if an error occurs
;
@wexcntrl.inc  ; Access to WEX stack and recursion control variables
@wexdisp.inc   ; Access to display control command sequences
@wexpath.inc   ; Access to global help file name and path
;
;Begin,
;
    Force = KeyWord_Set( Again )
    Quit  = 0
;
;   Check to see if the WEX System has be initialized yet:
;
    Not_Setup_Yet = (N_Elements( NWEXs ) eq 0)
;
    If ((Not_Setup_Yet) or Force) then begin
;
;      Initialize WEX stack and recursion variables:
;
       NWEXs       = 0
       CollapseWEX = ""
       WEXStack    = ""
;
;      Set-up display control command sequences:
;
       TextGrph
;
;      Set-up default file path and extensions:
;
       Case StrLowcase( !cgis_os ) of
          'vms' : OSPaths = [ "#CGIS_IDL#.*]" ]
          'windows' : OSPaths = [ " ","C:\RSI\IDL_5\CGIS\QHELP\",$
            "C:\RSI\IDL_5\CGIS\QUIDL\","C:\RSI\IDL_5\CGIS\QUIMAGE"]
          'unix' : OSPaths = [ " ", "$csrc/*/", $
                                   "$csrc/doc/*/" ]
       EndCase
;
       If (N_Elements( WEXPaths ) le 0) then WEXPaths = [ OSPaths ]
;
;      Set-up default mail facility:
;
       DefSysV, '!VMSMail',  'mail'
       DefSysV, '!UNIXMail', 'vmail'
       DefSysV, '!PCMail',   'mail'
;
;      Set-up default editors:
;
       DefSysV, '!VMSEdit',  'edit'
       DefSysV, '!UNIXEdit', 'emacs'
       DefSysV, '!PCEdit',   'edit'
;
;      Other WEX system constants:
;
       HlpExt     = '.hlp'
       ProExt     = '.pro'
       MainHelp   = 'main'+HlpExt
       ScrollHelp = 'scroll'+HlpExt
;
       If (N_Elements( Keyfiles ) le 0) then $
       Keyfiles = [ "oneliner.hlp" ]
;
;------Keyfiles = [ "oneliner.pro", "userlib.hlp" ]
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


