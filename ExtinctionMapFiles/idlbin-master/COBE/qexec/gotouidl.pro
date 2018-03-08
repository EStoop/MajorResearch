 PRO GotoUIDL, Quit
;+
;   Description
;   -----------
;   GotoUIDL simply writes instruction to the user as to how to
;   return to the CGIS menu system as it exits to UIDL.  This
;   procedure is invoked from the CGIS main menu.  The actual
;   exit is performed by the WEX system since the global WEX
;   variable Quit is set to 1.
;
;#
;   Called by
;   ---------
;   CGISMAIN.MNU (CGIS main menu parameter file)
;
;   Routines Called
;   ---------------
;   InitWEX
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  August 1992
;
;   11/23/92  KRT  Documented and integrated into the WEX system (SPR 10216).
;
;.Title
;   Procedure GotoUIDL
;-
 On_Error,2     ; Return to caller if an error occurs
;
@wexdisp.inc    ; Access to WEX global variable for display device control
;
;Begin,
    InitWEX
;
    Print, EraseScrn+CursrHome
    Print, " Type "+BoldText+"EXIT"+NormVideo+" to return to VMS or"
    Print, " type "+BoldText+"CGIS"+NormVideo+" to return to the main menu..."
    Print, " "
    Quit = 1
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


