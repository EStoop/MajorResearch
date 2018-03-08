 PRO CGIS
;+
;   Description
;   -----------
;   This is the CGIS main menu driver.  It initializes the WEX system
;   for the entire CGIS executive user interface and invokes the
;   main menu with the procedure Action.
;
;#
;   Called by
;   ---------
;   MAIN (UIDL)
;   CGISSTART (for startup invocation)
;
;   Routines Called
;   ---------------
;   Action
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  August 1992
;
;   11/20/92  KRT  Removed call to InitWEX (SPR 10216)
;   12/21/92  KRT  Removed include file WEXCntrl.Inc and
;                  add InitWEX with keyword /Again to force
;                  reinitialization.
;
;.Title
;   Procedure CGIS
;-
@wexcntrl.inc
;
 On_Error, 2   ; Return to calling program on error.
;
    SaveQuiet   = !Quiet
    !Quiet      = 1
;
    Quit = 0
    While (not Quit) do begin
       InitWEX, /Again
       Action, "CGISMAIN.MNU"
    EndWhile
;
    !Quiet      = SaveQuiet
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


