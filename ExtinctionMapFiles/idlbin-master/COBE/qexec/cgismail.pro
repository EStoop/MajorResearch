 PRO CGISMail
;+
;   Description
;   -----------
;   CGISMail uses the IDL Spawn procedure to initiate a mail facility.
;#
;   Called by
;   ---------
;   CGIS
;
;   Routines Called
;   ---------------
;   None
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  May 1992
;   SPR 11003  Jun 02 93   Change !version.os to !cgis_os.  J Newmark
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;
;.Title
;   Procedure CGISMail
;-
 On_Error, 2    ; Return to calling routine if an error occurs.
@wexdisp.inc    ; Access to WEX global variables for display device control.
@wexpath.inc    ; Access to WEX file path and utility names.
;
;Begin,
    InitWEX  ; Initialize WEX if it isn't.
;
    If ((StrLowCase( !D.Name ) ne 'x') and $
        (StrLowCase( !D.Name ) ne 'win')) then Print, EraseScrn+CursrHome
;
    Case StrLowCase( !cgis_os ) of
       'vms' : Spawn, !VMSMail
       'windows' : begin
                  Print, CursrHome, EraseScrn
                  Print, 'No mail available on PCs at this time.'
                  Print, 'Press any key to continue...'
                  Ans=Get_Kbrd( 1 )
               end
        'unix' : Spawn, !UNIXMail, /NoShell
    EndCase
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


