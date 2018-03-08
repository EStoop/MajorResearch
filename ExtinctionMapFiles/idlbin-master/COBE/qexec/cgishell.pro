 PRO CGIShell
;+
;   Description
;   -----------
;   CGIShell uses the IDL Spawn procedure to initiate a subprocess from
;   a WEX user interface.
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
;   Created by K. Turpie,  General Sciences Corporation,  March 1992
;   SPR 11003  Jun 02 93   Change !version.os to !cgis_os.  J Newmark
;
;.Title
;   Procedure CGIShell
;-
 On_Error, 2 ; return to calling procedure if an error occurs
;
@wexdisp.inc
;
;Begin
    Print, " "
    Case StrLowcase( !cgis_os ) of
       'vms' : Print, " Type "+BoldText+"LOGOUT"+NormVideo+" to return to CGIS"
        else : Print, " Type "+BoldText+"EXIT"  +NormVideo+" to return to CGIS"
    EndCase
    Print, " "
;
    Spawn
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


