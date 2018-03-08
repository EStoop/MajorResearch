 PRO TransTab, InString, OutString
;+
;   Description
;   -----------
;   This routine replaces tabs found in the input text array, InString,
;   with a fixed number of blanks.
;   (current the number has been "hard-wired" to 6)
;
;#
;   Called by
;   ---------
;   TScroll
;
;   Routines Called
;   ---------------
;   <none>
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  March 23, 1992
;
;   11/17/92  KRT  Removed include file (it didn't need it).  Documentation
;                  extended per SPR 10216.
;   02/22/93  KRT  Rewrote prolog description (SPR 10521)
;
;.Title
;   Procedure TransTab
;-
 On_Error, 2       ; Return to calling routine upon an error
;
;Begin,
;
    Tab       = String( 09b )
    Spaces    = "      "
;
    OutLen    = StrLen( InString )
    OutString = InString
;
    Len       = StrLen( Spaces )
    Mark      = StrPos( OutString, Tab )
;
    While (Mark ne -1) do begin
;
       LeftPart  = StrMid( OutString, 0, Mark )
       RightPart = StrMid( OutString, Mark+1, OutLen )
;
       OutLen    = OutLen + Len - 1
       OutString = LeftPart + Spaces + RightPart
;
       Mark      = StrPos( OutString, Tab )
;
    EndWhile
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


