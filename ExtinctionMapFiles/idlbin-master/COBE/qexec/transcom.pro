 PRO TransCom, InString, OutString, OutLen, InvLen, Plain=Plain
;+
;   Description
;   -----------
;   This routine replaces special text commands with escape sequences
;   to display text with special attributes (e.g., bold, underline, etc.)
;
;#
;   Called by
;   ---------
;   TScroll
;
;   Routines Called
;   ---------------
;   InitWEX
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  March 23, 1992
;
;   10/26/92  KRT  Added {T} command
;   11/17/92  KRT  Integrated into the WEX system and documented per
;                  SPR 10216.
;   11/19/92  PKS  Added /Plain to strip commands without replacement
;             KRT 
;   11/20/92  KRT  Added call to InitWEX
;
;.Title
;   Procedure TransCom
;-
;
 On_Error, 2       ; Return to calling routine upon an error
;
@wexdisp.inc       ; Access to display control escape sequences
;
;Begin,
;
    InitWEX
;
    Marker    = [ '{BT}',   '{N}',     '{R}',    '{U}',     '{G}',     '{T}'  ]
    SpcChr    = [ BoldText, NormVideo, RevVideo, Underline, GrphicsON, TextOn ]
;
    VarStat   = Size( Marker )
    NMarkers  = VarStat(1)
;
    OutString = InString
    OutLen    = StrLen( OutString )
    RemLen    = 0
    InvLen    = 0
;
    For i = 0, NMarkers-1 do begin
;
       Mark   = StrPos( OutString, Marker(i) )
       MrkLen = StrLen( Marker(i) )
       SpcLen = StrLen( SpcChr(i) )
;
       While (Mark ne -1) do begin
;
          LeftPart  = StrMid( OutString, 0, Mark )
          RightPart = StrMid( OutString, Mark+MrkLen, OutLen )
;
	  If (KeyWord_Set( Plain )) Then Begin
             OutString = LeftPart + RightPart
	  EndIf Else Begin
             OutString = LeftPart + SpcChr(i) + RightPart
             InvLen    = InvLen + SpcLen
	  EndElse

          RemLen    = RemLen + MrkLen
          Mark      = StrPos( OutString, Marker(i) )
;
       EndWhile
;
    EndFor
;
    OutLen = OutLen - RemLen
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


