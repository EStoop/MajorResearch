 Pro GotoOS
;+
;   Description
;   -----------
;   GotoOS is used to verify with the user that they wish exit the
;   CGIS (and UIDL) system and go to the operating system.  They are
;   warned that any data (in memory) they did not save would be lost.
;
;#
;   Called by
;   ---------
;   CGISMAIN.MNU (Parameter file for CGIS Main Menu)
;   
;   Routines Called
;   ---------------
;   InitWEX       DispText
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  August 1992
;
;   11/23/92  KRT  Documented and integrated into the WEX system (SPR 10216).
;
;.Title
;   Procedure GotoOS
;-
 On_Error,2     ; Return to caller if an error occurs
;
@wexdisp.inc    ; Access to WEX global variables for display device control
@wexcntrl.inc   ; Access to WEX stack, recursion control, and WEX components
;
;Begin
;
    InitWEX  ; Initialize WEX if it isn't.
;
    Print, EraseScrn+CursrHome
;
;        1234567890123456789012345678901234567890123456789
    a = "                                                 "
    b = " All data will be removed from memory upon exit. "
    c = "                                                 "
    d = " Be sure you have written to disk any images or  "
    e = " arrays you have created and wish to save.       "
    f = "                                                 "
;
    CurrRow    =  3
    CurrCol    = 15
;
    Width      = 49
    Top        = 0
    Bottom     = 5
;
    WEXTitle   = "**** CGIS EXIT ****"
    WEXMessage = "**** CGIS EXIT ****"
    Text       = [ a, b, c, d, e, f ]
;
    DispText, Text, Width, Top, Bottom
    Print, TxtCursrON  ; Enable text cursor
;
    Print, " "
    Print, " "
    Print, " "
    Print, " Exit CGIS? (Y/N)", Format = '($,15X,A)'
;
    Answer = ''
    Read, Answer
;
    If ((Answer eq 'Y') or (Answer eq 'y')) then Exit
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


