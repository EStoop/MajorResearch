 Pro TextGrph
;+
;   Description
;   -----------
;   TextGrph set-ups display control commands for common access
;
;#
;   Called by
;   ---------
;   InitWEX
;
;   Routines Called
;   ---------------
;   <none>
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  March 23, 1992
;
;   10/27/92  KRT  Add Cntl-B per SPR 10168.
;   11/17/92  KRT  Documented and integrated into the WEX system (SPR 10216).
;   12/22/92  KRT  Put OS case statement to set CR properly.
;   SPR 11003  Jun 02 93  Changed !version.os to !cgis_os. J Newmark
;   SPR 11127 Oct 18 93 IDL for Windows Compatability. J. Newmark
;
;.Title
;   Procedure TextGrph
;-
 On_Error, 2   ; Return to the calling procedure if an error occurs.
;
@wexdisp.inc   ; Access to screen display command sequences
;
;Begin,
;   Define special characters :
;
    CntlB = String(   2b )
    CntlF = String(   6b )
    CntlR = String(  19b )
    CntlW = String(  23b )
;
    Case StrLowcase( !cgis_os ) of
       "windows"  :  CR = String( 13b )
       "vms"  :  CR = String( 13b )
       "unix"  :  CR = String( 10b )
    EndCase
;
    Tab   = String(   9b )
    Bell  = String(   7b )
    ESC   = String(  27b )
    CSI   = String( 155b )
    SEQ   = ESC+'['
;
;   Define escape sequences for highlighting, bolding, line drawing, etc. :
;
    BoldText    = SEQ+'1m'   ; Display bold characters
    Underline   = SEQ+'4m'   ; Display underlined characters
    RevVideo    = SEQ+'7m'   ; Display reverse video text
    NormVideo   = SEQ+'m'    ; Set all text attributes to normal
;
    EraseToEoL  = SEQ+'?K'   ; Erase line from cursor to end of line
    EraseScrn   = SEQ+'2J'   ; Completely erase screen
;
    GrphicsON   = ESC+'(0'   ; Graphics text on
    TextON      = ESC+'(B'   ; Norm text graphics on
    TxtCursrON  = SEQ+'?25h' ; Enable text cursor
    TxtCursrOFF = SEQ+'?25l' ; Turn text cursor off
;
;   Define cursor movements commands :
;
    CursrUp     = SEQ+'A'    ; Move cursor up
    CursrDwn    = SEQ+'B'    ; Move cursor down
    CursrLft    = SEQ+'C'    ; Move cursor left
    CursrRgt    = SEQ+'D'    ; Move cursor right
;
    CursrHome   = SEQ+'H'    ; Move cursor to home
    Return
;
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


