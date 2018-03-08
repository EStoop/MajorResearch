 Function MovToCol, Col
;+
;   Description
;   -----------
;   MovToCol returns the escape sequence to move the cursor to column Col.
;   This movement is relative to the current cursor position.  If used
;   before the result of MovToRow, the MovToRow command will move the
;   cursor back to the first column.
;#
;   Called by
;   ---------
;   DispText
;   TErase
;
;   Routines Called
;   ---------------
;   <none>
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  February 1992
;
;   11/23/92  KRT  Documentation extended per SPR 10216.
;   01/25/93  KRT  Reduced to a one-liner.
;                  Used StrN instead of String.
;
;.Title
;   Function MovToCol
;-
 On_Error, 2 ; Return to calling procedure if there an error.
;
;Begin,
    If (Col-1 le 0) then Return, String( 27b ) + '[' + 'D' $
                    else Return, String( 27b ) + '[' + StrN( (Col-1) ) + 'C'
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


