 PRO TErase, XSize, YSize, XPos=XPos, YPos=YPos
;+
;   Purpose
;   -------
;   TErase will erase a rectangular area of text on the screen
;   that is XSize by YSize characters.  The position of the upper
;   left-hand corner of the rectangle is given by the global WEX
;   variables CurrCol and CurrRow found in the include file
;   WEXCntrl.Inc.
;
;   TErase, XSize, YSize, XPos=XPos, YPos=YPos
;
;   where,  XSize = the number of column to be erased (in characters)
;           YSize = the number of rows to be erased (in lines)
;
;           XPos  = the column location of the area's left side; this
;                   defaults to the global variable CurrCol
;           YPos  = the row location of the area's top side; this
;                   defaults to the global variable CurrRow
;
;   Note: TErase only works on terminals in VT220 mode or higher
;#
;   Called by
;   ---------
;   TScroll
;
;   Routines Called
;   ---------------
;   InitWEX     MovToCol     MovToRow
;
;   History
;   -------
;   Created by K. Turpie,  General Sciences Corporation,  November 1992
;
;   01/14/93  KRT  Cleaned-up algorithm per SPR 10455.  Added keywords
;                  XPos and YPos.  Improved the documentation.
;
;.Title
;   Procedure TErase
;-
;
 On_Error,2    ; Return to caller if an error occurs
;
@wexdisp.inc   ; Access to screen display command sequences
@wexcntrl.inc  ; Access to window attributes
;
;Begin,
;
    InitWEX    ; Initialize WEX if necessary
;
    If (KeyWord_Set( YPos )) then Row = YPos $
                             else Row = CurrRow
;
    If (KeyWord_Set( XPos )) then Col = XPos $
                             else Col = CurrCol
;
;   Define top and bottom of window :
    Top    = Row
    Bottom = Row+YSize
;
;   Erase window bottom-up and return :
    EraseNum = String( XSize, Format='(I2.2)' )
    LeftMarg = MovToCol( Col )
;
    For i= Bottom, Top, -1 do Print, MovToRow( i )  +       $
                                     LeftMarg       +       $
                                     SEQ + EraseNum + 'X',  $
                                     Format='(A,$)'
;
    Print, CursrHome, Format="(A,$)"
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


