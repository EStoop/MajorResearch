  PRO DocDisp, DocFile, DocName
;+
;   Description
;   -----------
;   DocDisp displays the document file specified in the input variable
;   DocFile using the WEX system routine Action.  A message is issued
;   to the user that the document is being retrieved and refers to
;   the document's name or title (which is specified by the input var-
;   iable DocName).
;
;   Calling Sequence:
;
;        DocDisp, DocFile, DocName
;
;   where,
;
;        DocFile - is an input string variable containing the document's
;                  file name.
;
;        DocName - is an input string variable containing the document's
;                  title.
;#
;   Called by
;   ---------
;   Action
;
;   Routines Called
;   ---------------
;   Action
;
;   History
;   -------
;   Created by C. Groden,  Universities Space Research Association,  Dec 1992
;
;.Title
;   Procedure DocDisp
;-
;
  On_Error, 2    ; Return to called on error
;
; Begin
;
  Print, " Please wait while "+DocName+" is being retrieved..."
  Action, DocFile
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


