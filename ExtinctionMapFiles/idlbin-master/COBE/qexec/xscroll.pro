PRO xscroll_ev,event
;+NAME/ONE LINE DESCRIPTION:
;    XSCROLL scrolls text in a scrollable widget (X windows only).
;
;DESCRIPTION:
;    This routine puts up text in a scrollable widget, in which the
;    user can employ a scroll bar to move through the text.
;
;CALLING SEQUENCE:
;    xscroll, text, xsize=xsize, ysize=ysize, title=title
;
;ARGUMENTS (I=input, o=output, []=optional): 
;    text       I   1-D arr, str        Array of text.
;    xsize     [I]  scalar,  int        X-size of scroll widget.
;    ysize     [I]  scalar,  int        Y-size of scroll widget.
;    title     [I]           str        Title to appear on the widget.
;
;WARNINGS: none
;
;EXAMPLE:
;    xscroll, text, xsize=xsize, ysize=ysize, title=title
;#
;COMMON BLOCKS: The common block is used to pass information from the
;               main routine to the event handler only. It's form is
;               extr,outtext,nlines
;
;PROCEDURE:
;    Eliminate any special character sequences such as "{B}" from the
;    text.  Calculate an appropriate X size of the widget.  Call the
;    appropriate WIDGET* routines.
;    The event handler XSCROLL_EV must be compiled first and so is
;    the first routine in this file.  It takes as its calling sequence 
;    the event structure associated with the mouse.  Almost all of the
;    real action (reading the cursor for done and file dump, and
;    outputting the text) occurs in the event handler routine.  The 
;    XSCROLL routine itself follows it.  It takes the calling
;    sequence as described in the prologue above.
;
;LIBRARY CALLS:  None.
;
;REVISION HISTORY:
;    Written by John Ewing,  ARC,  September 1992.
;    Eliminate substrings bounded by braces.  J Ewing  Nov 20, 1992.
;  SPR 10663  Mar 09 93  XSCROLL needs an adjustment for ULTRIX.  J Ewing
;  SPR 10333  Mar 10 93  Choose fixed width font so tables are
;                        displayed correctly.                 D. Bazell
;  SPR 11127  06 Jul 93  IDL for Windows compatability. J. Newmark
;  SPR 11381  Oct 16 93  Add file dump capability. J. Newmark
;
;.TITLE
;Routine XSCROLL
;-
COMMON extr,outtext,nlines
WIDGET_CONTROL, event.id, GET_UVALUE = eventval		;find the user value
CASE eventval OF
  "Print" : BEGIN
             Print, Format = '($,"Enter output file name")'
             FileName = ""
             Read, FileName
             If (StrTrim( FileName, 2 ) ne "") then begin
                Get_LUN, LUN
                OpenW, LUN, FileName, Error=IOS
                If (IOS eq 0) then begin
                   For i = 0, nlines-1 do begin
                      TransTab, outtext(i),    TabFreeStr
                      TransCom, TabFreeStr, PlainText, Len, SpcLen, /Plain
                      PrintF, LUN, PlainText
                   EndFor
                   Print, " "
                   Print, " Text was written to file "+FileName+"."
                 EndIf Else begin
                   Print, " "
                   Print, " File " + FileName + " cannot be opened."
                   Print, " Sorry, but no text was written."
                 EndElse
                 Close, LUN
                 Free_LUN, LUN
             EndIf
         END

  "Done": WIDGET_CONTROL, event.top, /DESTROY		
ENDCASE
END 

PRO xscroll, text, xsize=xsize, ysize=ysize, title=title , group=group
;
;  This is the main routine as described in the prologue to the event
;  handler, above.  The calling sequence, etc., is described in detail 
;  there.
;
COMMON extr,outtext,nlines
  IF((!D.NAME NE 'X') AND (!D.NAME NE 'WIN')) THEN BEGIN
    MESSAGE, 'This routine is only usable on X windows.', /CONT
    RETURN
  ENDIF
  IF(NOT KEYWORD_SET(ysize)) THEN ysize = 40
  IF(NOT KEYWORD_SET(title)) THEN title = 'Scrollable text'
  nlines = N_ELEMENTS(text)
  IF(ysize GT nlines) THEN ysize = nlines
;
;  Eliminate special character sequences bounded by braces.
;  --------------------------------------------------------
  FOR i=0,nlines-1 DO BEGIN
    ic = STRPOS(text(i), '{')
    WHILE(ic NE -1) DO BEGIN
      ic2 = STRPOS(text(i), '}')
      text(i) = STRMID(text(i),0,ic) + $
                STRMID(text(i),ic2+1,STRLEN(text(i)))
      ic = STRPOS(text(i), '{')
    ENDWHILE
  ENDFOR
  outtext=text
;
;  Calculate an appropriate X-size.
;  --------------------------------
  maxlen = 0
  FOR i=0,nlines-1 DO $
    IF(STRLEN(text(i)) GT maxlen) THEN maxlen = STRLEN(text(i))
  IF(NOT KEYWORD_SET(xsize)) THEN xsize = maxlen
;
;  Call the appropriate WIDGET* routines.
;  --------------------------------------
  base = WIDGET_BASE(TITLE = title)
  w1 = WIDGET_BUTTON(base, VALUE = "Done",uvalue='Done')
  w2 = WIDGET_BUTTON(base, xoffset=50,VALUE = "File Dump",uvalue='Print')
  w3 = WIDGET_TEXT(base, XSIZE = xsize, YSIZE = ysize, /SCROLL, VALUE = text,$
    YOFFSET = 26)
  WIDGET_CONTROL, base, /REALIZE
  XManager, "xscroll", base, EVENT_HANDLER = "xscroll_ev", GROUP_LEADER = GROUP
END

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


