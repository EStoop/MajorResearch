PRO uscroll, text, title = title, xsize=xsize, ysize=ysize, $
		xpos=xpos, ypos=ypos, quit=quit, retain=retain, init=init
;
;+NAME/ONE LINE DESCRIPTION:
;    USCROLL lets the user scroll through text.
;
;DESCRIPTION:
;    For non-X environments, and for X environments without widgets,
;    this routine allows text-examination by calling Kevin Turpie's
;    TSCROLL.  For X-environments with widgets, the text-examination
;    is performed by calling John Ewing's XSCROLL.
;
;CALLING SEQUENCE:
;    uscroll, text, title=title, xsize=xsize, ysize=ysize
;
;ARGUMENTS (I=input, o=output, []=optional):
;    text       I   1-D arr, string     Text to be examined.
;    title     [I]  scalar,  string     A title for the scrolling window.
;    xsize     [I]  scalar,  int        X-size of scrolling window.
;    ysize     [I]  scalar,  int        Y-size of scrolling window.
;    xpos      [I]  scalar,  int	X-position, column (0-80)
;    ypos      [I]  scalar,  int	Y-position, row (0-24)
;    quit      [O]  scalar,  int	0=Tab, or 1=Q (retained screen)
;    retain    [I]  scalar,  int        If set, TScroll screen remains after
;                                       exit.
;
;WARNINGS:
;    None.
;
;EXAMPLE:
;    uscroll, text
;#
;COMMON BLOCKS: none
;
;PROCEDURE:
;    If !D.NAME equals 'X' and widgets are available, then call XSCROLL,
;    else call TSCROLL.
;
;LIBRARY CALLS:  None.
;
;REVISION HISTORY:
;    Written by John Ewing (ARC)   16 September 1992
;
;    Mod	Pete KS	(GSC)	   18 November 1992 to match tscroll args
;						    optional args do not
;						    change calling seq
;    Mod        K.Turpie (GSC)     21 December 1992 added RETAIN and INIT
;                                                   keywords.
;                                                   Removed default settings
;                                                   and command to clear screen
;                                                   before call to TScroll.
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;.TITLE
;Routine USCROLL
;-
  old_quiet = !QUIET
  !QUIET = 1
  IF((!D.NAME EQ 'X') OR (!D.NAME EQ 'WIN')) THEN BEGIN
    xscroll, text, title=title, xsize=xsize, ysize=ysize
  ENDIF ELSE BEGIN
    tscroll, text, title=title, xsize=xsize, ysize=ysize, quit=quit, $
	xpos=xpos, ypos=ypos, retain=retain, init=init
  ENDELSE
  !QUIET = old_quiet
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


