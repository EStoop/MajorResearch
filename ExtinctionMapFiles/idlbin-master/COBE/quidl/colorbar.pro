        pro colorbar, winnum, label=label, min=min, max=max, rev=rev
;+NAME/ONE-LINE DESCRIPTION:
;     COLORBAR puts a color bar of any size, anywhere with the mouse
;
;DESCRIPTION:
;     The user uses mouse clicks to specify the size and location
;     of a colorbar in the specified window: three mouse clicks with
;     the left button specify the bar's center, left edge, and top.  
;     Through the use of keywords, the user can also:  specify the 
;     range of colors (default [0,255]); reverse the bar so that 
;     the bright edge is on the left or at the bottom; and activate 
;     the GUI labeling routine to place labels around the bar.
;
;CALLING SEQUENCE:
;     COLORBAR, WINNUM, [/LABEL], [/REVERSE], [MIN=<minimum>], [MAX=<max>]
;
;ARGUMENTS:
;     WINNUM       byte           Window number.
;     [/LABEL]     switch         Activates WINLABEL routine, which
;                                   uses a GUI to allow the user to
;                                   place whatever text (s)he wants
;                                   on or around the bar. 
;                                   Default = no label.
;     [/REVERSE]   switch         Causes dark part to be at the
;                                   right (or at the top for
;                                   vertically-oriented bars)
;                                   Default = dark left or bottom.
;     [MIN]        byte           Value between 0 and 254 to 
;                                   specify darkest color on bar.
;                                   Default = 0.
;     [MAX]        byte           Value between 1 and 255 to
;                                   specify lightest color on bar.
;                                   Default = 255. 
;
;WARNINGS:
;     1.  MIN must be less than MAX; both must be in [0,255].
;     2.  Routine uses the current active window.
;
;EXAMPLES:
;     The simplest effective invocation is just COLORBAR.  For a
;     rather grey-looking bar that uses only the greyscale range
;     [100,200] say COLORBAR,MIN=100,MAX=200.
;#
;COMMON BLOCKS:  none
;
;PROCEDURE AND OTHER PROGRAMMING NOTES:  nothing special
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    WINLABEL is called when the /LABEL switch is used.
;  
;MODIFICATION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  14 April 1993  SPR 10501
;    SPR 10907  May 07 93  Add new call parameter, WINNUM. J Newmark
;
;.title
;Routine COLORBAR
;-
;   Set default scaling range or quit if specified range is bad.
;
	ON_ERROR,2
	IF (N_PARAMS() LT 1) THEN MESSAGE, $
	    'You must specify a window number, e.g. COLORBAR, 3'
	IF (STRNUMBER(winnum) NE 1) THEN MESSAGE, $
            'Invalid window ID!  You must specify a numeric value.'
	IF NOT KEYWORD_SET(max) THEN max = 255
	IF NOT KEYWORD_SET(min) THEN min = 0
        IF (min LT 0) OR (max GT 255) OR (min GE max) THEN $
	     MESSAGE, 'Illegal min and max range; must be [0,255]'
;
;  Locate the center and edges of the bar
;
        WSET, winnum
        PRINT, 'Click cursor at center of color bar.'
        CURSOR,xctr,yctr,/DEVICE,/DOWN	
        PRINT, 'Click cursor at left edge of color bar.'
        CURSOR,xl,yl,/DEVICE,/DOWN
        PRINT, 'Click cursor at top of color bar.'
        CURSOR,xt,yt,/DEVICE,/DOWN
;
;   Figure out the dimensions, and, if it's higher than it is wide,
;   create the array "sideways".  The array itself is filled with
;   the numbers 0 - 255 by BYTSCLing the index number of every point.
;
        width = 2 * ABS(xctr - xl)
        height = 2 * ABS(yt - yctr)
	wedge = BYTSCL (FINDGEN(width<height, width>height))
;
;   Offset and scale to the min and max part of the color range.
;
        wedge = wedge * (max - min)/255. + min
;
;   Flip the wedge 90 degrees if necessary, and invert the grey scale
;   if the /REVERSE switch was set.
;
	IF (width GE height) THEN wedge = TRANSPOSE(wedge)
	IF KEYWORD_SET(reverse) THEN wedge = 255 - wedge
;
;   Put it on the screen and call the labelling routine if desired.
;
	TV,wedge,xl,yctr-height/2.
	IF KEYWORD_SET(label) THEN winlabel,winnum
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


