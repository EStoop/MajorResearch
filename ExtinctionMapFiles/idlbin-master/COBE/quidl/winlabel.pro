	PRO winlabel_event, ev
;
;+NAME/ONE-LINE DESCRIPTION:
;    WINLABEL uses a GUI to position labels in a window with the cursor.
;
;DESCRIPTION
;    This routine uses a GUI to allow the user to specify the text,
;    size, thickness, orientation, and justification of a string.  It 
;    also allows the user to see the string in a test window before
;    writing it into the window of interest, and to erase a string 
;    that has just been placed into the target window.  The string is
;    positioned in the window using the cursor. As many strings as
;    desired can be placed before exiting the routine.
;
;    The user must specify the desired target window number on the
;    command line.  The GUI that appears consists of the following:
;
;    Text box:  User clicks on here, then enters text.  TEXT MUST END
;               WITH A CARRIAGE RETURN!
;
;    Sliders:   Select the text size and thickness in units of
;               0.1 x IDL's default (i.e. default=10)
;
;               Select the orientation in degrees CCW from horizontal
;
;               Select the color, where the default 255 is brightest
;
;    Menu:      Select text justification with respect to cursor
;               position.  Default = 'Center'
;
;    Buttons:   Test:  Causes a horizontal, centered version of the
;                      string to appear in a test window.  If you
;                      don't like it, change the settings or text and
;                      click 'Test' again.
;
;               OK:    Will activate the cursor.  Move it to the desired    
;                      location in the window of interest, then click
;                      the left mouse button.  The text will appear with
;                      the appropriate orientation, aligned to the left,
;                      center, or right of the cursor as selected.
;
;               Erase: Will overwrite the most recently-placed string with
;                      the background color.
;
;               Done:  Exits WINLABEL, causing the widget box and test 
;                      window to disappear     
;                  
;CALLING SEQUENCE:
;    WINLABEL, window_number
;
;ARGUMENTS (I = input):
;    window_number      I       integer      Target window number
;
;WARNINGS:
;    1.  A valid window number must be supplied on the command line.
;    2.  The text string MUST end in a <CR>, or no text will appear. 
;    3.  Buffering idiosyncrasies when running under MS Windows on PCs
;          cause text output and erasing not to be completed until the
;          next mouse event happens.  Don't panic:  all buffers are
;          flushed and output completed when you leave the routine.
;
;EXAMPLE:  WINLABEL, 3
;#
;COMMON BLOCKS:
;    A common block is required in order for the main routine to 
;    preserve the text and other attributes during repeated calls to the
;    event handler.  The block is as follows:
;
;    COMMON TEXTINFO, size, color, thick, orient, $ ; font attributes
;                     align, xpos, ypos, $          ; cursor align and posn
;                     outstr,            $          ; output string
;                     testwin, firsttest,$          ; text window and flag 
;                     markmap            $          ; =1 if OK has happened
;
;PROCEDURE AND PROGRAMMING NOTES:
;    The event handler WINLABEL_EVENT must be compiled first and so is
;    the first routine in this file.  It takes as its calling sequence 
;    the event structure associated with the mouse.  Almost all of the
;    real action (extracting the attirbutes, reading the cursor, and
;    outputting the text) occurs in the event handler routine.  The 
;    WINLABEL routine itself follows it.  It takes the calling
;    sequence as described in the prologue above but otherwise does
;    nothing but establish the common block and set up the GUI widgets.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS:
;    STRNUMBER is called in the main routine to verify that the input 
;    parameter is a valid number.
;
;MODIFICATION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  26 Jan 1993
;    Modified by RBI to use STRNUMBER to test valid input  28 Jan 93
;    SPR 11980  8-Nov8-Nov-1994   IDL v3.5 changes   J. Newmark
;-
	COMMON textinfo, size, color, thick, orient, align, xpos, ypos, $
			 outstr, testwin, firsttest, markmap
	type = TAG_NAMES(ev, /STRUCTURE)
;
;  The <CR> at the end of the text string triggers the event handler to
;  get the string.
;
	IF (strmid(type,0,11) EQ 'WIDGET_TEXT') THEN WIDGET_CONTROL, ev.id, $
						GET_VALUE=outstr
;
;  Determine whether it was the size, orientation, thickness, or color
;  slider.
;
	IF (type EQ 'WIDGET_SLIDER') THEN BEGIN 
	    WIDGET_CONTROL, ev.id, GET_VALUE=value, GET_UVALUE=which
	    CASE which OF 's':  size = value
			  'c':  color= value
			  't':  thick= value
			  'o':  orient=value
	    ENDCASE
	ENDIF
;
;  Both the text alignment buttons and Test/Done/etc buttons send WIDGET_
;  BUTTON events, so we can treat them all in the same CASE statement.
;
	IF (type EQ 'WIDGET_BUTTON') THEN BEGIN
	    WIDGET_CONTROL, ev.id, GET_VALUE=value
	    CASE value OF 
		  'Left': align = 0.
		'Center': align = 0.5
		 'Right': align = 1.
 ;
 ;        The Test button looks to see whether a test window already exists.
 ;        If so, it deletes it, then establishes a new one with the same
 ;        label to display the test string.  The test string ignores
 ;        the orientation and alignment in the XYOUTS command, centering 
 ;        the string in the window instead.
 ;
		  'Test': BEGIN
			  IF (firsttest EQ 0) THEN WDELETE, testwin
			  WINDOW, /FREE, XSIZE=500, YSIZE=100, $
				  TITLE='Test Map Label (Horizontal)'
			  XYOUTS, 250, 50, outstr, COLOR=color, $
					CHARSIZE=size/10., $
					CHARTHICK=thick/10., $
					ALIGN=0.5, /DEVICE
			  testwin = !D.WINDOW
			  firsttest = 0             ; a test window exists
			  END
;
;         The OK button puts the string in the specified window with all
;         the appropriate slider settings, and at the position given by  
;         the cursor click.  The "markmap" flag is then
;         set so that the string can be erased if necessary.
;
		   'OK':  BEGIN
			  WIDGET_CONTROL, ev.top, GET_UVALUE=winnum
			  WSET, winnum
			  CURSOR,xpos,ypos,/DEVICE,/WAIT
			  XYOUTS, xpos, ypos, outstr, COLOR=color, $
				CHARSIZE=size/10., ORIENT=orient, $
				CHARTHICK=thick/10., $ 
				ALIGN=align, /DEVICE
			  markmap = 1
			  END
;
;         The erase operation only makes sense if a string has been
;         placed with an "OK", i.e. the "markmap" flag has been set.  
;         Erasing is done by overwriting in the background color.
;
		 'Erase':  BEGIN
			  IF markmap THEN XYOUTS, xpos, ypos, outstr, $
				COLOR=0, ALIGN=align, $
				CHARSIZE=size/10., ORIENT=orient, $
				CHARTHICK=thick/10., /DEVICE
			  END
;
;          If done, get rid of the widget and delete the test window if
;          the latter exists.
;
		 'Done':  BEGIN
			  WIDGET_CONTROL, /DESTROY, ev.top
			  IF (firsttest EQ 0) THEN WDELETE, testwin
			  END
	   ENDCASE
	ENDIF
	END


	PRO winlabel, winnum
;  This is the main routine as described in the prologue to the event
;  handler, above.  The calling sequence, etc., is described in detail 
;  there.
;
	COMMON textinfo, size, color, thick, orient, align, xpos, ypos, $
			 outstr, testwin, firsttest, markmap

	ON_ERROR, 2
	IF (N_PARAMS() LT 1) THEN MESSAGE, $
	    'You must specify a window number, e.g. WINLABEL, 3'
	IF (STRNUMBER(winnum) NE 1) THEN MESSAGE, $
            'Invalid window ID!  You must specify a numeric value.'
;
;  Initialize all of the flags and string characteristics.  This is           
;  necessary because if a slider is left in its initial position then
;  no event is sent and the value is not picked up.
;
	size = 10
	thick = 10
	color = 255
	orient = 0
	align = 0.5
	outstr = ' '
	firsttest = 1
	markmap = 0
;
;  Build the widget box with sliders, butons, etc.
;
	base = WIDGET_BASE (TITLE='Map Label', /COLUMN, UVALUE=winnum)
	
	textbox = WIDGET_BASE (base, /FRAME, /ROW)
	textspace = WIDGET_LABEL (textbox, VALUE='Enter text then <cr>: ')
	textspace = WIDGET_TEXT (textbox, /EDITABLE, XSIZE=45, YSIZE=1)

	newbase = WIDGET_BASE (base, /ROW, SPACE=70)
	lcol = WIDGET_BASE (newbase, /COLUMN,xsize=330)
	mcol = WIDGET_BASE (newbase, /COLUMN)
	rcol = WIDGET_BASE (newbase, /COLUMN)

	sizeslide = WIDGET_SLIDER (lcol, TITLE='Font Size (units of 0.1)',$
			 MIN=1, VALUE=10, UVALUE='s',ysize=75)

	thickslide = WIDGET_SLIDER (lcol, TITLE='Font Thickness (units of 0.1)',$
			 MIN=1, VALUE=10, UVALUE='t',ysize=75)

	orientslide = WIDGET_SLIDER (lcol, TITLE='Orientation (degrees)', $
			 MIN=0, MAX=359, VALUE=0, UVALUE='o',ysize=75)

	colorslide = WIDGET_SLIDER (lcol, TITLE='Label Color',$
			 MIN=0, MAX=255, VALUE=255, UVALUE='c',ysize=75)

	align_lbl = WIDGET_LABEL (mcol, VALUE='Text Align')
	XMENU, ['Left', 'Center', 'Right'], mcol, /COLUMN, $
	       /EXCLUSIVE, /FRAME

	donebox = WIDGET_BASE (rcol, /FRAME, /COLUMN)
	test = WIDGET_BUTTON (donebox, VALUE="Test")
	ok = WIDGET_BUTTON (donebox, VALUE="OK")
	erase = WIDGET_BUTTON (donebox, VALUE="Erase")
	done = WIDGET_BUTTON (donebox, VALUE="Done")
 ;
 ;  Bring the widget box onto the screen and pass control to the event
 ;  handler.
 ;
	WIDGET_CONTROL, base, /REALIZE
	XMANAGER,'winlabel',base
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


