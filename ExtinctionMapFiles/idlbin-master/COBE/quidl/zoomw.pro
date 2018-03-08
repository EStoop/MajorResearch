;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     ZOOMW zooms window and alters color tables on windowing devices
;
;DESCRIPTION:
;     Displays part of an image (or graphics) from the current window
;     expanded in another window. This window is retained for future 
;     use after the procedure finishes. The cursor is used to mark the
;     center of the zoom, and make selections from menus to alter the
;     color table. Instructions are provided in the session text 
;     window.
;
;CALLING SEQUENCE:
;     ZOOMW [,cur_w] [,FACT=fact] $ 
;               [,XSIZE=xsize] [,YSIZE=ysize] [,/CONTINUOUS]
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     cur_w        [I]  int        window to zoom
;     fact         [I]  flt        zoom expansion factor, default = 4
;     xsize        [I]  int        X size of new window, default = 500
;     ysize        [I]  int        Y size of new window, default = 500
;     /continuous  [I]  key        If specified, the zoom window 
;                                  tracks the mouse.
;
;WARNINGS:
;     1.  The display window to be zoomed should exist.
;     2.  Only works with color windowing systems.
;     3.  A window is created and retained after procedure finishes.
;     4.  If the option to save color tables as disk files is used,
;         then the directory to which the files are written must have 
;         appropriate write protection and sufficient free space 
;         (>5 blocks/file should do).
;
;EXAMPLE:
;     After displaying an image in a window using TV or TVSCL, suppose
;     you want to get a better look at fine detail in the image. Type:
;     
;       zoomw
;
;     A crosshair cursor will appear in the window. Center the cursor 
;     at the region of interest, and click the left mouse button (MB1)
;     to create an enlarged view in a new window. Move the cursor and
;     click again to zoom a different region. Clicking the center 
;     mouse button (MB2) will present a menu which allows you to reset 
;     the size of the zoom window and the magnification factor. This 
;     menu also leads to another menu providing access to LOADCT, 
;     ADJCT, and other procedures for altering, saving, and restoring 
;     color tables. Clicking the right mouse button (MB3) will exit 
;     from the procedure.
;     
;     The more complete command:
;
;       zoomw,2,fact=5,xsize=200,ysize=250,/continuous
;
;     first switches the active display window to window 2 (i.e. it 
;     performs WSET,2), and then provides a new 200x250 pixel window
;     in which a portion of window 2 is displayed at 5x magnification.
;     The size and magnification can easily be altered (click the 
;     middle mouse button, MB2), so it is seldom necessary to set
;     these in the command line. The /continuous keyword causes the 
;     magnified region to track the movements of the cursor without 
;     repeated clicks of the left mouse button (MB1). This only works 
;     well on fast computers which can recalculate the zoom window 
;     display as fast as the cursor moves. The COBECL computers do not 
;     normal qualify as fast.
;#
;COMMON BLOCKS:
;     zoomw1, zoomw2, zoomw3, ascii, box, box_save, colors, 
;     color_edit_zoomw, hist_ct
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     ZOOMW.PRO contains the procedure zoomw and 15 other procedures
;     and functions which were written to support zoomw, but could
;     be used for other procedures as well. The first four subroutines 
;     deal with drawing and erasing boxes on the TV screen. The next
;     pair deal with finding files and directories. Next are two for
;     processing text strings, and one for printing text to a window.
;     Five more subroutines provide some of the various aspects of
;     color table manipulation which are not found in the standard 
;     IDL user library. The last subroutine is used for read in 
;     changes in the zoom factor.
;
;PERTINANT ALGORITHMS, LIBRARY CALLS, ETC.:
;     adjct, loadct, palette, tvrdc
;
;MODIFICATION HISTORY:
;     Frank Varosi NASA/GSFC 1989: This version evolved from zoom.pro 
;     in IDL library with modifications to keep window, added menu 
;     items, show zoom box, etc. 
;     Rick Arendt ARC Nov 1991: Made minor modifications to run under
;     IDL v.2 on the COBECL 
;     Rick Arendt ARC Feb 1992: combined all Varosi subroutines into 
;     zoomw.pro.
;
; SPR 9616
;  SPR 11127 06 Jul 1993  IDL for Windows compatability. J. Newmark
;  SPR 11757 17 May 1994  Update for DEC Alpha. J. Newmark
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;.TITLE
;Routine ZOOMW
;-
;	SUBPROCEDURE # 1
;
function box_create, xLo,yLo, xHi,yHi
;+
; NAME:	
;	BOX_CREATE
; PURPOSE:
;	Interactive function to draw a box on the display.
; CATEGORY:
;	Display.
; CALLING SEQUENCE: 
;	result = BOX_CREATE(xLo, yLo, xHi, yHi)
; INPUTS:
;	None.
; OUTPUTS:
;	result = 2 if center button is pressed, 4 if right button is pressed
;	xLo, yLo, xHi, yHi are set to the corners of the final box
; COMMON BLOCKS:
;	box
; SIDE EFFECTS:
; RESTRICTIONS:
;	For use with ZOOMW.
;	Intended for TV devices.
; PROCEDURE:
;	Straightforward.
; MODIFICATION HISTORY:	
;Frank Varosi NASA/GSFC 1989
;-

  common box, Lox,Loy,Hix,Hiy, sav_Horiz_B, sav_Horiz_T, sav_Vert_L, sav_Vert_R

	cursor,/DEV, x1,y1

	if (!err GE 2) then return, -!err

LOOP:
	cursor,/DEV, x2,y2, 2

	x2 = x2+1
	y2 = y2+1	; so that box is not obscured by cursor.

	if (x2 LE 0) OR (y2 LE 0) OR	$
	   (x2 GE !D.x_vsize) OR (y2 GE !D.y_vsize) then begin
		wait,.3
		goto, LOOP
	   endif

	box_erase

	box_draw, x1,y1,x2,y2

	CASE !err OF

	2: BEGIN
		xLo = Lox
		yLo = Loy
		xHi = Hix
		yHi = Hiy
		return, 2
	     END

	4: BEGIN
		box_erase
		xLo = Lox
		yLo = Loy
		xHi = Hix
		yHi = Hiy
		return, 4
	     END

	else:	goto, LOOP

	ENDCASE	
end
;
;	SUBPROCEDURE # 2
;
pro box_draw, x1,y1, x2,y2, color
;+
; NAME:	
;	BOX_DRAW
; PURPOSE:
;	Draws a box on the display.
; CATEGORY:
;	Display.
; CALLING SEQUENCE:
;	box_draw, x1, y1, x2, y2 [, color]
; INPUTS:
;	x1,y1 = lower left corner of the box
;	x2,y2 = upper right corner of the box
;	optional:
;	color = color used to draw box
; OUTPUTS:
;	none
; COMMON BLOCKS:
;	box
; SIDE EFFECTS:
; RESTRICTIONS:
;	For use with ZOOMW.
;	Intended for TV devices.
; PROCEDURE:
;	Straightforward.
; MODIFICATION HISTORY:	
;Frank Varosi NASA/GSFC 1989
;-

  common box, Lox,Loy,Hix,Hiy, sav_Horiz_B, sav_Horiz_T, sav_Vert_L, sav_Vert_R

	Lox = (x1 < x2)
	Loy = (y1 < y2)

	Hix = (x1 > x2)
	Hiy = (y1 > y2)

	xsiz = Hix - Lox +1
	ysiz = Hiy - Loy +1

	sav_Horiz_B = tvrd( Lox,Loy, xsiz,1 )
	sav_Horiz_T = tvrd( Lox,Hiy, xsiz,1 )

	sav_Vert_L = tvrd( Lox,Loy, 1,ysiz )
	sav_Vert_R = tvrd( Hix,Loy, 1,ysiz )

	if N_elements( color ) NE 1 then  color = !D.n_colors-1
	color = byte( color )

	horiz = replicate( color, xsiz, 1 )
	vert = replicate( color, 1, ysiz )

	tv, horiz, Lox,Loy
	tv, vert, Lox,Loy
	tv, horiz, Lox,Hiy
	tv, vert, Hix,Loy

return
end
;
;	SUBPROCEDURE # 3
;
pro box_erase
;+
; NAME:	
;	BOX_ERASE
; PURPOSE:
;	Erases boxes drawn by BOX_DRAW
; CATEGORY:
;	Display.
; CALLING SEQUENCE:
;	BOX_ERASE
; INPUTS:
;	None.
; OUTPUTS:
;	None.
; COMMON BLOCKS:
;	box
; SIDE EFFECTS:
; RESTRICTIONS:
;	For use with ZOOMW.
;	Intended for TV devices.
; PROCEDURE:
;	Corners of box and original data are obtained from common block: box
; MODIFICATION HISTORY:
;Frank Varosi NASA/GSFC 1989
;-

  common box, Lox,Loy,Hix,Hiy, sav_Horiz_B, sav_Horiz_T, sav_Vert_L, sav_Vert_R

	if N_elements( sav_Horiz_B ) LE 1 then return

	tv, sav_Horiz_B, Lox,Loy
	tv, sav_Vert_L, Lox,Loy
	tv, sav_Horiz_T, Lox,Hiy
	tv, sav_vert_R, Hix,Loy

	sav_Horiz_B = 0
return
end
;
;	SUBPROCEDURE # 4
;
pro box_save, RESTORE=restore
;+
; NAME:	
;	BOX_SAVE
; PURPOSE:
;	Saves or restores values in the common block: box
; CATEGORY:
;	Display.
; CALLING SEQUENCE:
;	BOX_SAVE [,/restore]
; INPUTS:
;	restore -- if this keyword is set then values from prevoius 
;	BOX_SAVE are restored, otherwise the values of the common block: box
;	are save in the common block: box_save
; OUTPUTS:
;	None.
; COMMON BLOCKS:
;	box, box_save
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
;	For use with ZOOMW.
;	Straightforwards.
; MODIFICATION HISTORY:	
;Frank Varosi NASA/GSFC 1989.
;-

  common box, Lox,Loy,Hix,Hiy, sav_Horiz_B, sav_Horiz_T, sav_Vert_L, sav_Vert_R

  common box_save, xLo,yLo,xHi,yHi, saveHB, saveHT, saveVL, saveVR

	if keyword_set( restore ) AND N_elements( saveHB ) GT 0 then begin

		Lox = xLo
		Loy = yLo
		Hix = xHi
		Hiy = yHi
		sav_Horiz_B = saveHB
		sav_Horiz_T = saveHT
		sav_Vert_L = saveVL
		sav_Vert_R = saveVR

	 endif else if N_elements( sav_Horiz_B ) GT 0 then begin

		xLo = Lox
		yLo = Loy
		xHi = Hix
		yHi = Hiy
		saveHB = sav_Horiz_B
		saveHT = sav_Horiz_T
		saveVL = sav_Vert_L
		saveVR = sav_Vert_R

	  endif

return
end
;
;	SUBPROCEDURE # 5
;
function find_dir, dir
;+
; NAME:	
;	FIND_DIR
; PURPOSE:
;	Checks for the existence of a directory.
; CATEGORY:
; CALLING SEQUENCE:
;	dirname = FIND_DIR(dir)
; INPUTS:
;	dir = string containing the name of directory to find
; OUTPUTS:
;	dirname = system-dependently formatted string of the directory
;		null string if directory is not found
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
;	For use with ZOOMW.
; MODIFICATION HISTORY:
;	 Frank Varosi NASA/GSFC 1989.
;-
	f = findfile()
	w = where( f EQ dir, nf )

	if (nf GT 0) then begin

		version = !version.arch + " " + !cgis_os
                if (!cgis_os EQ 'windows') then dirnam = dir + "\"
                if (!cgis_os EQ 'unix') then dirnam = dir + "/"
		if (version EQ "vax vms") or (version EQ "alpha vms") $
                       then dirnam = "[." + dir + "]"	
		return, dirnam

	  endif else return,""
end
;
;	SUBPROCEDURE # 6
;
pro find_files, filext, files, filnams, DIR=dir
;+
; NAME:	
; 	FIND_FILES
; PURPOSE:
; 	search for files with file-extension in current directory and
; 	and specified subdirectory (or all subdirectories if not specified).
; CATEGORY:
; CALLING SEQUENCE:
;	FIND_FILES, filext, files, filenams [,DIR=dir]
; INPUTS:
;	filext = the file-extension (string after .) to look for (e.g. "dat")
;	DIR = optional specification of subdirectory to check in addition to 
;		current directory. All subdirectories should be checked by 
;		default.
; OUTPUTS:
;	files = complete file names found
;	filnams = file names without the directory in name.
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
;
; MODIFICATION HISTORY:	
; 	Frank Varosi NASA/Goddard 1989
;-

	if N_elements( dir ) NE 1 then dir = "*"
	fs = "*." + filext

	version = !version.arch + ' ' + !cgis_os
	if (version EQ "vax vms") OR (version EQ "alpha vms") then begin

		fsd = "[." + dir + "]" + fs
		dir_delim = "]"
                ver_delim = ";"

	  endif else begin

		fsd = dir + "/" + fs
		dir_delim = "/"
	   endelse
         if (!cgis_os EQ 'windows') then begin
		fsd = dir + "\" + fs
		dir_delim = "\"
          ENDIF

	files = [ findfile( fsd ), findfile( fs ) ]
	filnams = files

	w = where( strlen( files ) , Nfil )
	if (Nfil LE 0) then return

	files = files(w)
	filnams = files
	pos = strpos( filnams, dir_delim ) + 1
	if (version EQ "vax vms") OR (version EQ "alpha vms") then begin
        Len = strpos( filnams, ver_delim ) - strlen( filext ) - 1 - pos
	endif else begin
        Len = strlen( filnams ) - strlen( filext ) - 1 - pos
        endelse
	Nfil = N_elements(filnams)

	for i=0,Nfil-1 do  filnams(i) = strmid( filnams(i), pos(i), Len(i) )
return
end
;
;	SUBPROCEDURE # 7
;
function get_words, text
;+
; NAME:	
;	GET_WORDS
; PURPOSE:
;	Breaks a text string into a strarr of words
; CATEGORY:
; CALLING SEQUENCE:
;	words = GET_WORDS(text)
; INPUTS:
;	text = string(s) with words delimited by blanks or commas.
; OUTPUTS:
;	words = strarr of single words
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
;	For use with ZOOMW.
;	Straightforwards.
; MODIFICATION HISTORY:
; 	Frank Varosi NASA/Goddard 1989
;-

	if N_elements( text ) LE 0 then return,""

	if N_elements( text ) GT 1 then begin
		Lens = strlen( text )
		Lenmax = max( Lens, maxi )
		text(maxi) = text(maxi) + " "
	   endif

	textb = byte( text )
	Len = N_elements( textb )
	if (Len LE 1) then return, string( textb )

	bb = byte( " " )
	bb = bb(0)

	cb = byte( "," )
	cb = cb(0)
	compos = where( textb EQ cb, ncomma )
	if (ncomma GT 0) then textb(compos) = bb

	tb = byte( 9 )
	tabpos = where( textb EQ tb, ntab )
	if (ntab GT 0) then textb(tabpos) = bb

	zpos = where( textb EQ 0, nzero )
	if (nzero GT 0) then textb(zpos) = bb

	bpos = where( textb EQ bb, nblank )
	if (nblank LE 0) then return, string( textb )

	words = strarr( nblank+1 )
	tend = [ bpos, Len-1 ]
	bpos = [ 0, bpos ]

	for i=0,nblank do  words(i) = string( textb( bpos(i):tend(i) ) )

	words = strtrim( words, 2 )
	w = where( words NE "" )

return, words(w)
end
;
;	SUBPROCEDURE # 8
;
function next_word, text
;+
; NAME:	
;	NEXT_WORD
; PURPOSE:
;	Extracts the first word from a text string
; CATEGORY:
; CALLING SEQUENCE:
;	word = NEXT_WORD(text)
; INPUTS:
;	text = string with words delimited by blanks.
; OUTPUTS:
;	word = the first word of the input string
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
;	For use with ZOOMW.
; 	Extract the first word of text string and return it,
;  	the input string text is also returned without the first word.
; MODIFICATION HISTORY:
;	Frank Varosi NASA/Goddard 1989
;-

	Len = strlen( text )
	if (Len LE 0) then return,""
AGAIN:
	blank_pos = strpos( text, " " )
	if (blank_pos LT 0) then  blank_pos = Len

	if (blank_pos EQ 0) then begin

		text = strmid( text, 1, Len )
		goto,AGAIN
	   endif

	word = strmid( text, 0, blank_pos )

	text = strmid( text, blank_pos+1, Len )

return, word
end
;
;	SUBPROCEDURE # 9
;
pro printw, strings, LINE=Line, ERASE=erase, WINDOW=window, XOFFSET=xoffset
;+
; NAME:	
;	PRINTW
; PURPOSE:
; 	treat a window as a non-scrolling text screen and print one 
;	string per line.
; CATEGORY:
;	Display.
; CALLING SEQUENCE:
;	PRINTW, strings [,LINE=line] [,ERASE=erase] [,WINDOW=window] 
;		[,XOFFSET=xoffset]
; INPUTS:
; 	Line specifies starting point as # lines from top of window if negative,
;      		# lines from bottom of window if Line is non-negative.
; 		(Line=-1 is at top and Line=0 is at bottom, default value is -1)
; 	/ERASE will erase the region first,
;	WINDOW specifies the window which will be modified
; 	XOFFSET will horizontally shift by # characters to right each line 
;		(deflt=1).
; OUTPUTS:
;	None.
; COMMON BLOCKS:
; SIDE EFFECTS:
;	Text is written to a window.
; RESTRICTIONS:
; PROCEDURE:
;	Used by ZOOMW
; MODIFICATION HISTORY:
;	Frank Varosi NASA/GSFC 1989.
;-

	Nline = N_elements( strings )
	if (Nline LE 0) then return
	if (Nline EQ 1) then strings = [strings]
	if N_elements( Line ) NE 1 then Line = -1

	if N_elements( window ) EQ 1 then begin
		wset, window
		wshow, window
	   endif

	if (Line GE 0)  then  ystart = (Nline-1+Line) * !D.y_ch_size + 6   $
			else  ystart = Line * !D.y_ch_size + !D.y_vsize

	Nxoff = N_elements( xoffset )
	if (Nxoff LE 0) then xoff = replicate( 1, Nline ) else $
	 if (Nxoff EQ 1) then xoff = replicate( xoffset, Nline ) else $
	  if (Nxoff LT Nline) then $
	      xoff = [ xoffset, replicate( xoffset(Nxoff-1), Nline-Nxoff ) ]

	xoff = (( xoff * !D.x_ch_size ) > 0 ) < !D.x_vsize
	yline = (( ystart - indgen( Nline ) * !D.y_ch_size ) > 0 ) < !D.y_vsize

	if keyword_set( erase ) then begin

		if (Nline GT 1) then xLeft = min( xoff ) else xLeft = xoff(0)
		zero = bytarr( !D.x_vsize - xLeft, Nline * !D.y_ch_size )
		ybot = ( ystart - (Nline-1) * !D.y_ch_size ) > 0
		tv, zero, xLeft, ybot < (!D.y_vsize-1)
	   endif

	for i=0,Nline-1 do  xyouts, /DEV, xoff(i), yline(i), strings(i), FONT=0
	empty
return
end
;
;	SUBPROCEDURE # 10
;
pro hist_equal_ct_zoomw, vimage, WINDOW=image_window, INFO_WINDOW=info_window
;+
; NAME:
;	HIST_EQUAL_CT_ZOOMW
; PURPOSE:
;	Histogram equalize the color tables of an image or a region of display.
; CATEGORY:
;	Image processing.
; CALLING SEQUENCE:
;	HIST_EQUAL_CT_ZOOMW, VImage	... to histogram equalize from an image.
;	HIST_EQUAL_CT_ZOOMW		... to histogram equalize from a region
; INPUTS:
;      VImage = image whose histogram is to be used in determining
;		the new color tables.  If omitted, the user is prompted
;		to mark the diagonal corners of a region of the display.
;		The Image MUST be a byte image, scaled the same way as
;		the image loaded to the display.
;	WINDOW = keyword specifying the window for image region to equalize.
; OUTPUTS:
;	No explicit outputs.  The result is applied to the current color
;	tables.
; COMMON BLOCKS:
;	colors -- the color tables.
; SIDE EFFECTS:
;	Color tables are updated.
; RESTRICTIONS:
;	If a parameter is supplied, it is assumed to be an image that
;	was just loaded.
; PROCEDURE:
;	For use with ZOOMW.
;	Either the image parameter or the region of the display marked by
;	the user is used to obtain a pixel distribution histogram.  The
;	cumulative integral is taken and scaled.  This function is applied
;	to the current color tables.
; MODIFICATION HISTORY:
;	DMS, March, 1988, written.
;	modified by Frank Varosi NASA/GSFC 1989  to Loop and use rubber-box
;	(box_create, box_draw, and box_erase routines required).
;-

common colors, r,g,b, cur_red, cur_green, cur_blue
common ascii, BELL, LF
common hist_ct, h

nc = !d.n_colors	;# of colors in device

if n_elements(r) le 0 then begin		;color tables defined?
	r=indgen(nc) & g=r & b=r & endif

if n_elements(vimage) gt 0 then begin

	h = histogram(vimage)
	for i=1,n_elements(h)-1 do h(i) = h(i)+h(i-1)
	h = long(bytscl(h, top = nc-1))

	cur_red = r(h)
	cur_green = g(h)
	cur_blue = b(h)

	tvlct,cur_red, cur_green, cur_blue

 endif else begin

	if N_elements( image_window ) NE 1 then  image_window = !D.window

	if N_elements( info_window ) EQ 1 then begin

		wset, info_window
		save_window = tvrd( 0,0, !D.x_vsize, !D.y_vsize )
		erase
		func = [ "function = HIST_EQUAL_CT_ZOOMW", " " ]
		printw, func, LINE=-3
		instructions = ["LEFT & MIDDLE button:"			$
				+ " mark FIRST & SECOND corner of AREA,"   ,$
			"MIDDLE button: exit and keep new color table, "   ,$
			"RIGHT button: exit and restore orignal colors"," " ]
		printw, instructions, LINE=-5
		wshow, info_window

	  endif else begin

		print," LEFT button: mark FIRST corner of AREA," $
			+ " then SECOND corner with MIDDLE button"
		print," MIDDLE button: exit and keep new color table"
		print," RIGHT button: exit and restore orignal colors"
	    endelse

	wset, image_window
	wshow, image_window
	tvcrs, .3, .3, /NORM
	box_save	;in case there is an existing box on display, save it.

DEFBOX:
	button = box_create( x0,y0, x1,y1 )

	CASE button OF

	-2: BEGIN
		box_erase
		box_save, /RESTORE
		goto,EXIT
	     END

	-4: BEGIN
		tvlct, r, g, b
		box_erase
		box_save, /RESTORE
		goto,EXIT
	     END

	else:
	ENDCASE

	h = histogram( tvrd( x0,y0,x1-x0+1, y1-y0+1 ) )

	for i=1,n_elements(h)-1 do h(i) = h(i)+h(i-1)
	h = Long( bytscl( h, top = nc-1 ) )

	cur_red = r(h)
	cur_green = g(h)
	cur_blue = b(h)

	tvlct, cur_red, cur_green, cur_blue

	goto,DEFBOX
EXIT:
	if N_elements( info_window ) EQ 1 then begin

		wset, info_window
		tv, save_window
	   endif
  endelse

return
end
;
;	SUBPROCEDURE # 11
;
pro color_edit_bg		
;+
; NAME:
;	COLOR_EDIT_BG
; PURPOSE:
; 	create background for:  color_edit_zoomw
; CATAGORY:
; CALLING SEQUENCE:
;	COLOR_EDIT_BG
; INPUTS:
;	None.
; OUTPUTS:
;	None.
; COMMON BLOCKS:
;	color_edit_zoomw
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
;	Used by COLOR_EDIT_ZOOMW as part of ZOOMW
; MODIFICATION HISTORY:
;	Frank Varosi NASA/GSFC 1989?
;-

common color_edit_zoomw, nc, nc1, nc2, wxsize, wysize,r0,cx,cy,bar_ht, $
	colors, plot_xs, plot_ys, order, $
	bar_wid, bar_x0, bar_y, pv_y, pv_wid, pv_x0, names

ramp = bytscl(indgen(512),top=nc1)
for i=wysize-60,wysize-30 do tv,ramp,wxsize/2-256,i

r0 = 120  &  cx = wxsize/4.  &  cy = wysize - 90. - r0
angle = !dtor * 10 * findgen(37)
sina = sin(angle)
cosa = cos(angle)
plots,[cx,cx],[cy,cy],/dev	;Center of circle
polyfill,/dev,r0*cosa+cx,r0*sina+cy,col=nc2
for i=0.5, 1.0, 0.5 do plots,r0*i*cosa+cx, r0*i*sina+cy,/dev
c_name = ['Red','Yellow','Green','Cyan','Blue','Magenta']
align = [0,0,1,1,1,0]
for i=0,5 do begin
	a = i*60 * !dtor
	x = r0*cos(a) + cx
	y = r0 * sin(a) + cy
	plots,[cx,x],[cy,y],/dev
	xyouts,r0*1.1*cos(a)+cx,r0*1.1*sin(a)+cy,c_name(i),/dev,ali = align(i)
	end

bar_ht = 20		;Draw the brightness bar
bar_wid = wxsize/3
bar_x0 = wxsize/4-bar_wid/2
bar_y = cy - r0 - 40 - bar_ht
xyouts,wxsize/4, bar_y+bar_ht+3, names(order(2)), ali=0.5,/dev
xyouts,bar_x0,bar_y+bar_ht+3,ali=0,/dev,'0'
xyouts,bar_x0+bar_wid,bar_y+bar_ht+3,ali=1,/dev,'1.0'
polyfill,[bar_x0, bar_x0, bar_x0+bar_wid, bar_x0+bar_wid],$
	[ bar_y, bar_y+bar_ht, bar_y+bar_ht, bar_y],/dev

pv_y = bar_y - bar_ht - 30	;Pixel value bar
pv_wid = wxsize/3.
pv_x0 = wxsize/4-pv_wid/2
x = bytscl(findgen(pv_wid-1),top=nc-1)
for i=pv_y, pv_y+bar_ht do tv,x,pv_x0,i ;Another ramp
xyouts,wxsize/4, pv_y+bar_ht+3,'Pixel Value',ali=0.5,/dev
xyouts,pv_x0, pv_y+bar_ht+3,ali=0,/dev,'0'
xyouts,pv_x0+pv_wid, pv_y+bar_ht+3,ali=1,/dev,strtrim(nc1,2)
plots,[pv_x0, pv_x0, pv_x0+pv_wid, pv_x0+pv_wid,pv_x0],$
	[ pv_y, pv_y+bar_ht, pv_y+bar_ht, pv_y,pv_y],/dev

plot_xst = .6
plot_xend = .9
plot_ht = 0.2
yr = [360.,1.,1.]
plot_position = fltarr(4,3)
plot_xs = fltarr(2,3)
plot_ys = plot_xs
for i=0,2 do begin
	y = i/3.5+0.1	;Y of bottom
	plot_position(0,i) = [plot_xst,y, plot_xend, y+plot_ht]
	plot,/noer,colors(*,order(i)),yrange=[0,yr(i)],title=names(i),$
		pos=plot_position(*,i),ystyle=2,xstyle=3, tickl = -0.02
	plot_xs(0,i) = !x.s 
	plot_ys(0,i) = !y.s
	endfor
end
;
;	SUBPROCEDURE # 12
;
pro color_edit_zoomw, colors_out, hsv = hsv, hls = hls
;+
; NAME:		COLOR_EDIT_ZOOMW
; PURPOSE:	Interactive creation of color tables based on
;		the HLS or the HSV color systems using the mouse
;		and a color wheel.
; CATEGORY:	Color tables.
; CALLING SEQUENCE:
;	COLOR_EDIT_ZOOMW
; INPUTS:
;	None.
; KEYWORD PARAMETERS:
;	HLS = non zero to use the Hue Lightness Saturation system.
;	HSV = non zero to use the Hue Saturation Value system. (The Default)
; OUTPUTS:
;	colors_out = optional output parameter.  Will contain
;	  (Number_colors, 3) color values of final color table.
; COMMON BLOCKS:
;	COLORS - contains the current RGB color tables.
; SIDE EFFECTS:
;	Color tables are modified, values in COLORS common are changed.
;	A temporary window is used.
; RESTRICTIONS:
;	Only works with window systems.
; PROCEDURE:
;	A window is created with:
;	1) a color bar centered at the top.
;	2) A color wheel for the hue (azimuth of mouse from the center) 
;		and saturation (distance from the center).
;	3) Two "slider bars", the top one for the Value (HSV) or Lightness
;		(HLS), and the other for the pixel value.
;	4) 3 graphs showing the current values of the three parameters
;		versus pixel value.
;
;	The method is:	The left mouse button is used to mark values in the
;		left half of the window.  The middle button is used to erase
;		marked pixel values in the Pixel Value slider.  The right
;		button exits the procedure with the color tables updated.
;	To use: move the mouse into the circle or middle slider and depress
;		the left button to select a (hue, saturation) or a Value
;		or Lightness.  You can move the mouse with this button
;		depressed to select a color (shown in the circle).
;		  Next, move the mouse to the bottom slider, and select
;		a pixel value.  The three color parameters are interpolated
;		between pixel values that have been marked (called tie points).
;		Use the middle button on the Pixel Value slider to delete 
;		the nearest tie point.
;	Note that with the HSV system, a Value of 1.0 is maximum brightness
;		of the selected hue.  In the HLS system, a Lightness of 0.5
;		is maximum brightness of a chromatic hue, 0.0 is black,
;		and 1.0 is bright white.  In the HLS system, which models a 
;		double ended cone, the Saturation has no effect at the
;		extreme ends of the cone (lightness = 0 or 1).
;	You can access the new color tables by declaring the common block
;		colors, both the original and current arrays are set:
;   common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr
; MODIFICATION HISTORY:
;	DMS, July, 1988.
;	FV, 1989, load orig/curr colors if they exist,
;		color_edit_bg  and  color_interp  are in separate files.
;-

common color_edit_zoomw,nc, nc1, nc2, wxsize, wysize,r0,cx,cy,bar_ht, $
	colors, plot_xs, plot_ys, order, $
	bar_wid, bar_x0, bar_y, pv_y, pv_wid, pv_x0, names

common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

psave = !p		;Save !p
device,get_write=old_mask, set_write=255 ;Enable all bits
nc = !d.n_colors	;# of colors avail
nc1 = nc -1
nc2 = nc1-1		;Current color
!p.noclip = 1		;No clipping
!p.color = nc1		;Foreground color
!p.font = 0		;Hdw font
old_window = !d.window	;Previous window

if n_elements(hsv) eq 0 then hsv = 0
if n_elements(hls) eq 0 then hls = 0
if (hsv eq 0) and (hls eq 0) then hsv = 1	;default system

if hsv then begin
	names = ['Hue','Saturation','Value']
	order = [0,1,2]
endif else begin
	names = ['Hue','Lightness','Saturation']
	order = [0,2,1]
endelse

colors = [[fltarr(nc)],[replicate(1.0,nc)],[findgen(nc)/nc]]
if (windows = (!d.flags and 256) ne 0) then begin  ;Windows?
	wxsize = 640
	wysize = 600
	window,xs=wxsize, ys=wysize, title='Intensity transformation',/free
  endif else begin
	wxsize = !d.x_vsize
	wysize = !d.y_vsize
  endelse

color_edit_bg		;create background for color_edit_zoomw in window
tvcrs,.5,.5,/norm

if N_elements( r_curr ) GT 0 then tvlct,r_curr,g_curr,b_curr else	$
if N_elements( r_orig ) GT 0 then tvlct,r_orig,g_orig,b_orig else	$
	begin
		colors(nc1,*) = [0,0,1.]
	tvlct,colors(*,order(0)),colors(*,order(1)),colors(*,order(2)), $
	hsv = hsv, hls = hls
	endelse

npts = 2		;tie points
pts = indgen(nc)	;x values of tie points
pts(1) = nc1		;init values
h = 0.
s = 1.0
v = 1.0
pxl_ind = 0
old_v = 0.
;		*** Main loop ***

next:
tvrdc,x,y,/dev		;read mouse with wait
next1:
if !err eq 1 then begin
next2:
	tvrdc,x,y,0,/dev
	if !err ne 1 then goto,next
	if (y ge pv_y) and (y le pv_y + bar_ht) then goto, mark_pixel
	reload = 0
	d = sqrt((x - cx)^2 + (y-cy)^2) ;In circle?
	if d le (1.05 * r0) then begin	;New hue, sat & lightness
		h = atan(x-cx,y-cy) * !radeg ;hue
		s = d/r0 < 1.0		;saturation
		reload = 1
		endif
	if (y ge bar_y) and (y le bar_y+bar_ht) then begin ;lightness
		v = (x - bar_x0) / float(bar_wid) ;fract of brightness
		v = v > 0.0 < 1.0
		xyouts, bar_x0, bar_y,string(old_v,format='(f4.2,1x)'),$
			ali=1.,/dev,col=0
		xyouts, bar_x0, bar_y,string(v,format='(f4.2,1x)'),ali=1.,/dev
		old_v = v
		reload = 1
		endif
	if reload then if hls then tvlct,h,v,s,/hls,nc2 else $
				tvlct,h,s,v,/hsv,nc2
	!err = 0
	goto,next2

mark_pixel:		;Get a hit in the pixel value Y values
	x = x > pv_x0 < (pv_x0 + pv_wid)
	pxl_ind = fix(nc * float(x - pv_x0) / pv_wid) < nc1 > 0
	x = pv_x0 + float(pxl_ind) * pv_wid/nc
	plots,[x,x],[pv_y-8,pv_y],/dev
	p = where(pxl_ind eq pts(0:npts-1), n)
	if n eq 0 then begin	;already there?
		pts(npts) = pxl_ind
		pts(0) = pts(sort(pts(0:npts))) ;re sort
		npts = npts + 1
	endif
	while h lt 0 do h = h + 360.	;always positive

interp_it:
	 for i=0,2 do begin	;erase old plots
		!x.s = plot_xs(*,i)
		!y.s = plot_ys(*,i)
		oplot,colors(*,order(i)),col=0
		endfor
	colors(pxl_ind,*) = [h,s,v]
	color_interp, pts, npts, colors		;interpolate colors
	tvlct,colors(1:nc2-1,order(0)),colors(1:nc2-1,order(1)), $
		colors(1:nc2-1,order(2)),1, hls = hls, hsv = hsv
	for i=0,2 do begin	;draw new plots
		!x.s = plot_xs(*,i)
		!y.s = plot_ys(*,i)
		oplot,colors(*,order(i))
		endfor
endif $		;!err eq 1
;		here we delete a point:

else if (!err eq 2) and  (y ge pv_y) and (y le pv_y + bar_ht) then  begin
	x = x > pv_x0 < (pv_x0 + pv_wid)
	pxl_ind = fix(nc * float(x - pv_x0) / pv_wid) < nc1 > 0
	j = min(abs(pts(0:npts-1) - pxl_ind),i)	;get index of closest point
	x = pv_x0 + float(pts(i)) * pv_wid / nc
	plots,[x,x],[pv_y-8,pv_y-1],/dev,col=0 ;erase tick
		 ;never	delete first or last points
	if (i eq 0) or (i eq nc1) then goto, next
	pxl_ind = pts(i)	;old pixel index
	pts = [pts(0:i-1), pts(i+1:*)] ;remove it
	npts = npts - 1
	goto, interp_it
endif $	
else if !err eq 4 then goto,done		;all done

goto,next


done:	tvlct,r_orig, g_orig, b_orig,/get	;Read rgb, save in common
	r_curr = r_orig & g_curr = g_orig & b_curr = b_orig
	if n_params() ge 1 then colors_out = [r_orig, g_orig, b_orig]

	if (windows) then begin	;Using windows?
		wdelete			;kill window
		if old_window ge 0 then begin	;restore window?
			tvcrs,0.5,0.5,/norm	;show the table
			empty
			tvcrs			;hide cursor
			endif
	endif		;windows

!p = psave			;restore !p
device,set_write=old_mask	;Restore write mask

end
;
;	SUBPROCEDURE # 13
;
pro color_interp, pts, npts, colors	
					
;+
; NAME:	
;	COLOR_INTERP
; PURPOSE:
;	interpolates between selected colors in a color table
; CATEGORY:
;	Display.
; CALLING SEQUENCE:
;	COLOR_INTERP, pts, npts, colors
; INPUTS:
;	pts = array of subscripts of tie points.
;	npts = # of elements in pts
;	colors = (n,3) array of colors.  Interpolated between tie points.
; OUTPUTS:
;	None. 
; COMMON BLOCKS:
; SIDE EFFECTS:
;	The color table is modified.
; RESTRICTIONS:
; PROCEDURE:
; 	called by color_edit_zoomw; to be used with ZOOMW
; MODIFICATION HISTORY:
;	Frank Varosi NASA/GSFC 1989.
;-
;
	for i=0,npts-2 do begin			;interpolate
	   i0 = pts(i) & i1 = pts(i+1)
	   kc = i1 - i0 			;# of colors to interp+1
	   if (kc gt 1) then begin		;do it?
		c = colors(i0,*)
		dc = colors(i1,*) - c			;delta clockwise
		dc0 = dc(0)				;hue dist clockwise
		while (dc0 lt (-180.)) do dc0 = dc0 + 360.
		while (dc0 gt 180.) do  dc0 = dc0 - 360.
		dc1 = -dc(0)				;delta ccw
		while (dc1 lt (-180.)) do dc1 = dc1 + 360.
		while (dc1 gt 180.) do dc1 = dc1 - 360.
		if abs(dc1) lt abs(dc0) then begin	;Use closest
			dc(0) = dc1
		  endif else begin
			dc(0) = dc0
		  endelse
		dc = dc / kc
		colors(i0+1,0) = ((findgen(kc-1)+1) # dc) + $  ;Interpolate
				 (replicate(1,kc-1) # c)
	     endif					;kc gt 1
	  endfor					;i loop

	colors(0,0) = (colors(*,0) + 360.) mod 360.	;wrap the hue
end
;
;	SUBPROCEDURE # 14
;
pro color_tables, caller, $
	IMAGE_WINDOW=image_win, MENU_WINDOW=menu_win, INFO_WINDOW=info_window
;+
; NAME:	
;	COLOR_TABLES
; PURPOSE:
;	Allows interactive manipulation of color tables in many different modes 
; CATEGORY:
;	Display.
; CALLING SEQUENCE:
;	COLOR_TABLES [,caller] [,IMAGE_WINDOW=image_win] $
;		[,MENU_WINDOW=menu_win] [,INFO_WINDOW=info_window]
; INPUTS:
;	(All optional)
;	caller = name of the procedure calling COLOR_TABLES
;	IMAGE_WINDOW = window containing image of interest
;	MENU_WINDOW = window for display of menu of options
;	INFO_WINDOW = window for display of informational messages
; OUTPUTS:
;	None.
; COMMON BLOCKS:
;	colors
; SIDE EFFECTS:
;	The color table is modified.
; RESTRICTIONS:
; PROCEDURE:
;	Called by ZOOMW.
; MODIFICATION HISTORY:	
;	Frank Varosi NASA/GSFC 1989.
;-
   common colors, ro,go,bo, rc,gc,bc

	if N_elements( caller ) NE 1 then caller = "calling routine"
	if N_elements( image_win ) EQ 1 then  image_window = image_win $
					else  image_window = !D.window
	if N_elements( menu_win ) EQ 1  then  menu_window = menu_win   $
					else  menu_window = !D.window
	if (menu_window LT 0) then menu_window = !D.window

	menu = ["Color Table Options:"		,$
		"Adjust color table"		,$
		"load IDL color table"		,$
		"restore RGB colors"		,$
		"save RGB colors"		,$
		"histogram equalize color table",$
		"edit (create) colors"		,$
		"palette edit RGB"		,$
		"pause  (click mouse to resume)",$
		"return to " + caller		]
	sel = 1
	BELL = string( byte( 7 ) )
	LF = string( byte( 10 ) )

	version = !version.arch + ' ' + !cgis_os
	if (version EQ "vax vms") OR (version EQ "alpha vms")then  $
			idl_col_tbl = !DIR + "[000000]colors1.tbl"  $
		else    idl_col_tbl = !DIR + "/colors1.tbl"
        if (!cgis_os EQ 'windows') then idl_col_tbl = !DIR + "\colors1.tbl"

	get_lun, lun
	openr,lun, idl_col_tbl			;Get color table names
	ctrec = assoc(lun, bytarr(32,16))
	ctrec2 = string(ctrec(0))
	menuct = ["Select color table:", strtrim( indgen(16),2 ) + " : "]
	for i=1,16 do  menuct(i) =  menuct(i) + ctrec2(i-1)
	free_Lun,Lun
	ct = 12
	menuct = [menuct, "adjust color table"	  ,$
			  "return to options menu",$
			  "return to " + caller	  ]

MENU:	wset, menu_window
	wshow, menu_window
	wait,.1
	sel = wmenu( menu, INIT=sel, TITLE=0 )
	if (sel LE 0) then return

	task = next_word( menu(sel) )

	CASE task OF

	"Adjust":	adjct

	"load": BEGIN
			ctm = wmenu( menuct, INIT=ct+1, TITLE=0 )
			while (ctm GT 0) do begin
				request = menuct(ctm)
				Load = next_word( request )
				CASE Load OF
				"adjust":	adjct
				"return": BEGIN
						ret = get_words( request )
						CASE ret(1) OF
						"options":	goto,MENU
						     else:	return
						ENDCASE
					    END
				else: BEGIN
					ct = fix( Load )
					Loadct,ct
					END
				ENDCASE
				ctm = wmenu( menuct, INIT=ct+1, TITLE=0 )
			  endwhile
		  END

	"edit":		color_edit_zoomw

	"palette":	palette

	"histogram": BEGIN
			if (image_window LT 0) then begin

				print,LF+" image window is not available"+BELL
				wait,1

			  endif else   hist_equal_ct_zoomw, $ 
				WINDOW=image_window, INFO=info_window
		     END

	"pause": BEGIN
			wset, menu_window
			wshow, menu_window
			tvcrs, 8,8, /DEV
			cursor, x,y, /DEV
		   END

	"save": BEGIN
			filnam = ""
			read," enter file name for saving RGB colors " + $
                          "(omit extension): ",filnam
			if (filnam EQ "") then goto,MENU
			if N_elements( dir_colors ) NE 1 then $
					dir_colors = find_dir( "colors" )
			filnam = dir_colors + filnam + ".rgb"
			save,rc,gc,bc,file=filnam
			print," saved Red Green Blue components to: ",filnam
			filrgb = 0
		  END

	"restore": BEGIN
			if N_elements( filrgb ) LE 1 then begin
				find_files, "rgb", filrgb, filnams, DIR="colors"
				Nfil = N_elements( filrgb )
				sf = 0
				menu_RGB = ["Select file:", filnams	,$
					    "adjust color table"	,$
					    "return to options menu"	,$
					    "return to " + caller	]
			   endif
			sfm = wmenu( menu_RGB, INIT=sf+1, TITLE=0 )
			while (sfm GT 0) do begin
				request = menu_RGB(sfm)
				Load = next_word( request )
				CASE Load OF
				"adjust":	adjct
				"return": BEGIN
						ret = get_words( request )
						CASE ret(1) OF
						"options":	goto,MENU
						     else:	return
						ENDCASE
					    END
				else: BEGIN
					sf = sfm-1
					filnam = filrgb(sf)
					restore,filnam
					tvlct,rc,gc,bc
					ro = rc
					go = gc
					bo = bc
					END
				ENDCASE
				sfm = wmenu( menu_RGB, INIT=sf+1, TITLE=0 )
			  endwhile
		     END

	"return":	return

	ENDCASE

   goto,MENU
end
;
;	SUBPROCEDURE # 15
;
function set_zoom, maxz, menu_window, INIT=init
;+
; NAME:	
;	SET_ZOOM
; PURPOSE: 
;	Gets a new zoom factor for ZOOMW
; CATEGORY:
;	Display.
; CALLING SEQUENCE:
;	zoom_factor = SET_ZOOM(maxz [,menu_window] [,INIT= init]
; INPUTS:
;	maxz = maximum zoom factor to list in menu
;	optional:
;	menu_window = index of window to set
;	init = initially selected choice of zoom factor
; OUTPUTS:
;	zoom_factor = new zoom_factor to be used by ZOOMW
; COMMON BLOCKS:
; SIDE EFFECTS:
; RESTRICTIONS:
; PROCEDURE:
;	For use with ZOOMW.
; MODIFICATION HISTORY:
;Frank Varosi NASA/GSFC 1989.
;-

	if N_elements(menu_window) EQ 1 then wset,menu_window
	if N_elements(init) EQ 1 then initz=fix( init ) else initz=2

	menu = ['Zoom Factors:', strtrim( indgen(maxz)+1, 2 ) ]

	zoom_factor = wmenu( menu, init=initz, title=0 )

	if (zoom_factor LT 1) then zoom_factor = 1

return, zoom_factor
end

;-----	MAIN PROCEDURE
;
;-----  Descriptive comments are found at the beginnning of this file.
;
pro zoomw, cur_w, XSIZE=xs, YSIZE=ys, FACT=fact, CONTINUOUS=cont, INFO=info_w
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     ZOOMW zooms window and alters color tables on windowing devices
;
;DESCRIPTION:
;     Displays part of an image (or graphics) from the current window
;     expanded in another window. This window is retained for future 
;     use after the procedure finishes. The cursor is used to mark the
;     center of the zoom, and make selections from menus to alter the
;     color table. Instructions are provided in the session text 
;     window.
;
;CALLING SEQUENCE:
;     ZOOMW [,cur_w] [,FACT=fact] $ 
;               [,XSIZE=xsize] [,YSIZE=ysize] [,/CONTINUOUS]
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     cur_w        [I]  int        window to zoom
;     fact         [I]  flt        zoom expansion factor, default = 4
;     xsize        [I]  int        X size of new window, default = 500
;     ysize        [I]  int        Y size of new window, default = 500
;     /continuous  [I]  key        If specified, the zoom window 
;                                  tracks the mouse.
;
;WARNINGS:
;     1.  The display window to be zoomed should exist.
;     2.  Only works with color windowing systems.
;     3.  A window is created and retained after procedure finishes.
;     4.  If the option to save color tables as disk files is used,
;         then the directory to which the files are written must have 
;         appropriate write protection and sufficient free space 
;         (>5 blocks/file should do).
;
;EXAMPLE:
;     After displaying an image in a window using TV or TVSCL, suppose
;     you want to get a better look at fine detail in the image. Type:
;     
;       zoomw
;
;     A crosshair cursor will appear in the window. Center the cursor 
;     at the region of interest, and click the left mouse button (MB1)
;     to create an enlarged view in a new window. Move the cursor and
;     click again to zoom a different region. Clicking the center 
;     mouse button (MB2) will present a menu which allows you to reset 
;     the size of the zoom window and the magnification factor. This 
;     menu also leads to another menu providing access to LOADCT, 
;     ADJCT, and other procedures for altering, saving, and restoring 
;     color tables. Clicking the right mouse button (MB3) will exit 
;     from the procedure.
;     
;     The more complete command:
;
;       zoomw,2,fact=5,xsize=200,ysize=250,/continuous
;
;     first switches the active display window to window 2 (i.e. it 
;     performs WSET,2), and then provides a new 200x250 pixel window
;     in which a portion of window 2 is displayed at 5x magnification.
;     The size and magnification can easily be altered (click the 
;     middle mouse button, MB2), so it is seldom necessary to set
;     these in the command line. The /continuous keyword causes the 
;     magnified region to track the movements of the cursor without 
;     repeated clicks of the left mouse button (MB1). This only works 
;     well on fast computers which can recalculate the zoom window 
;     display as fast as the cursor moves. The COBECL computers do not 
;     normal qualify as fast.
;#
;COMMON BLOCKS:
;     zoomw1, zoomw2, zoomw3, ascii, box, box_save, colors, 
;     color_edit_zoomw, hist_ct
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     ZOOMW.PRO contains the procedure zoomw and 15 other procedures
;     and functions which were written to support zoomw, but could
;     be used for other procedures as well. The first four subroutines 
;     deal with drawing and erasing boxes on the TV screen. The next
;     pair deal with finding files and directories. Next are two for
;     processing text strings, and one for printing text to a window.
;     Five more subroutines provide some of the various aspects of
;     color table manipulation which are not found in the standard 
;     IDL user library. The last subroutine is used for read in 
;     changes in the zoom factor.
;
;PERTINANT ALGORITHMS, LIBRARY CALLS, ETC.:
;     adjct, loadct, palette, tvrdc
;
;MODIFICATION HISTORY:
;     Frank Varosi NASA/GSFC 1989: This version evolved from zoom.pro 
;     in IDL library with modifications to keep window, added menu 
;     items, show zoom box, etc. 
;     Rick Arendt ARC Nov 1991: Made minor modifications to run under
;     IDL v.2 on the COBECL 
;     Rick Arendt ARC Feb 1992: combined all Varosi subroutines into 
;     zoomw.pro.
;
; SPR 9616
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;.TITLE
;Routine ZOOMW
;-
;
   common zoomw1, zoom_w, zfact, zwxs, zwys, old_w, zx,zy
   common zoomw2, xwsizes, ywsizes, zoomwin
   common zoomw3, windows, winames
   common ascii, BELL, LF

if n_elements(xs) EQ 1 then zwxs = xs
if n_elements(ys) EQ 1 then zwys = ys
if n_elements(fact) EQ 1 then zfact=fact
if keyword_set(cont) then waitflg = 2 else waitflg = 1

if n_elements(zfact) le 0 then zfact=4
if N_elements( xwsizes ) LE 0 then xwsizes = [300,500,300,500,800,500,800]
if N_elements( ywsizes ) LE 0 then ywsizes = [300,300,500,500,500,800,700]
if n_elements( zwxs ) NE 1 then zwxs = xwsizes(3)
if n_elements( zwys ) NE 1 then zwys = ywsizes(3)

menu = ["Zoom Options:"		,$
	"Set Zoom Factor"	,$
	"Resize Zoom Window"	,$
	"Color Tables"		,$
	"return to main menu"	]
task = ""
menu_size = ["Select Size:",strtrim( xwsizes,2) + " x " + strtrim( ywsizes,2)]

if N_elements( cur_w ) EQ 1 then begin

	if (cur_w LT 0) then begin
		menuw = ["Select window to zoom:", winames]
		sel = wmenu( menuw, INIT=1, TIT=0 ) -1
		if (sel GE 0) then  cur_w = windows(sel)  else  return
		func = [ "function = ZOOM of " + winames(sel), " " ]
	   endif

  endif else  cur_w = !D.window

if N_elements( info_w ) EQ 1 then begin

	printw, func, LINE=-3, /ERASE, WIN=info_w
	instructions = ["LEFT button for zoom CENTER, "	,$
			"MIDDLE button for zoom MENU, "	,$
			"RIGHT button to quit zoom"," "	]
	printw, instructions, LINE=-5, /ERASE

  endif else print," LEFT for zoom CENTER, MIDDLE for zoom MENU, RIGHT to quit"

wset, cur_w
wshow, cur_w
if N_elements( old_w ) NE 1 then begin
	zx = !D.x_vsize/2
	zy = !D.y_vsize/2
 endif else if (cur_w NE old_w) then begin
	zx = !D.x_vsize/2
	zy = !D.y_vsize/2
  endif

tvcrs, zx, zy, /DEV
old_w = cur_w

AGAIN:  cursor,/DEV, x,y, waitflg	;Wait for change (main Loop)

CASE !err OF

4:   goto,EXIT

2:   BEGIN

	sel = wmenu( menu, INIT=1, TITLE=0 )
	task = next_word( menu(sel>0) )

	CASE task OF

	"Set":	zfact = set_zoom( 14, INIT=zfact )

	"Resize": BEGIN

		sel = wmenu( menu_size, INIT=2, TITLE=0 ) -1

		if (sel GE 0) then begin

		   zwxs = xwsizes(sel)
		   zwys = ywsizes(sel)

		   if N_elements( zoom_w ) GT 0 then wdelete,zoom_w

		   window, /FREE, XSIZ=zwxs, YSIZ=zwys, TIT="Zoom Window"
		   zoom_w = !D.window
		   zoomwin = zoom_w

		 endif
		END

	"Color":    color_tables, "zoom", IMAGE=cur_w, MENU=zoom_w, INFO=info_w

	"return":	goto,EXIT

	else:
	ENDCASE

	wset,cur_w
	wshow,cur_w
	tvcrs,x,y,/DEV	;Restore cursor
     END

else:	BEGIN

	if N_elements( x0 ) EQ 1 then  box_erase
	zx = x					;save zoom center.
	zy = y
	x0 = (x-zwxs/(zfact*2)) > 0 		;left edge from center
	y0 = (y-zwys/(zfact*2)) > 0 		;bottom
	nx = (zwxs/zfact) < (!D.x_vsize - x0)	;Size of new image
	ny = (zwys/zfact) < (!D.y_vsize - y0)

	area = tvrd( x0,y0, nx,ny )		;Read image
	box_draw, x0,y0, x0+nx-1, y0+ny-1

	if N_elements( zoom_w ) LE 0 then begin       ;make new window if none

		window, /FREE, XSIZE=zwxs, YSIZE=zwys, TITLE='Zoom Window'
		zoom_w = !D.window
		zoomwin = zoom_w

	 endif else if ((zoom_w LT 0) or (execute('wset,zoom_w') NE 1)) $
		 then begin       ;make window if deleted

		window, /FREE, XSIZE=zwxs, YSIZE=zwys, TITLE='Zoom Window'
		zoom_w = !D.window 
		zoomwin = zoom_w

	  endif else begin

		wset,zoom_w
		wshow,zoom_w
	   endelse

	xss = nx * zfact	;Make integer rebin factors
	yss = ny * zfact
	erase
	tv, rebin( area, xss,yss, /SAMPLE )
	area = 0
	wset, cur_w
     END

ENDCASE

goto,AGAIN

EXIT:
	if N_elements( x0 ) EQ 1 then  box_erase

	if N_elements( info_w ) EQ 1 then begin
		printw, [" "," "," "," "], LINE=-5, /ERASE, WIN=info_w
		wset, cur_w
	   endif
return
end
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


