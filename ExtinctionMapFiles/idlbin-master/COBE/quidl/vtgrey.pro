PRO vtgrey,vimage,title=title,tty=tty,termx=termx,termy=termy,$
     position=position,norescale=norescale
;+NAME/ONE LINE DESCRIPTION:
;    VTGREY displays grey-scale images on graphics terminals (e.g. VT240).
;
;DESCRIPTION:
;    This routine will put a grey-scale image on a graphics screen.
;    It is an alternative to contour plots for users of 
;    monochrome terminals.
;
;CALLING SEQUENCE:
;    vtgrey,vimage,[title=...],[tty=....],[termx=...],[termy=...],$
;           [position=...],[norescale=...]
;
;ARGUMENTS (I=input, o=output, []=optional): 
;    vimage     I   2-D arr  flt/int    Input image.
;    title     [I]  keyword  str        Title string.
;    tty       [I]  keyword  str        string specifying the terminal
;			                type, and thus setting termx 
;                                       and termy. Valid types are 
;                                       (case-insensitive): 
;				        'VT240' - Regis graphics 
;                                          (default)
;				        'Microterm' - Regis graphics 
;                                          (same as VT240)
;				        'GraphOn' - TEK graphics
;				        'PC' - NCSA Telnet 
;                                          TEK Grahpics
;				        'Mac' - VersaTerm PRO - TEK 
;                                          on Mac Plus, SE, etc.
;				        'LC' - VersaTerm TEK Graphics
;                                          on Mac LC
;    termx     [I]  keyword  int        Terminal X-resolution. 
;                                          (overrides tty)
;    termy     [I]  keyword  int        Terminal Y-resolution. 
;                                          (overrides tty)
;    position  [I]  1-D arr  int        A 1-D array with four elements
;                                       which tell where the plot
;                                       should be placed in an X
;                                       window.
;    norescale [I]  keyword             setting this keyword inhibits 
;                                       the default rescaling of IMAGE
;				        to fill as much of the screen  
;				        as possible.
;
;WARNINGS:
;     1. Can be very slow for images where most pixels are bright.
;     2. termx and termy may need to be set to produce acceptable 
;        output on non-standard terminals.
;
;EXAMPLES:
;     For a VT240 or Microterm in Regis graphics mode use the default.
;     IDL> vtgrey,band1a       
;     will plot a greyscale map of band1a on your terminal
;
;     For a Graphon in TEK graphics mode,
;     IDL> vtgrey,band1a,term='graphon'       
;     will plot a greyscale map of band1a on your terminal
;
;     For weird terminals or emulators the default or one of the 
;     defined terminal types may work well. If not, (meaning you get 
;     an awful plaid pattern for vtgrey,findgen(200,200) ) then you 
;     can play with termx and termy to find the right numbers for your
;     terminal. E.g.
;     IDL> greyterm,band1a,211,198  
;     might be what you need to get the proper spacing of pixels in 
;     the greyscale plot.
;#
;COMMON BLOCKS: None.
;
;PROCEDURE:
;    Divide the range of pixel values into 6 categories.  Plot a
;    symbol for each pixel, according to which category it belongs to.
;
;REVISION HISTORY:
;    Written by Rick Arendt.
;    TITLE keyword added by S.R.K. VIDYA SAGAR (ARC) MAR 92.
;    TERMX, TERMY , & POSITION keywords added by J.Ewing, May 92.
;    TTY options and auto rescaling added by Rick Arendt (ARC) Nov 1992
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;
;.TITLE
;Routine VTGREY
;-
  ON_ERROR, 1
  IF(NOT KEYWORD_SET(title)) THEN title = ''
;
; First check terminal type to find out how much room there is...
;
if (n_elements(termy) eq 0) then begin
  if (n_elements(tty) ne 0) then tty = strlowcase(tty) else tty='else'
  case tty of 
   'vt240'	: begin			; VT240 in REGIS mode
                    termx = 319
                    termy = 186
                  end
   'microterm'	: begin			; Mictoterm in REGIS mode
                    termx = 319
                    termy = 186
                  end
   'graphon'	: begin			; GraphOn in TEK mode
                    termx = 429
                    termy = 326
                  end
   'mac'	: begin			; Mac Plus, SE, or Classic running
                    termx = 354		; Versaterm Pro in TEK mode
                    termy = 261
                  end
   'lc'		: begin			; Mac LC (12" monitor) running
                    termx = 321		; Versaterm in TEK mode
                    termy = 240
                  end
   'pc'		: begin			; DIRBE PC running NCSA Telnet
                    termx = 429		; in TEK mode
                    termy = 292
                  end
   else		: begin
                    termx = 319
                    termy = 186
                  end
  endcase
endif
;
; If you want to rescale the image...
;
if not(keyword_set(norescale)) then begin
  imsiz = float(size(vimage))
  mag = (termx/imsiz(1)) < (termy/imsiz(2))
  im = congrid(vimage,imsiz(1)*mag,imsiz(2)*mag)
endif else begin
  mag = 1.
  im = vimage
endelse
;
; Now set up for the plotting...
;
  !x.style = 1 & !y.style = 1
  !x.range = [0,termx] / mag
  !y.range = [0,termy] / mag
  !p.psym = 3
  imsiz = SIZE(im)
  nx = imsiz(1)
  ny = imsiz(2)
  a = nx-(nx mod 2)-1
  b = ny-(ny mod 2)-1
;
; FLAG will indicate the grey level to use for each 2x2 pixel patch.
; The rebinning is necessary in order to display distinct grey levels.
;
  flag = REBIN(im(0:a,0:b),nx/2,ny/2)
  flag = flag-MIN(flag)
  flag = fix(5.0*flag/MAX(flag))
  sz = SIZE(position)
  IF(sz(0)+sz(1) EQ 0) THEN BEGIN
    PLOT,[0,termy],[0,termy],/nodata,title=title
  ENDIF ELSE BEGIN
    PLOT,[0,termy],[0,termy],/nodata,position=position,/dev,$
      xstyle=5,ystyle=5
  ENDELSE
  FOR j= 0,ny/2-1 DO BEGIN
    FOR i= 0,nx/2-1 DO BEGIN
      CASE flag(i,j) OF			
	0:	
	1: PLOTS,[i*2]/mag,[j*2]/mag,psy=3
	2: PLOTS,[i*2+1,i*2]/mag,[j*2+1,j*2]/mag
	3: PLOTS,[i*2,i*2,i*2+1]/mag,[j*2+1,j*2,j*2]/mag
	4: PLOTS,[i*2,i*2,i*2+1,i*2+1]/mag,[j*2,j*2+1,j*2+1,j*2]/mag
	5: PLOTS,[i*2,i*2,i*2+1,i*2+1,i*2]/mag,[j*2,j*2+1,j*2+1,j*2,j*2]/mag
	ELSE:
      ENDCASE
    ENDFOR
  ENDFOR
;
; cleanup and quit
;
  !x.style = 0
  !y.style = 0
  !x.range = [0,0]
  !y.range = [0,0]
  !p.psym = 0
  RETURN
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


