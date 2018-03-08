PRO BAR,MIN=MIN,MAX=MAX,TYPE=TYPE,TITLE=TITLE,REV=REV,THICK=THICK
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     BAR draws a color bar on the display.  
;
;CALLING SEQUENCE:  
;     BAR [,MIN = Min] [,MAX = Max] [,TYPE = Type] $
;         [,TITLE = Title]  [,THICK=Thick] [,/REV]
;
;ARGUMENTS: (I = input, O = output, [] = optional)
;     Min    [I]   byt   Minimum color number in bar; default is 0
;     Max    [I]   byt   Maximum color number in bar; default is 255
;     Type   [I]   chr   Positioning indicator.
;                         TYPE='XT': X bar and label at top (default).
;                         TYPE='XB': X bar and label at bottom.
;                         TYPE='YR': Y bar and label at right.
;                         TYPE='YL': Y bar and label at left.
;     Title  [I]   chr   Bar title.
;     Rev    [I]   key   If specified, bar is reversed.
;
;WARNINGS:
;     Works only with devices that allow TV operations.
;     Only values between 0 and 255 are allowed for Min and Max;
;        Min must be less than Max if both are specified.
;
;EXAMPLES:
;     1.  Put a color bar from right to left at top of picture.
;                BAR
;     2.  Put a reversed color bar of width 20 pixels of colors 0 
;            through 150 at top of picture.
;                Bar,/rev, max = 150
;     3.  Put a reversed color bar at bottom of picture with title
;            'REVERSED'
;                Bar, /rev, type = 'xb', title='REVERSED'
;     4.  Put a color bar to left of picture 30 pixels wide with 
;             title 'STRAIGHT'
;                Bar, type = 'yl', title='STRAIGHT', thick=30
;     5.  Put a color bar to right of picture showing only colors 
;             50 through 200
;                Bar, type = 'yr', min = 50, max = 200
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;     Uses system variables to determine pixel size of window; scales
;     size of title, if any, to correspond to thickness of bar.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     NONE
;
;MODIFICATION HISTORY:
;   R. Sterner. 12 NOV, 1987.
;      Johns Hopkins University Applied Physics Laboratory.
;   RES  27 Jan, 1988 --- upgraded to Y bar and allowed labels.
;   Upgraded to Y bar; allowed labels; used keywords -- IDL Version II
;           Alice R. Trenholme, General Sciences Corporation 
;.TITLE
;Routine BAR
;-
if !d.name eq 'TEK' or !d.name eq 'REGIS' then begin
        print, '    BAR will work only with X-windows or PostScript displays.'
endif else begin

    Bar_types = ['XT','XB','YL','YR']
    Bar_index = 0

    if n_elements(type) gt 0 then $
       for i=0,3 do if strupcase(type) eq bar_types(i) then bar_index = i

    if bar_index lt 2 then xlen = !d.x_vsize * .8 else $
          xlen = !d.y_vsize * .8

    IF N_ELEMENTS(thick) eq 1 THEN ylen = thick else ylen = 10

;      !d.x_vsize, !d.y_vsize = size of visible display area
    x_vis = !d.x_vsize
    y_vis = !d.y_vsize

;         Positions for the top, bottom, left, right bars
    x = [x_vis * .1,x_vis * .1,0,x_vis-ylen]
    y = [y_vis - ylen,0,y_vis * .1,y_vis * .1]
    
    if keyword_set(rev) then rotation_ix=[1,1,0,0] $
                         else rotation_ix=[4,4,2,2]
    bar_mat = bytscl(findgen(ylen,xlen))

    if ((n_elements(min) gt 0) or (n_elements(max) gt 0)) then begin
       if n_elements(min) eq 0 then min = 0
       if n_elements(max) eq 0 then max = 255
       if ((max lt min) or (((min lt 0) or (min gt 255)) or $
           ((max lt 0) or (max gt 255)))) then begin
                print,'  Disallowed values of minimum and/or maximum'
                return
       endif
       bar_mat = byte(.5 + min + (float(max-min)/255.)*bar_mat)
    endif

    tv,rotate(bar_mat,rotation_ix(bar_index)),x(bar_index),y(bar_index)

    if n_elements(title) ne 0 then begin
;         Positions for the top, bottom, left, right title
        xt = [x_vis*.5,x_vis*.5,2*ylen,x_vis-ylen-3]
        yt = [y_vis-2*ylen,ylen+3,y_vis*.5,y_vis*.5]
        ort = [0,0,90,90]

        xyouts, xt(bar_index),yt(bar_index),title,alignment=.5,  $
                orientation=ort(bar_index),/device,size = ylen/10.
    endif

endelse
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


