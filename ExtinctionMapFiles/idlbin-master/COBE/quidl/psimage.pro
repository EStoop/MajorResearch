pro psimage,vimage,min,max,output=output,scale=scale,color=color,$
    rcolor=rcolor,gcolor=gcolor,bcolor=bcolor,hsv=hsv,hls=hls,$
    seiko=seiko
;
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:    
;     PSIMAGE sets optimum size and displays (tvscl) image in postscript
;
;DESCRIPTION: 
;     Sets up an optimised display size for post_script output, 
;     proportional to the ratio of the x/y axes of the image; then 
;     displays the input array in this window, using the IDL display 
;     procedure 'tvscl'. The 'optimised' window selects the largest 
;     possible display size for letter size post-script printout, 
;     keeping the x/y ratio of the image the same as the original. The 
;     image is rotated if y is larger than x, before display. 'tvscl' 
;     which is used to display the image,defaults to the maxima and 
;     the minima of the image, unless otherwise specified by the user. 
;     Optional output file name can be put in if desired filename is 
;     other than 'idl.ps'. An optional scaling parameter can be used 
;     if the output should be less than full letter size. The output 
;     file can be sent directly to a postscript printer.
;
;CALLING SEQUENCE: 
;     PSIMAGE, vimage [,min] [,max] [,OUTPUT=output] [,SCALE=scale] 
;              [,COLOR=color] [,/HSV] [,/HLS] [,/SEIKO]
;
;ARGUMENTS (I = input, O = output, [] = optional): 
;     VIMAGE    I   flt [arr]  Name of input array to be displayed.    
;     MIN      [I]  flt        User specified minima for tvscl 
;                              display; Default is minima of array.  
;     MAX      [I]  flt        User specified maxima for tvscl 
;                              display; Default is maxima of array.   
;     OUTPUT   [O]  str        If specified, 'output' is the name for  
;                              postscript format output. Default is 
;                              'idl.ps'. The postscript format output
;                              file can be sent to a postscript 
;                              printer directly.
;     SCALE    [I]  flt        If specified (has to be less than 1.0), 
;                              the image covers less than the full size 
;                              letter page. Default is SCALE = 1.0.
;     COLOR    [I]  int        Number of IDL color table to use for
;                              color output. PS file must be sent to
;                              a color printer. Set to -1 to use pre-
;                              defined PS color table.
;     R(G,B)COLOR [I]  intarr  RGB colors to be used for defining output
;                              color table in TVLCT. PS file must be sent to 
;                              a color printer. Not compatible with COLOR
;                              keyword.
;     HSV      [I]  int        Set this keyword if using HSV color system
;                              for input to TVLCT.
;     HLS      [I]  int        Set this keyword if using HLS color system
;                              for input to TVLCT.
;     SEIKO   [I}  int        Set this keyword for proper margins if using
;                              a SEIKO color printer.
;
;WARNINGS:
;
;     1.   Scaling Factor has to be less than 1.0, as this routine is
;          set for maximum display size. 
;     2.   The default procedure with no optional keywords will produce
;          a black and white plot.
;     3.   The procedure sets the device type back to the previous one
;          if the device type was other than post-script, ex. 'X', 
;          'tek' or 'regis'. But if the previous device was also 'ps',
;          then the changed x/y sizes and offsets are retained. The 
;          user should remember to reset 'ps' device parameters to 
;          those desired, unless psimage is used again to optimize
;          the display scale.          
;     4.   The output file is created in the default or current 
;          directory. As postscript files may be quite large, the 
;          user should make sure there is enough available space.
;     5.   The postscript output uses 8 bits_per_pixel resolution.
;     6.   The default size of the output has been optimized for a
;          TEK2SD color printer. A keyword is available for optimization
;          on a Seiko color printer.
;
;EXAMPLE: 
;     1. To produce a balck and white display for image image1, letting 
;        the default maxima and minima of the array be the display scale, 
;        and let the default 'idl.ps' file be the output :
;
;        psimage, image1			     
;
;     2. To produce a  display for image array image1, scaled to 0.8 
;        times the maximum possible letter-size output display, using 
;        0 and 20 as minima and maxima of the display scale, and 
;        output the post-script file to a file called 'out.ps':
;     
;        psimage,image1,0,20, output = 'out.ps', scale= 0.8
;
;     3. To produce a display for image image1, letting the default 
;        maxima and minima of the array be the display scale, using
;        IDL color table 3 (Red Temperature), and
;        let the default 'idl.ps' file be the output :
;
;        psimage, image1, color=3
;
;     4. To produce a display for image image1, letting the default
;        maxima and minima of the array be the display scale, using
;        user defined colors (hsv) in the HSV color system, and send 
;        output to the file image1.ps:
;
;        psimage, image1, rcolor=h, gcolor=s, bcolor=v, /hsv
;
;#
;COMMON BLOCKS:
;     Colors
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES)
;     Uses IDL intrinsic routine 'tvscl' to display the input array in
;     a post-script format file. The post-script file is created in 
;     the default directory.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     none
;     
;MODIFICATION HISTORY:
;     Written by Mila Mitra, General Sciences Corp. 17  March 1992
;     SPR 11250  26-Aug-1993  Enable color output. J. Newmark
;     SPR 11906  07-Sep-1994  Change scaling from SEIKO to TEK2SD. J Newmark
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;
;.TITLE
;Routine PSIMAGE
;-
COMMON colors,r_orig,g_orig,b_orig,r_curr,g_curr,b_curr
s=size(vimage)
old_dev=!d.name
set_plot,'ps'
device,/portrait
;SET UP OUTPUT FILE, OPTIONAL COLORS AND OPTIONAL MIN, MAX AND SCALING 
; PARAMETERS:
;
if keyword_set(output) eq 0 then output='idl.ps'
IF (keyword_set(rcolor) and keyword_set(gcolor) and keyword_set(bcolor)) $
 THEN BEGIN
  device,/color,filename=output,bits_per_pixel=8 
  tvlct,rcolor,gcolor,bcolor
  IF keyword_set(hsv) then tvlct,rcolor,gcolor,bcolor,/hsv
  IF keyword_set(hls) then tvlct,rcolor,gcolor,bcolor,/hls
ENDIF ELSE IF (keyword_set(color)) THEN BEGIN
  device,/color,filename=output,bits_per_pixel=8 
  IF (color NE -1) THEN loadct,color
ENDIF ELSE device,filename=output,color=0,bits_per_pixel=8
if n_elements(scale) eq 0 then scale=1
if n_elements(min) eq 0 then min=min(vimage)
if n_elements(max) eq 0 then max=max(vimage)
if keyword_set(seiko) then begin
  rlim=0.74
  ys=0.75
  mn=9.5
endif else begin
  rlim=7/8.5
  ys=1.25
  mn=8.5
endelse
;
;IF YSIZE OF IMAGE IS GREATER THAN X SIZE, ROTATES IMAGE
;
r=float(s(1))/float(s(2))
if (r gt 1.0) then begin
imrot=rotate(vimage,1)
end
;
;CASE 1: iF X LESS THAN Y, AND X/Y < 0.74 (or 7/8.5 for tek2sd):
;
case 1 of
r le rlim : begin
n=mn
m=mn*r
device,/inches,xoff=0.75,yoff=ys,scale=scale,xsize=m,ysize=n
tvscl,vimage <max >min
end
;
;CASE 2: IF X LESS THAN Y, AND 0.74(or 7/8.5 for tek2sd) < (X/Y) < 1.0 :
;
(r gt rlim) and (r le 1.0):  begin
m=7.0
n=m/r
device,/inches,xoff=0.75,yoff=ys,scale=scale,xsize=m,ysize=n
tvscl,vimage <max >min
end
;
;CASE 3: IF X LESS THAN Y, AND 1.0 < (X/Y) < 1.35(or 8.5/7 for tek2sd), 
;         ROTATES IMAGE :
;
(r gt 1.0) and (r le 1./rlim): begin
m=7.0
n=m*r
device,/inches,xoff=0.75,yoff=ys,scale=scale,xsize=m,ysize=n
tvscl,imrot <max >min
end
;
;CASE 4: IF X LESS THAN Y, AND (X/Y) > 1.35(or 8.5/7 for tek2sd), 
;         ROTATES IMAGE :
;
(r gt 1./rlim) : begin
n=mn
m=n/r
device,/inches,xoff=0.75,yoff=ys,scale=scale,xsize=m,ysize=n
tvscl,imrot <max >min
end
;
else:print,'ratio of x/y out of range'
endcase
;
device,/close
;
;SETS DEVICE BACK TO PREVIOUS DEVICE
;
set_plot,old_dev
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


