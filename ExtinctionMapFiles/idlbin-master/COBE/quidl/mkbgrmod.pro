function mkbgrmod,map,dsrcmap,fit=fit,cut=cut,fsize=fsize,filter=filter, $
                  noplot=noplot,nocorn=nocorn,badval=bad_data,srcwidth=srcwidth
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    MKBGRMOD is a general 2-D background modeling tool.
;
;DESCRIPTION:                                   
;    IDL function to generate a 2-D background model from a rasterized
;    map.  The type of background model returned is determined by the
;    combination of input keyword parameters.  Some possibilities 
;    include:
;    
;    Model Type					Keywords/Values
;
;  o Quintic Interpolation To Pixels Specified	None
;    Interactively Using the Mouse.
;
;  o As Above But Including An nth-Order	Fit=n, where n>0.
;    Polynomial Surface Fit to the Interpolated
;    Background.
;
;  o Lower Envelope Determined Using An NxN	Fsize=N, where N>1.
;    Spatial Minimum Filter.
;
;  o As Above But Including An nth-Order	Fsize=N, where N>1, 
;    Polynomial Surface Fit to the Lower	Fit=n, where n>0.
;    Envelope.
;
;  o Median Filtered Using An NxN Spatial	Fsize=N, where N>1, 
;    Median Filter.				Filter='median'.
;
;  o Iterative Surface Fitting and Source	Cut=M, where M>0, 
;    Subtraction.  At each iteration, the       Fit=n, where n>0. 
;    routine first calculates an intermediate
;    background model using a surface fit.
;    The background model is then subtracted 
;    from the original map, and the standard 
;    deviation of the residuals is determined.
;    All pixels greater than M standard         
;    deviations above the background are 
;    located, a SRCWIDTH x SRCWIDTH patch 
;    centered on each deviant pixel is
;    sentinelized, and the next iteration 
;    begins.  Iteration stops when no more 
;    sources are found.
;
;  o As Above But Surface Fit is Done to the 	Cut=M, where M>0, 
;    NxN Lower Envelope.  			Fit=n, where n>0,
;						Fsize=N, where N>1.
;
;  o As Above But Using An Averaging Filter	Cut=M, where M>0, 
;    With No Surface Fitting. 		 	Fsize=N, where N>1,
;						Filter='avg'.    		  				
;
;CALLING SEQUENCE:
;    bgr = MKBGRMOD(map,dsrcmap,fit=fit,cut=cut,fsize=fsize,
;                   filter=filter,srcwidth=srcwidth,/nocorn,/noplot,
;                   badval=badval)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    map         I     arr         2-D array containing the input map
;    fit        [I]    scl         Order of polynomial surface fit 
;                                  (per dimension)
;    fsize      [I]    scl         Size of Spatial Filter in Pixels
;    filter     [I]    str         Name of filter function. Default=
;                                  'min' (i.e., the lower envelope), 
;                                  but any function which takes an
;			           array as input and returns a 
;                                  scalar is valid.
;    cut        [I]    scl         Number of Standard Deviations Used
;                                  to ID a Source
;    srcwidth   [I]    scl         Size of pixel patch to be removed 
;                                  when a source is found. (Default=3,
;                                  the central pixel and its nearest 
;                                  neighbors)
;    badval     [I]    scl         All data at this value will be 
;                                  ignored. (Default=0.0)
;    nocorn     [I]    scl         When selecting background points 
;                                  with the mouse, MKBGRMOD 
;                                  automatically includes the map 
;                                  corners. Specifying /nocorn will 
;                                  supress this.
;    noplot     [I]    scl         By default, MKBGRMOD will plot the
;                                  original map and the background 
;                                  model.  Specifying /NOPLOT will 
;                                  suppress the plotting.  The option 
;                                  is not valid for interactive 
;                                  modeling.     
;    bgr         O     arr         2-D Array containing the final 
;                                  background model.
;    dsrcmap    [O]    arr         If the CUT option is specified, 
;                                  DSRCMAP will contain the final 
;                                  source subtracted map.
;
;WARNINGS:
;    Plotting requires a TV terminal, and the input map must be 
;    smaller than 512x512. The iterative procedure may become 
;    unstable, particularly with high-order surface fits or low 
;    sigma cuts.  Specifying the interpolated model with /NOCORN and
;    a surface fit is not recommended.  The surface fit requires that
;    all pixels are sampled, but without the corners MKBGRMOD may  
;    have to extrapolate to fill the edges.
;
;EXAMPLES:
;    To specify background points with the mouse:
;      bgr = MKBGRMOD(map)
;      The input map will be displayed on th screen.  Press the left
;      mouse button to select points.  Select the last point with 
;      the right mouse button.  
;
;    To generate a 2nd-Order polynomial fit to the 4x4 lower envelope
;    while ignoring any data equal to -16000:
;      bgr = MKBGRMOD(map,fsize=4,fit=2,badval=-16000)
;
;    To generate an iteratively source subtracted model which uses a
;    4th-order polynomial surface fit to the 3x3 lower envelope and a
;    source cut of 5 sigma:
;      bgr = MKBGRMOD(map,dsrcmap,cut=5,fit=4,fsize=3)
;
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Interactive options use TRIANGULATE and TRIGRID to generate a 
;    uniformly gridded surface.  The lower surface envelope is 
;    generated using SPFILTER.pro with the MIN function.  Surface 
;    fits use the SURFACE_FIT procedure.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines Called:  RBUSTBGR, FILLSURF, BLANKPIX, SPFILTER, 
;                         MKSURF, UPKWHERE, LSPLINE. 
;
;MODIFICATION HISTORY
;    Written by B.A. Franz, Applied Research Corp.,   Feb 1992.
; SPR 12060  19-Jan-1995 Variable name change due to IDL3.6. J. Newmark
;
;------------------------------------------------------------------------------------
;
on_error,2
;
; Set defaults
;
if (n_elements(bad_data) eq 0) then bad_data=0.0
if (n_elements(cut)      eq 0) then cut=0
if (n_elements(fit)      eq 0) then fit=0
if (n_elements(fsize)    eq 0) then fsize=0
if (n_elements(filter)   eq 0) then filter='min'
if (n_elements(srcwidth) eq 0) then srcwidth=3
if (n_elements(noplot)   eq 0) then noplot=0
if (n_elements(nocorn)   eq 0) then nocorn=0
bgr = -1

if (cut      lt 0) then cut=0
if (fit      lt 0) then fit=0
if (fsize    lt 0) then fsize=0
if (srcwidth lt 1) then srcwidth=1
if (noplot   lt 0) then noplot=0
if (nocorn   lt 0) then nocorn=0

;
; Echo Input Parameters and Prompt for Acceptance
;
print,''
print,'----------------------------------------'
print,' Summary of MKBGRMOD Input Parameters'
print,''
print,'    cut      = ',cut
print,'    srcwidth = ',srcwidth
print,'    fit      = ',fit
print,'    fsize    = ',fsize
print,'    filter   = ',filter
print,'    badval   = ',bad_data
print,'    noplot   = ',noplot
print,'    nocorn   = ',nocorn
print,'----------------------------------------'
print,''
print,' Press <Return> to continue or <E> to exit.'
resp = ''
read,resp
if ((resp eq 'E') or (resp eq 'e')) then goto,endrun

;
; Get map size
;
s = size(map)
if (s(0) ne 2) then begin
    message,'MKBGRMOD requires a 2-D array.',/continue
    return,-1
endif else begin
    nx   = s(1)
    ny   = s(2)
endelse
;
; Display original map with pixels enlarged to fill 512x512 window
;
if (not keyword_set(noplot)) then begin
    ;
    if ((!d.name eq 'WIN') or (!d.name eq 'X')) then begin
        window,/free,xsize=550,ysize=600,title='MKBGRMOD: Original Map'
    endif
    ;
    ; Get rebin size
    ;
    f = min([fix(512/nx),fix(512/ny)])
    if (f eq 0) then begin
        message,'Maximum array dimension must be < 512 for plotting.',/continue
        return,-1
    endif
    xfact = f*nx
    yfact = f*ny
    ;
    ; Get scaling
    ;
    select = where(map ne bad_data)
    if (select(0) eq -1) then begin
        message,'No valid data points in original array.',/continue
        print,' *** BADVAL=',bad_data
        return,-1
    endif
    mx = max(map(select))
    mn = min(map(select))
    mx = mx + 0.05*(mx-mn)
    mn = mn - 0.05*(mx-mn)

    ;
    ; Display Original Map
    ;
    vimage = rebin(map,xfact,yfact,/sample)
    tv,bytscl(vimage,mn,mx),10,10
    bar
    xyouts, 40,560,strtrim(string(mn),2),color=240,/dev
    xyouts,440,560,strtrim(string(mx),2),color=240,/dev

endif

;
; Choose type of modeling desired
;
if ((cut le 0) and (fsize le 0)) then begin
    ;
    ; Interactive background modeling
    ;
    if (keyword_set(noplot)) then begin
        message,'Can not specify /noplot with interactive modeling.',/continue
        return,-1
    endif
    ;
    ; Get x,y Grid Maps
    ;
    xmap = lindgen(nx, ny) mod nx
    ymap = lindgen(nx, ny)  /  nx	
    ;
    ; Rebin x,y Grid Maps to Match Image
    ;
    ximage = rebin(xmap,xfact,yfact,/sample)
    yimage = rebin(ymap,xfact,yfact,/sample)

    ;
    ; Use corner points if user does not object
    ;
    if (not keyword_set(nocorn)) then begin
        xbgr = [0,0,   nx-1,nx-1]
        ybgr = [0,ny-1,ny-1,0   ]
        zbgr = map(xbgr,ybgr)
    endif    

    ;
    ; Begin point-and-click background specification
    ;
    print,''
    print,'Select background points with LEFT MOUSE BUTTON.'
    print,'Select last point with RIGHT MOUSE BUTTON to exit.'
    print,''
    ;
    while (!err ne 4) do begin  ; While right mouse button not pressed
      cursor,x,y,/dev,/down     ;   Get cursor position
      i = x-10                  ;   Subtract image offset
      j = y-10
      ;
      ; If position is within image boundaries, add it to background points list
      ;
      if ( (i ge 0) and (i lt xfact) and (j ge 0) and (j lt yfact) ) then begin
          if (n_elements(xbgr) eq 0) then begin
              ;
              ; First Point
              ;
              xbgr = ximage(i,j)
              ybgr = yimage(i,j)
              zbgr =  vimage(i,j)
          endif else begin
              xbgr = [xbgr,ximage(i,j)]
              ybgr = [ybgr,yimage(i,j)]
              zbgr = [zbgr, vimage(i,j)]         
          endelse
          ;
          ; Mark image with selected point
          ;
          plots,[x,x],[y,y],psym=1,color=240,/dev
      endif

    endwhile ; End point-and-click mode
    ;
    ; Eliminate Bad Data Values
    ;
    select = where(zbgr ne bad_data)
    if (select(0) eq -1) then begin
        message,'All selected points are BADVALs.',/continue
        print,' *** BADVAL=',bad_data
        return,-1
    endif
    xbgr = xbgr(select)
    ybgr = ybgr(select)
    zbgr = zbgr(select)
    ;
    ; Interpolate background points to fill (if possible) the original map grid
    ;
    if (n_elements(xbgr) gt 2) then begin
        message,'Interpolating Surface.',/continue
        bgr = mksurf(xbgr,ybgr,zbgr,ll=[0,0],ur=[nx-1,ny-1],dx=1,dy=1,badval=bad_data)
    endif else begin
        message,'Insufficient Number of Points for Interpolation.',/continue
        return,-1
    endelse
    ;
    ; Surface fit to smooth background if requested
    ;
    if (fit gt 0) then begin
        message,'Fitting Surface.',/continue
        bgr = fillsurf(bgr,bad_data)  ; Interpolate to fill sentinels
        bgr = surface_fit(bgr,fit)    ; Fit surface
    endif

endif else if ((cut le 0) and (fsize gt 0)) then begin
    ;
    ; Lower envelope fit
    ;
    bgr = spfilter(map,fsize,filter,badval=bad_data)
    ;
    ; Surface fit to smooth background if requested
    ;
    if (fit gt 0) then begin
        message,'Fitting Surface.',/continue
        bgr = fillsurf(bgr,bad_data)  ; Interpolate to fill sentinels
        bgr = surface_fit(bgr,fit)    ; Fit surface
    endif

endif else begin
    ;
    ; Iterative source subtracted background fit
    ;
    message,'Calculating Background, Please Be Patient.',/continue
    bgr = rbustbgr(map,cut,dsrcmap,fit=fit,fsize=fsize,filter=filter, $
                   badval=bad_data,srcwidth=srcwidth)

endelse

;
; Display background model
;
if (not keyword_set(noplot)) then begin
    ;
    if ((!d.name eq 'WIN') or (!d.name eq 'X')) then begin
        window,/free,xsize=550,ysize=600,title='MKBGRMOD: Background Model'
    endif
    ;
    ; Get scaling
    ;
    select = where(bgr ne bad_data)
    if (select(0) eq -1) then begin
        message,'All BGR points are BADVALs.',/continue
        print,' *** BADVAL=',bad_data
        return,-1
    endif
    mx = max(bgr(select))
    mn = min(bgr(select))
    mx = mx + 0.05*(mx-mn)
    mn = mn - 0.05*(mx-mn)

    ;
    tv,bytscl(rebin(bgr,xfact,yfact),mn,mx),10,10
    bar
    xyouts, 40,560,strtrim(string(mn),2),color=240,/dev
    xyouts,440,560,strtrim(string(mx),2),color=240,/dev

endif

endrun:

return,bgr
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


