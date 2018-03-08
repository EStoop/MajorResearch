;--------------------------------------------------------------------
; Function RBUSTBGR
;
; IDL function to calculate a background model from a rasterized map.  The
; background model is an nth-order polynomial surface fit to a "desourced" version
; of the input map.  The map is descourced by first identifying all pixels which are
; greater than CUT standard deviations above the background model, and then setting 
; those pixels and all neighboring pixels to sentinels.  The procedure iterates until 
; no sources are found.
;
; Written By: BA Franz, ARC, 1/93
;
; Inputs:
;     inmap - 2-dimensional array containing the input map
;     cut   - sigma limit used to identify sources
;     fit   - if set to zero, a simple mean background model will be 
;             used.  If set to n > 0, an nth order polynomial surface
;             fit will be used.
;    fsize  - fsize > 1, then the SPFILTER function will be applied 
;             first, and any surface fit will be done to the filtered data.
;    filter - string containing the name of the filter function (Default
;             is the IDL MIN function. See SPFILTER.pro)  
;    srcwidth - number of pixels to be blanked when a source is found.  A
;             srcwidth x srcwidth pixel region centered on the source peak
;             will be blanked. (Default=3)
;    maxiter - maximum number of iterations allowed (Default=100).
;    badval  - data at this value will be ignored (default=0.0).
;
; Outputs:
;     bgr    - the final background model
;     map    - the final desourced map (returned as the function)
;
;--------------------------------------------------------------------

function rbustbgr,inmap,cut,map,fit=fit,fsize=fsize,filter=filter,badval=bad_data, $
                  maxiter=maxiter,srcwidth=srcwidth
;
on_error,2
;

if (n_elements(fit)      eq 0) then fit=0
if (n_elements(fsize)    eq 0) then fsize=0
if (n_elements(filter)   eq 0) then filter='min'
if (n_elements(bad_data) eq 0) then bad_data=0.
if (n_elements(maxiter)  eq 0) then maxiter=100
if (n_elements(srcwidth) eq 0) then srcwidth=3

map  = inmap
bgr0 = inmap*0.0 
s    = size(map)
nx   = s(1)
ny   = s(2)

still_sources = 1   ; Source detection flag
cnt           = 0   ; Iteration count
nsources      = 0   ; Number of sources found
;
; While there are still sources detected: calculate the background 
; model, subtract background from map, calculate stdev of subtracted 
; map, search through map for any pixels greater than 5 stdevs and blank
; them and there nearest neighbors.
;
while (still_sources) and (cnt le maxiter) do begin

    cnt = cnt+1
    still_sources = 0

    ;
    ; Select all non-sentinel data
    ;
    select = where(map ne bad_data)
    if (select(0) ne -1) then begin
        ;
        ; Calculate background model
        ;
        if (fit gt 0) then begin
            ;
            ; Surface fit models
            ;
            if (fsize gt 1) then begin
                ;
                ; Form filtered map and interpolate to fill sentinels
                ;
                bgr = spfilter(map,fsize,filter,badval=bad_data)
                bgr = fillsurf(bgr,bad_data) 
            endif else begin
                ;
                ; Interpolate to fill sentinels
                ;
                bgr = fillsurf(map,bad_data)
            endelse
            ;
            ; Calculate polynomial surface fit
            ;
            bgr = surface_fit(bgr,fit)

        endif else begin                       
            if (fsize gt 1) then begin
                ;
                ; Form filtered map and interpolate to fill sentinels
                ;
                bgr = spfilter(map,fsize,filter,badval=bad_data)
                bgr = fillsurf(bgr,bad_data) 
            endif else begin
                ;
                ; Simple average model
                ;
                bgr = bgr0 + avg(map(select))
            endelse

        endelse

        ;
        ; Calculate background subtracted map, retain existing sentinels
        ;
        submap = map
        submap(select) = map(select) - bgr(select)
        ;
        ; Calculate standard deviation and source strength
        ;
        sigma  = stdev(submap(select))
        source = sigma*cut

        ;
        ; Locate and blank all sources based on current iteration 
        ; background model.
        ;
        still_peak = 1
        while (still_peak) do begin
            mpeak = max(submap)
            if (mpeak ge source) then begin
                ;
                ; A source was found
                ;
                still_sources = 1
                nsources = nsources+1
                ipeak = where(submap eq mpeak)
                upkwhere,submap,ipeak,i,j  &  i=i(0)  &  j=j(0)
                blankpix,submap,i,j,srcwidth,bad_data
                blankpix,map,i,j,srcwidth,bad_data
            endif else begin
                still_peak = 0
            endelse

        endwhile
 
    endif  ; Map contains non-sentinel data

endwhile  ; Still sources found

if (cnt gt maxiter) then message,'Iteration Limit Reached.'/continue

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


