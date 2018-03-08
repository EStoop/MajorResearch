;-----------------------------------------------------------------------------
; Function FILLSURF
;
; IDL function to fill bad data values in a rasterized map with quintic interpolated
; values, linearly extrapolated values, or (as a last resort) linear extrapolations
; of extrapolated values (user beware).
;
; Written By: BA Franz, ARC, 2/93
;
;-----------------------------------------------------------------------------
;
function fillsurf,surface,bad_data

on_error,2

s = size(surface)
nx = s(1)
ny = s(2)
coord_x = indgen(nx,ny) mod nx	;X coords.
coord_y = indgen(nx,ny) / nx	; & y.
surf = surface

select = where(surface eq bad_data)
if (select(0) ne -1) then begin
    ;
    ; Quintic interpolate over bad_data using triangulated surface
    ;
    ss = surface(*)
    xx = coord_x(*)
    yy = coord_y(*)
    select = where(ss ne bad_data)
    if (select(0) eq -1) then begin
        message,'Not enough valid points to fill surface.',/continue
    endif
    triangulate,xx(select),yy(select),tr,b
    surf = trigrid(xx(select),yy(select),ss(select),tr,/quintic, $
               [1,1],[0,0,nx-1,ny-1],missing=bad_data)

    ;
    ; Linearly interpolate to fill-in points outside of triangulation.
    ; Note: quintic extrapolation option of trigrid is too unstable.
    ;
    select = where(surf eq bad_data) 
    if (select(0) ne -1) then begin
        xspline = findgen(nx)
        yspline = findgen(ny)
        xsurf   = surf
        ysurf   = surf
        ;
        ; Spline interpolate along constant y
        ;
        for j=0,ny-1 do begin
            select = where(surf(*,j) ne bad_data)
            if (n_elements(select) gt 3) then $
                xsurf(*,j) = lspline(xspline(select),surf(select,j),xspline)
        endfor
        ;
        ; Spline interpolate along constant x
        ;
        for i=0,nx-1 do begin
            select = where(surf(i,*) ne bad_data)
            if (n_elements(select) gt 3) then $
                ysurf(i,*) = lspline(yspline(select),ysurf(i,select),yspline)
        endfor
        ;
        ; Average xspline and yspline maps 
        ;
        select = where((xsurf ne bad_data) and (ysurf eq bad_data))
        if (select(0) ne -1) then surf(select) = xsurf(select)
        select = where((xsurf eq bad_data) and (ysurf ne bad_data))
        if (select(0) ne -1) then surf(select) = ysurf(select)
        select = where((xsurf ne bad_data) and (ysurf ne bad_data))
        if (select(0) ne -1) then surf(select) = (xsurf(select)+ysurf(select))/2.0

    endif

    ;
    ; If original surface does not touch two edges, there may still be holes in corners.
    ; So, extrapolate the extrapolated values.
    ;
    select = where(surf eq bad_data) 
    if (select(0) ne -1) then begin
        xspline = findgen(nx)
        yspline = findgen(ny)
        xsurf   = surf
        ysurf   = surf
        ;
        ; Spline interpolate along constant y
        ;
        for j=0,ny-1 do begin
            select = where(surf(*,j) ne bad_data)
            if (n_elements(select) gt 3) then $
                xsurf(*,j) = lspline(xspline(select),surf(select,j),xspline)
        endfor
        ;
        ; Spline interpolate along constant x
        ;
        for i=0,nx-1 do begin
            select = where(surf(i,*) ne bad_data)
            if (n_elements(select) gt 3) then $
                ysurf(i,*) = lspline(yspline(select),ysurf(i,select),yspline)
        endfor
        ;
        ; Average xspline and yspline maps 
        ;
        select = where((xsurf ne bad_data) and (ysurf eq bad_data))
        if (select(0) ne -1) then surf(select) = xsurf(select)
        select = where((xsurf eq bad_data) and (ysurf ne bad_data))
        if (select(0) ne -1) then surf(select) = ysurf(select)
        select = where((xsurf ne bad_data) and (ysurf ne bad_data))
        if (select(0) ne -1) then surf(select) = (xsurf(select)+ysurf(select))/2.0

    endif

endif

return,surf
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


