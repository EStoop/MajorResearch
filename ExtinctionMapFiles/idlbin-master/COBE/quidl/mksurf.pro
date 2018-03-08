function mksurf,xin,yin,zin,xgrid,ygrid,ll=ll,ur=ur,dx=dx,dy=dy,badval=bad_data
;------------------------------------------------------------------------------------
; Function to generate a uniformly gridded surface from a random sampling of data.
;
; Written By:  BA Franz, ARC, 1/93
;
;
;------------------------------------------------------------------------------------
on_error,2
;
; Make inputs x,y,z into vectors
;
x = xin(*)
y = yin(*)
z = zin(*)

;
; Check inputs
;
sx = size(x)
sy = size(y)
sz = size(z)

status = ( (sx(3) eq sy(3)) and (sx(3) eq sz(3)) )
if (not status) then begin
    message,'Input vectors x,y,z are of incompatible size.',/continue
    return,-1
endif

status = ( (sx(2) ge 1) and (sx(2) le 5) and  $
           (sy(2) ge 1) and (sy(2) le 5) and  $
           (sz(2) ge 1) and (sz(2) le 5) )
if (not status) then begin
    message,'Input vectors x,y,z are of invalid type.',/continue
    return,-1
endif

;
; Set default grid parameters if not specified by user
;
if (n_elements(ll) ne 2) then ll = [min(x),min(y)]
if (n_elements(ur) ne 2) then ur = [max(x),max(y)]
if (n_elements(dx) eq 0) then dx = (max(x)-min(x))/50.0
if (n_elements(dy) eq 0) then dy = (max(y)-min(y))/50.0
if (n_elements(bad_data) eq 0) then bad_data = 0

;
; Make grid
;
xgrid = ll(0) + findgen(long( (ur(0)-ll(0))/dx + 1))*dx
ygrid = ll(1) + findgen(long( (ur(1)-ll(1))/dy + 1))*dy
triangulate,x,y,tr,b
grid = trigrid(x,y,z,tr,[dx,dy],[ll(0),ll(1),ur(0),ur(1)],/quintic,missing=bad_data)

return,grid
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


