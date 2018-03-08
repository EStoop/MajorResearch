;------------------------------------------------------------------------
; Procedure LSPLINE
;
; IDL Procedure to spline-fit a set of points using linear interpolation
; and/or extrapolation.
;
; Programmer: Bryan A. Franz, ARC
; Date Written: 2/93
;
; Inputs:
;   xin   - input x-vector (sorted)
;   yin   - input y-vector
;   xout  - output x-vector (sorted)
;
; Outputs:
;   yout  - output y-array
;
;------------------------------------------------------------------------
function lspline,xxin,yyin,xout
;
nin = n_elements(xxin)
xin = 1.*xxin
yin = 1.*yyin

nout = n_elements(xout)
yout = fltarr(nout)
;
; Calculate slopes (forward difference)
;
i   = lindgen(nin-1)
ip1 = i+1
m   = 1.*(yin(ip1)-yin(i))/(xin(ip1) - xin(i))
m   = [m,m(nin-2)]
;
; Forward extrapolation
;
nforw = 0
select = where(xout ge xin(nin-1))
if (select(0) ne -1) then begin
    nforw = n_elements(select)
    yout(select) = yin(nin-1) + m(nin-1)*(xout(select)-xin(nin-1))
endif
;
; Backward extrapolation
;
nback = 0
select = where(xout le xin(0))
if (select(0) ne -1) then begin
    nback = n_elements(select)
    yout(select) = yin(0) + m(0)*(xout(select)-xin(0))
endif
;
; Interpolation
;
j = long(nback)
i = 0L
while (j lt nout-nforw) do begin
    while (xin(i) le xout(j) ) do i=i+1
    i = i-1
    yout(j) = yin(i) + m(i)*(xout(j)-xin(i))
    j = j+1
endwhile     

return,yout
end

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


