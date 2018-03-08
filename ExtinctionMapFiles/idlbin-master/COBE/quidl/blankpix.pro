;---------------------------------------------------------------------
; IDL function to replace one or more pixels in a rasterized map with
; a sentinel value.  No input checking is done.
;
; Written By: BA Franz, ARC, 11/92
;
; Inputs:
;    map - the rasterized map
;    i,j - the indices of the central pixel to blank
;    w   - the width of the field to blank (wxw centered on i,j)
;    bad_data - blanking value (Default = 0)
;
;-----------------------------------------------------------------------
pro blankpix,map,i,j,w,bad_data

on_error,2

if (n_elements(bad_data) eq 0) then bad_data = 0

s = size(map)
nx = s(1)
ny = s(2)

;
; Find indices of region to blank...avoid edges
;
delta  = fix(w/2)

;
; Select range in x for region
;
xrange = [0,w-1] - delta + i
xrange(0) = max([0,   xrange(0)])
xrange(1) = min([nx-1,xrange(1)])

;
; Select range in y for region
;
yrange = [0,w-1] - delta + j
yrange(0) = max([0,   yrange(0)])
yrange(1) = min([ny-1,yrange(1)])

;
; Blank pixels in region
; 
map(xrange(0):xrange(1),yrange(0):yrange(1)) = bad_data 

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


