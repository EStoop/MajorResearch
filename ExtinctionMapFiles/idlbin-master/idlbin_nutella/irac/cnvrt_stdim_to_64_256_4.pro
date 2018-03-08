pro cnvrt_stdim_to_64_256_4,imin,imout
; converts an image in the usual col,row format
;   to (col,row,array_output)
; imin= array(256,256)
; imout = array(64,256,4)
;
type=size(imin,/type)
imsize=size(imin)

if type eq 2 then imout=intarr(imsize(1) / 4,imsize(1),4)
if type eq 3 then imout=lonarr(imsize(1) / 4,imsize(1),4)
if type eq 4 then imout=fltarr(imsize(1) / 4,imsize(1),4)
if type eq 5 then imout=dblarr(imsize(1) / 4,imsize(1),4)

for io=0,3 do begin
for ic=0,imsize(1) / 4 - 1 do begin
imout[ic,*,io]=imin[io+4*ic,*]
endfor
endfor
end
