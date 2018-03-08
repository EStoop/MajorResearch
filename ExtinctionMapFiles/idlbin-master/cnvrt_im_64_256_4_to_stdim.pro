pro cnvrt_im_64_256_4_to_stdim,imin,imout
; converts an image in (col,row,output) format 
;   to the usual (col,row) format
; imin = array(64,256,4)
; imout= array(256,256)
;
type=size(imin,/type)
imsize=size(imin)

if type eq 2 then imout=intarr(imsize(2),imsize(2))
if type eq 3 then imout=lonarr(imsize(2),imsize(2))
if type eq 4 then imout=fltarr(imsize(2),imsize(2))
if type eq 5 then imout=dblarr(imsize(2),imsize(2))

for io=0,3 do begin
for ic=0,imsize(2) / 4 -1 do begin
imout[io+4*ic,*]=imin[ic,*,io]
endfor
endfor
end
