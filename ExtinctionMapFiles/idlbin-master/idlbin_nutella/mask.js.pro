; reads in a list of files and mask images, then nans blank pixels
; in the masks, and write them to an output file.
; Note - the masks were originally IRAF .pl files, converted to fits
;  with imcopy.
PRO mask,filelist,masklist,proclist

readcol,filelist,mosaics,format="A"
readcol,masklist,masks,FORMAT="A"
readcol,proclist,procs,FORMAT="A"


nfiles=n_elements(mosaics)
nan=sqrt(-1)

for i=0,nfiles-1 do begin

mosaic=readfits(mosaics[i],hd)
exptime=sxpar(hd,'EXPTIME')

mask=readfits(masks[i])

proc=mosaic;/exptime
proc(where(mask NE 0))=nan

print, proc(172,983)

;sxaddpar,hd,'EXPTIME',1
writefits,procs[i],proc,hd



endfor

end
