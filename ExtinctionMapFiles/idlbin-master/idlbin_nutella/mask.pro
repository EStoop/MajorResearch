PRO mask

nan=sqrt(-1)
readcol,'/Users/jkrick/hst/raw/filelist',filename,format="A"
for n = 0, n_elements(filename) - 1 do begin
   print, filename(n)
   mosaic=readfits(filename(n),hd)
   check_fits, mosaic, hd, dimen

   for x = 0, dimen(0)-1, 1 do begin
      for y = 0, dimen(1)-1, 1 do begin

         if mosaic(x,y) EQ 0 then begin
            mosaic(x,y) = nan
         endif
         
      endfor
   endfor 
   
   writefits,filename(n),mosaic,hd
endfor


end
