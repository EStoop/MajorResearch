pro list
filename = "hello"

OPENR, lun,'/Users/jkrick/palomar/LFC/fixpix.cl', /GET_LUN  
openw, lunout, '/Users/jkrick/palomar/LFC/fixpix2.cl', /GET_LUN

WHILE (NOT EOF(lun)) DO BEGIN
   readf, lun, filename
   for i =1, 6, 1 do begin
      newname = filename+"["+string(i)+"]"
      printf,lunout, "fixpix ", strcompress(newname,/remove_all)
   endfor
endwhile

close, lun
free_lun, lun
close, lunout
free_lun, lunout
end
