pro display_raw_sub

  dir = '/Users/jkrick/irac_warm/snapshots/hd7924/r44605184/'
  CD, dir                       ; change directories to the correct AOR directory
  command  =  " ls ch2/raw/*_dce.fits > /Users/jkrick/irac_warm/snapshots/hd7924/rawlist.txt"
  spawn, command
  readcol,'/Users/jkrick/irac_warm/snapshots/hd7924/rawlist.txt',rawname, format = 'A', /silent
  
  for i = 0, n_elements(rawname) - 1 do begin
     fits_read, rawname(i), data, header
     mc = reform(data, 32, 32, 64)
;tvscl, mc(*,*,0)
     
     img = mc
; convert to real
     img=img*1.
     
; fix the problem with unsigned ints
     fix=where((img LT 0),count)
     if (count GT 0) then img[fix]=img[fix]+65536.
     
; flip the InSb
     img=65535.-img
     
;fits_write, '/Users/jkrick/irac_warm/snapshots/r44497152/ch2/raw/test_0023.fits', img, header
;then use ds9 to display
     if  img[15,16] gt 18000.0 then print, 'got a bright one', img[15,16], rawname(i)
     
;for i = 0 , 63 do begin
;   arr = img(*,*,i)
;   print, 'i', i, arr[15,16]
;endfor
 
  endfor
    
end
  
