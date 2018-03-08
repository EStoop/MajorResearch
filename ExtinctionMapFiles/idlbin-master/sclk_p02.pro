pro sclk_p02

;print out sclk values for the 0.02s subarray staring mode AORs
;so far these are 55cnc AORs

  aorname = ['r39524608', 'r43981056', 'r43981312','r43981568', 'r43981824'] ;ch2
  for a =0,   n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = '/Users/jkrick/irac_warm/pcrs_planets/55cnc/'+ string(aorname(a) ) 
     CD, dir                    ; change directories to the correct AOR directory
     command  =  "find ch2/bcd -name '*_bcd.fits' > /Users/jkrick/irac_warm/pcrs_planets/55cnc/bcdlist.txt"
     spawn, command
     readcol,'/Users/jkrick/irac_warm/pcrs_planets/55cnc/bcdlist.txt',fitsname, format = 'A', /silent
 
     outname = 'sclk_'+ string(aorname(a))+'.txt'
     openw, outlun, outname,/get_lun

     for i =0.D, n_elements(fitsname) - 1 do begin
        header = headfits(fitsname(i)) ;
        printf, outlun, sxpar(header, 'SCLK_OBS')
     endfor
  endfor



end
