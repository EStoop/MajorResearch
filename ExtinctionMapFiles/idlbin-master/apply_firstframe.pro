pro apply_firstframe
  
  restore,  '/Users/jkrick/virgo/irac/ffeffect_pc15.sav'
  at19 =  result(0) + result(1)*alog10(19)

  dirloc = '/Users/jkrick/virgo/irac/'

 aorname =['r35322368/', 'r35325952/', 'r35321856/', 'r35324928/','r35326208/','r35322112/','r35325440/'] ;pc15

 

  for ch =  0, 0 do begin
     
     for a = 0, n_elements(aorname) - 1 do begin
        
        cd, strcompress(dirloc + aorname(a) + 'ch' + string(ch + 1) + '/bcd',/remove_all)
        spawn, 'pwd'
        if ch eq 0 then command =  "find . -name 'corr_SPITZER_I1*_ndark.fits' > ndark_ch1.txt"
        if ch eq 1 then command =  "find . -name 'corr_SPITZER_I2*_ndark.fits' > ndark_ch2.txt"
        spawn, command
        
       command = 'rm -f SPITZER_I1*ndarkndark_ff.fits'
        spawn, command

                                ;read in a list of all images which
                                ;have been already dark subtracted
        readcol,strcompress( 'ndark_ch'+ string((ch+1))+'.txt',/remove_all), fitsname, format = 'A' ;, /silent
       
        for i =1, n_elements(fitsname) - 1 do begin
           fits_read, fitsname(i), data, header
           
                                ;ff effect
                                ;NORMALIZE the effect to delay time of 19 = no correction.
                                ;so I want to add (19) - (x) to the value of x to get the corrected value
           xdelay = sxpar(header, 'FRAMEDLY')
           ffcorr = at19 - ( result(0) + result(1)*alog10(xdelay))
           
           data = data + ffcorr
           
           fits_write, strcompress(strmid(fitsname(i), 0, 37) + 'ndark_ff.fits',/remove_all), data, header
        endfor                  ;for each fits image
        
     endfor                     ;for each AOR
  endfor ; for each channel
  

;-------------------------------------------
;-------------------------------------------

  restore,  '/Users/jkrick/virgo/irac/ffeffect_pc16.sav'
  at19 =  result(0) + result(1)*alog10(19)

  dirloc = '/Users/jkrick/virgo/irac/'
 aorname = ['r35320064/','r35320320/','r35320576/','r35320832/','r35321088/','r35321344/','r35321600/','r35322624/','r35322880/','r35323136/','r35323392/','r35323648/','r35323904/','r35324160/','r35324416/','r35324672/','r35325184/','r35325696/'] ;pc16

  for ch =  0, 0 do begin
     
     for a = 0, n_elements(aorname) - 1 do begin
        
        cd, strcompress(dirloc + aorname(a) + 'ch' + string(ch + 1) + '/bcd',/remove_all)
        spawn, 'pwd'

        command = 'rm -f SPITZER_I1*ndarkndark_ff.fits'
        spawn, command

        if ch eq 0 then command =  "find . -name 'corr_SPITZER_I1*_ndark.fits' > ndark_ch1.txt"
        if ch eq 1 then command =  "find . -name 'corr_SPITZER_I2*_ndark.fits' > ndark_ch2.txt"
        spawn, command
        
                                ;read in a list of all images which
                                ;have been already dark subtracted
        readcol,strcompress( 'ndark_ch'+ string((ch+1))+'.txt',/remove_all), fitsname, format = 'A' ;, /silent
       
        for i =1, n_elements(fitsname) - 1 do begin
           fits_read, fitsname(i), data, header
           
                                ;ff effect
                                ;NORMALIZE the effect to delay time of 19 = no correction.
                                ;so I want to add (19) - (x) to the value of x to get the corrected value
           xdelay = sxpar(header, 'FRAMEDLY')
           ffcorr = at19 - ( result(0) + result(1)*alog10(xdelay))
           
           data = data + ffcorr
           
           fits_write, strcompress(strmid(fitsname(i), 0, 37) + 'ndark_ff.fits',/remove_all), data, header
        endfor                  ;for each fits image
        
     endfor                     ;for each AOR
  endfor ; for each channel

end
