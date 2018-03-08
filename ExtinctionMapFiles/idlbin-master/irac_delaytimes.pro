pro irac_delaytimes
!P.multi = [0,1,2]
ps_open, filename='/Users/jkrick/virgo/irac/delaytimes.ps',/portrait,/square,/color


  dirloc = '/Users/jkrick/virgo/irac/'
  aorname = ['r35320064/','r35320320/','r35320576/','r35320832/','r35321088/','r35321344/','r35321600/','r35321856/','r35322112/','r35322368/','r35322624/','r35322880/','r35323136/','r35323392/','r35323648/','r35323904/','r35324160/','r35324416/','r35324672/','r35324928/','r35325184/','r35325440/','r35325696/','r35325952/','r35326208/']
  
  for ch = 0, 1 do begin
     ;generate a master list of all irac bcd files
     cd, dirloc
     if ch eq 0 then command =  "find . -name 'SPITZER_I1*_bcd.fits' > /Users/jkrick/virgo/irac/bcdlist_ch1.txt"
     if ch eq 1 then command =  "find . -name 'SPITZER_I2*_bcd.fits' > /Users/jkrick/virgo/irac/bcdlist_ch2.txt"
     spawn, command

     readcol,strcompress(dirloc + 'bcdlist_ch'+ string((ch+1))+'.txt',/remove_all), fitsname, format = 'A', /silent
        
     interdelayarr = fltarr(n_elements(fitsname))
     delayarr = fltarr(n_elements(fitsname))
     
     for i =1, n_elements(fitsname) - 1 do begin
        fits_read,fitsname(i), bcddata, bcdheader
        delayarr(i)  = sxpar(bcdheader, 'FRAMEDLY')
        interdelayarr(i) = sxpar(bcdheader, 'INTRFDLY')
     endfor                     ;for each bcd fits file
     
     plothist, delayarr, xhist, yhist, bin = 0.1, xrange = [0,20], xtitle = 'framedly (s)'
    ; plothist, interdelayarr, xhist, yhist, bin = 0.1, xrange= [0,20], xtitle = 'intrfdly (s)'


  endfor  ;for each channel
ps_close, /noprint,/noid

end
