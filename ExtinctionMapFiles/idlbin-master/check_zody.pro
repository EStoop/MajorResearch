pro check_zody

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)

  dirloc = '/Users/jkrick/virgo/irac/'
  aorname = ['r35320064/','r35320320/','r35320576/','r35320832/','r35321088/','r35321344/','r35321600/','r35321856/','r35322112/','r35322368/','r35322624/','r35322880/','r35323136/','r35323392/','r35323648/','r35323904/','r35324160/','r35324416/','r35324672/','r35324928/','r35325184/','r35325440/','r35325696/','r35325952/','r35326208/']

for ch = 0, 0 do begin
   zodiarr = fltarr(5000)
   timearr = fltarr(5000)
   delayarr = fltarr(5000)
   d = 0

   for a = 0, n_elements(aorname) - 1 do begin
      
      cd, strcompress(dirloc + aorname(a) + 'ch' + string(ch + 1) + '/bcd',/remove_all)
      spawn, 'pwd'
      if ch eq 0 then command =  "find . -name 'SPITZER_I1*0_ndark_ff.fits' > ndark_ff_ch1.txt"
      if ch eq 1 then command =  "find . -name 'SPITZER_I2*0_ndark_ff.fits' > ndark_ff_ch2.txt"
      spawn, command
      
      ;read in a list of all ndark and ff corrected images
      readcol,strcompress( 'ndark_ff_ch'+ string((ch+1))+'.txt',/remove_all), fitsname, format = 'A';, /silent
           
      for i =0, n_elements(fitsname) - 1 do begin
         ;fits_read, fitsname(i), data, header
         header = headfits(fitsname(i)) ;

        ;keep track of mean background levels, and time of
        ;the exposures.  
        timearr(d)  = sxpar(header, 'SCLK_OBS')
        delayarr(d)  = sxpar(header, 'FRAMEDLY')
        zodiarr(d) = sxpar(header, 'ZODY_EST')
        d = d + 1

      endfor   ;for each fits image
   endfor
timearr = timearr[0:d-1]
delayarr = delayarr[0:d-1]
zodiarr = zodiarr[0:d-1]

sorttime = timearr[sort(timearr)]
sortzodi = zodiarr[sort(timearr)]

x = findgen(n_elements(sorttime))
plot, x, sortzodi, psym = 2, yrange = [0.1,0.12]

endfor                          ;for each channel

save, /variables, filename ='/Users/jkrick/virgo/irac/zodi.sav'
end
