Pro run_bias_stability

!P.multi = [0,4,2]
!P.charsize = 1.5

;dir_name = [  '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158400/','/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158912/', '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159424/', '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159936/']
dir_name = [  '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159424/','/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159680/', '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159936/','/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037160192/']
command1 = 'ls IRAC.1*bcd_fp.fits >  bcd_ch1list.txt'
command2 = 'ls IRAC.2*bcd_fp.fits >  bcd_ch2list.txt'


gain = 3.82
for ch = 0, 1 do begin
   for j = 0, n_elements(dir_name) - 1 do begin
      print, 'working on di', dir_name(j)
      cd, dir_name(j)
      spawn, command1
      spawn, command2
      print, strcompress(dir_name(j) + 'bcd_ch'+string(ch+1) +'list.txt',/remove_all)
      readcol, strcompress(dir_name(j) + 'bcd_ch'+string(ch+1) +'list.txt',/remove_all), imagename, format = 'A'
      
      medianarr = fltarr(n_elements(imagename))
      med_4 = fltarr(4, n_elements(imagename))

      for i = 2,21 do begin
         
         fits_read, imagename(i), bcddata, bcdheader
         
                                ;divide into 4 readouts
         cnvrt_stdim_to_64_256_4, bcddata, bcddata_4
      
         for read = 0, 3 do begin
            median_read = median(bcddata_4[*,*,read])
            bcddata_4[*,*,read] = bcddata_4[*,*,read] - median_read
            med_4(read,i) = median_read
         endfor
         
         
      endfor
      
      plot, findgen(n_elements(imagename)), med_4(0,*)/ max(med_4(0,*)),  psym = 1, xtitle = 'frame number', $
            ytitle = 'median of background regions (MJy/sr)',xrange=[0,20],$
            title = "channel "+string( ch+1)+ "    exptime 6s ",yrange = [0.90,1.05]
      oplot, findgen(n_elements(imagename)), med_4(1,*)/ max(med_4(1,*)),  psym = 2
      oplot, findgen(n_elements(imagename)), med_4(2,*)/ max(med_4(2,*)),  psym = 4
      oplot, findgen(n_elements(imagename)), med_4(3,*)/ max(med_4(3,*)),  psym = 5
      
      
      print, 'x', findgen(n_elements(imagename))
      print, 'y', med_4(0,*)/ max(med_4(0,*))
   endfor
endfor



end
