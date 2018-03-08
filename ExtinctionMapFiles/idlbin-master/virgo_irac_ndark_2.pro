pro virgo_irac_ndark
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
!P.multi = [0,1,1]
;restore,  '/Users/jkrick/virgo/irac/ffeffect.sav'
;at19 =  result(0) + result(1)*alog10(19)

  dirloc = '/Users/jkrick/virgo/irac/'
  aorname = ['r35320064/','r35320320/','r35320576/','r35320832/','r35321088/','r35321344/','r35321600/','r35321856/','r35322112/','r35322368/','r35322624/','r35322880/','r35323136/','r35323392/','r35323648/','r35323904/','r35324160/','r35324416/','r35324672/','r35324928/','r35325184/','r35325440/','r35325696/','r35325952/','r35326208/']

 aorname =['r35322368/', 'r35325952/', 'r35321856/', 'r35324928/','r35326208/','r35322112/','r35325440/'] ;pc15

; aorname = ['r35320064/','r35320320/','r35320576/','r35320832/','r35321088/','r35321344/','r35321600/','r35322624/','r35322880/','r35323136/','r35323392/','r35323648/','r35323904/','r35324160/','r35324416/','r35324672/','r35325184/','r35325696/'] ;pc16

fext = [0.91, 0.94]
for ch = 0, 1 do begin
   meanarr = fltarr(5000)
   sigmaarr = fltarr(5000)
   timearr = fltarr(5000)
   delayarr = fltarr(5000)
   d = 0

   for a = 0, n_elements(aorname) - 1 do begin
      
      cd, strcompress(dirloc + aorname(a) + 'ch' + string(ch + 1) + '/bcd',/remove_all)
      spawn, 'pwd'
      if ch eq 0 then command =  "find . -name 'corr_SPITZER_I1*_bcd.fits' > bcdlist_ch1.txt"
      if ch eq 1 then command =  "find . -name 'corr_SPITZER_I2*_bcd.fits' > bcdlist_ch2.txt"
      spawn, command
      
      ;read in a list of all cbcds
      readcol,strcompress( 'bcdlist_ch'+ string((ch+1))+'.txt',/remove_all), fitsname, format = 'A';, /silent
           
      bigim = fltarr(256, 256, n_elements(fitsname))
      c = 0

      for i =1, n_elements(fitsname) - 1 do begin  ;don't use the zero frame
         ;run SExtractor
         cmd = 'sex ' + fitsname(i) + ' -c ../../../default.sex'
         spawn, cmd
         fits_read, fitsname(i), data, header
         fits_read, 'segmentation.fits', segdata, segheader

         ;write out image with all objects turned to NaN's
         m = where(segdata gt 0)
         data(m) = alog10(-1)   ;set the objects to NaN's.
         newname = strcompress(strmid(fitsname(i), 0, 37) + 'back.fits',/remove_all)
;;         fits_write, newname, data, segheader

        ;keep all images together to make median image for self-dark
        bigim(0,0,c) = data
        c = c + 1

        ;keep track of mean background levels, and time of
        ;the exposures.  
        timearr(d)  = sxpar(header, 'SCLK_OBS')
        zody  = sxpar(header, 'ZODY_EST')
        meanarr(d) = mean(data,/nan) - zody/ fext(ch)   ;correcting for extended source
        sigmaarr(d) = stddev(data, /nan)
        delayarr(d)  = sxpar(header, 'FRAMEDLY')
        d = d + 1

      endfor   ;for each fits image

      ;generate median image per AOR
      medarr, bigim, meddark

      ;only want to subtract structure in the dark, not the levels.
      ;so subtract the mean from the images so mean = 0
      print, "mean", mean(meddark,/nan)
      meddark = meddark - mean(meddark,/nan)

;;      fits_write, strcompress( 'med_ch' + string(ch + 1) + '.fits', /remove_all), meddark, imheader

      ;now go back and dark subtract all the cbcd's.
      ;for use in the mosaicing.
      ;also subtract the zody
      ;also correct for the first frame effect 
;taken from plot first frame.pro

      for i =1, n_elements(fitsname) - 1 do begin
         fits_read, fitsname(i), data, header
         data = data - meddark

         zody  = sxpar(header, 'ZODY_EST')
         data = data - zody/fext(ch)

         ;ff effect
         ;NORMALIZE the effect to delay time of 19 = no correction.
         ;so I want to add (19) - (x) to the value of x to get the corrected value
         ;xdelay = sxpar(header, 'FRAMEDLY')
         ;ffcorr = at19 - ( result(0) + result(1)*alog10(xdelay))

         ;if ch eq 0 then data = data + ffcorr
         
         ;now add some overall value so that the levels are above zero.
         data = data + 0.3

;;         fits_write, strcompress(strmid(fitsname(i), 0, 37) + 'ndark.fits',/remove_all), data, header
      endfor   ;for each fits image

   endfor ;for each AOR

   ;and plot mean background vs. time of all AORs.
   meanarr = meanarr[0:d-1]
   timearr = timearr[0:d-1]
   sigmaarr = sigmaarr[0:d-1]
;   if ch eq 0 then plot, (timearr - timearr(0)) / 3600, meanarr, psym = 2, xtitle = 'time(hrs)', ytitle = 'mean levels(Mjy/sr)', yrange = [-0.2, 0.1]
;   if ch eq 1 then oplot, (timearr - timearr(0)) / 3600, meanarr, psym = 2, color = redcolor

   ;try plotting first frame effect

   delayarr = delayarr[0:d-1]
   if ch eq 0 then plot, delayarr, meanarr, psym = 2, xtitle = 'delay time (s)', ytitle = 'mean levels(Mjy/sr)', yrange = [-0.2, 0.1], xrange = [0,20]
   if ch eq 1 then oplot, delayarr, meanarr, psym = 2, color = redcolor

   ;save these to work with them later
   ;instead of re-running SExtractor

   save, /variables, filename =strcompress( '/Users/jkrick/virgo/irac/ndark_ch' + string(ch + 1)+'_pc15.sav',/remove_all)

endfor  ;for each channel

end
