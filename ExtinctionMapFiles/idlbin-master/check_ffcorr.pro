pro check_ffcorr
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
!P.multi = [0,1,1]
;restore,  '/Volumes/irac_drive/virgo/irac/ffeffect.sav'
restore,  '/Users/jkrick/virgo/irac/ffeffect.sav'
at19 =  result(0) + result(1)*alog10(19)

  dirloc = '/Users/jkrick/virgo/irac/'
  aorname = ['r35320064/','r35320320/','r35320576/','r35320832/','r35321088/','r35321344/','r35321600/','r35321856/','r35322112/','r35322368/','r35322624/','r35322880/','r35323136/','r35323392/','r35323648/','r35323904/','r35324160/','r35324416/','r35324672/','r35324928/','r35325184/','r35325440/','r35325696/','r35325952/','r35326208/']

for ch = 0, 0 do begin
   meanarr = fltarr(5000)
   timearr = fltarr(5000)
   delayarr = fltarr(5000)
   d = 0

   oldmeanarr = fltarr(5000)
   oldtimearr = fltarr(5000)
   olddelayarr = fltarr(5000)
   od = 0

   for a = 0, n_elements(aorname) - 1 do begin
      
      cd, strcompress(dirloc + aorname(a) + 'ch' + string(ch + 1) + '/bcd',/remove_all)
      spawn, 'pwd'
      if ch eq 0 then command =  "find . -name 'corr_SPITZER_I1*0_ndark_ff.fits' > ndark_ff_ch1.txt"
      if ch eq 1 then command =  "find . -name 'corr_SPITZER_I2*0_ndark_ff.fits' > ndark_ff_ch2.txt"
      spawn, command
      
      ;read in a list of all ndark and ff corrected images
      readcol,strcompress( 'ndark_ff_ch'+ string((ch+1))+'.txt',/remove_all), fitsname, format = 'A';, /silent
           
      for i =0, n_elements(fitsname) - 1 do begin
         print, 'working on ', fitsname(i)
         ;run SExtractor
         ;spawn, 'which sex'
         cmd = 'sex ' + fitsname(i) + ' -c ../../../default.sex'
         spawn, cmd
         fits_read, fitsname(i), data, header
         fits_read, 'segmentation.fits', segdata, segheader

         ;write out image with all objects turned to NaN's
         m = where(segdata gt 0)
         data(m) = alog10(-1)   ;set the objects to NaN's.
;         newname = strcompress(strmid(fitsname(i), 0, 32) + 'back.fits',/remove_all)
;         fits_write, newname, data, segheader

        ;keep track of mean background levels, and time of
        ;the exposures.  
        timearr(d)  = sxpar(header, 'SCLK_OBS')
        meanarr(d) = mean(data,/nan)
        delayarr(d)  = sxpar(header, 'FRAMEDLY')

        d = d + 1

      endfor   ;for each fits image
;---------------------------------
;do the same for the old not first frame corrected images, so I can compare.

;      if ch eq 0 then command =  "find . -name 'corr_SPITZER_I1*0_ndark.fits' > ndark_ch1.txt"
;      if ch eq 1 then command =  "find . -name 'corr_SPITZER_I2*0_ndark.fits' > ndark_ch2.txt"
;      spawn, command
;      
;      ;read in a list of all ndark and ff corrected images
;      readcol,strcompress( 'ndark_ch'+ string((ch+1))+'.txt',/remove_all), fitsname, format = 'A';, /silent
;           
;      for i =0, n_elements(fitsname) - 1 do begin
;         ;run SExtractor
;         ;spawn, 'which sex'
;         cmd = 'sex ' + fitsname(i) + ' -c ../../../default.sex'
;         spawn, cmd
;         fits_read, fitsname(i), data, header
;         fits_read, 'segmentation.fits', segdata, segheader
;
;         ;write out image with all objects turned to NaN's
;         m = where(segdata gt 0)
;         data(m) = alog10(-1)   ;set the objects to NaN's.
;;         newname = strcompress(strmid(fitsname(i), 0, 32) + 'back.fits',/remove_all)
;;;         fits_write, newname, data, segheader
;
;        ;keep track of mean background levels, and time of
;        ;the exposures.  
;        oldtimearr(od)  = sxpar(header, 'SCLK_OBS')
;        oldmeanarr(od) = mean(data,/nan)
;        olddelayarr(od)  = sxpar(header, 'FRAMEDLY')
;
;        od = od + 1
;
;      endfor   ;for each fits image
;
   endfor ;for each AOR

   ;and plot mean background vs. time of all AORs.
;   oldmeanarr = oldmeanarr[0:d-1]
;   oldtimearr = oldtimearr[0:d-1]
;   if ch eq 0 then plot, (oldtimearr - oldtimearr(0)) / 3600, oldmeanarr, psym = 2, xtitle = 'time(hrs)', ytitle = 'mean levels(Mjy/sr)', yrange = [-0.2, 0.1]
;   if ch eq 1 then oplot, (timearr - timearr(0)) / 3600, meanarr, psym = 2, color = redcolor

   ;try plotting first frame effect

;   olddelayarr = olddelayarr[0:d-1]
;   if ch eq 0 then plot, delayarr, meanarr, psym = 2, xtitle = 'delay time (s)', ytitle = 'mean levels(Mjy/sr)', yrange = [-0.2, 0.1], xrange = [0,20]
;   if ch eq 1 then oplot, delayarr, meanarr, psym = 2, color = redcolor

   

   ;save these to work with them later
   ;instead of re-running SExtractor

   save, /variables, filename =strcompress( '/Users/jkrick/virgo/irac/ndark_ff_check_ch' + string(ch + 1)+'.sav',/remove_all)

endfor  ;for each channel

end
