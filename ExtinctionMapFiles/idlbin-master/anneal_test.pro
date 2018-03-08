pro anneal_test
greencolor = FSC_COLOR("Green", !D.Table_Size-4)

!P.charthick = 1
!P.thick = 1
!X.thick = 1
!Y.thick = 1
!P.charsize = 1
;junk = run_mopex('/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158400',  '/Applications/mopex/mopex-script-env.csh', 'mosaic_ch1.nl', 'mosaic_ch2.nl')
;print, 'finished mopex 0037158400'
;junk = run_mopex('/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158912',  '/Applications/mopex/mopex-script-env.csh', 'mosaic_ch1.nl', 'mosaic_ch2.nl')
;print, 'finished mopex 0037158912'
;junk = run_mopex('/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159424',  '/Applications/mopex/mopex-script-env.csh', 'mosaic_ch1.nl', 'mosaic_ch2.nl')
;print, 'finished mopex 0037159424'
;junk = run_mopex('/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159936',  '/Applications/mopex/mopex-script-env.csh', 'mosaic_ch1.nl', 'mosaic_ch2.nl')
;print, 'finished mopex 0037159936'
;junk = run_mopex('/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159680',  '/Applications/mopex/mopex-script-env.csh', 'mosaic_ch1.nl', 'mosaic_ch2.nl')
;print, 'finished mopex 0037159680'
;junk = run_mopex('/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158656',  '/Applications/mopex/mopex-script-env.csh', 'mosaic_ch1.nl', 'mosaic_ch2.nl')
;print, 'finished mopex 0037158656'
;junk = run_mopex('/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037160192',  '/Applications/mopex/mopex-script-env.csh', 'mosaic_ch1.nl', 'mosaic_ch2.nl')
;print, 'finished mopex 0037160192'
;junk = run_mopex('/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159168',  '/Applications/mopex/mopex-script-env.csh', 'mosaic_ch1.nl', 'mosaic_ch2.nl')
;print, 'finished mopex 0037159168'





;----------------------------------------------------------------------------------------
;bias stability

;run_bias_stability

;----------------------------------------------------------------------------------------
;make a histogram of the noise

;run_noise_histogram

;----------------------------------------------------------------------------------------
;how many hot pixels are there?
dir_name =['/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159424', '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159680','/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159936', '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037160192','/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158400', '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158656', '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158912','/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159168']

flux_conv = [ 0.1229,0.1469]    ;adopted from Mark's preliminary numbers
gain = 3.82
noiselimit = [6.3, 12.8]        ;e-/s of the zodiacal light
for ch = 0, 1 do begin
for i = 0, 0 do begin; n_elements(dir_name) - 1 do begin
   print, 'working on dir_name', dir_name(i)

;need to go through all the individual images in a directory
   command1 = 'ls IRAC.1*bcd_fp.fits >  bcd_ch1list.txt'
   command2 = 'ls IRAC.2*bcd_fp.fits >  bcd_ch2list.txt'
   cd, dir_name(i)
   spawn, command1
   spawn, command2

   openw, outlunred, dir_name(i) + '/noisypixmap.reg', /get_lun

   readcol, dir_name(i) + '/bcd_ch2list.txt', imagename, format = 'A'

;make smoothed images to work with.
   edataptr  = [Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256))] ;in e-/s



   for j = 0,21 do begin

      fits_read, imagename(j), bcddata, bcdheader
      
      edata = (bcddata/ flux_conv(ch))*gain
      smoothe =  filter_image(edata, median=5, /all)
      edata = edata - smoothe    
      *edataptr[j] = edata
;      fits_write, strcompress(dir_name +'/e'+strmid(imagename(j),2),/remove_all), *edataptr[j], eheader
   endfor


   hcount = fix(0)
   for x = 0,255 do begin
      for y = 0,255 do begin
;      print, 'x, y', x, y
         hotpixel = 0
         hpix = 0
         fin = 0
         for j = 0, n_elements(edataptr) - 1 do begin
                                ;need result of gaussian fit here so I can't do this in the first loop
            if finite((*edataptr[j])[x,y]) gt 0 then begin
               fin = fin + 1
       
            if (*edataptr[j])[x,y] gt noiselimit[ch] then begin
               hpix = hpix + 1
            endif            
         endif        
      endfor    
                                ;if that same pixel is hot on all of the images where finite then it must be
                                ;a hot pixel
         nfin = 18  
         if fin gt nfin and hpix gt fin - 1 then begin
            hcount = hcount + 1
            printf, outlunred, 'box( ', x, y, ' 3,3,0)'
            
         endif     
      endfor
   endfor

   print, 'ch, i, number of noisy pix', ch, i, hcount
close, outlunred
free_lun, outlunred

endfor
endfor

;----------------------------------------------
;put this information into some sort of plot form
!P.multi = [0,1,1]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
whitecolor = FSC_COLOR("white", !D.Table_Size-9)


mean_ch1_6 = [0.025012558,0.027296263,0.024949444]

mean_ch1_100 =[0.006016277,0.00598022,0.0076553]

mean_ch2_6 = [-0.053739655,-0.05536681,-0.046391015]

mean_ch2_100 = [-0.010824377,-0.007101799,0.006019855]
;---
sigma_ch1_6 = [0.017931831,0.019728481,0.015814558]

sigma_ch1_100 =[0.001437838,0.006116216,-0.003954054]

sigma_ch2_6 = [0.011250643,0.053628905,0.071616335]

sigma_ch2_100 = [-0.000735246,0.027851053,0.040994558]
;---
hot_ch1_6 = [-0.028225806,0.004032258,0.0524193554]

hot_ch1_100 =[0.02247191,0.112359551,-0.039325843]

hot_ch2_6 = [0.023952096,0.047904192,0.053892216]

hot_ch2_100 = [0.007874016,0.141732283,0.015748031]

x = findgen(3)

plot, x, mean_ch1_6, yrange=[-.1, .1], psym = 2, xrange=[-1,3], ytitle = 'Fractional change from pre-anneal values', charsize =1.5, charthick = 2
oplot, x, mean_ch1_6, linestyle = 2
oplot, x, mean_ch1_100, psym = 2
oplot, x, mean_ch1_100
oplot, x, mean_ch2_6, psym = 2, color = redcolor
oplot, x, mean_ch2_6,  color = redcolor, linestyle = 2
oplot, x, mean_ch2_100, psym = 2, color = redcolor
oplot, x, mean_ch2_100,color = redcolor
oplot, x, sigma_ch1_6 , psym = 4
oplot, x, sigma_ch1_6 , linestyle = 2
oplot, x, sigma_ch1_100 , psym = 4
oplot, x, sigma_ch1_100 
oplot, x, sigma_ch2_6  , psym = 4, color = redcolor
oplot, x, sigma_ch2_6,  color = redcolor, linestyle = 2
oplot, x, sigma_ch2_100 , psym = 4, color = redcolor
oplot, x, sigma_ch2_100 , color = redcolor
oplot, x, hot_ch1_6 , psym = 6
oplot, x, hot_ch1_6, linestyle = 2
oplot, x, hot_ch1_100, psym = 6
oplot, x, hot_ch1_100
oplot, x, hot_ch2_6, psym = 6, color = redcolor
oplot, x, hot_ch2_6, color = redcolor, linestyle = 2
oplot, x, hot_ch2_100, psym = 6, color = redcolor
oplot, x, hot_ch2_100, color = redcolor

legend, ['mean noise','sigma noise', 'hot pixels'], psym = [2,4,6], /top, /left, charsize =1.5, charthick = 2
legend, ['ch1', 'ch2'], linestyle = [0, 0], color = [whitecolor, redcolor], /top, /right, charsize =1.5, charthick = 2
legend, ['6s', '100s'], linestyle = [2,0], color = [redcolor, redcolor], /bottom, /right, charsize =1.5, charthick = 2

oplot, [-5,5], [0,0]
end
