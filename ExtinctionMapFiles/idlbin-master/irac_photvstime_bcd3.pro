pro irac_photvstime_bcd

ps_open, filename='/Users/jkrick/irac_warm/darks/compare_pctoc.ps',/portrait,/square,/color
!P.multi = [0,2,1]
  redcolor = FSC_COLOR("Red", !D.Table_Size-2)
  bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
  greencolor = FSC_COLOR("Green", !D.Table_Size-4)
  yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
  cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
  orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
  purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
  whitecolor = FSC_COLOR("white", !D.Table_Size-9)
  colorarr = [ redcolor ,bluecolor,greencolor,yellowcolor,cyancolor,purplecolor,orangecolor,whitecolor]

;time divide between cryo and warm
time_div = 9.2E8
;for array dependant photometric correction
fits_read, '/Users/jkrick/iwic/photcorr_cryo/ch1_photcorr_rj.fits', photcor_cryo_ch1, photcorhead_cryo_ch1
fits_read, '/Users/jkrick/iwic/ch1_photcorr_gauss.fits', photcor_warm_ch1, photcorhead_warm_ch1

fits_read, '/Users/jkrick/iwic/photcorr_cryo/ch2_photcorr_rj.fits', photcor_cryo_ch2, photcorhead_cryo_ch2
fits_read, '/Users/jkrick/iwic/ch2_photcorr_gauss.fits', photcor_warm_ch2, photcorhead_warm_ch2


;for the dark field
;ch1
ra_dark_ch1 = [264.95325D,265.04252D,265.04414D,264.96921D,264.9321D,264.9915D,264.89734D] ;265.00149D,264.91231D,265.0377D,264.88646D,
dec_dark_ch1 = [68.965557D,68.969728D,68.980173D,68.961311D,68.986669D,68.945809D,68.948206D] ;68.995166D,69.006478D,68.987122D,68.975255D,

;ch2
ra_dark_ch2 = [265.1829D, 265.05484D,265.15846D,265.24624D,265.12223D,265.1564D,265.04336D];265.01324D,265.21455D,,265.02307D
dec_dark_ch2 = [69.052068D, 69.050512D,69.032584D,69.07276D,69.05108D,69.043062D,69.068594D];69.06369D,69.085221D,,69.076821D

;dirname_dark = [ 'r24049152','r24046336','r23848448','r23845632','r27990784','r27993600', 'r37658880','r37812224','r37827840','r37831680','r37840640','r37923328','r38013696','r38013952','r38014208','r38054656']
dirname_dark = ['r23845632','r23848448','r24046336','r24049152','r27990784','r27993600','r38013952','r38014208','r38054656','r38013696']
dirname_dark_ch1 = '/Users/jkrick/irac_warm/darks/'+ dirname_dark+'/ch1/bcd'
dirname_dark_ch2 = '/Users/jkrick/irac_warm/darks/'+ dirname_dark+'/ch2/bcd'



for ch = 0, 1 do begin
   print, 'working on ch', ch+1
   chtitle = strcompress('ch '+string(ch + 1),/remove_all)

   percent_drop_arr= fltarr(n_elements(ra_dark_ch1) )
   electron_arr = fltarr(n_elements(ra_dark_ch1) )
   mjy_arr= fltarr(n_elements(ra_dark_ch1) )
   clock_arr= fltarr(n_elements(ra_dark_ch1) )

   n = 0

   if  ch eq 0 then begin
      ra =ra_dark_ch1
      dec =  dec_dark_ch1
      dirname = dirname_dark_ch1
   endif 
   if  ch eq 1 then begin
      ra =ra_dark_ch2
      dec =  dec_dark_ch2
      dirname = dirname_dark_ch2
   endif
      
   for r = 0,  n_elements(ra) - 1 do begin
      print, 'working on ra', ra(r), dec(r)
         
      mean_clock = fltarr(n_elements(dirname))
      mean_flux =fltarr(n_elements(dirname))
      k = 0
      
      for d = 0, n_elements(dirname) - 1 do begin
         print, 'working on ' ,dirname(d)
         j=0
         bcdarr = fltarr(200)
         clockarr = fltarr(200)
         
         command = 'ls '+dirname(d) +'/*_cbcd.fits > /Users/jkrick/irac_warm/ch1list.txt'
         spawn, command
         command = 'ls '+dirname(d) +'/*bcd_fp.fits >> /Users/jkrick/irac_warm/ch1list.txt'
         spawn, command
            
;read in the list of fits files
         readcol,'/Users/jkrick/irac_warm/ch1list.txt', fitsname, format = 'A', /silent
         
         for i =1, n_elements(fitsname) - 1 do begin
            fits_read,fitsname(i), bcddata, bcdheader
;               print, 'working on fitsname', fitsname(i)
            clock = sxpar(bcdheader, 'SCLK_OBS')
               
            adxy, bcdheader, ra(r), dec(r), bcdx, bcdy
                        
            if bcdx gt 20 and bcdy gt 20 and bcdy lt 240 and bcdx lt 240 then begin
             
               get_centroids, fitsname(i), t, dt, bcdxcen, bcdycen, bcdflux, xs, ys, fs, b, $
                              /WARM, /APER, RA=ra(r), DEC=dec(r),/silent
                                ;if cntrd fails, it returns positions
                                ;-1 -1, catch this error
               if bcdxcen lt 0 and bcdycen lt 0 then begin
                  bcdxcen = 25  ;then this aperture photometry will get thrown out in the 2sigma clipping
                  bcdycen = 25
               endif

              ;make a correction for pixel phase 
                                 ;cryo
               if clock lt time_div  then begin
                  if ch eq 0 then begin
                     xphase = (bcdxcen mod 1) - NINT(bcdxcen mod 1)
                     yphase = (bcdycen mod 1) - NINT(bcdycen mod 1)
                     p = sqrt(xphase^2 + yphase^2)                     
                     correction = 1 + 0.0535*(1/sqrt(2*!Pi) - p)
                     corrected_flux = bcdflux / correction
                  endif    
                  if ch eq 1 then corrected_flux = bcdflux ; there is no cryo correction for ch2
               endif

                                ;warm
               if clock gt time_div then begin
                  ;corrected_flux = correct_pixel_phase(ch+1,bcdxcen,bcdycen,bcdflux)
                  corrected_flux = pixel_phase_correct_gauss(bcdflux,bcdxcen,bcdycen,ch+1)

               endif
               
               ;apply array dependent correction
                                ;cold
               if clock lt time_div then begin
                  photcor_ch1 = photcor_cryo_ch1(bcdxcen, bcdycen)
                  photcor_ch2 = photcor_cryo_ch2(bcdxcen, bcdycen)
                  if ch eq 0 then corrected_flux = corrected_flux * photcor_ch1
                  if ch eq 1 then corrected_flux = corrected_flux * photcor_ch2
               endif
                                ;warm
               if clock gt time_div then begin
                  photcor_ch1 = photcor_warm_ch1(bcdxcen, bcdycen)
                  photcor_ch2 = photcor_warm_ch2(bcdxcen, bcdycen)
                  if ch eq 0 then corrected_flux = corrected_flux * photcor_ch1
                  if ch eq 1 then corrected_flux = corrected_flux * photcor_ch2
               endif
                
               bcdarr[j] = corrected_flux*1000.
               clockarr[j]= clock
               j = j +1
                  
            endif               ;the star is on the frame
               
         endfor                 ;end for each fits image

         if j gt 2 then begin
            bcdarr = bcdarr[0:j-1]
            clockarr = clockarr[0:j-1]
            meanclip, bcdarr, meanbcd, sigmabcd, clipsig = 2
            print, 'mean', meanbcd
            print, 'mean clock', mean(clockarr)
            mean_clock[k] = mean(clockarr)
            mean_flux[k] = meanbcd
            k = k + 1
         endif
 
      endfor                    ; end for each directory = AOR
      mean_clock = mean_clock[0:k-1]
      mean_flux = mean_flux[0:k-1]

         ;sort them
      mean_clock = mean_clock[sort(mean_clock)]
      mean_flux = mean_flux[sort(mean_clock)]
      print, 'mean_clock, mean_flux', mean_clock, mean_flux

      if r eq 0 then  begin
         plot, mean_clock, mean_flux, psym = 2, title = chtitle ,xrange = [8.6E8,9.5E8],$
               xtitle = 'sclk', ytitle = 'mean flux in mJy', charsize = 1,$
               yrange = [0,1.5], ystyle = 1
                                ;          yrange = [mean(mean_flux) -
                                ;          0.5*mean(mean_flux),
                                ;          mean(mean_flux) +
                                ;          0.5*mean(mean_flux)]
      endif else begin
         oplot, mean_clock, mean_flux, psym = 2, color = colorarr(r-1)
      endelse

                   

                                ;now try to track the percent change
                                ;in the mean from cryo to warm

      cryo = where(mean_clock lt time_div,cc)
      warm = where(mean_clock gt time_div,wc)
      print, 'counts', wc, cc

      if cc gt 0 then begin
         meanclip, mean_flux(cryo),cryo_mean, cryosigma, clipsig = 2
      endif

      if wc gt 0 then begin
         meanclip, mean_flux(warm),warm_mean, warmsigma, clipsig = 2
         
         percent_drop_arr[n] = (cryo_mean - warm_mean) / cryo_mean*100.
;      electron_arr[n] = mean_flux(0)/34.9837*sbtoe
         mjy_arr[n] = cryo_mean
         n = n + 1
      endif
      
   endfor                       ; end for each ra, for each star
   

   print, 'percentdroparr', percent_drop_arr
;   plot, electron_arr/122100.*100., percent_drop_arr, psym = 2, title = 'ch'+string(ch+1), xtitle = 'fraction of full well', ytitle = 'percent drop'
;   plot, mjy_arr, percent_drop_arr, /xlog, psym = 2, title = 'ch'+string(ch+1), xtitle = 'mJy', ytitle = 'percent drop', xrange = [1E1,1E6]
   
endfor                          ; for each channel
ps_close, /noprint,/noid

end
