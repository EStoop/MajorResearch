pro irac_photvstime_bcd

ps_open, filename='/Users/jkrick/irac_warm/test.ps',/portrait,/square,/color
!P.multi = [0,2,2]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)

;time divide between cryo and warm
time_div = 9.2E8
;for array dependant photometric correction
fits_read, '/Users/jkrick/iwic/photcorr/ch1_photcorr_rj.fits', photcor_cryo_ch1, photcorhead_cryo_ch1
fits_read, '/Users/jkrick/iwic/ch1_photcorr_warm.fits', photcor_warm_ch1, photcorhead_warm_ch1

fits_read, '/Users/jkrick/iwic/photcorr/ch2_photcorr_rj.fits', photcor_cryo_ch2, photcorhead_cryo_ch2
fits_read, '/Users/jkrick/iwic/ch2_photcorr_warm.fits', photcor_warm_ch2, photcorhead_warm_ch2


;for the dark field
;ch1
ra_dark_ch1 = [265.00149D,264.88646D,264.91231D,264.95325D,265.04252D,265.0377D,265.04414D,264.96921D,264.9321D,264.9915D,264.89734D]
dec_dark_ch1 = [68.995166D,68.975255D,69.006478D,68.965557D,68.969728D,68.987122D,68.980173D,68.961311D,68.986669D,68.945809D,68.948206D]

;ch2
ra_dark_ch2 = [265.1829D, 265.05484D,265.15846D,265.24624D,265.12223D,265.1564D,265.21455D,265.04336D,265.02307D];265.01324D,
dec_dark_ch2 = [69.052068D, 69.050512D,69.032584D,69.07276D,69.05108D,69.043062D,69.085221D,69.068594D,69.076821D];69.06369D,

dirname_dark = [ 'r24049152','r24046336','r23848448','r23845632','r27990784','r27993600', 'r37658880','r37812224','r37827840','r37831680','r37840640','r37923328','r38013696','r38013952','r38014208','r38054656']
dirname_dark_ch1 = '/Users/jkrick/irac_warm/'+ dirname_dark+'/ch1/bcd'
dirname_dark_ch2 = '/Users/jkrick/irac_warm/'+ dirname_dark+'/ch2/bcd'

;for the calstars
name = ['1812095_cal', 'HD165459', 'NPM1p60']
ra_cal = [273.04000000D,270.628089D , 261.2178D]
dec_cal = [63.49508333D, 58.627264D, 60.430817D]



;HD156896	38048512	
;HD156896	38048768	
;HD156896	38049024	
;HD156896	38049280	
;HD156896	38049536	
;HD156896	38049792	
;HD156896	38048256	
dirname_cal = ['r23525120','r23525376','r23526656','r23527424','r23527680','r23528960','r23840512','r23840768','r23842048','r23842816','r23843072','r23844352','r24040704','r24040960','r24042240','r24043008','r24043264','r24044544','r28003328','r28003584','r28004352','r28006400','r28006656','r28007424','r37651456','r37651968','r37652224','r37815552','r37816064','r37816320','r37851136','r37851648','r37851904','r37997056','r37993728','r37993984','r37913856','r37994496','r37996800','r37914112','r37994240','r37913344','r38052608', 'r38052352','r38051840']
dirname_cal_ch1 = '/Users/jkrick/irac_warm/calstars/'+ dirname_cal+'/ch1/bcd'
dirname_cal_ch2 = '/Users/jkrick/irac_warm/calstars/'+ dirname_cal+'/ch2/bcd'



for ch = 0, 1 do begin
   print, 'working on ch', ch+1
   chtitle = strcompress('ch '+string(ch + 1),/remove_all)

   percent_drop_arr= fltarr(n_elements(ra_dark_ch1) + n_elements(ra_cal))
   electron_arr = fltarr(n_elements(ra_dark_ch1) + n_elements(ra_cal))
   mjy_arr= fltarr(n_elements(ra_dark_ch1) + n_elements(ra_cal))
   n = 0

   for star = 0, 1 do begin

      if star eq 0 and ch eq 0 then begin
         ra =ra_dark_ch1
         dec =  dec_dark_ch1
         dirname = dirname_dark_ch1
      endif 
       if star eq 0 and ch eq 1 then begin
         ra =ra_dark_ch2
         dec =  dec_dark_ch2
         dirname = dirname_dark_ch2
      endif 
       if star eq 1 and ch eq 0 then begin
         ra = ra_cal
         dec = dec_cal
         dirname = dirname_cal_ch1
      endif
       if star eq 1 and ch eq 1 then begin
         ra = ra_cal
         dec = dec_cal
         dirname = dirname_cal_ch2
      endif
      
      
      for r = 0,  n_elements(ra) - 1 do begin
         print, 'working on ra', ra(r), dec(r)
         
            mean_clock = fltarr(n_elements(dirname))
            mean_flux =fltarr(n_elements(dirname))
            k = 0
         
         for d = 0, n_elements(dirname) - 1 do begin
;            print, 'working on ' ,dirname(d)
            j=0
            bcdarr = fltarr(200)
            clockarr = fltarr(200)
            
            command = 'ls '+dirname(d) +'/*_cbcd.fits > /Users/jkrick/irac_warm/ch1list.txt'
            spawn, command
            
;read in the list of fits files
            readcol,'/Users/jkrick/irac_warm/ch1list.txt', fitsname, format = 'A', /silent
 
            for i =1, n_elements(fitsname) - 1 do begin
               fits_read,fitsname(i), bcddata, bcdheader
;               print, 'working on fitsname', fitsname(i)
              clock = sxpar(bcdheader, 'SCLK_OBS')
               
               adxy, bcdheader, ra(r), dec(r), bcdx, bcdy
                        
               if bcdx gt 20 and bcdy gt 20 and bcdy lt 240 and bcdx lt 240 then begin
                  gcntrd, bcddata, bcdx, bcdy, bcdxcen, bcdycen, 3

                                ;if cntrd fails, it returns positions
                                ;-1 -1, catch this error

                  if bcdxcen lt 0 and bcdycen lt 0 then begin
                     bcdxcen = 25  ;then this aperture photometry will get thrown out in the 2sigma clipping
                     bcdycen = 25
                  endif

                 ;conver the data from MJy/sr to electrons
                  gain = 3.7
                  exptime =  sxpar(bcdheader, 'EXPTIME') 
                  fluxconv =  sxpar(bcdheader, 'FLUXCONV') 
                  sbtoe = gain*exptime/fluxconv
                  bcddata = bcddata*sbtoe

                                ;correct for array location dependence
              ; if clock lt time_div then bcddata = bcddata * photcor_cryo_ch1 ;cryo
               ;if clock gt time_div then bcddata = bcddata *photcor_warm_ch1  ;warm
               
                                ;do the photometry
                  aper, bcddata, bcdxcen, bcdycen, bcdflux, errap, sky, skyerr, 1, [2], [3,7],/NAN,/flux,/silent
 
                                ;convert back to Mjy/sr
                  bcdflux = bcdflux /sbtoe

              ;aperture correction
                  if ch eq 0 then bcdflux = bcdflux *1.124 ;for ch1
                  if ch eq 1 then bcdflux = bcdflux *1.127 ;for ch2

                                ;make a correction for pixel phase 
 
                  xphase = (bcdxcen mod 1) - NINT(bcdxcen mod 1)
                  yphase = (bcdycen mod 1) - NINT(bcdycen mod 1)
                  p = sqrt(xphase^2 + yphase^2)
                                ;cryo
                  if clock lt time_div and ch eq 0 then begin
                     correction = 1 + 0.0535*(1/sqrt(2*!Pi) - p)
                     corrected_flux = bcdflux / correction
                  endif
                  ; there is no cryo correction for ch2
                   if clock lt time_div and ch eq 1 then begin
                     corrected_flux = bcdflux 
                  endif
                 
                                ;warm
                   if ch eq 0 then restore, '/Users/jkrick/iwic/pixel_phase_img_ch1.sav'   
                   if ch eq 1 then restore, '/Users/jkrick/iwic/pixel_phase_img_ch2.sav'   

                  if clock gt time_div then begin
                     interpolated_relative_flux = interp2d(relative_flux,x_center,y_center,xphase,yphase,/regular)
                     corrected_flux = bcdflux / interpolated_relative_flux
                  endif
               
                     ;apply array dependent correction
               ;cold
                  if clock lt time_div then begin
                     photcor_ch1 = photcor_cryo_ch1(bcdxcen, bcdycen)
                     photcor_ch2 = photcor_cryo_ch2(bcdxcen, bcdycen)
;                  print, 'photcor', photcor_ch1, photcor_ch2
                     if ch eq 0 then corrected_flux = corrected_flux * photcor_ch1
                     if ch eq 1 then corrected_flux = corrected_flux * photcor_ch2
                  endif
                                ;warm
                  if clock gt time_div then begin
                     photcor_ch1 = photcor_warm_ch1(bcdxcen, bcdycen)
                     photcor_ch2 = photcor_warm_ch2(bcdxcen, bcdycen)
;                  print, 'photcor', photcor_ch1, photcor_ch2
                     if ch eq 0 then corrected_flux = corrected_flux * photcor_ch1
                     if ch eq 1 then corrected_flux = corrected_flux * photcor_ch2
                  endif
                  
                  
                                ;put into microjy
                  corrected_flux = corrected_flux * 34.9837
                  bcdarr[j] = corrected_flux
                  clockarr[j]= clock
                  j = j +1
;                  print, 'bcdflux', corrected_flux
                  
               endif   ;the star is on the frame
               
            endfor              ;end for each fits image

            if j gt 2 then begin
               bcdarr = bcdarr[0:j-1]
               clockarr = clockarr[0:j-1]
;               print, 'bcdarr',  bcdarr
               meanclip, bcdarr, meanbcd, sigmabcd, clipsig = 2
               print, 'mean', meanbcd
;               print, 'clockarr',clockarr
               print, 'mean', mean(clockarr)
               mean_clock[k] = mean(clockarr)
               mean_flux[k] = meanbcd
;               print, 'k', k
               k = k + 1
            endif
 
         endfor                 ; end for each directory = AOR
         mean_clock = mean_clock[0:k-1]
         mean_flux = mean_flux[0:k-1]

         ;sort them
         mean_clock = mean_clock[sort(mean_clock)]
         mean_flux = mean_flux[sort(mean_clock)]

         plot, mean_clock, mean_flux, psym = 2, title = chtitle + string(mean_flux(0)/34.9837*sbtoe),$
               xtitle = 'sclk', ytitle = 'mean flux in mJy', yrange = [mean_flux(0) - 0.15*mean_flux(0), mean_flux(0) + 0.15*mean_flux(0)]
                   

                                ;now try to track the percent change
                                ;in the mean from cryo to warm

         cryo = where(mean_clock lt time_div)
         warm = where(mean_clock gt time_div)
         meanclip, mean_flux(cryo),cryo_mean, cryosigma, clipsig = 2
         meanclip, mean_flux(warm),warm_mean, warmsigma, clipsig = 2

 ;        cryo_mean = mean(mean_flux(cryo))
;         warm_mean = mean(mean_flux(warm))
         print, 'cryomean, warmmean', cryo_mean, warm_mean
         percent_drop_arr[n] = (cryo_mean - warm_mean) / cryo_mean*100.
         electron_arr[n] = mean_flux(0)/34.9837*sbtoe
         mjy_arr[n] = cryo_mean
         print, 'percent drop', (cryo_mean - warm_mean) / cryo_mean*100.,  mean_flux(0)/34.9837*sbtoe/122100.0, cryo_mean
         n = n + 1
      endfor                    ; end for each ra, for each star
      
   endfor                       ; for calstars and dark field stars
   print, 'percentdroparr', percent_drop_arr
   plot, electron_arr/122100.*100., percent_drop_arr, psym = 2, title = 'ch'+string(ch+1), xtitle = 'fraction of full well', ytitle = 'percent drop'
   plot, mjy_arr, percent_drop_arr, /xlog, psym = 2, title = 'ch'+string(ch+1), xtitle = 'mJy', ytitle = 'percent drop', xrange = [1E1,1E6]

endfor                          ; for each channel
ps_close, /noprint,/noid

end
