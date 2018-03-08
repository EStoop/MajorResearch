pro irac_photvstime_bcd

ps_open, filename='/Users/jkrick/irac_warm/test.ps',/portrait,/square,/color
!P.multi = [0,2,2]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)


;for array dependant photometric correction
fits_read, '/Users/jkrick/iwic/photcorr/ch1_photcorr_rj.fits', photcor_cryo_ch1, photcorhead_cryo_ch1
fits_read, '/Users/jkrick/iwic/ch1_photcorr_warm_gauss_recal.fits', photcor_warm_ch1, photcorhead_warm_ch1

fits_read, '/Users/jkrick/iwic/photcorr/ch2_photcorr_rj.fits', photcor_cryo_ch2, photcorhead_cryo_ch2
fits_read, '/Users/jkrick/iwic/ch2_photcorr_warm_gauss_recal.fits', photcor_warm_ch2, photcorhead_warm_ch2

;for warm pixel phase correction
restore, '/Users/jkrick/iwic/pixel_phase_img_ch1.sav'   
;restore, '/Users/jkrick/iwic/pixel_phase_img_ch2.sav'   

;for the dark field
;ch1
ra_dark_ch1 = [265.00149D,264.88646D,264.91231D,264.95325D]
dec_dark_ch1 = [68.995166D,68.975255D,69.006478D,68.965557D]

;ch2
ra_dark_ch2 = [265.18437D, 265.01324D,265.05484D,265.15846D]
dec_dark_ch2 = [69.052516D, 69.06369D,69.050512D,69.032584D]

dirname_dark = [ '/Users/jkrick/irac_warm/r24049152/ch1/bcd','/Users/jkrick/irac_warm/r24046336/ch1/bcd','/Users/jkrick/irac_warm/r23848448/ch1/bcd','/Users/jkrick/irac_warm/r23845632/ch1/bcd','/Users/jkrick/irac_warm/r27990784/ch1/bcd','/Users/jkrick/irac_warm/r27993600/ch1/bcd', '/Users/jkrick/irac_warm/r37658880/ch1/bcd','/Users/jkrick/irac_warm/r37812224/ch1/bcd','/Users/jkrick/irac_warm/r37827840/ch1/bcd','/Users/jkrick/irac_warm/r37831680/ch1/bcd','/Users/jkrick/irac_warm/r37840640/ch1/bcd','/Users/jkrick/irac_warm/r37923328/ch1/bcd','/Users/jkrick/irac_warm/r38013696/ch1/bcd','/Users/jkrick/irac_warm/r38013952/ch1/bcd','/Users/jkrick/irac_warm/r38014208/ch1/bcd']

;for the calstars
name = ['1812095_cal', 'HD165459', 'NPM1p60']
ra_cal = [273.04000000D,270.628089D , 261.2178D]
dec_cal = [63.49508333D, 58.627264D, 60.430817D]


dirname_cal = ['/Users/jkrick/irac_warm/calstars/r23525120/ch1/bcd','/Users/jkrick/irac_warm/calstars/r23525376/ch1/bcd','/Users/jkrick/irac_warm/calstars/r23526656/ch1/bcd','/Users/jkrick/irac_warm/calstars/r23527424/ch1/bcd','/Users/jkrick/irac_warm/calstars/r23527680/ch1/bcd','/Users/jkrick/irac_warm/calstars/r23528960/ch1/bcd','/Users/jkrick/irac_warm/calstars/r23840512/ch1/bcd','/Users/jkrick/irac_warm/calstars/r23840768/ch1/bcd','/Users/jkrick/irac_warm/calstars/r23842048/ch1/bcd','/Users/jkrick/irac_warm/calstars/r23842816/ch1/bcd','/Users/jkrick/irac_warm/calstars/r23843072/ch1/bcd','/Users/jkrick/irac_warm/calstars/r23844352/ch1/bcd','/Users/jkrick/irac_warm/calstars/r24040704/ch1/bcd','/Users/jkrick/irac_warm/calstars/r24040960/ch1/bcd','/Users/jkrick/irac_warm/calstars/r24042240/ch1/bcd','/Users/jkrick/irac_warm/calstars/r24043008/ch1/bcd','/Users/jkrick/irac_warm/calstars/r24043264/ch1/bcd','/Users/jkrick/irac_warm/calstars/r24044544/ch1/bcd','/Users/jkrick/irac_warm/calstars/r28003328/ch1/bcd','/Users/jkrick/irac_warm/calstars/r28003584/ch1/bcd','/Users/jkrick/irac_warm/calstars/r28004352/ch1/bcd','/Users/jkrick/irac_warm/calstars/r28006400/ch1/bcd','/Users/jkrick/irac_warm/calstars/r28006656/ch1/bcd','/Users/jkrick/irac_warm/calstars/r28007424/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37651456/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37651968/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37652224/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37815552/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37816064/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37816320/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37851136/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37851648/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37851904/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37997056/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37993728/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37993984/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37913856/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37994496/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37996800/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37914112/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37994240/ch1/bcd','/Users/jkrick/irac_warm/calstars/r37913344/ch1/bcd']

for ch = 0, 0 do begin
   print, 'working on ch', ch+1

   for star = 0, 1 do begin

      if star eq 0 and ch eq 0 then begin
         ra =ra_dark_ch1
         dec =  dec_dark_ch1
         dirname = dirname_dark
      endif 
       if star eq 0 and ch eq 1 then begin
         ra =ra_dark_ch2
         dec =  dec_dark_ch2
         dirname = dirname_dark
      endif 
       if star eq 1 then begin
         ra = ra_cal
         dec = dec_cal
         dirname = dirname_cal
      endif
      
      
      for r = 0, n_elements(ra) - 1 do begin
         print, 'working on ra', ra(r), dec(r)
         
         j=0
         bcdarr = fltarr(200)
         clockarr = fltarr(200)
         medianx = fltarr(n_elements(dirname))
         mediany =fltarr(n_elements(dirname))
         
         
         for d = 0, n_elements(dirname) - 1 do begin
            print, 'working on ' ,dirname(d)
            command = 'ls '+dirname(d) +'/*_cbcd.fits > /Users/jkrick/irac_warm/ch1list.txt'
            spawn, command
            
;read in the list of fits files
            readcol,'/Users/jkrick/irac_warm/ch1list.txt', fitsname, format = 'A', /silent
         
         
            thisflux = fltarr(n_elements(fitsname))
            k = 0
            for i =1, n_elements(fitsname) - 1 do begin
               fits_read,fitsname(i), bcddata, bcdheader
;               print, 'working on fitsname', fitsname(i)
              clock = sxpar(bcdheader, 'SCLK_OBS')
               
               adxy, bcdheader, ra(r), dec(r), bcdx, bcdy
                        
               if bcdx gt 20 and bcdy gt 20 and bcdy lt 240 and bcdx lt 240 then begin
                  cntrd, bcddata, bcdx, bcdy, bcdxcen, bcdycen, 3

                                ;correct for array location dependence
              ; if clock lt 9.3E8 then bcddata = bcddata * photcor_cryo_ch1 ;cryo
               ;if clock gt 9.3E8 then bcddata = bcddata *photcor_warm_ch1  ;warm
               
                                ;do the photometry
                  aper, bcddata, bcdxcen, bcdycen, bcdflux, errap, sky, skyerr, 1, [2], [3,7],/NAN,/flux,/silent
               
              ;aperture correction
                  if ch eq 0 then bcdflux = bcdflux *1.124 ;for ch1
                  if ch eq 1 then bcdflux = bcdflux *1.127 ;for ch2

                                ;make a correction for pixel phase 
                  xphase = (bcdxcen mod 1) - NINT(bcdxcen mod 1)
                  yphase = (bcdycen mod 1) - NINT(bcdycen mod 1)
                  p = sqrt(xphase^2 + yphase^2)
                                ;cryo
                  if clock lt 9.3E8 then begin
                     correction = 1 + 0.0535*(1/sqrt(2*!Pi) - p)
                     corrected_flux = bcdflux / correction
                  endif
                  
                                ;warm
                  if clock gt 9.3E8 then begin
                     interpolated_relative_flux = interp2d(relative_flux,x_center,y_center,xphase,yphase,/regular)
                     corrected_flux = bcdflux / interpolated_relative_flux
                  endif
               
                     ;apply array dependent correction
               ;cold
                  if clock lt 9.3E8 then begin
                     photcor_ch1 = photcor_cryo_ch1(bcdxcen, bcdycen)
                     photcor_ch2 = photcor_cryo_ch2(bcdxcen, bcdycen)
;                  print, 'photcor', photcor_ch1, photcor_ch2
                     if ch eq 0 then corrected_flux = corrected_flux * photcor_ch1
                     if ch eq 1 then corrected_flux = corrected_flux * photcor_ch2
                  endif
                                ;warm
                  if clock gt 9.3E8 then begin
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
                  print, 'bcdflux', corrected_flux
                  
                  thisflux[k] = corrected_flux
                  k = k + 1
               endif
               
            endfor              ;end for each fits image
            
                                ;what is the mean value for just this
                                ;AOR?
            if k gt 1 then begin
               thisflux = thisflux[0:k-1]
               medianx[d] = (clockarr[j - 1]- clockarr[0])/3600/24
               ; do a sigma clipped mean
               meanclip, thisflux, meanthisflux, sigma, clipsig = 2
               mediany[d] = meanthisflux; mean(thisflux)
            endif
            
         endfor                 ; end for each directory = AOR
         
         if j gt 1 then begin
            bcdarr = bcdarr[0:j-1]
            clockarr = clockarr[0:j-1]
;      print, 'bcdarr',  bcdarr/ mediany(r)
;      print, 'median', mediany(r)
;      print, 'clockarr', (clockarr- clockarr[0])/3600/24
            
            if star eq 0 then t = 0
            if star eq 1 then t = r
;            if star eq 0 then print, 'dark', clockarr, 'darl',  (clockarr- clockarr[0])/3600/24
;            if star eq 1 then print, 'cal', clockarr, 'cal',  (clockarr- clockarr[0])/3600/24
            
            plot, (clockarr- clockarr[0])/3600/24, bcdarr/ mediany(t), psym = 2, xtitle = 'time in days', ytitle = 'relative microJy', xrange=[-100, 900], yrange = [0.7,1.2], ystyle = 1 ;yrange = [median(bcdarr) - .3*median(bcdarr), median(bcdarr) + .3*median(bcdarr)]
            print, 'median', median(bcdarr), stddev(bcdarr)/ median(bcdarr)
            
            for f = 0, n_elements(dirname) - 1 do begin
               xyouts, medianx(f), mediany(f)/ mediany(t), 'M', color = redcolor, alignment = 0.5
               xyouts, -100, 0.8, string(nint(mediany(t)))+ 'microJy'
            endfor
         endif
         
      endfor                    ; end for each ra, for each star
      
   endfor                       ; for calstars and dark field stars
   
endfor                          ; for each channel
ps_close, /noprint,/noid

end
