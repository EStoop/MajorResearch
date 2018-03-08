pro plot_sd1043_ch1
  colorarr = ['black', 'red', 'orange', 'green', 'blue', 'deep_pink', 'magenta', 'medium_purple', 'hot_pink', 'light_pink', 'rosy_brown', 'chocolate', 'saddle_brown', 'maroon', 'orange_red', 'dark_orange', 'peach_puff', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'light_cyan', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green']


  restore, '/Users/jkrick/irac_warm/pcrs_planets/sd1043/sd1043_phot_ch1.sav'
  aorname_sd1043 = ['r45758208', 'r43336448'] ;ch1

  ;get rid of the zeros
  print, 'test', AORsd1043[0].xcen
  good = where(AORsd1043[0].xcen ne 0, ngood)
  print, 'n_elements', n_elements(AORsd1043[0].xcen), ngood
;  print, 'test again', AORsd1043[0].xcen(good), AORsd1043[0].timearr(good)
  


  ;utmjd_center = [double(55947.20505),double(55949.41467), double(55957.51661),double(55944.25889)]  
 ; utmjd_center = [double(55944.25889), double(55947.20505),double(55949.41467), double(55957.51661)]  
;  duration = 105.7              ;in minutes
;  duration = duration /60./24.  ; in days
;  period = 0.73654 ;days
  

;set time = 0.5 in the middle of the eclipse
;     bin_time = bin_time -tt
;  for a = 0, n_elements(aorname_sd1043) - 1 do  begin
;     test =  ((AORsd1043[a].bmjdarr - utmjd_center[a]) / period) + 0.5
     
;     ;if a eq 1 then begin
;        ;for t = 0, 100*63,63  do begin
;        ;   print, t/63, ' ', AORsd1043[1].bmjdarr[t], ' ',test[t], ' ', utmjd_center[a], format = '(F0,A, F0, A, F0, A, F0)'
;       ; endfor
;     ;endif
;
;     AORsd1043[a].bmjdarr = ((AORsd1043[a].bmjdarr - utmjd_center[a]) / period) + 0.5
;
;  endfor


; 
;  for a = 0, n_elements(aorname_sd1043) - 1 do begin
;     xcen500 = 500.* (AORsd1043[a].xcen - 14.5)
;     ycen500 = 500.* (AORsd1043[a].ycen - 14.5)
;     an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a])
; endfor

;-----
  
    for a = 0, n_elements(aorname_sd1043) - 1 do begin
     if a eq 0 then begin
         am = plot( (AORsd1043[a].timearr(good) - AORsd1043[a].timearr(0))/60./60., AORsd1043[a].xcen(good),'6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'X pix', yrange = [22.5, 23.5])
     endif else begin
        am = plot( (AORsd1043[a].timearr - AORsd1043[0].timearr(0))/60./60., AORsd1043[a].xcen,'6r1s', sym_size = 0.2,   sym_filled = 3, color = colorarr[a],/overplot)
     endelse

  endfor

;------

  for a = 0, n_elements(aorname_sd1043) - 1 do begin
     if a eq 0 then begin
         am = plot( (AORsd1043[a].timearr(good) - AORsd1043[a].timearr(0))/60./60., AORsd1043[a].ycen(good),'6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Y pix', yrange = [230.5, 231.5])
     endif else begin
        am = plot( (AORsd1043[a].timearr - AORsd1043[0].timearr(0))/60./60., AORsd1043[a].ycen,'6r1s', sym_size = 0.2,   sym_filled = 3, color = colorarr[a],/overplot)
     endelse

  endfor

;------
 for a = 0, n_elements(aorname_sd1043) - 1 do begin
     if a eq 0 then begin
         am = plot( (AORsd1043[a].timearr(good) - AORsd1043[a].timearr(0))/60./60., AORsd1043[a].bkgd(good),'6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Bkgd')
     endif else begin
        am = plot( (AORsd1043[a].timearr - AORsd1043[0].timearr(0))/60./60., AORsd1043[a].bkgd,'6r1s', sym_size = 0.2,   sym_filled = 3, color = colorarr[a],/overplot)
     endelse

  endfor

;------
for a = 0, n_elements(aorname_sd1043) - 1 do begin
     if a eq 0 then begin
         am = plot( (AORsd1043[a].timearr(good) - AORsd1043[a].timearr(0))/60./60., AORsd1043[a].flux(good),'6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Flux (Jy)', yrange = [0.0014, 0.0017])
     endif else begin
        am = plot( (AORsd1043[a].timearr - AORsd1043[0].timearr(0))/60./60., AORsd1043[a].flux,'6r1s', sym_size = 0.2,   sym_filled = 3, color = colorarr[a],/overplot)
     endelse

  endfor

;------

  for a = 0, n_elements(aorname_sd1043) - 1 do begin

      ;remove outliers
     good = where(AORsd1043[a].xcen lt mean(AORsd1043[a].xcen) + 2.5*stddev(AORsd1043[a].xcen) and AORsd1043[a].xcen gt mean(AORsd1043[a].xcen) -2.5*stddev(AORsd1043[a].xcen) and AORsd1043[a].ycen lt mean(AORsd1043[a].ycen) +3.0*stddev(AORsd1043[a].ycen) and AORsd1043[a].ycen gt mean(AORsd1043[a].ycen) - 3.0*stddev(AORsd1043[a].ycen),ngood_pmap, complement=bad) 
     print, 'bad ',n_elements(bad)

     xarr = AORsd1043[a].xcen[good]
     yarr = AORsd1043[a].ycen[good]
     timearr = AORsd1043[a].timearr[good]
     flux = AORsd1043[a].flux[good]
     corrflux = AORsd1043[a].corrflux[good]
     corrfluxerr = AORsd1043[a].corrfluxerr[good]
     bmjdarr = AORsd1043[a].bmjdarr[good]
     bkgdarr = AORsd1043[a].bkgd[good]
     bkgderrarr = AORsd1043[a].bkgderr[good]
     
; binning
     bin_level = 5L;*60
     nframes = bin_level
     nfits = long(n_elements(corrflux)) / nframes
     print, 'nframes, nfits', n_elements(corrflux), nframes, nfits
     bin_flux = dblarr(nfits)
     bin_corrflux= dblarr(nfits)
     bin_corrfluxerr= dblarr(nfits)
     bin_timearr = dblarr(nfits)
     bin_bmjdarr = dblarr(nfits)
     bin_bkgd = dblarr(nfits)
     bin_bkgderr = dblarr(nfits)
     bin_xcen = dblarr(nfits)
     bin_ycen = dblarr(nfits)

     for si = 0L, long(nfits) - 1L do begin
                                ;print, 'working on si', n_elements(AORsd1043[a].corrflux), si, si*nframes, si*nframes + (nframes - 1)
        idata = corrflux[si*nframes:si*nframes + (nframes - 1)]
        id = flux[si*nframes:si*nframes + (nframes - 1)]
        idataerr = corrfluxerr[si*nframes:si*nframes + (nframes - 1)]
        ib = bkgdarr[si*nframes:si*nframes + (nframes - 1)]
        iberr = bkgderrarr[si*nframes:si*nframes + (nframes - 1)]
       ; if si lt 2 then print, 'ib', ib
       ; print, 'mean', mean(ib,/nan)
        bin_flux[si] = mean(id, /nan)
        bin_corrflux[si] = mean(idata,/nan)
        bin_corrfluxerr[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
        bin_timearr[si]= mean(timearr[si*nframes:si*nframes + (nframes - 1)])
        bin_bmjdarr[si]= mean(bmjdarr[si*nframes:si*nframes + (nframes - 1)])
        bin_bkgd[si] = mean(ib,/nan)
        bin_bkgderr[si] = mean(iberr,/nan)
        bin_xcen[si]= mean(xarr[si*nframes:si*nframes + (nframes - 1)])
        bin_ycen[si]= mean(yarr[si*nframes:si*nframes + (nframes - 1)])

     endfor                     ;for each fits image
     
 ;try getting rid of flux outliers.
  ;do some running mean with clipping
     start = 0
     print, 'nflux', n_elements(bin_flux)
     for ni = 50, n_elements(bin_flux) -1, 50 do begin
        meanclip,flux[start:ni], m, s, subs = subs ;,/verbose
                                ;print, 'good', subs+start
                                ;now keep a list of the good ones
        if ni eq 50 then good_ni = subs+start else good_ni = [good_ni, subs+start]
        start = ni
     endfor
                                ;see if it worked
     ;xarr = xarr[good_ni]
     ;yarr = yarr[good_ni]
     bin_timearr = bin_timearr[good_ni]
     bin_bmjdarr = bin_bmjdarr[good_ni]
     bin_corrfluxerr = bin_corrfluxerr[good_ni]
     bin_flux = bin_flux[good_ni]
     bin_corrflux = bin_corrflux[good_ni]
     bin_bkgd = bin_bkgd[good_ni]
     bin_bkgderr = bin_bkgderr[good_ni]
     bin_xcen = bin_xcen[good_ni]
     bin_ycen = bin_ycen[good_ni]

     print, 'nflux', n_elements(bin_flux)
     
 xcen_coeff = [1.5,0.5,-0.1, -0.8]
 ycen_coeff = [1.0,0,0,-0.7]
     if a eq 0 then begin
;        am = plot( (AORsd1043[a].timearr - AORsd1043[a].timearr(0))/60./60., AORsd1043[a].flux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14], xtitle = 'Time(hrs)', ytitle = 'Flux', xrange = [0,6], yrange = [3,5])
;        am = plot( (AORsd1043[a].timearr - AORsd1043[a].timearr(0))/60./60., AORsd1043[a].corrflux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black', /overplot)
        
;        am = plot(bin_bmjdarr, bin_corrflux,  '1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14], xtitle = 'Phase', ytitle = 'Corrected Flux')

;        ab = plot(bin_bmjdarr  , bin_flux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14], xtitle = 'period', ytitle = 'Flux (Jy)')
;        ab = plot((bin_timearr - bin_timearr(0))/60./60., bin_xcen ,'6r1s', sym_size = 0.3,   sym_filled = 1, color = colorarr[a],  xtitle = 'Period', ytitle = 'Xcen');, xrange = [0.3,0.65])
  
;        ab = plot(bin_bmjdarr, bin_flux/mean(bin_flux), '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], title = 'aperture [3,11,15.5]', xtitle = 'Period', ytitle = 'Flux');, xrange = [0.3,0.65])
;        ab = plot(bin_bmjdarr, bin_bkgd, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], title = 'aperture [3,11,15.5]', xtitle = 'Period', ytitle = 'Bkgd');, xrange = [0.3,0.65])


   ;  print,'time', bin_bmjdarr - bin_bmjdarr(0) 
     endif else begin
     ; am = plot( (AORsd1043[a].timearr - AORsd1043[a].timearr(0))/60./60., AORsd1043[a].flux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14],/overplot)
; am = plot( (AORsd1043[a].timearr - AORsd1043[a].timearr(0))/60./60., AORsd1043[a].corrflux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black', /overplot)
;        am = plot(bin_bmjdarr, bin_corrflux- 1.5E-2*a,  '1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14],/overplot)

;        am = plot((bin_timearr - bin_timearr(0))/60./60., bin_xcen ,'1s', sym_size = 0.3,   sym_filled = 1, color = colorarr[a], /overplot)
;        am = plot(bin_bmjdarr, bin_bkgd-0.00005*(a),'1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], /overplot)

     endelse


  endfor

end
