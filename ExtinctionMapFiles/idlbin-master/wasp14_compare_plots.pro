pro wasp14_compare_plots

  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/wasp14_phot_ch2.sav'

;------
;first plot the bcd fluxes, uncorrected
;  for a = 0, n_elements(aorname) - 1 do begin
;     if a eq 0 then begin
;        p = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].flux,'61s', sym_size = 0.1,   sym_filled = 1, color ='black', xtitle = 'Time(hrs)', ytitle = 'Flux (Jy)', xrange = [0,40], yrange =[0.17,0.18])
;        p = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].corrflux,'61s', sym_size = 0.1,   sym_filled = 1, color = 'blue',/overplot)
;     endif else begin
;        p = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].flux,'61s', sym_size = 0.1,   sym_filled = 1, color ='black',/overplot)
;       p = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].corrflux,'61s', sym_size = 0.1,   sym_filled = 1, color ='blue',/overplot)
;     endelse

 ; endfor

restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_timearr.sav'
restore,  '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_flux.sav'
restore,  '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_corrflux.sav'
print, 'nele', n_elements(bin_timearr), n_elements(bin_flux), n_elements(bin_corrflux)
     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_flux/mean(bin_flux), '61s', sym_size = 0.2,   sym_filled = 12, xrange = [0,65], xtitle = 'Time (hrs)', ytitle = 'Normalized Flux', yrange = [0.97,1.01], title = 'Wasp 14 Ch2')
     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_corrflux/median(bin_corrflux) - 0.012, '6r1s', sym_size = 0.2,   sym_filled = 1,/overplot)

;------
;then plot the self calibrated corrected fluxes

  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/selfcal_timearr.sav'
  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/selfcal_flux.sav'
   restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_time.sav'
  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_sub.sav'
  p2 = plot(bin_time,((bin_sub+0.06) / mean(bin_sub +0.06)) - 0.024,'6bs1', color = 'light_sky_blue', sym_size = 0.2,   sym_filled = 1, /overplot)
  
  t1 = text(15, 0.978, 'Self-Calibrated', color ='light_sky_blue',/data)
  t1 = text(15, 0.99, 'Pmap Corrected', color ='red',/data)
  t1 = text(15, 1.003, 'Original', color ='black',/data)
  ;t1 = text([15, 15, 15], [0.978, 0.99, 1.003],/data,[ 'Self-Calibrated','Pmap Corrected', 'Original'], color = ['light_sky_blue', 'red','black'])




;----------------------------------------------------------------------
;CH1
;----------------------------------------------------------------------

restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/wasp14_phot_ch1.sav'

;------
;first plot the bcd fluxes, uncorrected
;  for a = 0, n_elements(aorname) - 1 do begin
;     if a eq 0 then begin
;        p = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].flux,'61s', sym_size = 0.1,   sym_filled = 1, color ='black', xtitle = 'Time(hrs)', ytitle = 'Flux (Jy)', xrange = [0,40], yrange =[0.17,0.18])
;        p = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].corrflux,'61s', sym_size = 0.1,   sym_filled = 1, color = 'blue',/overplot)
;     endif else begin
;        p = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].flux,'61s', sym_size = 0.1,   sym_filled = 1, color ='black',/overplot)
;       p = plot( (AORwasp14[a].timearr - AORwasp14[0].timearr(0))/60./60., AORwasp14[a].corrflux,'61s', sym_size = 0.1,   sym_filled = 1, color ='blue',/overplot)
;     endelse

 ; endfor

restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_timearr_ch1.sav'
restore,  '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_flux_ch1.sav'
restore,  '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_corrflux_ch1.sav'
print, 'nele', n_elements(bin_timearr), n_elements(bin_flux), n_elements(bin_corrflux)
     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_flux/mean(bin_flux), '61s', sym_size = 0.2,   sym_filled = 12, xrange = [0,65], xtitle = 'Time (hrs)', ytitle = 'Normalized Flux', yrange = [0.93,1.02], title = 'Wasp 14 Ch1')
     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_corrflux/median(bin_corrflux) - 0.035, '6r1s', sym_size = 0.2,   sym_filled = 1,/overplot)

;------
;then plot the self calibrated corrected fluxes

  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/selfcal_timearr_ch1.sav'
  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/selfcal_flux_ch1.sav'
   restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_time_ch1.sav'
  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_sub_ch1.sav'
  p2 = plot(bin_time,((bin_sub+0.06) / mean(bin_sub +0.06)) - 0.06,'6bs1', color = 'light_sky_blue', sym_size = 0.2,   sym_filled = 1, /overplot)
  
  t1 = text(15, 0.948, 'Self-Calibrated', color ='light_sky_blue',/data)
  t1 = text(15, 0.97, 'Pmap Corrected', color ='red',/data)
  t1 = text(15, 1.015, 'Original', color ='black',/data)
  ;t1 = text([15, 15, 15], [0.978, 0.99, 1.003],/data,[ 'Self-Calibrated','Pmap Corrected', 'Original'], color = ['light_sky_blue', 'red','black'])



 end
