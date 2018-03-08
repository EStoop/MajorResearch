pro wasp33_compare_plots

  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/wasp33_phot_ch2.sav'
  aorname = ['r45383424', 'r45384448', 'r45384704'] ;ch2

;------
;first plot the bcd fluxes, uncorrected
;  for a = 0, n_elements(aorname) - 1 do begin
;     if a eq 0 then begin
;        p = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].flux,'61s', sym_size = 0.1,   sym_filled = 1, color ='black', xtitle = 'Time(hrs)', ytitle = 'Flux (Jy)', xrange = [0,40], yrange =[0.17,0.18])
;        p = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].corrflux,'61s', sym_size = 0.1,   sym_filled = 1, color = 'blue',/overplot)
;     endif else begin
;        p = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].flux,'61s', sym_size = 0.1,   sym_filled = 1, color ='black',/overplot)
;       p = plot( (AORwasp33[a].timearr - AORwasp33[0].timearr(0))/60./60., AORwasp33[a].corrflux,'61s', sym_size = 0.1,   sym_filled = 1, color ='blue',/overplot)
;     endelse

 ; endfor

restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/bin_timearr.sav'
restore,  '/Users/jkrick/irac_warm/pcrs_planets/wasp33/bin_flux.sav'
restore,  '/Users/jkrick/irac_warm/pcrs_planets/wasp33/bin_corrflux.sav'

     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_flux, '61s', sym_size = 0.1,   sym_filled = 1, xrange = [0,40], yrange = [0.166, 0.178], xtitle = 'Time (hrs)', ytitle = 'Flux (Jy)', title = 'wasp33')
     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_corrflux-0.002, '6r1s', sym_size = 0.1,   sym_filled = 1,/overplot)

;------
;then plot the self calibrated corrected fluxes

  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/selfcal_timearr.sav'
  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/selfcal_flux.sav'
   restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/bin_time.sav'
  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/bin_sub.sav'

  p2 = plot(bin_time, bin_sub+0.168, '6rs1', sym_size = 0.1,   sym_filled = 1, color = 'light_sky_blue',/overplot)
   t1 = text(8, 0.1685, 'Self-Calibrated', color ='light_sky_blue',/data)
  t1 = text(8, 0.1725, 'Pmap Corrected', color ='red',/data)
  t1 = text(8, 0.177, 'Original', color ='black',/data)


  end
