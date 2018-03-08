pro plot_CR

;add CR to this plot from Joe's house keeping statistics
  readcol, '/Users/jkrick/irac_warm/calstars/pixel_stats_warm_2014.csv',time, MJD, date, hot_ch1, noisy_ch1, dead_ch1, cr_ch1, hot_ch2, noisy_ch2, dead_ch2, cr_ch2, skipline = 6
  
  pcr = plot((MJD - min(MJD))/365., cr_ch1 , '1o', sym_filled =1, color = 'green', xtitle = 'Time(years) of Warm Mission', ytitle = 'CR', yrange = [2, 8])
  pcr = plot((MJD- min(MJD))/365., cr_ch2, '1s',sym_filled = 1, color = 'blue', overplot = pcr)
  tt1 = text(1.0, 3.0,'Ch2', /data, color = 'blue', font_style = 'bold')
  tt1 = text(1.0, 2.6,'Ch1', /data, color = 'green', font_style = 'bold')
;  pcr = plot(cr_ch1, cr_ch2, '1s', color = 'black', xtitle = 'CR ch1', ytitle = 'CR ch2', xrange = [2, 8], yrange = [2,8])
end
