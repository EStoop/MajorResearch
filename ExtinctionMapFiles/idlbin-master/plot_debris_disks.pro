pro plot_debris_disks
;make plots for the senior review 2016
;----------------------------------------------------
 ;;make Eta Corvi plot from Marengo preliminary data
  readcol, '/Users/jkrick/external/irac_warm/senior_review_2016/marengo.txt', year, month, day, ch1excess, ch1err, ch2excess, ch2err, format = '(I10, I10, I10, F10, F10, F10, F10)'

  jd_arr = julday( month, day, year)
  print, 'jd_arr', jd_arr
  junk = LABEL_DATE(Date_format = ['%D-%M-%Y'])
;;  p1 = errorplot(jd_arr, ch1excess/ch1excess(0), ch1err/ch1excess(0), '1o',  ytitle = 'Percentage Increase in Disk Excess', color = 'blue', errorbar_color = 'blue', sym_filled = 1, yrange = [0.8, 1.5],   xrange = [jd_arr(0) - 20, jd_arr(n_elements(jd_arr)-1) + 20], XTICKFORMAT='(C(CDI,1x,CMoA,1x,CYI))', xtickunits = 'years',xminor = 11,name = '3.6micron') ;xtickformat = 'LABEL_DATE', xstyle = 1,xtext_orientation= 45,xtitle = 'Time',

;;  p2 = errorplot(jd_arr, ch2excess/ch2excess(0), ch2err/ch2excess(0), '1o', color = 'green', errorbar_color = 'green', sym_filled = 1, overplot = p1, name = '4.5micron')



  ;;or maybe plot only the times after 2015, with a line which
  ;;includes the average value before march 2013
  ch1normfactor = mean(ch1excess[0:2])
  ch2normfactor = mean(ch2excess[0:1])

  p1 = errorplot(jd_arr, ch1excess/ch1normfactor, ch1err/ch1normfactor, '1o',  ytitle = 'Percentage Increase in Disk Excess', color = 'blue', errorbar_color = 'blue', sym_filled = 1, yrange = [0.8, 1.5],   xrange = [jd_arr(3) - 10, jd_arr(n_elements(jd_arr)-1) + 10], XTICKFORMAT='(C(CDI,1x,CMoA,1x,CYI))', xtickunits = 'months',xminor = 10,name = '3.6micron')
  p2 = errorplot(jd_arr, ch2excess/ch2normfactor, ch2err/ch2normfactor, '1o', color = 'green', errorbar_color = 'green', sym_filled = 1, overplot = p1, name = '4.5micron')
  p3 = plot(jd_arr + 20., fltarr(n_elements(jd_arr)) + 1.0, thick = 2, overplot = p1)

  leg = legend(target = [p1, p2], position = [2457145, 0.95], /data, /auto_text_color)



;-----------------------------------------------------------
  ;;make plot from Meng et al.2014, 2015 for ID8 and P1121.

  readcol, '/Users/jkrick/external/irac_warm/senior_review_2016/P1121.phot.ch1.20150820.txt', aorkey, BMJD,  starflux, starflux_err, flux, flux_err, format = '(I10,D10, F10, F10, F10, F10)'
  jd = BMJD + 2400000.5
  q1 = errorplot(jd, flux/flux(0), flux_err, '1o', sym_filled = 1, XTICKFORMAT='(C(CDI,1x,CMoA,1x,CYI))', xtickunits = 'years',xminor = 11,name = 'P1121 3.6$\mu$m', ytitle = 'Normalized Disk Flux Density(mJy)', color=[166,206,227], errorbar_color = [166,206,227], errorbar_capsize = 0.1, xrange = [jd(0) -30, jd(n_elements(jd) - 1) + 30], sym_size = 0.5, yrange = [0, 3.4], xthick = 2, ythick = 2, font_style = 'bf')
  readcol, '/Users/jkrick/external/irac_warm/senior_review_2016/P1121.phot.ch2.20150820.txt', aorkey, BMJD,  starflux, starflux_err, flux, flux_err, format = '(I10,D10, F10, F10, F10, F10)'
  jd = BMJD + 2400000.5
  q2 = errorplot(jd, flux/flux(0), flux_err, '1o', sym_filled = 1, name = 'P1121 4.5$\mu$m',  color=[31,120,180], errorbar_color = [31,120,180], errorbar_capsize = 0.1, overplot = q1, sym_size = 0.5)

;----

    readcol, '/Users/jkrick/external/irac_warm/senior_review_2016/id8.phot.ch1.20150820.txt', aorkey, BMJD,  starflux, starflux_err, flux, flux_err, format = '(I10,D10, F10, F10, F10, F10)'
  jd = BMJD + 2400000.5
  q3 = errorplot(jd, flux/flux(0) , flux_err, '1o', sym_filled = 1, name = 'ID8 3.6$\mu$m',  color=[178,223,138], errorbar_color = [178,223,138],  errorbar_capsize = 0.1,overplot = q1, sym_size = 0.5)
  readcol, '/Users/jkrick/external/irac_warm/senior_review_2016/id8.phot.ch2.20150820.txt', aorkey, BMJD,  starflux, starflux_err, flux, flux_err, format = '(I10,D10, F10, F10, F10, F10)'
  jd = BMJD + 2400000.5
  q4 = errorplot(jd, flux/flux(0), flux_err, '1o', sym_filled = 1, name = 'ID8 4.5$\mu$m',  color=[51,160,44], errorbar_color = [51,160,44], errorbar_capsize = 0.1, overplot = q1, sym_size = 0.5)

    leg = legend(target = [q3,q4,q1,q2], position = [2456455, 3.2], /data, /auto_text_color)


end
