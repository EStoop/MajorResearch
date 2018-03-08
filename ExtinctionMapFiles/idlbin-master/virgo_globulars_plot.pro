pro virgo_globulars_plot
  restore, '/Users/jkrick/Virgo/IRAC/globulars/gclus_back.sav'


;plot a histogram of ch1 - ch2
;------------------------
  good = where(gclus.ch2_mag gt 0 and gclus.ch2_mag lt 90 and gclus.ch1_mag gt 0 and gclus.ch1_mag lt 90 and finite(gclus.ch2_mag) eq 1 and finite(gclus.ch1_mag) eq 1 ,countgood) ;and gclus.flag_confuse lt 1
  print, countgood
  plothist, gclus[good].ch1_mag - gclus[good].ch2_mag, xarr, yarr, bin = 0.05, /noplot
  d = barplot(xarr, yarr, xtitle='ch1 - ch2 (AB mag)', ytitle='Number', fill_color = 'blue', xrange = [-2, 2])

  start = [-2.5, 5, 30000]
  noise = fltarr(n_elements(xarr)) + 1.
  P= MPFITFUN('mygauss',xarr,yarr, noise, start) ;ICL
  d = plot( xarr, P(2)/sqrt(2.*!Pi) * exp(-0.5*((xarr - P(0))/P(1))^2.), '6r2', /overplot)


;what about a color mag diagram
;--------------------------
c = plot(gclus[good].ch1_mag - gclus[good].ch2_mag, gclus[good].ch1_mag,  '6r1.', xtitle ='ch1 - ch2 (AB mag)', ytitle='ch1 mag AB', title = 'Virgo',xrange = [-8, 5], yrange = [15, 27])

good  = where(finite(gclus.ch1_flux) eq 1 and finite(gclus.ch2_flux) eq 1)
;c2 = plot(gclus[good ].ch1_flux*1E6 - gclus[good ].ch2_flux*1E6, gclus[good ].ch1_flux*1E6,  '6r1.', /ylog, /xlog, xrange = [1E-4,1E4], yrange = [1E-4,1E4],xtitle ='ch1 - ch2 (flux)', ytitle='ch1 flux', title = 'Virgo')


c3 = plot(gclus[good].ch2_flux*1E6, gclus[good].ch1_flux*1E6,  '6r1.', /ylog, /xlog, xrange = [1E-4,1E4], yrange = [1E-4,1E4],xtitle ='ch2 (flux)', ytitle='ch1 flux', title = 'Virgo')

;try an optical to irac color mag diagram for optically bright sources
;-------------------------------
good = where(gclus.mag_auto[0] gt 0 and gclus.mag_auto[0] lt 28 and finite(gclus.mag_auto[0]) eq 1 and finite(gclus.ch1_mag) eq 1 and gclus.ch1_mag gt 0 and gclus.ch1_mag lt 90 and gclus.flag_confuse lt 1,countgood)

d = plot(gclus[good].mag_auto[0] - (gclus[good].ch2_mag ), gclus[good].mag_auto[0], '6r1.', xtitle ='r - ch2 (AB )', ytitle='r mag (AB)', yrange = [12, 28], xrange = [-4, 8], title = "Virgo")




end
