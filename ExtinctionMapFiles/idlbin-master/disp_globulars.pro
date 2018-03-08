pro disp_globulars

  restore,  '/Users/jkrick/Virgo/m87/gclus.sav'
  
a = where(gclus.ch1_mag gt 0 and gclus.ch1_mag lt 90 and gclus.ch2_mag gt 0 and gclus.ch2_mag lt 90 and gclus.ch1_magerr lt 0.4 and gclus.ch2_magerr lt 0.5)
color = gclus[a].ch1_mag - gclus[a].ch2_mag
plothist, color, xarr, yarr, bin = 0.1, xtitle = '[3.6 - 4.5] microns', ytitle = 'number'

p = plot(xarr, yarr,  xtitle = '[3.6 - 4.5] microns', ytitle = 'number')

;print, gclus[a].ch1_magerr
;print, gclus[a].ch2_magerr

end
