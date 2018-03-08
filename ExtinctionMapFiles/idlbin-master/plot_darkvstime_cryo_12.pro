pro plot_darkvstime
ps_open, filename='/Users/jkrick/irac_darks/darksvstime_cryo_12.ps',/portrait,/square,/color
;ps_start, filename='/Users/jkrick/irac_darks/darksvstime_3.ps'
!P.thick = 3
!P.charthick = 3
!P.multi = [0,1,4]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)

;zodi from header not removed from the data
restore, filename = '/Users/jkrick/irac_darks/cryodarks_12_1.sav'
cryo_nozodi_ch1_sclksort = ch1_sclksort
cryo_nozodi_ch1_meansort = ch1_meansort
cryo_nozodi_ch1_sigmasort = ch1_sigmasort
cryo_nozodi_ch2_sclksort = ch2_sclksort
cryo_nozodi_ch2_meansort = ch2_meansort
cryo_nozodi_ch2_sigmasort = ch2_sigmasort
cryo_nozodi_ch3_sclksort = ch3_sclksort
cryo_nozodi_ch3_meansort = ch3_meansort
cryo_nozodi_ch3_sigmasort = ch3_sigmasort
cryo_nozodi_ch4_sclksort = ch4_sclksort
cryo_nozodi_ch4_meansort = ch4_meansort
cryo_nozodi_ch4_sigmasort = ch4_sigmasort

;zodi in the header has been subtracted 
restore, filename = '/Users/jkrick/irac_darks/cryodarks_12_2.sav'
cryo_ch1_sclksort = ch1_sclksort
cryo_ch1_meansort = ch1_meansort
cryo_ch1_sigmasort = ch1_sigmasort
cryo_ch2_sclksort = ch2_sclksort
cryo_ch2_meansort = ch2_meansort
cryo_ch2_sigmasort = ch2_sigmasort
cryo_ch3_sclksort = ch3_sclksort
cryo_ch3_meansort = ch3_meansort
cryo_ch3_sigmasort = ch3_sigmasort
cryo_ch4_sclksort = ch4_sclksort
cryo_ch4_meansort = ch4_meansort
cryo_ch4_sigmasort = ch4_sigmasort

;the warm dark measurements
;zodi from header not removed from the data
restore, filename = '/Users/jkrick/irac_darks/alldarks_12_4.sav'
nozodi_ch1_sclksort = ch1_sclksort
nozodi_ch1_meansort = ch1_meansort
nozodi_ch1_sigmasort = ch1_sigmasort
nozodi_ch2_sclksort = ch2_sclksort
nozodi_ch2_meansort = ch2_meansort
nozodi_ch2_sigmasort = ch2_sigmasort

;zodi in the header has been subtracted 
restore, filename = '/Users/jkrick/irac_darks/alldarks_12_3.sav'
;-------------------------------------------------------------------
;need to combine cryo and warm arrays.
nozodi_ch1_sclksort_all = [cryo_nozodi_ch1_sclksort,nozodi_ch1_sclksort]
nozodi_ch1_meansort_all = [cryo_nozodi_ch1_meansort,nozodi_ch1_meansort]
nozodi_ch1_sigmasort_all = [cryo_nozodi_ch1_sigmasort,nozodi_ch1_sigmasort]
nozodi_ch2_sclksort_all = [cryo_nozodi_ch2_sclksort,nozodi_ch2_sclksort]
nozodi_ch2_meansort_all = [cryo_nozodi_ch2_meansort,nozodi_ch2_meansort]
nozodi_ch2_sigmasort_all = [cryo_nozodi_ch2_sigmasort,nozodi_ch2_sigmasort]

ch1_sclksort_all = [cryo_ch1_sclksort,ch1_sclksort]
ch1_meansort_all = [cryo_ch1_meansort,ch1_meansort]
ch1_sigmasort_all = [cryo_ch1_sigmasort,ch1_sigmasort]
ch2_sclksort_all = [cryo_ch2_sclksort,ch2_sclksort]
ch2_meansort_all = [cryo_ch2_meansort,ch2_meansort]
ch2_sigmasort_all = [cryo_ch2_sigmasort,ch2_sigmasort]

;-------------------------------------------------------------------
plot, (ch1_sclksort_all - ch1_sclksort_all(0)) / 86400 / 365, ch1_meansort_all, psym = 2, ytitle = 'Mean flux', xtitle = 'Time in years starting from Dec 01 2003', title = 'Channel 1 F12s', yrange = [0.01, 3.5], ystyle = 1

oplot, (nozodi_ch1_sclksort_all - nozodi_ch1_sclksort_all(0)) / 86400 / 365, nozodi_ch1_meansort_all, psym = 4

;-----
plot, (ch2_sclksort_all - ch2_sclksort_all(0)) / 86400 / 365, ch2_meansort_all ,  psym = 2, ytitle = 'Mean flux', xtitle = 'Time in years starting from Dec 01 2003', title = 'Channel 2 F12s', yrange = [-0.2,0.6], ystyle = 1

oplot, (nozodi_ch2_sclksort_all - nozodi_ch2_sclksort_all(0)) / 86400 / 365, nozodi_ch2_meansort_all ,  psym = 4
;-----

plot, (cryo_ch3_sclksort - cryo_ch3_sclksort(0)) / 86400 / 365, cryo_ch3_meansort ,  psym = 2, ytitle = 'Mean flux', xtitle = 'Time in years starting from Dec 01 2003', title = 'Channel 3 F12s', yrange = [-4.0,2.0], ystyle = 1

oplot, (cryo_nozodi_ch3_sclksort - cryo_nozodi_ch3_sclksort(0)) / 86400 / 365, cryo_nozodi_ch3_meansort ,  psym = 4
;-----

plot, (cryo_ch4_sclksort - cryo_ch4_sclksort(0)) / 86400 / 365, cryo_ch4_meansort ,  psym = 2, ytitle = 'Mean flux', xtitle = 'Time in years starting from Dec 01 2003', title = 'Channel 4 F12s', yrange = [-1.0,8.0], ystyle = 1

oplot, (cryo_nozodi_ch4_sclksort - cryo_nozodi_ch4_sclksort(0)) / 86400 / 365, cryo_nozodi_ch4_meansort ,  psym = 4

;-------------------------------------------------------
;plot the normalized versions
plot, (ch1_sclksort_all - ch1_sclksort_all(0)) / 86400 / 365, ch1_meansort_all/ ch1_meansort_all(1) , psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in years starting from Dec 01 2003', title = 'Channel 1 F12s', yrange = [0.0,3.0]
oplot, (nozodi_ch1_sclksort_all - nozodi_ch1_sclksort_all(0)) / 86400 / 365, nozodi_ch1_meansort_all/ nozodi_ch1_meansort_all(1) , psym = 4
a = where((ch1_sclksort_all - ch1_sclksort_all(0)) / 86400 / 365 gt 5.5)
oplot, (ch1_sclksort_all(a) - ch1_sclksort_all(0)) / 86400 / 365, ch1_meansort_all(a)/ ch1_meansort_all(a(0)) , psym = 2
oplot, (ch1_sclksort_all(a) - ch1_sclksort_all(0)) / 86400 / 365, nozodi_ch1_meansort_all(a)/ nozodi_ch1_meansort_all(a(0)) , psym = 4

plot, (ch2_sclksort_all - ch2_sclksort_all(0)) / 86400 / 365, ch2_meansort_all / ch2_meansort_all(1),  psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in years starting from Dec 01 2003', title = 'Channel 2 F12s', yrange = [0.0,2.4]
oplot, (nozodi_ch2_sclksort_all - nozodi_ch2_sclksort_all(0)) / 86400 / 365, nozodi_ch2_meansort_all / nozodi_ch2_meansort_all(1),  psym = 4
a = where((ch2_sclksort_all - ch2_sclksort_all(0)) / 86400 / 365 gt 5.5)
oplot, (ch2_sclksort_all(a) - ch2_sclksort_all(0)) / 86400 / 365, ch2_meansort_all(a)/ ch2_meansort_all(a(0)) , psym = 2
oplot, (ch2_sclksort_all(a) - ch2_sclksort_all(0)) / 86400 / 365, nozodi_ch2_meansort_all(a)/ nozodi_ch2_meansort_all(a(0)) , psym = 4

plot, (cryo_ch3_sclksort - cryo_ch3_sclksort(0)) / 86400 / 365, cryo_ch3_meansort / cryo_ch3_meansort(1),  psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in years starting from Dec 01 2003', title = 'Channel 3 F12s', yrange = [0.5,3.5], ystyle = 1
oplot, (cryo_nozodi_ch3_sclksort - cryo_nozodi_ch3_sclksort(0)) / 86400 / 365, cryo_nozodi_ch3_meansort / cryo_nozodi_ch3_meansort(1),  psym = 4

plot, (cryo_ch4_sclksort - cryo_ch4_sclksort(0)) / 86400 / 365, cryo_ch4_meansort / cryo_ch4_meansort(1),  psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in years starting from Dec 01 2003', title = 'Channel 4 F12s', yrange = [0.0,2.0]
oplot, (cryo_nozodi_ch4_sclksort - cryo_nozodi_ch4_sclksort(0)) / 86400 / 365, cryo_nozodi_ch4_meansort / cryo_nozodi_ch4_meansort(1),  psym = 4

;-------------------------------------------------------

;fit with sin curve

;noise = fltarr(n_elements(ch2_meansort_all))
noise = ch2_sigmasort_all / ch2_meansort_all(1)

start = [0.03,1.03,60.,3.1]
x =  (ch2_sclksort_all - ch2_sclksort_all(0)) / 86400 / 365
result= MPFITFUN('sin_func',x, ch2_meansort_all / ch2_meansort_all(1), noise, start)    ;ICL
;oplot, x,  result(0)*sin(X/result(2) + result(3)) + result(1), thick = 3;,psym = 2;, color = colors.green

period = 2 * !PI * result(2)
print, 'period', period

;xyouts, 50, 0.92, strcompress('Period' + string(round(period))+ ' days')
;xyouts, 350, 0.92, strcompress('Amplitude' + strmid(string(fix(1000*result(0)) / 10.),0,9)+ '%')
;--------------------------------------------------------------
;try plotting colors
;don't use the normalized means
!P.multi = [0,1,2]

ch1flux_zodi = ch1_meansort_all;/ ch1_meansort_all(1) 
ch1flux_nozodi = nozodi_ch1_meansort_all;/ nozodi_ch1_meansort_all(1)
ch2flux_zodi = ch2_meansort_all;/ ch2_meansort_all(1) 
ch2flux_nozodi = nozodi_ch2_meansort_all;/ nozodi_ch2_meansort_all(1)

plot,  (ch1_sclksort_all - ch1_sclksort_all(0)) / 86400 / 365, ch1flux_nozodi - ch2flux_nozodi, psym = 4, ytitle = 'Ch1 - Ch2', xtitle = 'Time in years starting from Dec 01 2003', title = ' F12s', yrange = [-0.4, 5.2]
oplot,  (ch2_sclksort_all - ch2_sclksort_all(0)) / 86400 / 365, ch1flux_zodi - ch2flux_zodi, psym = 2

plot,  (cryo_ch2_sclksort - cryo_ch2_sclksort(0)) / 86400 / 365, cryo_nozodi_ch2_meansort - cryo_nozodi_ch4_meansort, psym = 4, ytitle = 'Ch2 - Ch4', xtitle = 'Time in years starting from Dec 01 2003', title = ' F12s', yrange = [-7.5, 1.0], ystyle = 1
oplot, (cryo_ch2_sclksort - cryo_ch2_sclksort(0)) / 86400 / 365, cryo_ch2_meansort - cryo_ch4_meansort, psym = 2

;--------------------------------------------------------------

ps_close, /noprint,/noid
;ps_end;, /png



end
