pro plot_darkvstime
ps_open, filename='/Users/jkrick/irac_darks/darksvstime_cryo.ps',/portrait,/square,/color
;ps_start, filename='/Users/jkrick/irac_darks/darksvstime_3.ps'
!P.thick = 3
!P.charthick = 3
!P.multi = [0,1,4]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)

;zodi from header not removed from the data
restore, filename = '/Users/jkrick/irac_darks/cryodarks_1.sav'
nozodi_ch1_sclksort = ch1_sclksort
nozodi_ch1_meansort = ch1_meansort
nozodi_ch1_sigmasort = ch1_sigmasort
nozodi_ch2_sclksort = ch2_sclksort
nozodi_ch2_meansort = ch2_meansort
nozodi_ch2_sigmasort = ch2_sigmasort
nozodi_ch3_sclksort = ch3_sclksort
nozodi_ch3_meansort = ch3_meansort
nozodi_ch3_sigmasort = ch3_sigmasort
nozodi_ch4_sclksort = ch4_sclksort
nozodi_ch4_meansort = ch4_meansort
nozodi_ch4_sigmasort = ch4_sigmasort

;zodi in the header has been subtracted 
restore, filename = '/Users/jkrick/irac_darks/cryodarks_2.sav'

plot, (ch1_sclksort - ch1_sclksort(0)) / 86400, ch1_meansort, psym = 2, ytitle = 'Mean flux', xtitle = 'Time in days starting from Dec 01 2003', title = 'Channel 1 F30s', yrange = [0.0, 0.08], ystyle = 1
;ploterror, (ch1_sclksort - ch1_sclksort(0)) / 86400, ch1_meansort/ ch1_meansort(1) , ch1_sigmasort/ ch1_meansort(1), psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in days starting from pc5', title = 'Channel 1 F30s', yrange = [0.96,1.04]

oplot, (nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400, nozodi_ch1_meansort, psym = 4

;-----
plot, (ch2_sclksort - ch2_sclksort(0)) / 86400, ch2_meansort ,  psym = 2, ytitle = 'Mean flux', xtitle = 'Time in days starting from Dec 01 2003', title = 'Channel 2 F30s', yrange = [0.0,0.4], ystyle = 1
;ploterror, (ch2_sclksort - ch2_sclksort(0)) / 86400, ch2_meansort / ch2_meansort(1),  ch2_sigmasort/ ch2_meansort(1), psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in days starting from pc5', title = 'Channel 2 F30s', yrange = [0.9,1.2]

oplot, (nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400, nozodi_ch2_meansort ,  psym = 4
;-----

plot, (ch3_sclksort - ch3_sclksort(0)) / 86400, ch3_meansort ,  psym = 2, ytitle = 'Mean flux', xtitle = 'Time in days starting from Dec 01 2003', title = 'Channel 3 F30s', yrange = [-4.0,2.0], ystyle = 1

oplot, (nozodi_ch3_sclksort - nozodi_ch3_sclksort(0)) / 86400, nozodi_ch3_meansort ,  psym = 4
;-----

plot, (ch4_sclksort - ch4_sclksort(0)) / 86400, ch4_meansort ,  psym = 2, ytitle = 'Mean flux', xtitle = 'Time in days starting from Dec 01 2003', title = 'Channel 4 F30s', yrange = [-1.0,8.0], ystyle = 1

oplot, (nozodi_ch4_sclksort - nozodi_ch4_sclksort(0)) / 86400, nozodi_ch4_meansort ,  psym = 4

;-------------------------------------------------------
;plot the normalized versions
plot, (ch1_sclksort - ch1_sclksort(0)) / 86400, ch1_meansort/ ch1_meansort(1) , psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in days starting from Dec 01 2003', title = 'Channel 1 F30s', yrange = [0.5,2.5]
oplot, (nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400, nozodi_ch1_meansort/ nozodi_ch1_meansort(1) , psym = 4

plot, (ch2_sclksort - ch2_sclksort(0)) / 86400, ch2_meansort / ch2_meansort(1),  psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in days starting from Dec 01 2003', title = 'Channel 2 F30s', yrange = [0.8,1.4]
oplot, (nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400, nozodi_ch2_meansort / nozodi_ch2_meansort(1),  psym = 4

plot, (ch3_sclksort - ch3_sclksort(0)) / 86400, ch3_meansort / ch3_meansort(1),  psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in days starting from Dec 01 2003', title = 'Channel 3 F30s', yrange = [-0.5,2.5], ystyle = 1
oplot, (nozodi_ch3_sclksort - nozodi_ch3_sclksort(0)) / 86400, nozodi_ch3_meansort / nozodi_ch3_meansort(1),  psym = 4

plot, (ch4_sclksort - ch4_sclksort(0)) / 86400, ch4_meansort / ch4_meansort(1),  psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in days starting from Dec 01 2003', title = 'Channel 4 F30s', yrange = [0.0,2.0]
oplot, (nozodi_ch4_sclksort - nozodi_ch4_sclksort(0)) / 86400, nozodi_ch4_meansort / nozodi_ch4_meansort(1),  psym = 4

;-------------------------------------------------------

;fit with sin curve

;noise = fltarr(n_elements(ch2_meansort))
noise = ch2_sigmasort / ch2_meansort(1)

start = [0.03,1.03,60.,3.1]
x =  (ch2_sclksort - ch2_sclksort(0)) / 86400
result= MPFITFUN('sin_func',x, ch2_meansort / ch2_meansort(1), noise, start)    ;ICL
;oplot, x,  result(0)*sin(X/result(2) + result(3)) + result(1), thick = 3;,psym = 2;, color = colors.green

period = 2 * !PI * result(2)
print, 'period', period

;xyouts, 50, 0.92, strcompress('Period' + string(round(period))+ ' days')
;xyouts, 350, 0.92, strcompress('Amplitude' + strmid(string(fix(1000*result(0)) / 10.),0,9)+ '%')
;--------------------------------------------------------------
;try plotting colors
;don't use the normalized means
!P.multi = [0,1,2]

ch1flux_zodi = ch1_meansort;/ ch1_meansort(1) 
ch1flux_nozodi = nozodi_ch1_meansort;/ nozodi_ch1_meansort(1)
ch2flux_zodi = ch2_meansort;/ ch2_meansort(1) 
ch2flux_nozodi = nozodi_ch2_meansort;/ nozodi_ch2_meansort(1)

plot,  (ch1_sclksort - ch1_sclksort(0)) / 86400, ch1flux_nozodi - ch2flux_nozodi, psym = 4, ytitle = 'Ch1 - Ch2', xtitle = 'Time in days starting from Dec 01 2003', title = ' F30s';, yrange = [0.8, 1.2]
plot,  (ch2_sclksort - ch2_sclksort(0)) / 86400, ch1flux_zodi - ch2flux_zodi, psym = 2, ytitle = 'Ch1 - Ch2 with zodi subtraction', xtitle = 'Time in days starting from Dec 01 2003', title = ' F30s';, yrange = [1.10, 1.15]

;--------------------------------------------------------------

ps_close, /noprint,/noid
;ps_end;, /png



end
