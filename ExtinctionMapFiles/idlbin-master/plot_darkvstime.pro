pro plot_darkvstime
ps_open, filename='/Users/jkrick/irac_darks/darksvstime_7.ps',/portrait,/square,/color
;ps_start, filename='/Users/jkrick/irac_darks/darksvstime_3.ps'
!P.thick = 3
!P.charthick = 3
!P.multi = [0,1,2]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)

;zodi from header not removed from the data
restore, filename = '/Users/jkrick/irac_darks/alldarks_4.sav'
nozodi_ch1_sclksort = ch1_sclksort
nozodi_ch1_meansort = ch1_meansort
nozodi_ch1_sigmasort = ch1_sigmasort
nozodi_ch2_sclksort = ch2_sclksort
nozodi_ch2_meansort = ch2_meansort
nozodi_ch2_sigmasort = ch2_sigmasort

;zodi in the header has been subtracted 
restore, filename = '/Users/jkrick/irac_darks/alldarks_3.sav'



;remove 101 and 102 from the arays, these are not real darks, but got run through the dark pipeline so they show up in this directory.
;n = n_elements(ch2_sclksort) - 2
;test_sclksort = fltarr(n_elements(ch2_sclksort) - 2)
;test_meansort = fltarr(n_elements(ch2_meansort) - 2)
;test_sigmasort = fltarr(n_elements(ch2_sigmasort) - 2)
;print, 'meansort', ch2_meansort / ch2_meansort(1)
;test_meansort = [ch2_meansort[0:100], ch2_meansort[103:*]]
;test_sclksort =  [ch2_sclksort[0:100], ch2_sclksort[103:*]]
;test_sigmasort =  [ch2_sigmasort[0:100], ch2_sigmasort[103:*]]
;print, 'test', test_meansort / test_meansort(1)

plot, (ch1_sclksort - ch1_sclksort(0)) / 86400, ch1_meansort, psym = 2, ytitle = 'Mean flux', xtitle = 'Time in days starting from pc5', title = 'Channel 1 F30s', yrange = [1.20,1.40]
oploterror, (ch1_sclksort(120) - ch1_sclksort(0)) / 86400, ch1_meansort(120), ch1_sigmasort(120)
print, 'sigma', ch1_sigmasort(120)


oplot, (nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400, nozodi_ch1_meansort, psym = 4
oploterror, (nozodi_ch1_sclksort(130) - nozodi_ch1_sclksort(0)) / 86400, nozodi_ch1_meansort(130), nozodi_ch1_sigmasort(130)
print, 'nozodi sigma', nozodi_ch1_sigmasort(120)


plot, (ch2_sclksort - ch2_sclksort(0)) / 86400, ch2_meansort ,  psym = 2, ytitle = 'Mean flux', xtitle = 'Time in days starting from pc5', title = 'Channel 2 F30s', yrange = [0.1,0.4]
oploterror, (ch2_sclksort(120) - ch2_sclksort(0)) / 86400, ch2_meansort(120), ch2_sigmasort(120)
print, 'sigma', ch2_sigmasort(120)

oplot, (nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400, nozodi_ch2_meansort ,  psym = 4
oploterror, (nozodi_ch2_sclksort(130) - nozodi_ch2_sclksort(0)) / 86400, nozodi_ch2_meansort(130), nozodi_ch2_sigmasort(130)
print, 'nozodi sigma', nozodi_ch2_sigmasort(120)
;-------------------------------------------------------
;plot the normalized versions
plot, (ch1_sclksort - ch1_sclksort(0)) / 86400, ch1_meansort/ ch1_meansort(1) , psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in days starting from pc5', title = 'Channel 1 F30s', yrange = [0.96,1.02]
oploterror, (ch1_sclksort(120) - ch1_sclksort(0)) / 86400, ch1_meansort(120)/ ch1_meansort(1) , ch1_sigmasort(120)/ ch1_meansort(1) 
oplot, (nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400, nozodi_ch1_meansort/ nozodi_ch1_meansort(1) , psym = 4
oploterror, (nozodi_ch1_sclksort(130) - nozodi_ch1_sclksort(0)) / 86400, nozodi_ch1_meansort(130)/ nozodi_ch1_meansort(1) , nozodi_ch1_sigmasort(130)/ nozodi_ch1_meansort(1) 

plot, (ch2_sclksort - ch2_sclksort(0)) / 86400, ch2_meansort / ch2_meansort(1),  psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in days starting from pc5', title = 'Channel 2 F30s', yrange = [0.9,1.2]
oploterror, (ch2_sclksort(120) - ch2_sclksort(0)) / 86400, ch2_meansort(120)/ ch2_meansort(1) , ch2_sigmasort(120)/ ch2_meansort(1) 
oplot, (nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400, nozodi_ch2_meansort / nozodi_ch2_meansort(1),  psym = 4
oploterror, (nozodi_ch2_sclksort(130) - nozodi_ch2_sclksort(0)) / 86400, nozodi_ch2_meansort(130)/ nozodi_ch2_meansort(1) , nozodi_ch2_sigmasort(130)/ nozodi_ch2_meansort(1) 


;-------------------------------------------------------

;fit with sin curve

;noise = fltarr(n_elements(ch2_meansort))
noise = ch2_sigmasort / ch2_meansort(1)

start = [0.03,1.03,60.,3.1]
x =  (ch2_sclksort - ch2_sclksort(0)) / 86400
result= MPFITFUN('sin_func',x, ch2_meansort / ch2_meansort(1), noise, start)    ;ICL
oplot, x,  result(0)*sin(X/result(2) + result(3)) + result(1), thick = 3;,psym = 2;, color = colors.green
;oplot, x, 0.03*sin(x/57 + 3.1) + 1.03

period = 2 * !PI * result(2)
print, 'period', period

xyouts, 50, 0.92, strcompress('Period' + string(round(period))+ ' days')
xyouts, 350, 0.92, strcompress('Amplitude' + strmid(string(fix(1000*result(0)) / 10.),0,9)+ '%')
;--------------------------------------------------------------
;try plotting colors
;don't use the normalized means

ch1flux_zodi = ch1_meansort;/ ch1_meansort(1) 
ch1flux_nozodi = nozodi_ch1_meansort;/ nozodi_ch1_meansort(1)
ch2flux_zodi = ch2_meansort;/ ch2_meansort(1) 
ch2flux_nozodi = nozodi_ch2_meansort;/ nozodi_ch2_meansort(1)

plot,  (ch1_sclksort - ch1_sclksort(0)) / 86400, ch1flux_nozodi - ch2flux_nozodi, psym = 4, ytitle = 'Ch1 - Ch2', xtitle = 'Time in days starting from pc5', title = ' F30s', yrange = [0.8, 1.2]
oplot,  (ch2_sclksort - ch2_sclksort(0)) / 86400, ch1flux_zodi - ch2flux_zodi, psym = 2

;--------------------------------------------------------------

ps_close, /noprint,/noid
;ps_end;, /png



end
