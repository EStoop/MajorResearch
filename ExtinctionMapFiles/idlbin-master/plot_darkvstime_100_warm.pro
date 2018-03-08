pro plot_darkvstime_100_warm

;the warm dark measurements
;zodi from header not removed from the data
restore, filename = '/Users/jkrick/irac_darks/warmdarks_100_nozodisub.sav'
nozodi_ch1_sclksort = ch1_sclksort
nozodi_ch1_meansort = ch1_meansort
nozodi_ch1_sigmasort = ch1_sigmasort
nozodi_ch2_sclksort = ch2_sclksort
nozodi_ch2_meansort = ch2_meansort
nozodi_ch2_sigmasort = ch2_sigmasort

;zodi in the header has been subtracted 
restore, filename = '/Users/jkrick/irac_darks/warmdarks_100_zodisub.sav'
;-------------------------------------------------------------------
;-------------------------------------------------------------------
;begin the plotting
z = pp_multiplot(multi_layout=[1, 2], global_xtitle='Time in Years Starting from Oct 15, 2009',global_ytitle='Flux (MJy/sr)')

;plot the mean levels
mp=z.plot( (ch1_sclksort - ch1_sclksort(0)) / 86400 / 365, ch1_meansort, '1ro', sym_filled = 1, sym_size = 0.6, yrange = [0.65, 0.75],multi_index = 0, yminor = 0, ymajor = 5)
mp=z.errorplot( (ch1_sclksort - ch1_sclksort(0)) / 86400 / 365, ch1_meansort,ch1_sigmasort,errorbar_color = 'red', linestyle = 6,multi_index = 0,/overplot)
mp = z.plot((nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400 / 365, nozodi_ch1_meansort,'1bs', sym_filled = 1, sym_size = 0.6,/overplot, multi_index = 0)
mp = z.errorplot((nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400 / 365, nozodi_ch1_meansort, nozodi_ch1_sigmasort,errorbar_color = 'blue', linestyle = 6,/overplot, multi_index = 0)

;fit with mean values with a sin curve

;need to get rid of outliers
a = [74,75,76,77,78,79]
ch1_meansort[a] = alog10(-1)
ch1_sigmasort[a] = alog10(-1)

noise = ch1_sigmasort ;/ ch1_meansort(1)

start = [0.01,mean(ch1_meansort,/nan),1.0/(2*!Pi),3.1]
x =  (ch1_sclksort - ch1_sclksort(0)) / 86400 / 365
print, 'x'


result= MPFITFUN('sin_func',x, ch1_meansort , noise, start,/nan)    ;ICL
mp= z.plot( x,  result(0)*sin(x/result(2) + result(3)) + result(1),/overplot, multi_index = 0)

period = 2 * !PI * result(2)
print, 'period', period

;t1 = text( 0.7, 0.68, strcompress('Period' + string(period)+ ' Years'),/data)
;t2 = text( 1.5, 0.68, strcompress('Amplitude ' + string(result(0))+ 'MJy/sr'),/data)

;-----
mp = z.plot( (ch2_sclksort - ch2_sclksort(0)) / 86400 / 365, ch2_meansort , '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 1,  yrange = [0.1,0.4], yminor = 0, ymajor = 5)
mp = z.errorplot( (ch2_sclksort - ch2_sclksort(0)) / 86400 / 365, ch2_meansort , ch2_sigmasort,errorbar_color='red', linestyle = 6, multi_index = 1, /overplot)

mp = z.plot( (nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400 / 365, nozodi_ch2_meansort, '1bs', sym_filled = 1, sym_size = 0.6, /overplot, multi_index = 1)
mp = z.errorplot( (nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400 / 365, nozodi_ch2_meansort, nozodi_ch2_sigmasort, errorbar_color = 'blue', linestyle = 6, /overplot, multi_index = 1)


;-------------------------------------------------------

;fit with mean values with a sin curve

noise = ch2_sigmasort ;/ ch2_meansort(1)

start = [0.01,mean(ch2_meansort),1.0/(2*!Pi),3.1]
x =  (ch2_sclksort - ch2_sclksort(0)) / 86400 / 365
print, 'x'
result= MPFITFUN('sin_func',x, ch2_meansort , noise, start)    ;ICL
mp= z.plot( x,  result(0)*sin(x/result(2) + result(3)) + result(1),/overplot, multi_index = 1)

period = 2 * !PI * result(2)
print, 'period', period

;t1 = text( 0.7, 0.175, strcompress('Period' + string(period)+ ' Years'),/data)
;t2 = text( 1.5, 0.175, strcompress('Amplitude ' + string(result(0))+ 'MJy/sr'),/data)

;-------------------------------------------------------------------
;-------------------------------------------------------------------
;plot the mean levels with means subtracted.
z2 = pp_multiplot(multi_layout=[1, 2], global_xtitle='Time in Years Starting from Oct 15, 2009',global_ytitle='Mean Subtracted Flux (MJy/sr)')
print, 'max time in years', max((ch1_sclksort - ch1_sclksort(0)) / 86400 / 365)

sp = z2.plot( (ch1_sclksort - ch1_sclksort(0)) / 86400 / 365, ch1_meansort - mean(ch1_meansort),'1ro', sym_filled = 1, sym_size = 0.6, multi_index = 0,yrange = [-0.02, 0.02],   ymajor = 1)
sp = z2.errorplot( (ch1_sclksort - ch1_sclksort(0)) / 86400 / 365, ch1_meansort - mean(ch1_meansort),ch1_sigmasort,multi_index = 0, /overplot, errorbar_color = 'red', linestyle = 6)
sp = z2.plot( (nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400 / 365, nozodi_ch1_meansort - mean(nozodi_ch1_meansort),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 0,/overplot)
sp = z2.errorplot( (nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400 / 365, nozodi_ch1_meansort - mean(nozodi_ch1_meansort), nozodi_ch1_sigmasort, errorbar_color = 'blue', linestyle = 6, multi_index = 0,/overplot)
;-----
sp = z2.plot( (ch2_sclksort - ch2_sclksort(0)) / 86400 / 365, ch2_meansort - mean(ch2_meansort) , '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 1, yrange = [-0.05,0.05], yminor = 0, ymajor = 5)
sp = z2.errorplot( (ch2_sclksort - ch2_sclksort(0)) / 86400 / 365, ch2_meansort - mean(ch2_meansort) , ch2_sigmasort, errorbar_color = 'red', linestyle = 6,multi_index = 1,/overplot)
sp = z2.plot( (nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400 / 365, nozodi_ch2_meansort - mean(nozodi_ch2_meansort),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 1,/overplot)
sp = z2.errorplot( (nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400 / 365, nozodi_ch2_meansort - mean(nozodi_ch2_meansort),nozodi_ch2_sigmasort, errorbar_color = 'blue', linestyle = 6,multi_index = 1,/overplot)
;-----

;-------------------------------------------------------------------
;-------------------------------------------------------
;plot the normalized versions
z3 = pp_multiplot(multi_layout=[1, 2], global_xtitle='Time in Years Starting from Oct 15, 2009',global_ytitle='Normalized Flux')

np = z3.plot( (ch1_sclksort - ch1_sclksort(0)) / 86400 / 365, ch1_meansort/ ch1_meansort(0), '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 0,yrange = [0.95,1.05], yminor = 0, ymajor = 5)
np = z3.errorplot( (ch1_sclksort - ch1_sclksort(0)) / 86400 / 365, ch1_meansort/ ch1_meansort(0), ch1_sigmasort/ch1_meansort(0), errorbar_color = 'red', linestyle = 6,  multi_index = 0,/overplot)
np=z3.plot( (nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400 / 365, nozodi_ch1_meansort/ nozodi_ch1_meansort(0) ,'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 0,/overplot)
np=z3.errorplot( (nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400 / 365, nozodi_ch1_meansort/ nozodi_ch1_meansort(0) ,  nozodi_ch1_sigmasort/nozodi_ch1_meansort(0), errorbar_color = 'blue', linestyle = 6,multi_index = 0,/overplot)

np=z3.plot( (ch2_sclksort - ch2_sclksort(0)) / 86400 / 365, ch2_meansort / ch2_meansort(0),  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 1,yrange = [0.9,1.3], yminor = 0, ymajor = 5)
np=z3.errorplot( (ch2_sclksort - ch2_sclksort(0)) / 86400 / 365, ch2_meansort / ch2_meansort(0), ch2_sigmasort/ ch2_meansort(0), errorbar_color = 'red', linestyle = 6,  multi_index = 1,/overplot)
np=z3.plot( (nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400 / 365, nozodi_ch2_meansort / nozodi_ch2_meansort(0),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 1,/overplot)
np=z3.errorplot( (nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400 / 365, nozodi_ch2_meansort / nozodi_ch2_meansort(0), nozodi_ch2_sigmasort / nozodi_ch2_meansort(0), errorbar_color = 'blue', linestyle = 6, multi_index = 1,/overplot)


;--------------------------------------------------------------
;try plotting colors
;don't use the normalized means
!P.multi = [0,1,2]

ch1flux_zodi = ch1_meansort;/ ch1_meansort(1) 
ch1flux_nozodi = nozodi_ch1_meansort;/ nozodi_ch1_meansort(1)
ch2flux_zodi = ch2_meansort;/ ch2_meansort(1) 
ch2flux_nozodi = nozodi_ch2_meansort;/ nozodi_ch2_meansort(1)

;plot,  (ch1_sclksort - ch1_sclksort(0)) / 86400 / 365, ch1flux_nozodi - ch2flux_nozodi, psym = 4, ytitle = 'Ch1 - Ch2', xtitle = 'Time in years starting from Dec 01 2003', title = ' F100s', yrange = [-0.4, 0.1], ystyle = 1
;oplot,  (ch2_sclksort - ch2_sclksort(0)) / 86400 / 365, ch1flux_zodi - ch2flux_zodi, psym = 2

;plot,  (ch2_sclksort - ch2_sclksort(0)) / 86400 / 365, nozodi_ch2_meansort - nozodi_ch2_meansort, psym = 4, ytitle = 'Ch2 - Ch3', xtitle = 'Time in years starting from Dec 01 2003', title = ' F100s', yrange = [-0.1, 0.1], ystyle = 1, xrange = [0,8]
;oplot, (ch2_sclksort - ch2_sclksort(0)) / 86400 / 365, ch2_meansort - ch4_meansort, psym = 2

;--------------------------------------------------------------

;save the plots
; mp.Save , '/Users/jkrick/irac_darks/meanvstime_100.png'
;  np.Save ,  '/Users/jkrick/irac_darks/normalizedvstime_100.png'
;  sp.Save, '/Users/jkrick/irac_darks/meansubvstime_100.png'

end
