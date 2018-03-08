pro plot_darkvstime_cryo

;zodi from header not removed from the data
restore, filename = '/Users/jkrick/irac_darks/cryodarks_100_nozodisub.sav'
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
restore, filename = '/Users/jkrick/irac_darks/cryodarks_100_zodisub.sav'
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
z = pp_multiplot(multi_layout=[1, 4], global_xtitle='Time in Years Starting from Dec 01 2003',global_ytitle='Mean Flux (units?)')

;plot the mean levels
mp=z.plot( (cryo_ch1_sclksort - cryo_ch1_sclksort(0)) / 86400 / 365, cryo_ch1_meansort, '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 0,yrange = [0.03, 0.1], yminor = 0, ymajor = 5)
;mp = z.errorplot( (cryo_ch1_sclksort(80:81) - cryo_ch1_sclksort(0)) / 86400/365, cryo_ch1_meansort(80:81), cryo_ch1_sigmasort(80:81),/overplot, multi_index = 0)
;mp = z.errorplot( (cryo_ch1_sclksort(250:251) - cryo_ch1_sclksort(0)) / 86400/365, cryo_ch1_meansort(250:251), cryo_ch1_sigmasort(250:251),/overplot, multi_index = 0)

mp = z.plot((cryo_nozodi_ch1_sclksort - cryo_nozodi_ch1_sclksort(0)) / 86400 / 365, cryo_nozodi_ch1_meansort, '1bs', sym_filled = 1, sym_size = 0.6,/overplot, multi_index = 0)
;mp = z.errorplot(  (cryo_nozodi_ch1_sclksort(40:41) - cryo_nozodi_ch1_sclksort(0)) / 86400 / 365, cryo_nozodi_ch1_meansort(40:41), cryo_nozodi_ch1_sigmasort(40:41),/overplot, multi_index = 0)
;mp = z.errorplot( (cryo_nozodi_ch1_sclksort(250:251) - cryo_nozodi_ch1_sclksort(0)) / 86400 / 365, cryo_nozodi_ch1_meansort(250:251), cryo_nozodi_ch1_sigmasort(250:251),/overplot, multi_index = 0)

;-----
mp = z.plot( (cryo_ch2_sclksort - cryo_ch2_sclksort(0)) / 86400 / 365, cryo_ch2_meansort , '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 1,  yrange = [0.0,0.4], yminor = 0, ymajor = 5)
;mp = z.errorplot( (cryo_ch2_sclksort(80:81) - cryo_ch2_sclksort(0)) / 86400/365, cryo_ch2_meansort(80:81), cryo_ch2_sigmasort(80:81),/overplot, multi_index = 1)
;mp = z.errorplot( (cryo_ch2_sclksort(230:231) - cryo_ch2_sclksort(0)) / 86400/365, cryo_ch2_meansort(230:231), cryo_ch2_sigmasort(230:231),/overplot, multi_index = 1)

mp = z.plot( (cryo_nozodi_ch2_sclksort - cryo_nozodi_ch2_sclksort(0)) / 86400 / 365, cryo_nozodi_ch2_meansort, '1bs', sym_filled = 1, sym_size = 0.6, /overplot, multi_index = 1)
;mp = z.errorplot(  (cryo_nozodi_ch2_sclksort(40:41) - cryo_nozodi_ch2_sclksort(0)) / 86400 / 365, cryo_nozodi_ch2_meansort(40:41), cryo_nozodi_ch2_sigmasort(40:41),/overplot, multi_index = 1)
;mp = z.errorplot( (cryo_nozodi_ch2_sclksort(240:241) - cryo_nozodi_ch2_sclksort(0)) / 86400 / 365, cryo_nozodi_ch2_meansort(240:241), cryo_nozodi_ch2_sigmasort(240:241),/overplot, multi_index = 1)


;-------------------------------------------------------

;fit with mean values with a sin curve

noise = ch2_sigmasort ;/ ch2_meansort(1)

start = [0.03,1.03,60.,3.1]
x =  (ch2_sclksort - ch2_sclksort(0)) / 86400 / 365
;result= MPFITFUN('sin_func',x, ch2_meansort , noise, start)    ;ICL
;mp= z.plot( x,  result(0)*sin(X/result(2) + result(3)) + result(1),/overplot, multi_index = 1)

;period = 2 * !PI * result(2)
;print, 'period', period

;xyouts, 50, 0.92, strcompress('Period' + string(round(period))+ ' days')
;xyouts, 350, 0.92, strcompress('Amplitude' + strmid(string(fix(1000*result(0)) / 10.),0,9)+ '%')

;-----
mp = z.plot( (cryo_ch3_sclksort - cryo_ch3_sclksort(0)) / 86400 / 365, cryo_ch3_meansort ,  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 2, yrange = [-4.0,2.0], yminor = 0, ymajor = 5)
;mp = z.errorplot( (cryo_ch3_sclksort(80:81) - cryo_ch3_sclksort(0)) / 86400/365, cryo_ch3_meansort(80:81), cryo_ch3_sigmasort(80:81),/overplot, multi_index = 2)
mp = z.plot( (cryo_nozodi_ch3_sclksort - cryo_nozodi_ch3_sclksort(0)) / 86400 / 365, cryo_nozodi_ch3_meansort,'1bs', sym_filled = 1, sym_size = 0.6, /overplot , multi_index = 2)
;mp = z.errorplot( (cryo_nozodi_ch3_sclksort(40:41) - cryo_nozodi_ch3_sclksort(0)) / 86400/365, cryo_nozodi_ch3_meansort(40:41), cryo_nozodi_ch3_sigmasort(40:41),/overplot, multi_index = 2)
;-----
mp = z.plot( (cryo_ch4_sclksort - cryo_ch4_sclksort(0)) / 86400 / 365, cryo_ch4_meansort ,  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 3, yrange = [-1.0,9.0], yminor = 0, ymajor = 5)
;mp = z.errorplot( (cryo_ch4_sclksort(240:241) - cryo_ch4_sclksort(0)) / 86400/365, cryo_ch4_meansort(240:241), cryo_ch4_sigmasort(240:241),/overplot, multi_index = 3)
mp = z.plot( (cryo_nozodi_ch4_sclksort - cryo_nozodi_ch4_sclksort(0)) / 86400 / 365, cryo_nozodi_ch4_meansort,'1bs', sym_filled = 1, sym_size = 0.6, /overplot , multi_index = 3)
;mp = z.errorplot( (cryo_nozodi_ch4_sclksort(80:81) - cryo_nozodi_ch4_sclksort(0)) / 86400/365, cryo_nozodi_ch4_meansort(80:81), cryo_nozodi_ch4_sigmasort(80:81),/overplot, multi_index = 3)


;-------------------------------------------------------------------
;-------------------------------------------------------------------
;plot the mean levels with means subtracted.
z2 = pp_multiplot(multi_layout=[1, 4], global_xtitle='Time in Years Starting from Dec 01 2003',global_ytitle='Mean Subtracted Flux (units?)')
print, 'max', max((ch1_sclksort - ch1_sclksort(0)) / 86400 / 365)
sp = z2.plot( (cryo_ch1_sclksort - cryo_ch1_sclksort(0)) / 86400 / 365, cryo_ch1_meansort - mean(cryo_ch1_meansort),'1ro', sym_filled = 1, sym_size = 0.6, multi_index = 0,yrange = [-0.02, 0.02],  yminor = 0, ymajor = 5)
sp = z2.plot( (cryo_nozodi_ch1_sclksort - cryo_nozodi_ch1_sclksort(0)) / 86400 / 365, cryo_nozodi_ch1_meansort - mean(cryo_nozodi_ch1_meansort),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 0,/overplot)
;-----
sp = z2.plot( (cryo_ch2_sclksort - cryo_ch2_sclksort(0)) / 86400 / 365, cryo_ch2_meansort - mean(cryo_ch2_meansort) , '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 1, yrange = [-0.1,0.1], yminor = 0, ymajor = 5)
sp = z2.plot( (cryo_nozodi_ch2_sclksort - cryo_nozodi_ch2_sclksort(0)) / 86400 / 365, cryo_nozodi_ch2_meansort - mean(cryo_nozodi_ch2_meansort),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 1,/overplot)
;-----
sp = z2.plot( (cryo_ch3_sclksort - cryo_ch3_sclksort(0)) / 86400 / 365, cryo_ch3_meansort - mean(cryo_ch3_meansort), '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 2,  yrange = [-0.5,0.5], yminor = 0, ymajor = 5)
sp = z2.plot( (cryo_nozodi_ch3_sclksort - cryo_nozodi_ch3_sclksort(0)) / 86400 / 365, cryo_nozodi_ch3_meansort - mean(cryo_nozodi_ch3_meansort),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 2,/overplot)
;-----
sp = z2.plot( (cryo_ch4_sclksort - cryo_ch4_sclksort(0)) / 86400 / 365, cryo_ch4_meansort - mean(cryo_ch4_meansort),  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 3, yrange = [-1.5,1.5], yminor = 0, ymajor = 5)
sp = z2.plot( (cryo_nozodi_ch4_sclksort - cryo_nozodi_ch4_sclksort(0)) / 86400 / 365, cryo_nozodi_ch4_meansort - mean(cryo_nozodi_ch4_meansort) ,'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 3,/overplot)

;-------------------------------------------------------------------
;-------------------------------------------------------
;plot the normalized versions
z3 = pp_multiplot(multi_layout=[1, 4], global_xtitle='Time in Years Starting from Dec 01 2003',global_ytitle='Normalized Flux')

np = z3.plot( (cryo_ch1_sclksort - cryo_ch1_sclksort(0)) / 86400 / 365, cryo_ch1_meansort/ cryo_ch1_meansort(0), '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 0,yrange = [0.5,1.5], yminor = 0, ymajor = 5)
np=z3.plot( (cryo_nozodi_ch1_sclksort - cryo_nozodi_ch1_sclksort(0)) / 86400 / 365, cryo_nozodi_ch1_meansort/ cryo_nozodi_ch1_meansort(0) ,'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 0,/overplot)

np=z3.plot( (cryo_ch2_sclksort - cryo_ch2_sclksort(0)) / 86400 / 365, cryo_ch2_meansort / cryo_ch2_meansort(0),  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 1,yrange = [0.6,1.8], yminor = 0, ymajor = 5)
np=z3.plot( (cryo_nozodi_ch2_sclksort - cryo_nozodi_ch2_sclksort(0)) / 86400 / 365, cryo_nozodi_ch2_meansort / cryo_nozodi_ch2_meansort(0),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 1,/overplot)

np=z3.plot( (cryo_ch3_sclksort - cryo_ch3_sclksort(0)) / 86400 / 365, cryo_ch3_meansort / cryo_ch3_meansort(0),  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 2,yrange = [-2.0, 1.5], yminor = 0, ymajor = 5)
np=z3.plot( (cryo_nozodi_ch3_sclksort - cryo_nozodi_ch3_sclksort(0)) / 86400 / 365, cryo_nozodi_ch3_meansort / cryo_nozodi_ch3_meansort(0),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 2,/overplot)

np=z3.plot( (cryo_ch4_sclksort - cryo_ch4_sclksort(0)) / 86400 / 365, cryo_ch4_meansort / cryo_ch4_meansort(0),  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 3,yrange = [-2.0,2.0], yminor = 0, ymajor = 5)
np=z3.plot( (cryo_nozodi_ch4_sclksort - cryo_nozodi_ch4_sclksort(0)) / 86400 / 365, cryo_nozodi_ch4_meansort / cryo_nozodi_ch4_meansort(0),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 3,/overplot)

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

;plot,  (cryo_ch2_sclksort - cryo_ch2_sclksort(0)) / 86400 / 365, cryo_nozodi_ch2_meansort - cryo_nozodi_ch2_meansort, psym = 4, ytitle = 'Ch2 - Ch3', xtitle = 'Time in years starting from Dec 01 2003', title = ' F100s', yrange = [-0.1, 0.1], ystyle = 1, xrange = [0,8]
;oplot, (cryo_ch2_sclksort - cryo_ch2_sclksort(0)) / 86400 / 365, cryo_ch2_meansort - cryo_ch4_meansort, psym = 2

;--------------------------------------------------------------

;save the plots
; mp.Save , '/Users/jkrick/irac_darks/meanvstime_100.png'
;  np.Save ,  '/Users/jkrick/irac_darks/normalizedvstime_100.png'
;  sp.Save, '/Users/jkrick/irac_darks/meansubvstime_100.png'

end
