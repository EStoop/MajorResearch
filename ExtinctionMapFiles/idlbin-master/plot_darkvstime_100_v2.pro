pro plot_darkvstime

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
;need to combine cryo and warm arrays.

;make a real day array for warm data
warmdayarr = fltarr(n_elements(nozodi_ch1_sclksort))
for i = 0, n_elements(nozodi_ch1_sclksort)- 1 do begin
   year = (nozodi_ch1_sclksort[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   warmdayarr[i] = day
;   print, 'warm',nozodi_ch1_sclksort[i], day
endfor

ch2_warmdayarr = fltarr(n_elements(nozodi_ch2_sclksort))
for i = 0, n_elements(nozodi_ch2_sclksort)- 1 do begin
   year = (nozodi_ch2_sclksort[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   ch2_warmdayarr[i] = day
;   print, 'warm',nozodi_ch2_sclksort[i], day
endfor

;make a real day array for cryo data
cryodayarr = fltarr(n_elements(cryo_nozodi_ch1_sclksort))
for i = 0, n_elements(cryo_nozodi_ch1_sclksort)- 1 do begin
   year = (cryo_nozodi_ch1_sclksort[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   cryodayarr[i] = day
;   print, 'cryo',cryo_nozodi_ch1_sclksort[i], day
endfor

ch2_cryodayarr = fltarr(n_elements(cryo_nozodi_ch2_sclksort))
for i = 0, n_elements(cryo_nozodi_ch2_sclksort)- 1 do begin
   year = (cryo_nozodi_ch2_sclksort[i] / 60./60./24./365.) + 1980.
   day = year - fix(year)
   day = day *365.
   ch2_cryodayarr[i] = day
 ;  print, 'cryo',cryo_nozodi_ch2_sclksort[i], day
endfor

early = 278
late = 305

wd = where (warmdayarr gt early and warmdayarr lt late)
cd = where(cryodayarr gt early and cryodayarr lt late)

;what are the levels at the same time of year?
;print, 'warm sclks', mean(warmdayarr[wd])
;print, 'warm mean', mean(nozodi_ch1_meansort[wd])
;print, 'cryo sclks', mean(cryodayarr[cd])
;print, 'cryo mean', mean(cryo_nozodi_ch1_meansort[cd])
ch1_delta_nozodi = mean(nozodi_ch1_meansort[wd]) - mean(cryo_nozodi_ch1_meansort[cd])
ch1_delta = mean(ch1_meansort[wd]) - mean(cryo_ch1_meansort[cd])

;print, 'ch2_warm sclks',mean(ch2_warmdayarr[ch2_wd])
;print, 'ch2_warm mean', mean(nozodi_ch2_meansort[ch2_wd])
;print, 'ch2_cryo sclks', mean(ch2_cryodayarr[ch2_cd])
;print, 'ch2_cryo mean', mean(cryo_nozodi_ch2_meansort[ch2_cd])
ch2_delta_nozodi = mean(nozodi_ch2_meansort[wd]) - mean(cryo_nozodi_ch2_meansort[cd])
ch2_delta = mean(ch2_meansort[wd]) - mean(cryo_ch2_meansort[cd])

nozodi_ch1_sclksort_all = [cryo_nozodi_ch1_sclksort,nozodi_ch1_sclksort]
nozodi_ch1_meansort_all = [cryo_nozodi_ch1_meansort,nozodi_ch1_meansort - ch1_delta_nozodi]
nozodi_ch2_sclksort_all = [cryo_nozodi_ch2_sclksort,nozodi_ch2_sclksort]
nozodi_ch2_meansort_all = [cryo_nozodi_ch2_meansort,nozodi_ch2_meansort-ch2_delta_nozodi]

ch1_sclksort_all = [cryo_ch1_sclksort,ch1_sclksort]
ch1_meansort_all = [cryo_ch1_meansort,ch1_meansort- ch1_delta]
ch2_sclksort_all = [cryo_ch2_sclksort,ch2_sclksort]
ch2_meansort_all = [cryo_ch2_meansort,ch2_meansort-ch2_delta]

;make sure the sigmas are still the same percentage of the means after delta subtraction

ch1_sigma_percent = ch1_sigmasort / ch1_meansort
ch1_sigma_delta = ch1_sigma_percent * (ch1_meansort - ch1_delta)
ch1_sigmasort_all = [cryo_ch1_sigmasort,ch1_sigma_delta]

ch2_sigma_percent = ch2_sigmasort / ch2_meansort
ch2_sigma_delta = ch2_sigma_percent * (ch2_meansort - ch2_delta)
ch2_sigmasort_all = [cryo_ch2_sigmasort,ch2_sigma_delta]

nozodi_ch1_sigma_percent = nozodi_ch1_sigmasort / nozodi_ch1_meansort
nozodi_ch1_sigma_delta = nozodi_ch1_sigma_percent * (nozodi_ch1_meansort - ch1_delta_nozodi)
nozodi_ch1_sigmasort_all = [cryo_nozodi_ch1_sigmasort,nozodi_ch1_sigma_delta]

nozodi_ch2_sigma_percent = nozodi_ch2_sigmasort / nozodi_ch2_meansort
nozodi_ch2_sigma_delta = nozodi_ch2_sigma_percent * (nozodi_ch2_meansort - ch2_delta_nozodi)
nozodi_ch2_sigmasort_all = [cryo_nozodi_ch2_sigmasort,nozodi_ch2_sigma_delta]


;-------------------------------------------------------------------
;-------------------------------------------------------------------
;begin the plotting
z = pp_multiplot(multi_layout=[1, 4], global_xtitle='Time in Years Starting from Dec 01 2003',global_ytitle='Mean Flux (units?)')

;plot the mean levels
mp=z.plot( (ch1_sclksort_all - ch1_sclksort_all(0)) / 86400 / 365, ch1_meansort_all, '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 0,yrange = [0.03, 0.1], yminor = 0, ymajor = 4)
mp = z.errorplot( (ch1_sclksort_all(80:81) - ch1_sclksort_all(0)) / 86400/365, ch1_meansort_all(80:81), ch1_sigmasort_all(80:81),/overplot, multi_index = 0)
mp = z.errorplot( (ch1_sclksort_all(250:251) - ch1_sclksort_all(0)) / 86400/365, ch1_meansort_all(250:251), ch1_sigmasort_all(250:251),/overplot, multi_index = 0)

mp = z.plot((nozodi_ch1_sclksort_all - nozodi_ch1_sclksort_all(0)) / 86400 / 365, nozodi_ch1_meansort_all, '1bs', sym_filled = 1, sym_size = 0.6,/overplot, multi_index = 0)
mp = z.errorplot(  (nozodi_ch1_sclksort_all(40:41) - nozodi_ch1_sclksort_all(0)) / 86400 / 365, nozodi_ch1_meansort_all(40:41), nozodi_ch1_sigmasort_all(40:41),/overplot, multi_index = 0)
mp = z.errorplot( (nozodi_ch1_sclksort_all(250:251) - nozodi_ch1_sclksort_all(0)) / 86400 / 365, nozodi_ch1_meansort_all(250:251), nozodi_ch1_sigmasort_all(250:251),/overplot, multi_index = 0)

;-----
mp = z.plot( (ch2_sclksort_all - ch2_sclksort_all(0)) / 86400 / 365, ch2_meansort_all , '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 1,  yrange = [0.0,0.4], yminor = 0, ymajor = 4)
mp = z.errorplot( (ch2_sclksort_all(80:81) - ch2_sclksort_all(0)) / 86400/365, ch2_meansort_all(80:81), ch2_sigmasort_all(80:81),/overplot, multi_index = 1)
mp = z.errorplot( (ch2_sclksort_all(230:231) - ch2_sclksort_all(0)) / 86400/365, ch2_meansort_all(230:231), ch2_sigmasort_all(230:231),/overplot, multi_index = 1)
mp = z.plot( (nozodi_ch2_sclksort_all - nozodi_ch2_sclksort_all(0)) / 86400 / 365, nozodi_ch2_meansort_all, '1bs', sym_filled = 1, sym_size = 0.6, /overplot, multi_index = 1)
mp = z.errorplot(  (nozodi_ch2_sclksort_all(40:41) - nozodi_ch2_sclksort_all(0)) / 86400 / 365, nozodi_ch2_meansort_all(40:41), nozodi_ch2_sigmasort_all(40:41),/overplot, multi_index = 1)
mp = z.errorplot( (nozodi_ch2_sclksort_all(240:241) - nozodi_ch2_sclksort_all(0)) / 86400 / 365, nozodi_ch2_meansort_all(240:241), nozodi_ch2_sigmasort_all(240:241),/overplot, multi_index = 1)


;-------------------------------------------------------

;fit with mean values with a sin curve

noise = ch2_sigmasort_all ;/ ch2_meansort_all(1)

start = [0.03,1.03,60.,3.1]
x =  (ch2_sclksort_all - ch2_sclksort_all(0)) / 86400 / 365
;result= MPFITFUN('sin_func',x, ch2_meansort_all , noise, start)    ;ICL
;mp= z.plot( x,  result(0)*sin(X/result(2) + result(3)) + result(1),/overplot, multi_index = 1)

;period = 2 * !PI * result(2)
;print, 'period', period

;xyouts, 50, 0.92, strcompress('Period' + string(round(period))+ ' days')
;xyouts, 350, 0.92, strcompress('Amplitude' + strmid(string(fix(1000*result(0)) / 10.),0,9)+ '%')

;-----
mp = z.plot( (cryo_ch3_sclksort - cryo_ch3_sclksort(0)) / 86400 / 365, cryo_ch3_meansort ,  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 2, yrange = [-4.0,2.0], yminor = 0, ymajor = 4)
mp = z.errorplot( (cryo_ch3_sclksort(80:81) - cryo_ch3_sclksort(0)) / 86400/365, cryo_ch3_meansort(80:81), cryo_ch3_sigmasort(80:81),/overplot, multi_index = 2)
mp = z.plot( (cryo_nozodi_ch3_sclksort - cryo_nozodi_ch3_sclksort(0)) / 86400 / 365, cryo_nozodi_ch3_meansort,'1bs', sym_filled = 1, sym_size = 0.6, /overplot , multi_index = 2)
mp = z.errorplot( (cryo_nozodi_ch3_sclksort(40:41) - cryo_nozodi_ch3_sclksort(0)) / 86400/365, cryo_nozodi_ch3_meansort(40:41), cryo_nozodi_ch3_sigmasort(40:41),/overplot, multi_index = 2)
;-----
mp = z.plot( (cryo_ch4_sclksort - cryo_ch4_sclksort(0)) / 86400 / 365, cryo_ch4_meansort ,  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 3, yrange = [-1.0,9.0], yminor = 0, ymajor = 4)
mp = z.errorplot( (cryo_ch4_sclksort(240:241) - cryo_ch4_sclksort(0)) / 86400/365, cryo_ch4_meansort(240:241), cryo_ch4_sigmasort(240:241),/overplot, multi_index = 3)
mp = z.plot( (cryo_nozodi_ch4_sclksort - cryo_nozodi_ch4_sclksort(0)) / 86400 / 365, cryo_nozodi_ch4_meansort,'1bs', sym_filled = 1, sym_size = 0.6, /overplot , multi_index = 3)
mp = z.errorplot( (cryo_nozodi_ch4_sclksort(80:81) - cryo_nozodi_ch4_sclksort(0)) / 86400/365, cryo_nozodi_ch4_meansort(80:81), cryo_nozodi_ch4_sigmasort(80:81),/overplot, multi_index = 3)


;-------------------------------------------------------------------
;-------------------------------------------------------------------
;plot the mean levels with means subtracted.
z2 = pp_multiplot(multi_layout=[1, 4], global_xtitle='Time in Years Starting from Dec 01 2003',global_ytitle='Mean Subtracted Flux (units?)')

sp = z2.plot( (ch1_sclksort_all - ch1_sclksort_all(0)) / 86400 / 365, ch1_meansort_all - mean(ch1_meansort_all),'1ro', sym_filled = 1, sym_size = 0.6, multi_index = 0,yrange = [-0.02, 0.02], yminor = 0, ymajor = 5)
sp = z2.plot( (nozodi_ch1_sclksort_all - nozodi_ch1_sclksort_all(0)) / 86400 / 365, nozodi_ch1_meansort_all - mean(nozodi_ch1_meansort_all),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 0,/overplot)
;-----
sp = z2.plot( (ch2_sclksort_all - ch2_sclksort_all(0)) / 86400 / 365, ch2_meansort_all - mean(ch2_meansort_all) , '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 1, yrange = [-0.1,0.1], yminor = 0, ymajor = 5)
sp = z2.plot( (nozodi_ch2_sclksort_all - nozodi_ch2_sclksort_all(0)) / 86400 / 365, nozodi_ch2_meansort_all - mean(nozodi_ch2_meansort_all),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 1,/overplot)
;-----
sp = z2.plot( (cryo_ch3_sclksort - cryo_ch3_sclksort(0)) / 86400 / 365, cryo_ch3_meansort - mean(cryo_ch3_meansort), '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 2,  yrange = [-0.5,0.5], yminor = 0, ymajor = 5)
sp = z2.plot( (cryo_nozodi_ch3_sclksort - cryo_nozodi_ch3_sclksort(0)) / 86400 / 365, cryo_nozodi_ch3_meansort - mean(cryo_nozodi_ch3_meansort),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 2,/overplot)
;-----
sp = z2.plot( (cryo_ch4_sclksort - cryo_ch4_sclksort(0)) / 86400 / 365, cryo_ch4_meansort - mean(cryo_ch4_meansort),  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 3, yrange = [-1.5,1.5], yminor = 0, ymajor = 5)
sp = z2.plot( (cryo_nozodi_ch4_sclksort - cryo_nozodi_ch4_sclksort(0)) / 86400 / 365, cryo_nozodi_ch4_meansort - mean(cryo_nozodi_ch4_meansort) ,'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 3,/overplot)

;-------------------------------------------------------------------
;-------------------------------------------------------
;plot the normalized versions
;z3 = pp_multiplot(multi_layout=[1, 4], global_xtitle='Time in Years Starting from Dec 01 2003',global_ytitle='Normalized Flux')

;np = z3.plot( (ch1_sclksort_all - ch1_sclksort_all(0)) / 86400 / 365, ch1_meansort_all/ ch1_meansort_all(0), '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 0,yrange = [0.5,1.5], yminor = 0, ymajor = 4)
;np=z3.plot( (nozodi_ch1_sclksort_all - nozodi_ch1_sclksort_all(0)) / 86400 / 365, nozodi_ch1_meansort_all/ nozodi_ch1_meansort_all(0) ,'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 0,/overplot)

;np=z3.plot( (ch2_sclksort_all - ch2_sclksort_all(0)) / 86400 / 365, ch2_meansort_all / ch2_meansort_all(0),  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 1,yrange = [0.6,1.8], yminor = 0, ymajor = 4)
;np=z3.plot( (nozodi_ch2_sclksort_all - nozodi_ch2_sclksort_all(0)) / 86400 / 365, nozodi_ch2_meansort_all / nozodi_ch2_meansort_all(0),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 1,/overplot)

;np=z3.plot( (cryo_ch3_sclksort - cryo_ch3_sclksort(0)) / 86400 / 365, cryo_ch3_meansort / cryo_ch3_meansort(0),  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 2,yrange = [-2.0, 1.5], yminor = 0, ymajor = 4)
;np=z3.plot( (cryo_nozodi_ch3_sclksort - cryo_nozodi_ch3_sclksort(0)) / 86400 / 365, cryo_nozodi_ch3_meansort / cryo_nozodi_ch3_meansort(0),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 2,/overplot)

;np=z3.plot( (cryo_ch4_sclksort - cryo_ch4_sclksort(0)) / 86400 / 365, cryo_ch4_meansort / cryo_ch4_meansort(0),  '1ro', sym_filled = 1, sym_size = 0.6, multi_index = 3,yrange = [-2.0,2.0], yminor = 0, ymajor = 4)
;np=z3.plot( (cryo_nozodi_ch4_sclksort - cryo_nozodi_ch4_sclksort(0)) / 86400 / 365, cryo_nozodi_ch4_meansort / cryo_nozodi_ch4_meansort(0),'1bs', sym_filled = 1, sym_size = 0.6, multi_index = 3,/overplot)

;--------------------------------------------------------------
;try plotting colors
;don't use the normalized means
!P.multi = [0,1,2]

ch1flux_zodi = ch1_meansort_all;/ ch1_meansort_all(1) 
ch1flux_nozodi = nozodi_ch1_meansort_all;/ nozodi_ch1_meansort_all(1)
ch2flux_zodi = ch2_meansort_all;/ ch2_meansort_all(1) 
ch2flux_nozodi = nozodi_ch2_meansort_all;/ nozodi_ch2_meansort_all(1)

;plot,  (ch1_sclksort_all - ch1_sclksort_all(0)) / 86400 / 365, ch1flux_nozodi - ch2flux_nozodi, psym = 4, ytitle = 'Ch1 - Ch2', xtitle = 'Time in years starting from Dec 01 2003', title = ' F100s', yrange = [-0.4, 0.1], ystyle = 1
;oplot,  (ch2_sclksort_all - ch2_sclksort_all(0)) / 86400 / 365, ch1flux_zodi - ch2flux_zodi, psym = 2

;plot,  (cryo_ch2_sclksort - cryo_ch2_sclksort(0)) / 86400 / 365, cryo_nozodi_ch2_meansort - cryo_nozodi_ch2_meansort, psym = 4, ytitle = 'Ch2 - Ch3', xtitle = 'Time in years starting from Dec 01 2003', title = ' F100s', yrange = [-0.1, 0.1], ystyle = 1, xrange = [0,8]
;oplot, (cryo_ch2_sclksort - cryo_ch2_sclksort(0)) / 86400 / 365, cryo_ch2_meansort - cryo_ch4_meansort, psym = 2

;--------------------------------------------------------------

;save the plots
; mp.Save , '/Users/jkrick/irac_darks/meanvstime_100.png'
;  np.Save ,  '/Users/jkrick/irac_darks/normalizedvstime_100.png'
;  sp.Save, '/Users/jkrick/irac_darks/meansubvstime_100.png'

end
