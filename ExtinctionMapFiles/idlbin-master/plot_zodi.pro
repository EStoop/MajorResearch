pro plot_zodi, channel

;COBE DIRBE data
  restore, '/Users/jkrick/irac_darks/dirbe/dirbe.sav'
  timearr = timearr - timearr(0) ; in seconds
  timearr = timearr / 60/60/24.  ; in days
  timearr = timearr / 365. ; in years
  
  
  ;do a quick outlier check
  ;3------------------
  resistant_mean,phot3aarr,2, rmean, rsig
  rsigma = robust_sigma(phot3aarr)
  print, '3', rmean, rsigma
  g = where(phot3aarr gt (rmean - 5*rsigma) and phot3aarr lt (rmean + 5*rsigma), gcount)
  print, '3 gcount, orig', gcount, n_elements(phot3aarr)
  phot3aarr = phot3aarr[g]
  phot3barr = phot3barr[g]
  phot3carr = phot3carr[g]
  timearr3 = timearr[g]
  g2 = where(moon2losarr gt 10. and jup2los gt 1.5 and xsnoise ne 8)
  phot3aarr = phot3aarr[g2]
  phot3barr = phot3barr[g2]
  phot3carr = phot3carr[g2]
  timearr3 = timearr3[g2]
  phot3err = fltarr(n_elements(phot3aarr)) + .033

   ;4------------------
  resistant_mean,phot4arr,2, rmean, rsig
  rsigma = robust_sigma(phot4arr)
  print, '4', rmean, rsigma
  g = where(phot4arr gt (rmean - 5*rsigma) and phot4arr lt (rmean + 5*rsigma), gcount)
  print, '4 gcount, orig', gcount, n_elements(phot4arr)
  phot4arr = phot4arr[g]
  timearr4 = timearr[g]
  g2 = where(moon2losarr gt 10. and jup2los gt 1.5 and xsnoise ne 16)
  phot4arr = phot4arr[g2]
  timearr4 = timearr4[g2]
  phot4err = fltarr(n_elements(phot4arr)) + .027
  ;5------------------
  ;need to do something different here XXX
  resistant_mean,phot5arr,2, rmean, rsig
  rsigma = robust_sigma(phot5arr)
  print, '5', rmean, rsigma
  g = where(phot5arr gt (rmean - 5*rsigma) and phot5arr lt (rmean + 5*rsigma), gcount)
  print, '5 gcount, orig', gcount, n_elements(phot5arr)
  phot5arr = phot5arr[g]
  timearr5 = timearr[g]
  g2 = where(moon2losarr gt 10. and jup2los gt 1.5 and xsnoise ne 32 and phot5arr gt 0. and phot5arr lt 21.)
  phot5arr = phot5arr[g2]
  timearr5 = timearr5[g2]
  phot5err = fltarr(n_elements(phot5arr)) + .051


;------------------------------------------------------------------
;binning
;------------------------------------------------------------------
  bin_level = 15L
  nframes = bin_level


;3-------------
  nfits = long(n_elements(phot3aarr)) / nframes
  bin_flux3 = dblarr(nfits)
  bin_flux3err = dblarr(nfits)
  bin_time3= dblarr(nfits)
  for si = 0L, long(nfits) - 1L do begin
     idata = phot3aarr[si*nframes:si*nframes + (nframes - 1)]
     idataerr = phot3err[si*nframes:si*nframes + (nframes - 1)]
     bin_flux3err[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
     bin_flux3[si] = mean(idata,/nan)
     bin_time3[si] = mean(timearr3[si*nframes:si*nframes + (nframes - 1)])
  endfor                       
;4-------------
  nfits = long(n_elements(phot4arr)) / nframes
  bin_flux4 = dblarr(nfits)
  bin_flux4err = dblarr(nfits)
  bin_time4= dblarr(nfits)
  for si = 0L, long(nfits) - 1L do begin
     idata = phot4arr[si*nframes:si*nframes + (nframes - 1)]
     idataerr = phot4err[si*nframes:si*nframes + (nframes - 1)]
     bin_flux4err[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
     bin_flux4[si] = mean(idata,/nan)
     bin_time4[si] = mean(timearr4[si*nframes:si*nframes + (nframes - 1)])
  endfor                       
 ;5-------------
  nfits = long(n_elements(phot5arr)) / nframes
  bin_flux5 = dblarr(nfits)
  bin_flux5err = dblarr(nfits)
  bin_time5= dblarr(nfits)
  for si = 0L, long(nfits) - 1L do begin
     idata = phot5arr[si*nframes:si*nframes + (nframes - 1)]
     idataerr = phot5err[si*nframes:si*nframes + (nframes - 1)]
     bin_flux5err[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
     bin_flux5[si] = mean(idata,/nan)
     bin_time5[si] = mean(timearr5[si*nframes:si*nframes + (nframes - 1)])
  endfor                       
  
print, 'after binning', n_elements(bin_flux3), n_elements(bin_flux4), n_elements(bin_flux5)

;------------------------------------------------------------------
  ;fit with a sine curve
;------------------------------------------------------------------
 ;P(0)*sin(X/P(2) + P(3)) + P(1)
;P(0) = amplitude,
;P(1) = y-intercept
;P(2) = x amplitude = period
;p(3) = phase

   
;5-------------
 ; print, 'working on DIRBE 12'
 ; start5 = [6.,13.,1./(2.*!Pi),2.6]
  ;fix the period at a year
 ; parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(start5))
 ; parinfo[2].fixed = 1          
 ; result5= MPFITFUN('sin_func',bin_time5, bin_flux5 , bin_flux5err, start5, parinfo = parinfo) ;ICL
 ; period = 2 * !PI * result5(2)
;4-------------
;  print, 'working on DIRBE 4.5'
;  noise4 = fltarr(n_elements(bin_flux4)) + .027
;  start4 = [0.3,0.5,1./(2.*!Pi),result5(3)]
  ;fix the period at a year
 ; parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(start4))
;  parinfo[2].fixed = 1          
;  parinfo[3].fixed = 1
 ; result4= MPFITFUN('sin_func',bin_time4, bin_flux4 , bin_flux4err, start4, parinfo = parinfo) ;ICL
;  period = 2 * !PI * result4(2)
;3-------------
;  print, 'working on DIRBE 3.6'
;  noise3 = fltarr(n_elements(bin_flux3)) + .033
;  start3 = [0.1,0.1,1./(2.*!Pi),result5(3)]
 ; ;fix the period at a year
;  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(start3))
;  parinfo[2].fixed = 1          
;  parinfo[3].fixed = 1
 ; result3= MPFITFUN('sin_func',bin_time3, bin_flux3 , bin_flux3err, start3, parinfo = parinfo) ;ICL
 ; period = 2 * !PI * result3(2)
;


;-------------------------------------------------------------------
;IRAC data
;-------------------------------------------------------------------
;zodi from header not removed from the data
restore, filename = '/Users/jkrick/irac_darks/cryodarks_100_nozodisub.sav'
cryo_nozodi_ch1_sclksort = ch1_sclksort
cryo_nozodi_ch1_meansort = ch1_meansort ;/ median(ch1_meansort)
cryo_nozodi_ch1_sigmasort = ch1_sigmasort;/ median(ch1_meansort)
cryo_nozodi_ch1_zodisort = ch1_zodisort
cryo_nozodi_ch2_sclksort = ch2_sclksort
cryo_nozodi_ch2_meansort = ch2_meansort;/ median(ch2_meansort)
cryo_nozodi_ch2_sigmasort = ch2_sigmasort;/ median(ch2_meansort)
cryo_nozodi_ch2_zodisort = ch2_zodisort
cryo_nozodi_ch3_sclksort = ch3_sclksort
cryo_nozodi_ch3_meansort = ch3_meansort;/ median(ch3_meansort)
cryo_nozodi_ch3_sigmasort = ch3_sigmasort;/ median(ch3_meansort)
cryo_nozodi_ch3_zodisort = ch3_zodisort
cryo_nozodi_ch4_sclksort = ch4_sclksort
cryo_nozodi_ch4_meansort = ch4_meansort;/ median(ch4_meansort)
cryo_nozodi_ch4_sigmasort = ch4_sigmasort;/ median(ch4_meansort)
cryo_nozodi_ch4_zodisort = ch4_zodisort


;the warm dark measurements
;zodi from header not removed from the data
restore, filename = '/Users/jkrick/irac_darks/warmdarks_100_nozodisub.sav'
nozodi_ch1_sclksort = ch1_sclksort
nozodi_ch1_meansort = ch1_meansort ;/ median(ch1_meansort)
nozodi_ch1_sigmasort = ch1_sigmasort;/ median(ch1_meansort)
nozodi_ch1_zodisort = ch1_zodisort
nozodi_ch2_sclksort = ch2_sclksort
nozodi_ch2_meansort = ch2_meansort;/ median(ch2_meansort)
nozodi_ch2_sigmasort = ch2_sigmasort;/ median(ch2_meansort)
nozodi_ch2_zodisort = ch2_zodisort

;outlier check on IRAC data
;resistant_mean,nozodi_ch2_meansort,2, rmean, rsig
;rsigma = robust_sigma(nozodi_ch2_meansort)
;print, 'ch2 outliers', rmean, rsigma
;g = where(nozodi_ch2_meansort gt (rmean - 5*rsigma) and nozodi_ch2_meansort lt (rmean + 5*rsigma), gcount)
;print, 'ch2 gcount, orig', gcount, n_elements(nozodi_ch2_meansort)

;delta of measured vs. model flux is too large
;warm ch1
;plothist,nozodi_ch1_meansort - (nozodi_ch1_zodisort + 0.665048) , bin = 0.01
g = where(abs(nozodi_ch1_meansort - (nozodi_ch1_zodisort + 0.665048)) le 0.02, gcount)
print, 'ch1 gcount, orig', gcount, n_elements(nozodi_ch1_meansort)

nozodi_ch1_meansort = nozodi_ch1_meansort[g]
nozodi_ch1_sclksort = nozodi_ch1_sclksort[g]
nozodi_ch1_sigmasort = nozodi_ch1_sigmasort[g]
nozodi_ch1_zodisort = nozodi_ch1_zodisort[g]
;ch1_camparr = ch1_camparr[g]
;ch1_aorkey = ch1_aorkey[g]

;warm ch2
;plothist,nozodi_ch2_meansort - (nozodi_ch2_zodisort + 0.117269) , bin = 0.01
g = where(abs(nozodi_ch2_meansort - (nozodi_ch2_zodisort + 0.117269)) le 0.02, gcount)
print, 'ch2 gcount, orig', gcount, n_elements(nozodi_ch2_meansort)

nozodi_ch2_meansort = nozodi_ch2_meansort[g]
nozodi_ch2_sclksort = nozodi_ch2_sclksort[g]
nozodi_ch2_sigmasort = nozodi_ch2_sigmasort[g]
nozodi_ch2_zodisort = nozodi_ch2_zodisort[g]
ch2_camparr = ch2_camparr[g]
ch2_aorkey = ch2_aorkey[g]
;----------------------
;cryo ch1
;plothist,cryo_nozodi_ch1_meansort - (cryo_nozodi_ch1_zodisort + 0.0227086) , bin = 0.01
g = where(abs(cryo_nozodi_ch1_meansort - (cryo_nozodi_ch1_zodisort +0.0227086)) le 0.05, gcount)
print, 'ch1 gcount, orig', gcount, n_elements(cryo_nozodi_ch1_meansort)

cryo_nozodi_ch1_meansort = cryo_nozodi_ch1_meansort[g]
cryo_nozodi_ch1_sclksort = cryo_nozodi_ch1_sclksort[g]
cryo_nozodi_ch1_sigmasort = cryo_nozodi_ch1_sigmasort[g]
cryo_nozodi_ch1_zodisort = cryo_nozodi_ch1_zodisort[g]

;cryo ch2
;plothist,cryo_nozodi_ch2_meansort - (cryo_nozodi_ch2_zodisort + 0.0627375) , bin = 0.01
g = where(abs(cryo_nozodi_ch2_meansort - (cryo_nozodi_ch2_zodisort +0.0627375)) le 0.025, gcount)
print, 'ch2 gcount, orig', gcount, n_elements(cryo_nozodi_ch2_meansort)

cryo_nozodi_ch2_meansort = cryo_nozodi_ch2_meansort[g]
cryo_nozodi_ch2_sclksort = cryo_nozodi_ch2_sclksort[g]
cryo_nozodi_ch2_sigmasort = cryo_nozodi_ch2_sigmasort[g]
cryo_nozodi_ch2_zodisort = cryo_nozodi_ch2_zodisort[g]

;cryo ch3
;plothist,cryo_nozodi_ch3_meansort - (cryo_nozodi_ch3_zodisort -1.77718) , bin = 0.01
g = where(abs(cryo_nozodi_ch3_meansort - (cryo_nozodi_ch3_zodisort -1.77718)) le 0.3, gcount)
print, 'ch3 gcount, orig', gcount, n_elements(cryo_nozodi_ch3_meansort)

cryo_nozodi_ch3_meansort = cryo_nozodi_ch3_meansort[g]
cryo_nozodi_ch3_sclksort = cryo_nozodi_ch3_sclksort[g]
cryo_nozodi_ch3_sigmasort = cryo_nozodi_ch3_sigmasort[g]
cryo_nozodi_ch3_zodisort = cryo_nozodi_ch3_zodisort[g]

;cryo ch4
;plothist,cryo_nozodi_ch4_meansort - (cryo_nozodi_ch4_zodisort +1.85452) , bin = 0.05, xrange = [-1,1]
g = where(abs(cryo_nozodi_ch4_meansort - (cryo_nozodi_ch4_zodisort +1.85452)) le 0.5, gcount)
print, 'ch4 gcount, orig', gcount, n_elements(cryo_nozodi_ch4_meansort)

cryo_nozodi_ch4_meansort = cryo_nozodi_ch4_meansort[g]
cryo_nozodi_ch4_sclksort = cryo_nozodi_ch4_sclksort[g]
cryo_nozodi_ch4_sigmasort = cryo_nozodi_ch4_sigmasort[g]
cryo_nozodi_ch4_zodisort = cryo_nozodi_ch4_zodisort[g]
;-------------------------------------------------------------------
;-------------------------------------------------------------------
;begin the plotting

;ch1
if channel eq 1 then begin
   print, 'working on ch', channel
   deltatime = (nozodi_ch1_sclksort(0) - cryo_nozodi_ch1_sclksort(0)) / 86400/365
   dcryo = mean(cryo_nozodi_ch1_meansort)- mean( cryo_nozodi_ch1_zodisort) 
   print, 'ch1 dcryo', dcryo
  dwarm = mean(nozodi_ch1_meansort)- mean( nozodi_ch1_zodisort) 
   print, 'ch1 dwarm', dwarm

   p4cryo = errorplot(((cryo_nozodi_ch1_sclksort - cryo_nozodi_ch1_sclksort(0)) / 86400 / 365) , cryo_nozodi_ch1_meansort , cryo_nozodi_ch1_sigmasort,'1gs', sym_filled = 1, sym_size = 0.6,xtitle = 'Time in Years Starting from Dec 01 2003', ytitle = 'Surface Brightness (Mjy/sr)', title = '3.6 micron', xrange = [0, 8.5], yrange = [0.02, 0.1], name = 'IRAC cryo', position = [0.15, .35, .90, .95], xshowtext = 0)
   p4dirbe = errorplot(bin_time3+0.15, bin_flux3-0.12, bin_flux3err,'1s', sym_filled = 1, sym_size = 0.6, /overplot, name = 'DIRBE 3.5$\mu$m', color = 'black') ;/median(bin_flux4)
   p4modelc = plot(((cryo_nozodi_ch1_sclksort - cryo_nozodi_ch1_sclksort(0)) / 86400 / 365) , cryo_nozodi_ch1_zodisort + dcryo, color = 'black', name = 'cryo model',/overplot) ;/ median(cryo_nozodi_ch1_zodisort)

   yaxis = axis('Y', location = [8.5,0],  tickvalues = [0.66,0.68,0.7,.72,.74]-0.64, tickname = ['0.66','0.68','0.7','0.72','0.74'], target = p4cryo, textpos = 1.0, ticklayout = 1, ticklen = 0., text_color = 'blue')

   
   p4warm = errorplot(((nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400 / 365) + deltatime, nozodi_ch1_meansort, nozodi_ch1_sigmasort,'1bs', axis_style=0, sym_filled = 1, sym_size = 0.6, name = 'IRAC warm', xrange = [5.5, 8.5], position = [0.666, 0.35, 0.90, 0.95], /current, yrange = [0.66, 0.74])

   

    p4modelw = plot(((nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400 / 365) +deltatime, nozodi_ch1_zodisort+dwarm , color = 'black', name = 'Kelsall et al. 1998',axis_style = 0, /overplot) ;/median(nozodi_ch1_zodisort)

   l = legend(target = [p4cryo,p4warm,p4dirbe, p4modelw], position = [5, 0.045], /data)

;new plot with residuals
   residcryo  =  cryo_nozodi_ch1_meansort - (cryo_nozodi_ch1_zodisort + dcryo)
   residwarm = nozodi_ch1_meansort - (nozodi_ch1_zodisort + dwarm)

  print, 'residcryo', mean(residcryo), stddev(residcryo)
   print, 'residwarm', mean(residwarm), stddev(residwarm)

   r4c = plot(((cryo_nozodi_ch1_sclksort - cryo_nozodi_ch1_sclksort(0)) / 86400 / 365) , (residcryo/(cryo_nozodi_ch1_zodisort + dcryo))*100.,'1gs', sym_filled = 1, sym_size = 0.6,xtitle = 'Time in Years Starting from Dec 01 2003', ytitle = 'Residuals (%)', xrange = [0, 8.5], name = 'IRAC cryo', yrange = [-40,40], position = [0.15, .1, .90, 0.3], /current)

  yaxis = axis('Y', location = [8.5,0],  tickvalues = [-2,-1,0,1, 2]*20., tickname = ['-2','-1','0','1','2'], target = r4c, textpos = 1.0, ticklayout = 1, ticklen = 0., text_color = 'blue')

   r4w = plot(((nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400 / 365) + deltatime, (residwarm/(nozodi_ch1_zodisort+dwarm))*100.,'1bs', sym_filled = 1, sym_size = 0.6,name = 'IRAC warm', xrange = [5.5, 8.5], axis_style = 0,position = [0.666,0.1,0.90,0.3],/current, yrange = [-2,2])

;   r4c = plot(((cryo_nozodi_ch1_sclksort - cryo_nozodi_ch1_sclksort(0)) / 86400 / 365) , residcryo*1000.,'1gs', sym_filled = 1, sym_size = 0.6,xtitle = 'Time in Years Starting from Dec 01 2003', ytitle = 'Residuals (KJy/sr)', xrange = [0, 8.5], name = 'IRAC cryo', position = [0.15, .1, .95, 0.3], yrange = [-20,20],/current)
;   r4w = plot(((nozodi_ch1_sclksort - nozodi_ch1_sclksort(0)) / 86400 / 365) + deltatime, residwarm*1000,'1bs', sym_filled = 1, sym_size = 0.6,name = 'IRAC warm', /overplot)


 
endif

;---------
;ch2
if channel eq 2 then begin
   print, 'working on ch', channel

   p4cryo = errorplot(((cryo_nozodi_ch2_sclksort - cryo_nozodi_ch2_sclksort(0)) / 86400 / 365) , cryo_nozodi_ch2_meansort , cryo_nozodi_ch2_sigmasort,'1gs', sym_filled = 1, sym_size = 0.6,xtitle = 'Time in Years Starting from Dec 01 2003', ytitle = 'Surface Brightness (Mjy/sr)', title = '4.5 micron', xrange = [0, 8.5], name = 'IRAC cryo', position = [0.15, .35, .95, .95], xshowtext = 0)

   deltatime = (nozodi_ch2_sclksort(0) - cryo_nozodi_ch2_sclksort(0)) / 86400/365

   p4warm = errorplot(((nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400 / 365) + deltatime, nozodi_ch2_meansort, nozodi_ch2_sigmasort,'1bs', sym_filled = 1, sym_size = 0.6,/overplot, name = 'IRAC warm')

   p4dirbe = errorplot(bin_time4+0.15, bin_flux4-0.18, bin_flux4err,'1s', sym_filled = 1, sym_size = 0.6, /overplot, name = 'DIRBE 4.9$\mu$m', color = 'black') ;/median(bin_flux4)

   dcryo = mean(cryo_nozodi_ch2_meansort)- mean( cryo_nozodi_ch2_zodisort) 
   print, 'dcryo', dcryo

   p4modelc = plot(((cryo_nozodi_ch2_sclksort - cryo_nozodi_ch2_sclksort(0)) / 86400 / 365) , cryo_nozodi_ch2_zodisort + dcryo, color = 'black', name = 'cryo model',/overplot) ;/ median(cryo_nozodi_ch2_zodisort)

   dwarm = mean(nozodi_ch2_meansort)- mean( nozodi_ch2_zodisort) 
   print, 'dwarm', dwarm

   p4modelw = plot(((nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400 / 365) +deltatime, nozodi_ch2_zodisort+dwarm-0.005, color = 'black', name = 'Kelsall et al. 1998',/overplot) ;/median(nozodi_ch2_zodisort)
   l = legend(target = [p4cryo,p4warm,p4dirbe, p4modelw], position = [0.05, 0.38], /data)

;fit warm with a sine curve

  
 xw = ((nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400 / 365) + deltatime
 yw = nozodi_ch2_meansort
 noise = nozodi_ch2_sigmasort
 start4 = [0.3,0.05,1./(2.*!Pi),5.]

 for test = 0, n_elements(xw) - 1 do begin
    print, test, xw(test), yw(test)
    ;36 and 39 are the outliers
 endfor

xw = [xw[0:35],xw[37:38],xw[40:*]]
yw = [yw[0:35],yw[37:38],yw[40:*]]
noise = [noise[0:35],noise[37:38],noise[40:*]]
 resultw= MPFITFUN('sin_func',xw, yw , noise, start4) ;ICL
; r4wf = plot(xw, resultw(0)*sin(xw/resultw(2) + resultw(3)) + resultw(1),/overplot, color = 'red')

sub = yw - (resultw(0)*sin(xw/resultw(2) + resultw(3)) + resultw(1))

;new plot with residuals
   residcryo  =  cryo_nozodi_ch2_meansort - (cryo_nozodi_ch2_zodisort + dcryo)
   residwarm = nozodi_ch2_meansort - (nozodi_ch2_zodisort + dwarm ) 

;   r4c = plot(((cryo_nozodi_ch2_sclksort - cryo_nozodi_ch2_sclksort(0)) / 86400 / 365) , (residcryo/(cryo_nozodi_ch2_zodisort + dcryo))*100.,'1gs', sym_filled = 1, sym_size = 0.6,xtitle = 'Time in Years Starting from Dec 01 2003', ytitle = 'Residuals (%)', xrange = [0, 8.5], yrange = [-8,6], name = 'IRAC cryo', position = [0.15, .1, .95, 0.3], /current)
;   r4w = plot(((nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400 / 365) + deltatime, (residwarm/(nozodi_ch2_zodisort+dwarm))*100.,'1bs', sym_filled = 1, sym_size = 0.6,name = 'IRAC warm', /overplot)
;   subplot = plot(xw, (sub/(nozodi_ch2_zodisort+dwarm))*100., color = 'red', /overplot)

   r4c = plot(((cryo_nozodi_ch2_sclksort - cryo_nozodi_ch2_sclksort(0)) / 86400 / 365) , residcryo*1000.,'1gs', sym_filled = 1, sym_size = 0.6,xtitle = 'Time in Years Starting from Dec 01 2003', ytitle = 'Residuals (KJy/sr)', xrange = [0, 8.5], name = 'IRAC cryo', position = [0.15, .1, .95, 0.3], /current)
   r4w = plot(((nozodi_ch2_sclksort - nozodi_ch2_sclksort(0)) / 86400 / 365) + deltatime, residwarm*1000,'1bs', sym_filled = 1, sym_size = 0.6,name = 'IRAC warm', /overplot)


endif

;fit residwarm with a linear function, it will have a slope.
; xr = (cryo_nozodi_ch2_sclksort - cryo_nozodi_ch2_sclksort(0)) / 86400 / 365
; yr = (residcryo/(cryo_nozodi_ch2_zodisort + dcryo))*100.
; noise = fltarr(n_elements(yr)) + 1.0
; start = [0.5, 0]
; result= MPFITFUN('linear',xr, yr , noise, start) ;ICL
; r4cf =plot(xr, result(0)*xr +result(1),/overplot, color = 'green')
 

;---------
;ch3
if channel eq 3 then begin
   print, 'working on ch', channel
   p4cryo = errorplot(((cryo_nozodi_ch3_sclksort - cryo_nozodi_ch3_sclksort(0)) / 86400 / 365) , cryo_nozodi_ch3_meansort , cryo_nozodi_ch3_sigmasort,'1gs', sym_filled = 1, sym_size = 0.6,xtitle = 'Time in Years Starting from Dec 01 2003', ytitle = 'Surface Brightness (Mjy/sr)', title = '5.8 micron', xrange = [0, 8.5], name = 'IRAC cryo', position = [0.15, .35, .95, .95], xshowtext = 0)
   dcryo = mean(cryo_nozodi_ch3_meansort)- mean( cryo_nozodi_ch3_zodisort) 
   print, 'dcryo', dcryo
   p4modelc = plot(((cryo_nozodi_ch3_sclksort - cryo_nozodi_ch3_sclksort(0)) / 86400 / 365) , cryo_nozodi_ch3_zodisort + dcryo, color = 'black', name = 'Kelsall et al. 1998',/overplot) ;/ median(cryo_nozodi_ch3_zodisort)
   l = legend(target = [p4cryo, p4modelc], position = [5.1, 0.5], /data)

;new plot with residuals
   residcryo  =  cryo_nozodi_ch3_meansort - (cryo_nozodi_ch3_zodisort + dcryo)

 ;  r4c = plot(((cryo_nozodi_ch3_sclksort - cryo_nozodi_ch3_sclksort(0)) / 86400 / 365) , (residcryo/(cryo_nozodi_ch3_zodisort + dcryo))*100.,'1gs', sym_filled = 1, sym_size = 0.6,xtitle = 'Time in Years Starting from Dec 01 2003', ytitle = 'Residuals (%)', xrange = [0, 8.5], name = 'IRAC cryo', yrange = [-50, 50], position = [0.15, .1, .95, 0.3], /current)
 r4c = plot(((cryo_nozodi_ch3_sclksort - cryo_nozodi_ch3_sclksort(0)) / 86400 / 365) , residcryo*1000.,'1gs', sym_filled = 1, sym_size = 0.6,xtitle = 'Time in Years Starting from Dec 01 2003', ytitle = 'Residuals (kJy/sr)', xrange = [0, 8.5], name = 'IRAC cryo', position = [0.15, .1, .95, 0.3], /current)
endif
;---------

;ch4
if channel eq 4 then begin
   print, 'working on ch', channel
   p4cryo = errorplot(((cryo_nozodi_ch4_sclksort - cryo_nozodi_ch4_sclksort(0)) / 86400 / 365) , cryo_nozodi_ch4_meansort , cryo_nozodi_ch4_sigmasort,'1gs', sym_filled = 1, sym_size = 0.6,xtitle = 'Time in Years Starting from Dec 01 2003', ytitle = 'Surface Brightness (Mjy/sr)', title = '8.0 micron', xrange = [0, 8.5], name = 'IRAC cryo', position = [0.15, .35, .95, .95], xshowtext = 0)
   p4dirbe = errorplot(bin_time5+0.15, bin_flux5-4.9, bin_flux5err,'1s', sym_filled = 1, sym_size = 0.6, /overplot, name = 'DIRBE 12$\mu$m', color = 'black') ;/median(bin_flux4)  
   dcryo = mean(cryo_nozodi_ch4_meansort)- mean( cryo_nozodi_ch4_zodisort) 
   print, 'dcryo', dcryo
   p4modelc = plot(((cryo_nozodi_ch4_sclksort - cryo_nozodi_ch4_sclksort(0)) / 86400 / 365) , cryo_nozodi_ch4_zodisort + dcryo, color = 'black', name = 'Kelsall et al. 1998',/overplot) ;/ median(cryo_nozodi_ch4_zodisort)
   l = legend(target = [p4cryo, p4dirbe,p4modelc], position = [5.1, 8.4], /data)

;new plot with residuals
   residcryo  =  cryo_nozodi_ch4_meansort - (cryo_nozodi_ch4_zodisort + dcryo)

;   r4c = plot(((cryo_nozodi_ch4_sclksort - cryo_nozodi_ch4_sclksort(0)) / 86400 / 365) , (residcryo/(cryo_nozodi_ch4_zodisort + dcryo))*100.,'1gs', sym_filled = 1, sym_size = 0.6,xtitle = 'Time in Years Starting from Dec 01 2003', ytitle = 'Residuals (%)', xrange = [0, 8.5], name = 'IRAC cryo', yrange = [-5, 5], position = [0.15, .1, .95, 0.3], /current)
   r4c = plot(((cryo_nozodi_ch4_sclksort - cryo_nozodi_ch4_sclksort(0)) / 86400 / 365) , residcryo*1000.,'1gs', sym_filled = 1, sym_size = 0.6,xtitle = 'Time in Years Starting from Dec 01 2003', ytitle = 'Residuals (kJy/sr)', xrange = [0, 8.5], name = 'IRAC cryo', position = [0.15, .1, .95, 0.3], /current)
endif

end

