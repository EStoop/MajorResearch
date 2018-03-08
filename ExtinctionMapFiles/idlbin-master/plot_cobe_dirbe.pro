pro plot_cobe_dirbe
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
  


;------------------------------------------------------------------
  ;fit with a sine curve
;------------------------------------------------------------------
 ;P(0)*sin(X/P(2) + P(3)) + P(1)
;P(0) = amplitude,
;P(1) = y-intercept
;P(2) = x amplitude = period
;p(3) = phase

   
;5-------------
  start5 = [6.,13.,1./(2.*!Pi),2.6]
  ;fix the period at a year
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(start5))
  parinfo[2].fixed = 1          
  result5= MPFITFUN('sin_func',bin_time5, bin_flux5 , bin_flux5err, start5, parinfo = parinfo) ;ICL
  period = 2 * !PI * result5(2)
;4-------------
  noise4 = fltarr(n_elements(bin_flux4)) + .027
  start4 = [0.3,0.5,1./(2.*!Pi),result5(3)]
  ;fix the period at a year
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(start4))
  parinfo[2].fixed = 1          
  parinfo[3].fixed = 1
  result4= MPFITFUN('sin_func',bin_time4, bin_flux4 , bin_flux4err, start4, parinfo = parinfo) ;ICL
  period = 2 * !PI * result4(2)
;3-------------
  noise3 = fltarr(n_elements(bin_flux3)) + .033
  start3 = [0.1,0.1,1./(2.*!Pi),result5(3)]
  ;fix the period at a year
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(start3))
  parinfo[2].fixed = 1          
  parinfo[3].fixed = 1
  result3= MPFITFUN('sin_func',bin_time3, bin_flux3 , bin_flux3err, start3, parinfo = parinfo) ;ICL
  period = 2 * !PI * result3(2)


;
  samplex = findgen(100) / 100.
 p3 = errorplot(bin_time3, bin_flux3,bin_flux3err, '1*',  xtitle = 'Time (years)', ytitle = 'Intensity (MJy/sr)', title = 'DIRBE 3.5 micron')
    p3s = plot(samplex, result3(0)*sin(samplex/result3(2) + result3(3)) + result3(1), '2-',/overplot)

  p4 = errorplot(bin_time4, bin_flux4,bin_flux4err, '1*',  xtitle = 'Time (years)', ytitle = 'Intensity (MJy/sr)', title = 'DIRBE 4.9 micron')
  p4s = plot(samplex, result4(0)*sin(samplex/result4(2) + result4(3)) + result4(1), '2-',/overplot)

  p5 = errorplot(bin_time5, bin_flux5, bin_flux5err, '1*',  xtitle = 'Time (years)', ytitle = 'Intensity (MJy/sr)', title = 'DIRBE 12 micron');yrange = [12, 21],
  p5s = plot(samplex, result5(0)*sin(samplex/result5(2) + result5(3)) + result5(1), '2-',/overplot)




end
