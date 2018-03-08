pro plot_55cnc_ch2_39524608_v2
  colorarr = ['black', 'red', 'orange', 'green', 'blue', 'deep_pink', 'magenta', 'medium_purple', 'hot_pink', 'light_pink', 'rosy_brown', 'chocolate', 'saddle_brown', 'maroon', 'orange_red', 'dark_orange', 'peach_puff', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'light_cyan', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green']

   fits_read, '/Users/jkrick/idlbin/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_120124.fits', pmapdata, pmapheader

; c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', aspect_ratio = 1, xrange = [0,500], yrange = [0,500])


  restore, '/Users/jkrick/irac_warm/pcrs_planets/55cnc/55cnc_phot_ch2_39524608.sav'
  ;aorname_55cnc = ['r43981312','r43981568','r43981824','r43981056'] ;ch2
  aorname_55cnc  = ['r39524608'] ;ch2

  
  ;utmjd_center = [double(55947.20505),double(55949.41467), double(55957.51661),double(55944.25889)]  
 ; utmjd_center = [double(55944.25889), double(55947.20505),double(55949.41467), double(55957.51661)]  

  duration = 105.7              ;in minutes
  duration = duration /60./24.  ; in days
  period = 0.73654 ;days
  

;set time = 0.5 in the middle of the eclipse
;     bin_time = bin_time -tt
;  for a = 0, n_elements(aorname_55cnc) - 1 do  begin
;     test =  ((AOR55cnc[a].bmjdarr - utmjd_center[a]) / period) + 0.5
     
;     ;if a eq 1 then begin
;        ;for t = 0, 100*63,63  do begin
;        ;   print, t/63, ' ', AOR55cnc[1].bmjdarr[t], ' ',test[t], ' ', utmjd_center[a], format = '(F0,A, F0, A, F0, A, F0)'
;       ; endfor
;     ;endif

;     AOR55cnc[a].bmjdarr = ((AOR55cnc[a].bmjdarr - utmjd_center[a]) / period) + 0.5

;  endfor


; 
;  for a = 0, n_elements(aorname_55cnc) - 1 do begin
;     xcen500 = 500.* (AOR55cnc[a].xcen - 14.5)
;     ycen500 = 500.* (AOR55cnc[a].ycen - 14.5)
;     an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14],/overplot)
; endfor

;-----
  
 ;        am = plot( (AOR55cnc[0].timearr - AOR55cnc[0].timearr(0))/60./60., AOR55cnc[0].xcntrd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14], xtitle = 'Time(hrs)', ytitle = 'X CNTRD', yrange = [14.2, 15.4])
  
;------
;
;        am = plot( (AOR55cnc[0].timearr - AOR55cnc[0].timearr(0))/60./60., AOR55cnc[0].ycntrd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14], xtitle = 'Time(hrs)', ytitle = 'Y CNTRD',  yrange = [14.2, 15.4])
 ;-----
  
;        am = plot( (AOR55cnc[0].timearr - AOR55cnc[0].timearr(0))/60./60., AOR55cnc[0].xgcntrd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14], xtitle = 'Time(hrs)', ytitle = 'X CNTRD', yrange = [14.2, 15.4])

;------

;        am = plot( (AOR55cnc[0].timearr - AOR55cnc[0].timearr(0))/60./60., AOR55cnc[0].ygcntrd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14], xtitle = 'Time(hrs)', ytitle = 'Y CNTRD',  yrange = [14.2, 15.4])
 

;------
;        am = plot( (AOR55cnc[0].timearr - AOR55cnc[0].timearr(0))/60./60., AOR55cnc[0.bkgd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Background')
 
;        am = plot((AOR55cnc[0].timearr - AOR55cnc[0].timearr(0))/60./60., AOR55cnc[0].bkgd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black', xtitle = 'Time (hrs)', ytitle = 'Background (Mjy/sr)')
;        am = plot((AOR55cnc[0].timearr - AOR55cnc[0].timearr(0))/60./60., AOR55cnc[0].xcen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black', xtitle = 'Time (hrs)', ytitle = 'X', yrange = [14.4,14.9])
;  am = plot((AOR55cnc[0].timearr - AOR55cnc[0].timearr(0))/60./60., AOR55cnc[0].ycen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black', xtitle = 'Time (hrs)', ytitle = 'Y', yrange = [14.9,15.1])
; am = plot((AOR55cnc[0].timearr - AOR55cnc[0].timearr(0))/60./60., AOR55cnc[0].flux / mean(AOR55cnc[0].flux),'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black', xtitle = 'Time (hrs)', ytitle = 'Normalized BCD fluxes', yrange = [0.9, 1.1])
;------

  for a = 0, n_elements(aorname_55cnc) - 1 do begin
;sort by time
     s = sort(AOR55cnc[a].bmjdarr)
     sorttime = AOR55cnc[a].timearr[s]
     sortxcen= AOR55cnc[a].xcen[s]
     sortycen = AOR55cnc[a].ycen[s]
     sortflux = AOR55cnc[a].flux[s]
     sortcorrflux = AOR55cnc[a].corrflux[s]
     sortcorrfluxerr = AOR55cnc[a].corrfluxerr[s]
     sortbmjd = AOR55cnc[a].bmjdarr[s]
     sortutcs = AOR55cnc[a].utcsarr[s]
     sortbkgd = AOR55cnc[a].bkgd[s]
     sortbkgderr = AOR55cnc[a].bkgderr[s]
     sortxgcntrd = AOR55cnc[a].xgcntrd[s]
     sortygcntrd = AOR55cnc[a].ygcntrd[s]
     sortxcntrd = AOR55cnc[a].xcntrd[s]
     sortycntrd = AOR55cnc[a].ycntrd[s]

      ;remove outliers
     good = where(sortxcen lt mean(sortxcen) + 2.5*stddev(sortxcen) and sortxcen gt mean(sortxcen) -2.5*stddev(sortxcen) and sortycen lt mean(sortycen) +3.0*stddev(sortycen) and sortycen gt mean(sortycen) - 3.0*stddev(sortycen),ngood_pmap, complement=bad) 
     print, 'bad ',n_elements(bad)
     xarr = sortxcen[good]
     yarr = sortycen[good]
     timearr = sorttime[good]
     flux = sortflux[good]
     corrflux = sortcorrflux[good]
     corrfluxerr = sortcorrfluxerr[good]
     bmjdarr = sortbmjd[good]
     bkgdarr = sortbkgd[good]
     bkgderrarr = sortbkgderr[good]
     xgcntrd = sortxgcntrd[good]
     ygcntrd = sortygcntrd[good]
     xcntrd = sortxcntrd[good]
     ycntrd = sortycntrd[good]
;------------------------------------------------------------------
; binning
     nframes = bin_level
     nfits = long(n_elements(corrflux)) / nframes
     print, 'nframes, nfits', n_elements(corrflux), nframes, nfits
     bin_flux = dblarr(nfits)
     bin_corrflux= dblarr(nfits)
     bin_corrfluxerr= dblarr(nfits)
     bin_timearr = dblarr(nfits)
     bin_bmjdarr = dblarr(nfits)
     bin_bkgd = dblarr(nfits)
     bin_bkgderr = dblarr(nfits)
     bin_xcen = dblarr(nfits)
     bin_ycen = dblarr(nfits)
     bin_xgcntrd= dblarr(nfits)
     bin_ygcntrd= dblarr(nfits)
     bin_xcntrd= dblarr(nfits)
     bin_ycntrd= dblarr(nfits)

     for si = 0L, long(nfits) - 1L do begin
                                ;print, 'working on si', n_elements(AOR55cnc[a].corrflux), si, si*nframes, si*nframes + (nframes - 1)
        idata = corrflux[si*nframes:si*nframes + (nframes - 1)]
        id = flux[si*nframes:si*nframes + (nframes - 1)]
        idataerr = corrfluxerr[si*nframes:si*nframes + (nframes - 1)]
        ib = bkgdarr[si*nframes:si*nframes + (nframes - 1)]
        iberr = bkgderrarr[si*nframes:si*nframes + (nframes - 1)]
       ; if si lt 2 then print, 'ib', ib
       ; print, 'mean', mean(ib,/nan)
        bin_flux[si] = mean(id, /nan)
        bin_corrflux[si] = mean(idata,/nan)
        bin_corrfluxerr[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
        bin_timearr[si]= mean(timearr[si*nframes:si*nframes + (nframes - 1)])
        bin_bmjdarr[si]= mean(bmjdarr[si*nframes:si*nframes + (nframes - 1)])
        bin_bkgd[si] = mean(ib,/nan)
        bin_bkgderr[si] = mean(iberr,/nan)
        bin_xcen[si]= mean(xarr[si*nframes:si*nframes + (nframes - 1)])
        bin_ycen[si]= mean(yarr[si*nframes:si*nframes + (nframes - 1)])
        bin_xgcntrd[si] = mean(xgcntrd[si*nframes:si*nframes + (nframes - 1)])
        bin_ygcntrd[si] = mean(ygcntrd[si*nframes:si*nframes + (nframes - 1)])
        bin_xcntrd[si] = mean(xcntrd[si*nframes:si*nframes + (nframes - 1)])
        bin_ycntrd[si] = mean(ycntrd[si*nframes:si*nframes + (nframes - 1)])
     endfor                     ;for each fits image
     
 

;----------------------------
;try binning
;histogram numberarr by 64 as a technique to bin
     bin_level = 63L            ; 63L;*60
     h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
     print, 'omin', om, 'nh', n_elements(h)

;mean together the flux values in each phase bin
  binxcen = dblarr(n_elements(h))
  binnum = dblarr(n_elements(h))
  binycen= dblarr(n_elements(h))
  binbkgd =  dblarr(n_elements(h))
  binflux=  dblarr(n_elements(h))
  binback=  dblarr(n_elements(h))
  c = 0
  for j = 0L, n_elements(h) - 1 do begin

;get rid of the bins with no values and low numbers, meaning low overlap
     if (ri[j+1] gt ri[j])  then begin
        ;print, 'binning together', n_elements(numberarr[ri[ri[j]:ri[j+1]-1]])
        ;print, 'binning', numberarr[ri[ri[j]:ri[j+1]-1]]
        
        meanclip, xarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        binxcen[c] = meanx; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

        meanclip, yarr[ri[ri[j]:ri[j+1]-1]], meany, sigmay
        binycen[c] = meany; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

        meanclip, skyarr[ri[ri[j]:ri[j+1]-1]], meansky, sigmasky
        binbkgd[c] = meansky; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

        binnum[c] = mean(numberarr[ri[ri[j]:ri[j+1]-1]])

        meanclip, fluxarr[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
        binflux[c] = meanflux; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

        meanclip, backarr[ri[ri[j]:ri[j+1]-1]], meanback, sigmaback
        binback[c] = meanback; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

        c = c + 1
        ;print, 'testing', j, phasearr[ri[ri[j]:ri[j+1]-1]]
     endif
  endfor

  binxcen = binxcen[0:c-1]
  binycen = binycen[0:c-1]
  binnum = binnum[0:c-1]
  binbkgd = binbkgd[0:c-1]
  binflux = binflux[0:c-1]
  binback = binback[0:c-1]

;----------------------------------
;try getting rid of flux outliers.
  ;do some running mean with clipping
     start = 0
     print, 'nflux', n_elements(bin_flux)
     for ni = 50, n_elements(bin_flux) -1, 50 do begin
        meanclip,flux[start:ni], m, s, subs = subs ;,/verbose
                                ;print, 'good', subs+start
                                ;now keep a list of the good ones
        if ni eq 50 then good_ni = subs+start else good_ni = [good_ni, subs+start]
        start = ni
     endfor
                                ;see if it worked
     ;xarr = xarr[good_ni]
     ;yarr = yarr[good_ni]
     bin_timearr = bin_timearr[good_ni]
     bin_bmjdarr = bin_bmjdarr[good_ni]
     bin_corrfluxerr = bin_corrfluxerr[good_ni]
     bin_flux = bin_flux[good_ni]
     bin_corrflux = bin_corrflux[good_ni]
     bin_bkgd = bin_bkgd[good_ni]
     bin_bkgderr = bin_bkgderr[good_ni]
     bin_xcen = bin_xcen[good_ni]
     bin_ycen = bin_ycen[good_ni]
     bin_xgcntrd = bin_xgcntrd[good_ni]
     bin_ygcntrd = bin_ygcntrd[good_ni]
     bin_xcntrd = bin_xcntrd[good_ni]
     bin_ycntrd = bin_ycntrd[good_ni]

     print, 'nflux', n_elements(bin_flux)
     


 xcen_coeff = [1.5,0.5,-0.1, -0.8]
 ycen_coeff = [1.0,0,0,-0.7]
;     if a eq 0 then begin
;        am = plot( (AOR55cnc[a].timearr - AOR55cnc[a].timearr(0))/60./60., AOR55cnc[a].flux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14], xtitle = 'Time(hrs)', ytitle = 'Flux', xrange = [0,6], yrange = [3,5])
;        am = plot( (AOR55cnc[a].timearr - AOR55cnc[a].timearr(0))/60./60., AOR55cnc[a].corrflux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black', /overplot)
        
;        am = plot(bin_bmjdarr, bin_corrflux/mean(bin_corrflux,/nan),  '1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], xtitle = 'Phase', ytitle = 'Corrected Flux')

;        ab = plot(bin_bmjdarr  , bin_flux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14], xtitle = 'period', ytitle = 'Flux (Jy)')
;        ax = plot(bin_bmjdarr, bin_xcen ,'6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], xtitle = 'Period', ytitle = 'Xcen');, xrange = [0.3,0.65])
;          ay = plot(bin_bmjdarr, bin_ycen ,'6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], xtitle = 'Period', ytitle = 'Ycen');, xrange = [0.3,0.65])
;        axc = plot(bin_bmjdarr, bin_xcntrd ,'6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], xtitle = 'Period', ytitle = 'X CNTRD');, xrange = [0.3,0.65])
;          ayc = plot(bin_bmjdarr, bin_ycntrd ,'6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], xtitle = 'Period', ytitle = 'Y CNTRD');, xrange = [0.3,0.65])
;        axg = plot(bin_bmjdarr, bin_xgcntrd ,'6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], xtitle = 'Period', ytitle = 'X GCNTRD');, xrange = [0.3,0.65])
;          ayg = plot(bin_bmjdarr, bin_ygcntrd ,'6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], xtitle = 'Period', ytitle = 'Y GCNTRD');, xrange = [0.3,0.65])

        af = plot((bin_timearr - bin_timearr(0))/60./60., bin_flux/mean(bin_flux), '6r1s', sym_size = 0.2,   sym_filled = 1, color = 'black', xtitle = 'Time (hrs)', ytitle = 'Normalized Flux', yrange = [0.991, 1.006]);, xrange = [0.3,0.65])

        af = plot((bin_timearr - bin_timearr(0))/60./60., bin_bkgd, '6r1s', sym_size = 0.2,   sym_filled = 1, color = 'black', xtitle = 'Time (hrs)', ytitle = 'Bkgd (MJy/sr)');, xrange = [0.3,0.65])

        af = plot((bin_timearr - bin_timearr(0))/60./60., bin_xcen, '6r1s', sym_size = 0.2,   sym_filled = 1, color = 'black', xtitle = 'Time (hrs)', ytitle = 'X', yrange = [14.4,14.9]);, xrange = [0.3,0.65])

        af = plot((bin_timearr - bin_timearr(0))/60./60., bin_ycen, '6r1s', sym_size = 0.2,   sym_filled = 1, color = 'black', xtitle = 'Time (hrs)', ytitle = 'Y', yrange = [14.9,15.1]);, xrange = [0.3,0.65])

;        ab = plot(bin_timearr, bin_bkgd, '6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], title = 'aperture [3,11-15.5]', xtitle = 'Period', ytitle = 'Bkgd');, xrange = [0.3,0.65])


   ;  print,'time', bin_bmjdarr - bin_bmjdarr(0) 
 

     endfor


end
