pro plot_wasp39_ch2
  colorarr = ['blue_violet']
  
  fits_read, '/Users/jkrick/idlbin/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_120124.fits', pmapdata, pmapheader
  
  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp39/wasp39_phot_ch2.sav'
  aorname = 'r42804736'

; AORwasp39 = replicate({wasp39ob, ra:0D,dec:0D,xcen:dblarr(nfits),ycen:dblarr(nfits),flux:dblarr(nfits),fluxerr:dblarr(nfits), corrflux:dblarr(nfits), corrfluxerr:dblarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ', bmjdarr:dblarr(nfits), utcsarr:dblarr(nfits), bkgd:dblarr(nfits), bkgderr:dblarr(nfits)},n_elements(aorname))


;  utjd_center = 2455947.70505 ;2455948.07332 ; center of eclipse
;  jdtomjd = 2400000.5
;  utmjd_center = utjd_center - jdtomjd

 
;  duration = 105.7              ;in minutes
;  duration = duration /60./24.  ; in days
;  period = 0.73654 ;days
  

;set time = 0 in the middle of the transit
;     bin_time = bin_time -tt
 ;AOR55cnc.bmjdarr = AOR55cnc.bmjdarr - utmjd_center

;set time in orbital phase from 0 at eclipse to 1 at next transit
  ;    AOR55cnc.bmjdarr= AOR55cnc.bmjdarr/ (period)

 c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', aspect_ratio = 1, xrange = [0,500], yrange = [0,500])

 xcen500 = 500.* (AORwasp39[0].xcen - 14.5)
 ycen500 = 500.* (AORwasp39[0].ycen - 14.5)
 an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[0],/overplot)


;-----
 
 am = plot( (AORwasp39.timearr - AORwasp39.timearr(0))/60./60., AORwasp39.xcen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[0], xtitle = 'Time(hrs)', ytitle = 'X pix', yrange = [14.8, 15.3], xrange = [0,8])

;------
 
 am = plot( (AORwasp39.timearr - AORwasp39.timearr(0))/60./60., AORwasp39.ycen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[0], xtitle = 'Time(hrs)', ytitle = 'Y pix',  yrange = [14.8, 15.3], xrange = [0,8])

;------
 am = plot( (AORwasp39.timearr - AORwasp39.timearr(0))/60./60., AORwasp39.bkgd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[0], xtitle = 'Time(hrs)', ytitle = 'Background', xrange = [0,8])

;------

 
      ;remove outliers
 good = where(AORwasp39.xcen lt mean(AORwasp39.xcen) + 2.5*stddev(AORwasp39.xcen) and AORwasp39.xcen gt mean(AORwasp39.xcen) -2.5*stddev(AORwasp39.xcen) and AORwasp39.ycen lt mean(AORwasp39.ycen) +3.0*stddev(AORwasp39.ycen) and AORwasp39.ycen gt mean(AORwasp39.ycen) - 3.0*stddev(AORwasp39.ycen),ngood_pmap, complement=bad) 
 print, 'bad ',n_elements(bad)
 xarr = AORwasp39.xcen[good]
 yarr = AORwasp39.ycen[good]
 timearr = AORwasp39.timearr[good]
 flux = AORwasp39.flux[good]
 corrflux = AORwasp39.corrflux[good]
 corrfluxerr = AORwasp39.corrfluxerr[good]
 bmjdarr = AORwasp39.bmjdarr[good]
 bkgdarr = AORwasp39.bkgd[good]
     
     
; binning
 bin_level = 63L                ;*60
 nframes = bin_level
 nfits = long(n_elements(corrflux)) / nframes
 print, 'nframes, nfits', n_elements(corrflux), nframes, nfits
 bin_flux = dblarr(nfits)
 bin_corrflux= dblarr(nfits)
 bin_corrfluxerr= dblarr(nfits)
 bin_timearr = dblarr(nfits)
 bin_bmjdarr = dblarr(nfits)
 bin_bkgd = dblarr(nfits)

 for si = 0L, long(nfits) - 1L do begin
                                ;print, 'working on si', n_elements(AORwasp39[a].corrflux), si, si*nframes, si*nframes + (nframes - 1)
    idata = corrflux[si*nframes:si*nframes + (nframes - 1)]
    id = flux[si*nframes:si*nframes + (nframes - 1)]
    idataerr = corrfluxerr[si*nframes:si*nframes + (nframes - 1)]
    ib = bkgdarr[si*nframes:si*nframes + (nframes - 1)]
                                ; if si lt 2 then print, 'ib', ib
                                ; print, 'mean', mean(ib,/nan)
    bin_flux[si] = mean(id, /nan)
    bin_corrflux[si] = mean(idata,/nan)
    bin_corrfluxerr[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
    bin_timearr[si]= mean(timearr[si*nframes:si*nframes + (nframes - 1)])
    bin_bmjdarr[si]= mean(bmjdarr[si*nframes:si*nframes + (nframes - 1)])
    bin_bkgd[si] = mean(ib,/nan)

 endfor                         ;for each fits image
     
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
 
 print, 'nflux', n_elements(bin_flux)
;     print, 'ti', bin_timearr - bin_timearr(0)
 
 am = plot( (AORwasp39.timearr - AORwasp39[0].timearr(0))/60./60., AORwasp39.flux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[0], xtitle = 'Time(hrs)', ytitle = 'Flux', xrange = [0,8], yrange =[0.011, 0.015])
 am = plot( (AORwasp39.timearr - AORwasp39[0].timearr(0))/60./60., AORwasp39.corrflux-0.002,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black', /overplot)

ab = plot((bin_timearr - bin_timearr(0))/60./60., bin_flux, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[0], xtitle = 'Time(hrs)', ytitle = 'Flux', xrange = [0,8])

ab = plot((bin_timearr - bin_timearr(0))/60./60., bin_corrflux, '6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black',/overplot)

  

end
