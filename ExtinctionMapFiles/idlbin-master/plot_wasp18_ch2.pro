pro plot_wasp18_ch2
  colorarr = ['blue', 'red','deep_pink', 'magenta', 'medium_purple', 'hot_pink', 'light_pink', 'rosy_brown', 'chocolate', 'saddle_brown', 'maroon', 'orange_red', 'dark_orange', 'peach_puff', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'light_cyan', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green']

;  fits_read, '/Users/jkrick/idlbin/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_120124.fits', pmapdata, pmapheader

  restore,'/Users/jkrick/irac_warm/wasp18/wasp18_phot_ch2.sav'
  aorname = ['r40269056', 'r40269312'] ;ch2


;-------------------------------------------------------
;taken from nsted planet properties calculator
  utmjd_center = double(56086.74059); 2456087.24059

  duration = 226.4             ;transit duration in minutes
  duration = duration /60./24.  ; in days
  period = 2.519961 ;days
;set phase = 0.0 in the middle of the transit
;  for a = 0, n_elements(aorname) - 1 do  begin
;      AORwasp18[a].bmjdarr = ((AORwasp18[a].bmjdarr - utmjd_center) / period) 
;   endfor

;print, 'phase', AORwasp18[0].bmjdarr[0:100]


;  fits_read, '/Users/jkrick/idlbin/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_120124.fits', pmapdata, pmapheader
;  c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', aspect_ratio = 1, xrange = [0,500], yrange = [0,500])
;  for a = 0, n_elements(aorname) - 1 do begin
;     xcen500 = 500.* (AORwasp18[a].xcen - 14.5)
;     ycen500 = 500.* (AORwasp18[a].ycen - 14.5)
;     an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
; endfor

;-----
 
    for a = 0, n_elements(aorname) - 1 do begin
     if a eq 0 then begin
        am = plot( (AORwasp18[a].timearr - AORwasp18[0].timearr(0))/60./60., AORwasp18[a].xcen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'X pix',  title = 'wasp18 ch2', xrange = [0,30], yrange = [14,15])
     endif else begin
        am = plot( (AORwasp18[a].timearr - AORwasp18[0].timearr(0))/60./60., AORwasp18[a].xcen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
     endelse

  endfor

;------

  for a = 0, n_elements(aorname) - 1 do begin
     if a eq 0 then begin
        am = plot( (AORwasp18[a].timearr - AORwasp18[0].timearr(0))/60./60., AORwasp18[a].ycen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Y pix',  title = 'wasp18 ch2',  xrange = [0,30], yrange = [14.5,15.5])
     endif else begin
        am = plot( (AORwasp18[a].timearr - AORwasp18[0].timearr(0))/60./60., AORwasp18[a].ycen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
     endelse

  endfor

;------
;  for a = 0, n_elements(aorname) - 1 do begin
;     if a eq 0 then begin
;        am = plot( (AORwasp18[a].timearr - AORwasp18[0].timearr(0))/60./60., AORwasp18[a].bkgd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Background', xrange = [0,10])
;     endif else begin
;        am = plot( (AORwasp18[a].timearr - AORwasp18[0].timearr(0))/60./60., AORwasp18[a].bkgd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;     endelse

;  endfor

;------
  for a = 0, n_elements(aorname) - 1 do begin
     if a eq 0 then begin
        am = plot( (AORwasp18[a].timearr - AORwasp18[0].timearr(0))/60./60., AORwasp18[a].flux, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Flux',  title = 'wasp18 ch2',  xrange = [0,30], yrange = [0.09, 0.10])
;        am = plot( (AORwasp18[a].timearr - AORwasp18[0].timearr(0))/60./60., AORwasp18[a].corrflux-0.003,,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], /overplot)
     endif else begin
        am = plot( (AORwasp18[a].timearr - AORwasp18[0].timearr(0))/60./60., AORwasp18[a].flux,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], /overplot)
; am = plot( (AORwasp18[a].timearr - AORwasp18[0].timearr(0))/60./60., AORwasp18[a].corrflux-0.003,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
     endelse

  endfor

;------
;a = where(finite(AORwasp18[1].corrflux) gt 0 and finite(AORwasp18[1].flux) gt 0)

; print, (abs(AORwasp18[1].corrflux[a[0:100]] - AORwasp18[1].flux[a[0:100]])/ AORwasp18[1].flux[a[0:100]])*100
;plothist,(abs(AORwasp18[1].corrflux[a[0:100]] - AORwasp18[1].flux[a[0:100]])/ AORwasp18[1].flux[a[0:100]])*100, xhist, yhist, bin = 0.1,/fill


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;try this diffferently, try putting the AORs togetether as one array; then bin

timearr = [AORwasp18[0].timearr, AORwasp18[1].timearr]
fluxarr = [AORwasp18[0].flux, AORwasp18[1].flux]
fluxerrarr = [AORwasp18[0].fluxerr, AORwasp18[1].fluxerr]
corrfluxarr = [AORwasp18[0].corrflux, AORwasp18[1].corrflux]
corrfluxerrarr = [AORwasp18[0].corrfluxerr, AORwasp18[1].corrfluxerr]
xarr = [AORwasp18[0].xcen, AORwasp18[1].xcen]
yarr = [AORwasp18[0].ycen, AORwasp18[1].ycen]
bkgd = [AORwasp18[0].bkgd, AORwasp18[1].bkgd]
bmjd = [AORwasp18[0].bmjdarr, AORwasp18[1].bmjdarr]

print, 'n', n_elements(timearr), n_elements(AORwasp18[0].timearr), n_elements(AORwasp18[1].timearr)
;now for the binning

      ;remove outliers
good = where(xarr lt mean(xarr) + 2.5*stddev(xarr) and xarr gt mean(xarr) -2.5*stddev(xarr) and xarr lt mean(xarr) +3.0*stddev(yarr) and yarr gt mean(yarr) - 3.0*stddev(yarr),ngood_pmap, complement=bad) 
     print, 'bad ',n_elements(bad)
     xarr = xarr[good]
     yarr = yarr[good]
     timearr = timearr[good]
     flux = fluxarr[good]
     fluxerr = fluxerrarr[good]
     corrflux = corrfluxarr[good]
     corrfluxerr = corrfluxerrarr[good]
     bmjdarr = bmjd[good]
     bkgdarr = bkgd[good]
     
     
; binning
     bin_level = 63*1L;*60
     nframes = bin_level
     nfits = long(n_elements(corrflux)) / nframes
     print, 'nframes, nfits', n_elements(corrflux), nframes, nfits
     bin_flux = dblarr(nfits)
     bin_fluxerr = dblarr(nfits)
     bin_corrflux= dblarr(nfits)
     bin_corrfluxerr= dblarr(nfits)
     bin_timearr = dblarr(nfits)
     bin_bmjdarr = dblarr(nfits)
     bin_bkgd = dblarr(nfits)

     for si = 0L, long(nfits) - 1L do begin
                                ;print, 'working on si', n_elements(AORwasp18[a].corrflux), si, si*nframes, si*nframes + (nframes - 1)
        idata = corrflux[si*nframes:si*nframes + (nframes - 1)]
        id = flux[si*nframes:si*nframes + (nframes - 1)]
        idataerr = corrfluxerr[si*nframes:si*nframes + (nframes - 1)]
        ib = bkgdarr[si*nframes:si*nframes + (nframes - 1)]
       ; if si lt 2 then print, 'ib', ib
       ; print, 'mean', mean(ib,/nan)
        bin_flux[si] = mean(id, /nan)
        bin_fluxerr[si] = sqrt(total(id^2))/ (n_elements(id))
        bin_corrflux[si] = mean(idata,/nan)
        bin_corrfluxerr[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
        bin_timearr[si]= mean(timearr[si*nframes:si*nframes + (nframes - 1)])
        bin_bmjdarr[si]= mean(bmjdarr[si*nframes:si*nframes + (nframes - 1)])
        bin_bkgd[si] = mean(ib,/nan)

     endfor                     ;for each fits image
     
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
     bin_fluxerr = bin_fluxerr[good_ni]
     bin_corrflux = bin_corrflux[good_ni]
     bin_bkgd = bin_bkgd[good_ni]

     print, 'nflux', n_elements(bin_flux)
     
;now try plotting

     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_flux/median(bin_flux), '61s', sym_size = 0.1,   sym_filled = 1, xrange = [0,30],  yrange = [0.99, 1.02], xtitle = 'Time (hrs)', ytitle = 'Normalized Binned Flux',  title = 'wasp18 ch2')
     ;pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_corrflux/median(bin_corrflux)-0.005, '6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black',/overplot)

;     pl = errorplot(bin_bmjdarr, bin_flux/median(bin_flux) ,bin_fluxerr, '6r1s', sym_size = 0.1,   sym_filled = 1,   xtitle = 'Phase', ytitle = 'Normalized Flux', errorbar_capsize = 0., xrange = [-1,1], yrange = [0.99, 1.01])
;     pl = errorplot(bin_bmjdarr, (bin_corrflux/median(bin_corrflux))-0.003, bin_corrfluxerr,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black',/overplot, errorbar_capsize = 0.)

     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_bkgd/median(bin_bkgd), '61s', sym_size = 0.1,   sym_filled = 1,  xtitle = 'Time (hrs)', ytitle = 'Normalized binned Bkgd',  title = 'wasp18 ch2', xrange = [0,30], yrange = [0.8, 1.2])

end


