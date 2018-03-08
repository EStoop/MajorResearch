pro plot_secondaryeclipse_ch2, planetname, utmjd_center, transit_duration, period, n_aors
;example call
;plot_secondaryeclipse_ch2, 'wasp15', 56086.74059D, 226.4, 2.519961, 2
  
  colorarr = ['blue', 'red','deep_pink', 'magenta', 'medium_purple', 'hot_pink', 'light_pink', 'rosy_brown', 'chocolate', 'saddle_brown', 'maroon', 'orange_red', 'dark_orange', 'peach_puff', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'light_cyan', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green']

  savefilename = strcompress( '/Users/jkrick/irac_warm/pcrs_planets/'+planetname+'/'+planetname+'_phot_ch2.sav',/remove_all)
  restore, savefilename
 
;-------------------------------------------------------
  transit_duration = transit_duration /60./24.  ; in days
  for a = 0, n_aors - 1 do  begin
       planetob[a].bmjdarr = ((planetob[a].bmjdarr - utmjd_center) / period) 
      print, 'doing something'
   endfor



  fits_read, '/Users/jkrick/idlbin/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_120827.fits', pmapdata, pmapheader
  c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', aspect_ratio = 1, xrange = [0,500], yrange = [0,500])
  for a = 0, n_aors - 1 do begin
     xcen500 = 500.* (planetob[a].xcen - 14.5)
     ycen500 = 500.* (planetob[a].ycen - 14.5)
     an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
 endfor

;-----
 
  print, 'min', min((planetob[a].timearr - planetob[0].timearr(0))/60./60.)

    for a = 0, n_aors - 1 do begin
     if a eq 0 then begin
        am = plot( (planetob[a].timearr - planetob[0].timearr(0))/60./60., planetob[a].xcen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'X pix')
     endif else begin
        am = plot( (planetob[a].timearr - planetob[0].timearr(0))/60./60., planetob[a].xcen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
     endelse

  endfor

;------

  for a = 0, n_aors - 1 do begin
     if a eq 0 then begin
        am = plot( (planetob[a].timearr - planetob[0].timearr(0))/60./60., planetob[a].ycen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Y pix')
     endif else begin
        am = plot( (planetob[a].timearr - planetob[0].timearr(0))/60./60., planetob[a].ycen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
     endelse

  endfor

;------
;  for a = 0, n_aors - 1 do begin
;     if a eq 0 then begin
;        am = plot( (planetob[a].timearr - planetob[0].timearr(0))/60./60., planetob[a].bkgd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Background', xrange = [0,10])
;     endif else begin
;        am = plot( (planetob[a].timearr - planetob[0].timearr(0))/60./60., planetob[a].bkgd,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;     endelse

;  endfor

;------
;  for a = 0, n_aors - 1 do begin
;     if a eq 0 then begin
;        am = errorplot( (planetob[a].timearr - planetob[0].timearr(0))/60./60., planetob[a].flux,planetob[a].fluxerr, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'Flux', xrange = [0,10], yrange = [0.03,0.045], errorbar_capsize = 0.)
;        am = errorplot( (planetob[a].timearr - planetob[0].timearr(0))/60./60., planetob[a].corrflux-0.003,planetob[a].corrfluxerr,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], /overplot, errorbar_capsize = 0.)
;     endif else begin
;        am = errorplot( (planetob[a].timearr - planetob[0].timearr(0))/60./60., planetob[a].flux,planetob[a].fluxerr,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot, errorbar_capsize = 0.)
; am = errorplot( (planetob[a].timearr - planetob[0].timearr(0))/60./60., planetob[a].corrflux-0.003,planetob[a].corrfluxerr,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot, errorbar_capsize = 0.)
;     endelse

;  endfor

;------
;a = where(finite(planetob[1].corrflux) gt 0 and finite(planetob[1].flux) gt 0)

; print, (abs(planetob[1].corrflux[a[0:100]] - planetob[1].flux[a[0:100]])/ planetob[1].flux[a[0:100]])*100
;plothist,(abs(planetob[1].corrflux[a[0:100]] - planetob[1].flux[a[0:100]])/ planetob[1].flux[a[0:100]])*100, xhist, yhist, bin = 0.1,/fill


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;try this diffferently, try putting the AORs togetether as one array; then bin

;timearr = [planetob[0].timearr, planetob[1].timearr]
;fluxarr = [planetob[0].flux, planetob[1].flux]
;fluxerrarr = [planetob[0].fluxerr, planetob[1].fluxerr]
;corrfluxarr = [planetob[0].corrflux, planetob[1].corrflux]
;corrfluxerrarr = [planetob[0].corrfluxerr, planetob[1].corrfluxerr]
;xarr = [planetob[0].xcen, planetob[1].xcen]
;yarr = [planetob[0].ycen, planetob[1].ycen]
;bkgd = [planetob[0].bkgd, planetob[1].bkgd]
;bmjd = [planetob[0].bmjdarr, planetob[1].bmjdarr]

;actually don't include the initial 30min AOR.
timearr = [planetob[1].timearr]
fluxarr = [ planetob[1].flux]
fluxerrarr = [ planetob[1].fluxerr]
corrfluxarr = [ planetob[1].corrflux]
corrfluxerrarr = [planetob[1].corrfluxerr]
xarr = [ planetob[1].xcen]
yarr = [planetob[1].ycen]
bkgd = [ planetob[1].bkgd]
bmjd = [ planetob[1].bmjdarr]

;print, 'n', n_elements(timearr), n_elements(planetob[0].timearr), n_elements(planetob[1].timearr)
;now for the binning

      ;remove outliers
good = where(xarr lt mean(xarr) + 2.5*stddev(xarr) and xarr gt mean(xarr) -2.5*stddev(xarr) and xarr lt mean(xarr) +3.0*stddev(yarr) and yarr gt mean(yarr) - 3.0*stddev(yarr) and finite(corrfluxarr) gt 0 ,ngood_pmap, complement=bad) 
     print, 'bad ',n_elements(bad), n_elements(good)
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
     numberarr = findgen(n_elements(xarr))

     bin_level = 64L            ; 63L;*60
     h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
     print, 'omin', om, 'nh', n_elements(h)

;mean together the flux values in each phase bin
     bin_flux = dblarr(n_elements(h))
     bin_fluxerr = dblarr(n_elements(h))
     bin_corrflux= dblarr(n_elements(h))
     bin_corrfluxerr= dblarr(n_elements(h))
     bin_timearr = dblarr(n_elements(h))
     bin_bmjdarr = dblarr(n_elements(h))
     bin_bkgd = dblarr(n_elements(h))
     bin_bkgderr = dblarr(n_elements(h))
     bin_xcen = dblarr(n_elements(h))
     bin_ycen = dblarr(n_elements(h))
     bin_xgcntrd= dblarr(n_elements(h))
     bin_ygcntrd= dblarr(n_elements(h))
     bin_xcntrd= dblarr(n_elements(h))
     c = 0
     for j = 0L, n_elements(h) - 1 do begin

;get rid of the bins with no values and low numbers, meaning low overlap
        if (ri[j+1] gt ri[j])  then begin
;           print, 'binning together', n_elements(numberarr[ri[ri[j]:ri[j+1]-1]])
        ;print, 'binning', numberarr[ri[ri[j]:ri[j+1]-1]]
        
           meanclip, xarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_xcen[c] = meanx   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;           print, 'finished x'

           meanclip, yarr[ri[ri[j]:ri[j+1]-1]], meany, sigmay
           bin_ycen[c] = meany   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;           print, 'finished y'

           meanclip, bkgdarr[ri[ri[j]:ri[j+1]-1]], meansky, sigmasky
           bin_bkgd[c] = meansky ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;           print, 'finished bkgd'

           meanclip, flux[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
           bin_flux[c] = meanflux ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;           print, 'finished flux'

           meanclip, corrflux[ri[ri[j]:ri[j+1]-1]], meancorrflux, sigmacorrflux
           bin_corrflux[c] = meancorrflux
;           print, 'finished corrfluxx'

           meanclip, timearr[ri[ri[j]:ri[j+1]-1]], meantimearr, sigmatimearr
           bin_timearr[c]=meantimearr
;           print, 'finished time'

           meanclip, bmjdarr[ri[ri[j]:ri[j+1]-1]], meanbmjdarr, sigmabmjdarr
           bin_bmjdarr[c]= meanbmjdarr
;           print, 'finished first bmjd'


           ;xxxx fix this to divide out N
           meanclip, corrfluxerr[ri[ri[j]:ri[j+1]-1]], meancorrfluxerr, sigmacorrfluxerr
           bin_corrfluxerr[c] = meancorrfluxerr
           meanclip, fluxerr[ri[ri[j]:ri[j+1]-1]], meanfluxerr, sigmafluxerr
           bin_fluxerr[c] = meanfluxerr

;           meanclip, bkgderrarr[ri[ri[j]:ri[j+1]-1]], meanbkgderrarr, sigmabkgderrarr
 ;          bin_bkgderr[c] = meanbkgderrarr
;           print, 'finished last meanclips'
           c = c + 1
        ;print, 'testing', j, phasearr[ri[ri[j]:ri[j+1]-1]]
        endif
     endfor
     
  bin_xcen = bin_xcen[0:c-1]
  bin_ycen = bin_ycen[0:c-1]
  bin_bkgd = bin_bkgd[0:c-1]
  bin_flux = bin_flux[0:c-1]
  bin_fluxerr = bin_fluxerr[0:c-1]
  bin_corrflux = bin_corrflux[0:c-1]
  bin_timearr = bin_timearr[0:c-1]
  bin_bmjdarr = bin_bmjdarr[0:c-1]
  bin_corrfluxerr = bin_corrfluxerr[0:c-1]
;  bin_bkgderr = bin_bkgderr[0:c-1]


 ;try getting rid of flux outliers.
  ;do some running mean with clipping
     start = 0
     print, 'nflux', n_elements(bin_corrflux)
     for ni = 50, n_elements(bin_corrflux) -1, 50 do begin
        meanclip,bin_corrflux[start:ni], m, s, subs = subs ;,/verbose
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

;     pl = errorplot((bin_timearr - bin_timearr(0))/60./60., bin_flux/median(bin_flux),bin_fluxerr, '6r1s', sym_size = 0.1,   sym_filled = 1, xrange = [0,9],  yrange = [0.99, 1.02], xtitle = 'Time (hrs)', ytitle = 'Normalized Flux', errorbar_capsize = 0.)

      pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_flux/median(bin_flux), '6r1s', sym_size = 0.1,   sym_filled = 1,  xrange = [0,9],yrange = [0.98, 1.01], xtitle = 'Time (hrs)', ytitle = 'Normalized Flux', title = planetname)
 
;    pl = errorplot((bin_timearr - bin_timearr(0))/60./60., bin_corrflux/median(bin_corrflux)-0.005,bin_corrfluxerr, '6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black',/overplot, errorbar_capsize = 0.)

    pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_corrflux/median(bin_corrflux)-0.01, '1s', sym_size = 0.1,   sym_filled = 1, color = 'black',/overplot)

;     pl = errorplot(bin_bmjdarr, bin_flux/median(bin_flux) ,bin_fluxerr, '6r1s', sym_size = 0.1,   sym_filled = 1,   xtitle = 'Phase', ytitle = 'Normalized Flux', errorbar_capsize = 0., xrange = [-1,1], yrange = [0.99, 1.01])
;     pl = errorplot(bin_bmjdarr, (bin_corrflux/median(bin_corrflux))-0.003, bin_corrfluxerr,'6r1s', sym_size = 0.1,   sym_filled = 1, color = 'black',/overplot, errorbar_capsize = 0.)

;     pl = plot((bin_timearr - bin_timearr(0))/60./60., bin_bkgd, '6r1s', sym_size = 0.1,   sym_filled = 1, xrange = [0,9], xtitle = 'Time (hrs)', ytitle = 'Bkgd')

end


