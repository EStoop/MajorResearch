pro plot_calstar_stability, chname, binperstar = binperstar, fit_slope = fit_slope, montecarlo = montecarlo
  
  fits_read, '/Users/jkrick/irac_warm/calstars/arrayloccorr/ch1_photcorr_ap_5.fits', photcor_ch1, photcorhead_ch1
  fits_read, '/Users/jkrick/irac_warm/calstars/arrayloccorr/ch2_photcorr_ap_5.fits', photcor_ch2, photcorhead_ch2
  ;potential calstar names
;     names = [ '1812095', 'KF08T3_', 'KF06T2_', 'KF06T1_',  'NPM1p67','NPM1p60']
     names = [ '1812095', 'KF08T3_', 'KF06T2_', 'KF06T1_', 'KF09T1_',  'NPM1p67','NPM1p60'];, 'NPM1p68'] ;, 'HD16545']
     normvals = fltarr(n_elements(names))
  for chan = 1, 2 do begin
     chname = strtrim(chan,1)
  for m = 0, 1 do begin         ; for cryo and warm
     
     if m eq 1 then filename = strcompress('/Users/jkrick/irac_warm/calstars/allch'+chname +'phot.sav',/remove_all)
     if m eq 0 then filename = strcompress('/Users/jkrick/irac_warm/calstars/cryo/allch'+chname +'phot.sav',/remove_all)
     restore, filename
     
;binning_ch1 = [5L,5L, 5L,5L,5L,5L,10L, 10L]
;colors = ['black', 'grey', 'red', 'blue', 'chocolate', 'cyan', 'orange', 'maroon', 'medium purple']
;start with the first four, the rest were observed differently
;keyname = ['f12s', 'F12sc', 'F12sc', 'F12s', 'F2s', 'F2sc', 'F0p4s', 'F0p4sc', 'FS0p4sc']
     
;add a correction for pixel phase effect 
     corrflux = pixel_phase_correct_gauss(fluxarr,xcenarr,ycenarr,1, '3_3_7')
     
;for array dependant photometric correction warm
     if chname eq '1' then photcorr = photcor_ch1(xcenarr, ycenarr) else photcorr = photcor_ch2(xcenarr, ycenarr)
     corrflux= corrflux * photcorr
     
     
     
     for n = 0, n_elements(names) - 1 do begin
        
        a = where(starnamearr eq names(n), count)
        print, names(n), count
        timeb = timearr(a)
        fluxb = fluxarr(a)
        fluxerrb = fluxerrarr(a)
        starnameb = starnamearr(a)
        corrfluxb  = corrflux(a)

;   p1 = errorplot((timearr(a) - min(timearr))/60./60./ 24./365, fluxarr(a)/median(fluxarr(a)), fluxerrarr(a)/median(fluxarr(a)),$
;                  '1s', sym_size = 0.2,   sym_filled = 1, color = colors(n),errorbar_color = colors(n), yrange = [0.9, 1.1],$
;                   xtitle = 'Time(years)', ytitle = 'Normalized Flux', title = 'Ch 1 Primary Calibrators', overplot = p1)
        
;----------------------------------------------------
;need some binning on time, since there are either 5 or 6 observations per visit.
;XXXch2 has 6 frames, ch1 has 5 frames.  Would like to ignore the first frame on ch2 anyway.
;aka only works for ch1 right now
        if keyword_set(binperstar) then zero = binperstar(timeb, fluxb, corrfluxb, fluxerrb, starnameb, binning_ch1, n, colors, pb)
   
        
;----------------------------------------------------
;fill in multiple star array
;a bit messy because I don't know how big each array will be up front
;need to keep this going because I need to normalize each star before
;I bin together.

        if m eq 0 then normvals(n) = median(corrfluxb[n_elements(corrfluxb)/3.:*])

        if n eq 0 then begin
;        all_corrflux = corrfluxb/mean(corrfluxb[n_elements(corrfluxb)/3.:*],/nan)
;           all_corrflux = corrfluxb/median(corrfluxb[n_elements(corrfluxb)/3.:*])
           all_corrflux = corrfluxb/normvals(n)
           all_time = timeb
           all_starname = starnameb
           all_fluxerr = fluxerrb
        endif else begin
;        all_corrflux = [all_corrflux,corrfluxb/mean(corrfluxb[n_elements(corrfluxb)/3.:*],/nan)]
;           all_corrflux = [all_corrflux,corrfluxb/median(corrfluxb[n_elements(corrfluxb)/3.:*])]
           all_corrflux = [all_corrflux,corrfluxb/normvals(n)]
           all_time = [all_time, timeb]
           all_starname = [all_starname, starnameb]
           all_fluxerr = [all_fluxerr, fluxerrb]
        endelse
        
;     print, 'n, norm', n, normvals(n), median(corrfluxb[n_elements(corrfluxb)/3.:*])
        
     endfor                     ; for each star
     
;----------------------------------------------------
;bin all stars together as a function of 2 week PAO!
;turns out PAOs change length slightly so this isn't perfect at
;keeping equal numbers per bin
;print, 'n corrflux', n_elements(all_corrflux)
     if m eq 0 then sclk_time0 = min(all_time)
     time_years = (all_time - sclk_time0)/60./60./ 24./365. ; aka time in days
     s = sort(time_years)
     sort_time_years = time_years(s)
;print, 'sort tiem', sort_time_years[0:100]
     sort_corrflux = all_corrflux(s)
     sort_fluxerr = all_fluxerr(s)
     sort_starname = all_starname(s)
     
;for now need to get rid of the first few campaigns because of a
;different fluxconv.
;     clip = where(sort_time_years gt 0.9 ) ;and sort_time_years lt 5.0)
;     sort_time_years = sort_time_years(clip)
;     sort_corrflux = sort_corrflux(clip)
;     sort_fluxerr = sort_fluxerr(clip)
;     sort_starname = sort_starname(clip)

     h = histogram(sort_time_years, OMIN=om, binsize = 14.0/365., reverse_indices = ri)
     print,  n_elements(h)
     bin_all_corrflux = findgen(n_elements(h))
     bin_all_fluxerr = bin_all_corrflux
     bin_all_time = bin_all_corrflux
     c = 0
     for j = 0L, n_elements(h) - 1 do begin
        if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
           
           meanclip, sort_corrflux[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_all_corrflux[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           idataerr = sort_fluxerr[ri[ri[j]:ri[j+1]-1]]
           bin_all_fluxerr[c] =   sigmax/sqrt(n_elements(idataerr)) ; sqrt(total(idataerr^2))/ (n_elements(idataerr))
           
           meanclip, sort_time_years[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_all_time[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           c = c + 1
        endif
     endfor
     
     bin_all_corrflux = bin_all_corrflux[0:c-1]
     bin_all_fluxerr = bin_all_fluxerr[0:c-1]
     bin_all_time = bin_all_time[0:c-1]
     
 ;print out some output
     openw, outlun, strcompress('/Users/jkrick/irac_warm/calstars/stability_ch'+chname+ '.txt',/remove_all), /get_lun,/append
     for count = 0, n_elements(bin_all_time) - 1 do begin
        printf, outlun, bin_all_time(count), bin_all_corrflux(count) / median(bin_all_corrflux), bin_all_fluxerr(count) / median(bin_all_corrflux)
     endfor
     free_lun, outlun


     
;----------------------------------------------------
;and some binned plotting
                                ;get rid of negative outlier
;     good = where(bin_all_corrflux/median(bin_all_corrflux) gt 0.990)
;     bin_all_time = bin_all_time(good)
;     bin_all_corrflux = bin_all_corrflux(good)
;     bin_all_fluxerr = bin_all_fluxerr(good)
     print, 'chname', chan
     if chname eq '1' then begin
        print, 'plotting chname 1'
        help, chname
        pch1 = errorplot(bin_all_time, bin_all_corrflux, bin_all_fluxerr,$
                        '1s', sym_size = 0.2,   sym_filled = 1, color = 'black', yrange = [0.99, 1.02],$
                        xtitle = 'Time (Years)',$; ytitle = 'Normalized Flux', $ ;title = 'Ch '+ chname+ ' Primary Calibrators', $
                         overplot = pch1, dimensions = [700,400],$ ;MARGIN = [0.15, 0.15, 0.20, 0.15],
                        position = [0.2,0.5,0.9,0.8] ,   xshowtext = 0   )                       ;  axis_style = 1 
     endif else begin
        print, 'plotting chname 2'
        pch2 = errorplot(bin_all_time, bin_all_corrflux, bin_all_fluxerr, yrange = [0.99, 1.02],$
                        xtitle = 'Time Since Launch (Years)', $;ytitle = 'Normalized Flux',$
                        '1s', sym_size = 0.2,   sym_filled = 1, color = 'black',  overplot = pch2,$
                        position = [0.2, 0.2, 0.9, 0.5],/current) ;  axis_style = 1 
     endelse
     
     pt1 = text(0.3, 0.7, 'Ch1')
     pt2 = text(0.3, 0.4, 'Ch2')
     pty = text(0.1, 0.35, 'Normalized Flux', baseline = [0, 1.0, 0], updir = [-1.0,0,0])

;----------------------------------------------------
;and now what about fitting the best fit slope?
;fit with a liner function.
     if keyword_set(fit_slope) then begin
        
        start = [-1E-3,median(bin_all_corrflux)]
        startflat = [1.0]
        result= MPFITFUN('linear',bin_all_time, bin_all_corrflux/median(bin_all_corrflux), bin_all_fluxerr/median(bin_all_corrflux), $
                         start, perror = perror, bestnorm = bestnorm)    
        
        fp = plot(bin_all_time, result(0)*(bin_all_time) + result(1),  color = 'black', overplot = pb2, thick = 3, axis_style = 1)
        
        DOF     = N_ELEMENTS(bin_all_time) - N_ELEMENTS(result) ; deg of freedom
        PCERROR = PERROR * SQRT(BESTNORM / DOF)                 ; scaled uncertainties
        print, 'linear result', result, 'pcerror', pcerror, ' perror', perror, 'bestnorm', bestnorm, 'dof', DOF, sqrt(bestnorm/dof)
        t1 = text(0.7, 0.991, sigfig(result(0), 3, /scientific) + '    ' + sigfig(bestnorm/DOF, 5), /data, color = 'black', font_style = 'bold') ;pcerror(0)
        
        result2 = MPFITFUN('flatline',bin_all_time,bin_all_corrflux/median(bin_all_corrflux), bin_all_fluxerr/median(bin_all_corrflux), $
                           startflat, perror = perror, bestnorm = bestnorm)    
        fp = plot(bin_all_time, 0*bin_all_time + result2(0),  color = 'blue', overplot = pb2,thick = 3, linestyle = 2, axis_style = 1)
        tt1 = text(3.0, 0.991, '0    ' + sigfig(bestnorm/DOF, 5), /data, color = 'blue', font_style = 'bold')
     endif
     
;--------------
;add CR to this plot from Joe's house keeping statistics
     if m eq 0 then begin
        readcol, '/Users/jkrick/irac_warm/calstars/pixel_stats_nomops_2009.csv',month, day, year, hours, min, hot_ch1, noisy_ch1, dead_ch1, cr_ch1, hot_ch2, noisy_ch2, dead_ch2, cr_ch2, hot_ch3, noisy_ch3, dead_ch3, cr_ch3, hot_ch4, noisy_ch4, dead_ch4, cr_ch4, skipline = 18
                                ;arg this one only has time in a date
                                ;format, need to get that into sclk.
        JDCNV, year, month, day, hours, JD ; now in JD
        mjd = jd - 2400000.5
        crsclk = mjd_to_sclk(mjd)
        if chname eq '1' then cr = cr_ch1  else cr = cr_ch2
     endif else begin
        readcol, '/Users/jkrick/irac_warm/calstars/pixel_stats_warm_2014.csv',time, MJD, date, hot_ch1, noisy_ch1, dead_ch1, cr_ch1, hot_ch2, noisy_ch2, dead_ch2, cr_ch2, skipline = 6
;need to put time on same scale as the original plot
        crsclk =  [crsclk, mjd_to_sclk(MJD)]
     if chname eq '1' then cr = [cr, cr_ch1] else cr = [cr, cr_ch2]

     endelse

     crsclk_days = (crsclk -sclk_time0)/60./60./ 24./365. ; in years since time 0 of the calstar dataset
;     if m eq 0 then for cs = 0, n_elements(crsclk) - 1 do print, crsclk_days(cs), month(cs), year(cs)

     if m eq 1 then begin
;        pcr = plot(crsclk_days, cr , '1*', color = 'green',  axis_style = 0, /current, /data, $
;                   MARGIN = [0.15, 0.15, 0.20, 0.15], yrange = [-5,8], xrange = pb2.xrange)
        
;add other axes
;        a_psx = axis('y', target = pcr, LOCATION = [max(pb2.xrange),0,0], textpos = 1, title = 'CR', color = 'green')
;        xaxis = axis('x', target = pb2, LOCATION = [0, max(pb2.yrange)], showtext = 0, tickdir = 1)
     endif

;--------------
;add predicted solar cycle.
;  readcol, '/Users/jkrick/irac_warm/calstars/solarpredict.txt',year, mon, ninetyfive, fifty, five,  skipline = 3
     
;--------------
     
; make another plot to look at the correlation between CR and
; normalized fluxes.
;tricky part here is that the times are not binned on the same
;scale. not so trivial to do this.
;faking it for now
;print, n_elements(cr), n_elements(bin_all_corrflux)
;;want to make them have the same number by chopping the beginning of cr_ch1
;chop_cr = cr[n_elements(cr) - n_elements(bin_all_corrflux): *]
;;and now normalized chop_cr_ch1
;chop_cr = chop_cr / median(chop_cr)
;print, n_elements(cr), n_elements(bin_all_corrflux), n_elements(chop_cr)
;pcrf = plot(bin_all_corrflux, chop_cr, '1s', sym_filled = 1, xtitle = 'Normalized Fluxes', $
;            ytitle = 'CR', xrange = [.99, 1.01], yrange = [0.5,1.5]);
;
;lin = findgen(4)
;pcrf = plot(lin, lin, '1-', overplot = pcrf)

;----------------------------------------------------
;Monte Carlo
;randomize the time variables, re-measure the slopes.  
; What is the distribution of slopes measured from S such randomizations?
  
     if keyword_set(montecarlo) then begin
        MC = 1000               ; how many implementations should I use?
        slopearr = findgen(MC)
        for m = 0, MC - 1 do begin
;use Fisher-Yates shuffle to randomize the time array
           sbin_all_time = fisher_yates_shuffle(bin_all_time)
           
;re-fit
           l_result= MPFITFUN('linear',sbin_all_time, bin_all_corrflux/median(bin_all_corrflux), $
                              bin_all_fluxerr/median(bin_all_corrflux), $
                              start, perror = perror, bestnorm = bestnorm,/quiet)    
           f_result= MPFITFUN('flatline',sbin_all_time,bin_all_corrflux/median(bin_all_corrflux), $
                              bin_all_fluxerr/median(bin_all_corrflux), $
                              startflat, perror = perror, bestnorm = bestnorm,/quiet) 
           
;keep an array of slopes
           slopearr[m] = l_result(0)
        endfor
        
;now plot a distribution of the slopes
        plothist, slopearr, xhist, yhist,  bin=0.0001, /noprint,/noplot
        titlename = string(MC) + ' Monte Carlo realizations'
        ph = barplot(xhist, yhist,  xtitle = 'Slope', ytitle = 'Number', fill_color = 'sky blue' , title = titlename)
        ph = plot(intarr(200) + result(0), indgen(200), linestyle = 2, thick = 2, overplot = ph)
        
        
;XXXtake this one step further, what fraction of the slopes are closer to the mean?
;fit a gaussian
        start = [0.000, 0.0002, MC]
        junkerr = fltarr(n_elements(xhist)) + 1.0
        g_result= MPFITFUN('mygauss',xhist, yhist,junkerr, start, perror = perror, bestnorm = bestnorm,/quiet)    
        ph = plot(xhist, g_result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - g_result(0))/g_result(1))^2.), color = 'black', overplot = ph)
        
        g = where(slopearr lt -1*(result(0)) and slopearr gt result(0), gcount) ; this works because the mean is zero slope
        print, 'fraction with slopes closer to zero',100*(gcount / float(MC)), gcount
     endif
     
  endfor                        ; for each mission
endfor    ; for each channel
end



;---------------------------------------------------------------------------------------------------
function binperstar, timeb, fluxb, corrfluxb, fluxerrb, starnameb, binning_ch1, n, colors, pb
  
  timeb = (timeb - min(timeb))/60./60./ 24./365. ; aka time in days
  s2 = sort(timeb)
  timeb = timeb(s2)
;print, 'sort tiem', sort_time_years[0:100]
  corrfluxb = corrfluxb(s2)
  fluxerrb = fluxerrb(s2)
  starnameb = starnameb(s2)
  
;for now need to get rid of the first few campaigns because of a
;different fluxconv.
  clip2 = where(timeb gt 0.9 );and timeb lt 5.0)
  timeb = timeb(clip2)
  corrfluxb = corrfluxb(clip2)
  fluxerrb = fluxerrb(clip2)
  starnameb = starnameb(clip2)

  x = indgen(n_elements(timeb))
  h = histogram(x, OMIN=om, binsize = binning_ch1(n), reverse_indices = ri)
  print, n, n_elements(h)
  bin_flux = dblarr(n_elements(h))
  bin_corrflux = bin_flux
  bin_fluxerr = bin_flux
  bin_time = bin_flux
  bin_starname = strarr(n_elements(h))
  c = 0
  for j = 0L, n_elements(h) - 1 do begin
     if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
        
        meanclip, fluxb[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_flux[c] = meanx     ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        
        meanclip, corrfluxb[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_corrflux[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        
        idataerr = fluxerrb[ri[ri[j]:ri[j+1]-1]]
        bin_fluxerr[c] = sigmax/sqrt(n_elements(idataerr)) ; sqrt(total(idataerr^2))/ (n_elements(idataerr))     ; just the standard deviation in the bins
                                ;
        
        meanclip, timeb[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_time[c] = meanx     ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        
        bin_starname[c] = starnameb[ri[ri[j]]] ; is true unless the observations changed!
        c = c + 1
     endif
  endfor
  
  bin_flux = bin_flux[0:c-1]
  bin_corrflux = bin_corrflux[0:c-1]
  bin_fluxerr = bin_fluxerr[0:c-1]
  bin_time = bin_time[0:c-1]
  bin_starname = bin_starname[0:c-1]
  
 

;----------------------------------------------------
;and some binned plotting
  btime = (bin_time); - min(bin_time))
  pb = errorplot(btime, bin_corrflux/median(bin_corrflux), bin_fluxerr/median(bin_flux),$
                 '1s', sym_size = 0.2,   sym_filled = 1, color = colors(n),errorbar_color = colors(n), yrange = [0.96, 1.04],$
                 xtitle = 'Time(years)', ytitle = 'Normalized Flux', title = 'Primary Calibrators', overplot = pb)
                                ;and now what about fitting the best fit slope?
;fit with a liner function.
  
     ;get rid of negative outlier
  goodb = where(bin_corrflux / median(bin_corrflux) gt 0.980)
  btime = btime(goodb)
  bin_corrflux = bin_corrflux(goodb)
  bin_fluxerr = bin_fluxerr(goodb)
  start = [1E-6,1.0]
  startflat = [1.0]
  result= MPFITFUN('linear',btime, bin_corrflux/median(bin_corrflux), bin_fluxerr/median(bin_corrflux), $
                   start, perror = perror, bestnorm = bestnorm)    
  
  fp = plot(btime, result(0)*(btime) + result(1),  color = colors(n), overplot = pb, thick = 3)
  
  DOF     = N_ELEMENTS(time) - N_ELEMENTS(result) ; deg of freedom
  PCERROR = PERROR * SQRT(BESTNORM / DOF)         ; scaled uncertainties
  print, 'linear result', result, 'pcerror', pcerror, ' perror', perror, 'bestnorm', bestnorm, 'dof', DOF, sqrt(bestnorm/dof)
  t1 = text(0.1, 0.981- (0.003*n), sigfig(result(0), 3, /scientific) + '    ' + sigfig(bestnorm/DOF, 5) +'   ' + bin_starname[0], /data, color = colors(n), font_style = 'bold') ;pcerror(0)
  
  result2 = MPFITFUN('flatline',btime,bin_corrflux/median(bin_corrflux), bin_fluxerr/median(bin_corrflux), startflat, perror = perror, bestnorm = bestnorm)    
;  fp = plot(btime, 0*btime + result2(0),  color = colors(n), overplot = pb,thick = 3, linestyle = 2)
;  tt1 = text(6E7, 0.981- (0.003*n), '0    ' + sigfig(bestnorm/DOF, 5), /data, color =colors(n), font_style = 'bold')
  
  
  return, 0

end
