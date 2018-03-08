pro plot_wasp33_ch2
  planetname = 'WASP-33b'
  chname = '2'
  apradius = 2.25
  bin_level = 10L

  COMMON bin_block, aorname, planethash, bin_xcen, bin_ycen, bin_bkgd, bin_flux, bin_fluxerr,  bin_timearr, bin_phase, bin_ncorr,bin_np, bin_npcent, bin_xcenp, bin_ycenp, bin_bkgdp, bin_fluxp, bin_fluxerrp,  bin_corrfluxp,  bin_timearrp, bin_corrfluxerrp,  bin_phasep,  bin_ncorrp, bin_nparrp, bin_npcentarrp, bin_bmjdarr

  colorarr = ['blue', 'red','black','green']

  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  stareaor = planetinfo[planetname, 'stareaor']
  print, 'stareaor', stareaor
  dirname = strcompress(basedir + planetname +'/');+'/hybrid_pmap_nn/')
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'superdark.sav',/remove_all)
  print, 'restoring ', savefilename
  restore, savefilename
  print, 'aorname', aorname(0)


;for debugging: skip some AORs
startaor = 0;5
stopaor =  n_elements(aorname) - 1

;setup some arrays to hold binned values
  bin_corrfluxfinal = fltarr(n_elements(aorname))
  bin_phasefinal = fltarr(n_elements(aorname))
  bin_corrfluxerrfinal = fltarr(n_elements(aorname))

  mean_corrfluxarr = fltarr(stopaor - startaor + 1,/nozero)
  mean_dsweetarr = mean_corrfluxarr
  stddev_corrfluxarr = mean_corrfluxarr
  stddev_dsweetarr = mean_corrfluxarr
 
  for a = startaor,stopaor do begin
     print, '------------------------------------------------------'
     print, 'working on AOR', a, '   ', aorname(a)

     ;check if I should be using pmap corr or not
     ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
     ;if 20% of the values are correctable than go with the pmap corr 
     print, '0.2nflux, ncorr, ', 0.2*n_elements([planethash[aorname(a),'flux']]), corrcount
     if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0
     print, 'pmapcorr', pmapcorr
     junkpar = binning_function(a, bin_level, pmapcorr)


     if a eq 0 then total_timearr = bin_timearr else total_timearr = [total_timearr, bin_timearr]
     if a eq 0 then total_phase = bin_phase else total_phase = [total_phase, bin_phase]
     if a eq 0 then total_flux = bin_flux else total_flux = [total_flux, bin_flux]
     if a eq 0 then total_corrfluxp = bin_corrfluxp else total_corrfluxp = [total_corrfluxp, bin_corrfluxp]
  endfor



;------------------------------------------------------
;now try plotting
  setxrange = [0,40]            ;, [-0.5, 0.5]
  corrnormoffset =  -0.0015 ; 0.02
  corroffset = 0                         ;0.001
  setynormfluxrange = [0.98, 1.01]       ;[0.97, 1.005]
  
  
  plot_norm = median(total_flux)
  plot_corrnorm = mean(total_corrfluxp,/nan)
  pr = plot((total_timearr - total_timearr(0))/60./60., total_flux/plot_norm, '1s', sym_size = 0.2,   sym_filled = 1,  $
            color = 'black',  xtitle = 'Time (Hrs.)', ytitle = 'Normalized Flux', $
            title = planetname, yrange = setynormfluxrange, xrange = setxrange,/nodata) 
  pr = plot((total_timearr - total_timearr(0))/60./60., (total_corrfluxp/plot_corrnorm) -corrnormoffset ,/overplot, '1s', $
            sym_size = 0.2,   sym_filled = 1, color = 'blue')
;    t = text(26, 0.97, 'Krick', color = 'blue',/data)
;    t = text(26, 0.99, 'Raw', color = 'black',/data)

  
;-----------------------------------------------------------------------------------------------
  name = 'Hardegree-Ullman'
  readcol, '/Users/jkrick/conferences/IRACaas2014/data_challenge/Hardegree-UllmanWASPdata/WASP-33_lightcurve.dat', t, flux, err, comment = '#'
  
;;  pr = errorplot(t, flux - 0.05, err, '1s',  sym_size = 0.2,   sym_filled = 1, color = 'cyan', errorbar_color = 'cyan', /overplot)
;  pr = plot(t, flux - 0.04, '1s',  sym_size = 0.2,   sym_filled = 1, color = 'chocolate', errorbar_color = 'chocolate', /overplot)
;  t = text(26, 0.95, 'Hardegree-Ullman', color = 'chocolate',/data)

;-----------------------------------------------------------------------------------------------
  name = 'Lewis'
  readcol, '/Users/jkrick/conferences/IRACaas2014/data_challenge/lewis/lewis_wasp33_corrected_data.dat', t, flux, err, comment = '#'
;  t in BJD -- need to convert to hours starting at 0
  t = t - t(0)
  t = t * 24.

  ;need to bin this down.  ; use same binning level as my plot
  numberarr = findgen(n_elements(t))
     
  h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
  print, 'omin', om, 'nh', n_elements(h)


;mean together the flux values in each phase bin
     bin_corrflux= dblarr(n_elements(h))
     bin_corrfluxerr= bin_corrflux
     bin_timearr = bin_corrflux
 
     c = 0
     for j = 0L, n_elements(h) - 1 do begin

        if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
        
           meanclip, t[ri[ri[j]:ri[j+1]-1]], meant, sigmat
           bin_timearr[c] = meant  
 
           meanclip, flux[ri[ri[j]:ri[j+1]-1]], meanf, sigmaf
           bin_corrflux[c] = meanf  
           
           icorrdataerr = err[ri[ri[j]:ri[j+1]-1]]
           bin_corrfluxerr[c] =   sqrt(total(icorrdataerr^2))/ (n_elements(icorrdataerr))

           c = c + 1

        endif
     endfor
     bin_corrflux = bin_corrflux[0:c-1]
     bin_corrfluxerr = bin_corrfluxerr[0:c-1]
     bin_timearr = bin_timearr[0:c-1]
     

;;  pr = errorplot(t, flux - 0.09, err, '1s',  sym_size = 0.2,   sym_filled = 1, color = 'dark red',errorbar_color = 'dark red',/overplot)
;;  pr = plot(t, flux - 0.09, '1s',  sym_size = 0.2,   sym_filled = 1, color = 'dark red',errorbar_color = 'dark red',/overplot)
;  pr = plot(bin_timearr, bin_corrflux- 0.07, '1s',  sym_size = 0.2,   sym_filled = 1, color = 'dark red',errorbar_color = 'dark red',/overplot)
;  t = text(26, 0.92, 'Lewis', color = 'dark red',/data)

  ;and now the model 
 readcol, '/Users/jkrick/conferences/IRACaas2014/data_challenge/lewis/lewis_wasp33_model_lc.dat', t, flux, format='(D0, F0, F0)'
  ;  t in BJD -- need to convert to hours starting at 0
  t = t - t(0)
  t = t * 24.
;  pr = plot(t, flux - 0.07, '1s',  sym_size = 0.2,   sym_filled = 1, color = 'black', errorbar_color = 'black', /overplot)
  
  ;what does the Lewis model look like on my pmap corrected light curve?
  pr = plot(t+0.6, flux , '1s',  sym_size = 0.2,   sym_filled = 1, color = 'black', errorbar_color = 'black', /overplot, linestyle = 2)

;-----------------------------------------------------------------------------------------------
  name = 'Stevenson'
  readcol, '/Users/jkrick/conferences/IRACaas2014/data_challenge/stevenson/Stevenson-wasp33-CorrectedLC.txt', t, flux, err, format='(D0, F0, F0)'
  ;  t in BJD -- need to convert to hours starting at 0
  t = t - t(0)
  t = t * 24.
  
;also needs binning
  ;need to bin this down.  ; use same binning level as my plot
  numberarr = findgen(n_elements(t))
     
  h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
  print, 'omin', om, 'nh', n_elements(h)


;mean together the flux values in each phase bin
     bin_corrflux= dblarr(n_elements(h))
     bin_corrfluxerr= bin_corrflux
     bin_timearr = bin_corrflux
 
     c = 0
     for j = 0L, n_elements(h) - 1 do begin

        if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
        
           meanclip, t[ri[ri[j]:ri[j+1]-1]], meant, sigmat
           bin_timearr[c] = meant  
 
           meanclip, flux[ri[ri[j]:ri[j+1]-1]], meanf, sigmaf
           bin_corrflux[c] = meanf  
           
           icorrdataerr = err[ri[ri[j]:ri[j+1]-1]]
           bin_corrfluxerr[c] =   sqrt(total(icorrdataerr^2))/ (n_elements(icorrdataerr))

           c = c + 1

        endif
     endfor
     bin_corrflux = bin_corrflux[0:c-1]
     bin_corrfluxerr = bin_corrfluxerr[0:c-1]
     bin_timearr = bin_timearr[0:c-1]
     

;;  print, 't', t[0:100]
;;  pr = errorplot(t, flux - 0.05, err, '1s',  sym_size = 0.2,   sym_filled = 1, color = 'cyan', errorbar_color = 'cyan', /overplot)
;;  pr = plot(t, flux - 0.1, '1s',  sym_size = 0.2,   sym_filled = 1, color = 'purple', errorbar_color = 'purple', /overplot)
;  pr = plot(bin_timearr, bin_corrflux- 0.1, '1s',  sym_size = 0.2,   sym_filled = 1, color = 'purple',errorbar_color = 'purple',/overplot)
;    t = text(26, 0.885, 'Stevenson', color = 'purple',/data)

  ;and now the model 
 readcol, '/Users/jkrick/conferences/IRACaas2014/data_challenge/stevenson/Stevenson-wasp33-ModelLC.txt', t, flux, format='(D0, F0, F0)'
  ;  t in BJD -- need to convert to hours starting at 0
  t = t - t(0)
  t = t * 24.
pr = plot(t, flux - 0.1, '1s',  sym_size = 0.2,   sym_filled = 1, color = 'black', errorbar_color = 'black', /overplot)

;what does the stevenson model look like on my light curve
pr = plot(t, flux,  '1s',  sym_size = 0.2,   sym_filled = 1, color = 'black', errorbar_color = 'black', /overplot)

;-----------------------------------------------------------------------------------------------
;  pr.save, '/Users/jkrick/conferences/IRACaas2014/data_challenge/wasp33_results.png'
end
