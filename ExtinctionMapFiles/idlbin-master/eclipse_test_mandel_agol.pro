pro eclipse_test_mandel_agol, planetname, bin_level, apradius, chname,nobin = nobin
;planetname = 'XO3'

  COMMON bin_block, aorname, planethash, bin_xcen, bin_ycen, bin_bkgd, bin_flux, bin_fluxerr,  bin_timearr, bin_phase, bin_ncorr,bin_np, bin_npcent, bin_xcenp, bin_ycenp, bin_bkgdp, bin_fluxp, bin_fluxerrp,  bin_corrfluxp,  bin_timearrp, bin_corrfluxerrp,  bin_phasep,  bin_ncorrp, bin_nparrp, bin_npcentarrp, bin_bmjdarr, bin_xfwhm, bin_yfwhm,  bin_corrflux_dp, bin_bmjdarrp

 colorarr = ['burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']

aornum = 11

;run code to read in all the input planet parameters
  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  stareaor = planetinfo[planetname, 'stareaor']
  print, 'stareaor', stareaor
  exoplanet_data_file = '/Users/jkrick/idlbin/exoplanets.csv'
  exosystem = strmid(planetname, 0, 8 )+ ' b' ;'HD 209458 b' ;
  if planetname eq 'simul_XO3' then exosystem = 'XO-3 b'

  dirname = strcompress(basedir + planetname +'/')       
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_150723.sav',/remove_all) ;
  print, 'restoring ', savefilename
  restore, savefilename
  print, 'aorname', aorname(0)
  time_0 = (planethash[aorname(0),'timearr']) 
  time_0 = time_0(0)

  modelfilename = strcompress(dirname + planetname +'_model_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)

 ;pick the normalization
  phase0 = planethash[aorname(aornum),'phase']
  corrflux_d = planethash[aorname(aornum),'corrflux_d']
  se = where(phase0 gt 0.47 and phase0 lt 0.51)
  
  plot_corrnorm =  mean(corrflux_d(se),/nan)
  startaor = 1
  stopaor = 20
  for a = startaor,stopaor, 2 do begin
     print, '------------------------------------------------------'
     print, 'working on AOR', a, '   ', aorname(a), startaor
     ;unbinned stuff
     corrflux = planethash[aorname(a),'corrflux']
     corrfluxerr =  planethash[aorname(a),'corrfluxerr']
     bmjd =  planethash[aorname(a),'bmjdarr']

     ;check if I should be using pmap corr or not
     ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
     ;if 20% of the values are correctable than go with the pmap corr 
     ;;print, '0.2nflux, ncorr, ', 0.2*n_elements([planethash[aorname(a),'flux']]), corrcount
     if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0
     ;;print, 'pmapcorr', pmapcorr
     
     if keyword_set(nobin) then begin
        bin_corrfluxp = corrflux_d
        bin_phasep = phase0
        bin_corrfluxerrp = corrfluxerr
        print, n_elements(bin_corrfluxp), n_elements(bin_phasep), n_elements(bin_corrfluxerrp)
     endif else begin
        junkpar = binning_function(a, bin_level, pmapcorr,chname)
        ;;use time degraded corrfluxes
        bin_corrfluxp = bin_corrflux_dp
     endelse
     

     ;fake that this is a transit/ not eclipse
     bin_phasep = bin_phasep - 0.585; 0.585-0.685;

     ;pick the normalization
     se = where(bin_phasep gt -0.015 and bin_phasep lt 0.015)
     plot_corrnorm =  mean(bin_corrfluxp(se),/nan)
    ;; print, 'plot_corrnorm', plot_corrnorm, mean(bin_corrfluxp)

     pu = errorplot(bin_phasep, (bin_corrfluxp/plot_corrnorm), $
                             bin_corrfluxerrp/plot_corrnorm,  '1s', sym_size = 0.3, yrange = [0.97, 1.03], $
                              sym_filled = 1,color =colorarr[a] ,xtitle = 'Orbital Phase',$
                             errorbar_color =  colorarr[a], title = aorname(a), ytitle = 'Pmap Corrected Flux')

    
     test = function_fit_lightcurve(planetname, bin_phasep, bin_corrfluxp/plot_corrnorm,  bin_corrfluxerrp/plot_corrnorm, modelfilename)
     p2 = plot(findgen(20) - 5, fltarr(20) + 1.0 + 0.001875, overplot = pu , xrange = [-0.08, 0.08]) ; 0.00158
     plotname = strcompress(dirname + planetname +'_'+aorname(a)+'_'+string(apradius)+'_fit_150723.png',/remove_all) 
     p2.save, plotname

     ;write output to be used by TAP
     outfilename =  strcompress(dirname +'phot_ch'+chname+'_TAP_'+aorname(a)+'_bin10.ascii',/remove_all)
     openw, outlun, outfilename,/GET_LUN
     for te = 0, n_elements(bin_bmjdarrp) -1 do begin ;n_elements(bin_phasep) -1
        printf, outlun,  bin_bmjdarrp(te) , ' ', bin_corrfluxp(te)/plot_corrnorm,  bin_corrfluxerrp(te)/plot_corrnorm, format = '(F0, A,F0, A, F0)'
 ;;       if finite(corrflux(te)) then printf, outlun, bmjd(te), ' ', corrflux(te)/plot_corrnorm, ' ', corrfluxerr(te)/plot_corrnorm, format = '(F0, A,F0, A, F0)'
     endfor
     free_lun, outlun
     
  endfor



end
