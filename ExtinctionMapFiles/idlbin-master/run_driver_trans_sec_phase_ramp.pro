pro run_driver_trans_sec_phase_ramp, planetname, chname, apradius, pmapname, bin_level, error_bars = error_bars, unbin = unbin, shuffle=shuffle
  ;;run_driver_trans_sec_phase_ramp, 'WASP-14b', '2', 2.25, '150723',63L
  COMMON bin_block, aorname, planethash, bin_xcen, bin_ycen, bin_bkgd, bin_flux, bin_fluxerr,  bin_timearr, bin_phase, bin_ncorr,bin_np, bin_npcent, bin_xcenp, bin_ycenp, bin_bkgdp, bin_fluxp, bin_fluxerrp,  bin_corrfluxp,  bin_timearrp, bin_corrfluxerrp,  bin_phasep,  bin_ncorrp, bin_nparrp, bin_npcentarrp, bin_bmjdarr, bin_xfwhm, bin_yfwhm,  bin_corrflux_dp, bin_bmjdarrp

  ;;COMMON data, bjd_tot, flux_tot, nbr_ind, gw, err, time_tot, nonlin_c, err_tot
  COMMON data, bjd_tot, flux_tot, nbr_ind, gw, err_tot, time_tot, nonlin_c

  
  ;;---------------------------------------------------------------
  ;;make the appropriate save file as input to Nikole's code
  ;;---------------------------------------------------------------
  savename1 = strcompress('/Users/jkrick/irac_warm/pcrs_planets/' + planetname + '/'+ planetname +'_phot_ch'+chname+'_'+string(apradius)+'_' + pmapname + '_newtime.sav',/remove_all)
  ;;savename1 = strcompress('/Users/jkrick/irac_warm/pcrs_planets/' + planetname + '/'+ planetname +'_phot_ch2_2.25000_160126.sav',/remove_all)
  restore, savename1
  ;; aorname = ['r45428224', 'r45428480', 'r45428736','r45428992' ] ; start with stares only
  ;;first stare?:  'r45426688',

  aorname= ['r45840128','r45839616','r45839104','r45838592','r45847040','r45846784','r45846528','r45846272','r45846016','r45843200','r45842944','r45842688','r45842432','r45842176','r45841920','r45841664','r45841408','r45845504','r45840896','r45845760','r45845504','r45845248','r45844992','r45844736','r45844480','r45844224','r45843968','r45843712','r45843456','r45840640','r45840384','r45839872','r45839360','r45838848','r45838336','r48688384','r48688128','r48687872','r48687616','r48683776','r48683264','r48682752','r48682240','r48681472','r48681216','r48680704']

  colorarr = ['burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia'] ;'gray', 'gray','gray','gray','gray',

  T0 = 56042.687D0
  period = 2.2437651D0

  for a = 0, n_elements(aorname) - 1 do begin
     ;;check if I should be using pmap corr or not
     ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
     ;;if 20% of the values are correctable than go with the pmap corr 
     print,a, ' 0.2nflux, ncorr, ', 0.2*n_elements([planethash[aorname(a),'flux']]), corrcount
     if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0

     if keyword_set(unbin) then begin
        if a eq 0 then begin
           bjd_tot = planethash[aorname(a),'bmjdarr']
           flux_tot = [planethash[aorname(a),'corrflux_d']]
           err_tot = [planethash[aorname(a),'corrfluxerr']]
           phase_tot = [planethash[aorname(a),'phase']]
        endif else begin
           bjd_tot = [bjd_tot, planethash[aorname(a),'bmjdarr']]
           flux_tot = [flux_tot, planethash[aorname(a),'corrflux_d']]
           err_tot = [err_tot, planethash[aorname(a),'corrfluxerr']]
           phase_tot = [phase_tot, planethash[aorname(a),'phase']]
        endelse
        
     endif else begin
        ;;do some binning at least for the tests
        ;;junkpar = binning_function(a, bin_level, pmapcorr,chname)
        junkpar = binning_function(a, bin_level, pmapcorr,chname,/set_nbins, n_nbins = 15L)
        ;;print, 'bin_bmjdarr', bin_bmjdarr[0:10]
        ;;if a eq 0 then bin_corrflux_dp = bin_corrflux_dp + 0.0001
        print, 'pmapcorr', pmapcorr
        extraerr = 0;0.0001         ;  0.0001
        if pmapcorr gt 0 then begin
           if a eq 0 then begin
              bjd_tot =bin_bmjdarrp
              flux_tot = bin_corrflux_dp 
              err_tot = bin_corrfluxerrp + extraerr
              phase_tot = bin_phasep
             endif else begin
              if a le 37 and a ge 35 then begin
                 bin_corrfluxerrp = bin_corrfluxerrp + 0.0001
                 bin_corrflux_dp = bin_corrflux_dp - 0.00015
              endif
              bjd_tot = [bjd_tot, bin_bmjdarrp]
              flux_tot = [flux_tot, bin_corrflux_dp]
              
              err_tot = [err_tot, bin_corrfluxerrp + extraerr]
              phase_tot = [phase_tot, bin_phasep]
           endelse
        endif  ;;pmapcorr gt 0
     endelse  ;;binning

     bjd_plot = (bin_phasep*period) +T0 

     p1 = errorplot(bjd_plot, bin_corrflux_dp/0.057606741, bin_corrfluxerrp/0.057606741, $
                    xtitle = 'BJD', ytitle = 'Relative Flux',xrange = [5.60412E4, 5.60440E4],  $
                    overplot = p1, color = colorarr[a], errorbar_color = colorarr[a],$
                   yrange = [0.986, 1.007], xmajor = 4);

  endfor  ;for all AORs

  if keyword_set(unbin) then begin
     ;;get rid of Nan's
     nonan = where(finite(flux_tot) gt 0 and finite(err_tot) gt 0)
     flux_tot = flux_tot(nonan)
     bjd_tot = bjd_tot(nonan)
     err_tot = err_tot(nonan)
     phase_tot = phase_tot(nonan)
     
  endif
  
  ;;phase all bjds into one orbit for fitting purposes
  ;;first sort
  sortco = sort(phase_tot)
  phase_tot = phase_tot(sortco)
  flux_tot = flux_tot(sortco)
  err_tot = err_tot(sortco)
  bjd_tot = bjd_tot(sortco)

  bjd_single = (phase_tot*period) +T0
  bjd_tot = bjd_single

  ;;mess with timing to see what happens to the fits
  at = where(bjd_tot gt 5.60436E4 and bjd_tot lt 5.60437E4, nat)
 ;; print, 'n_elements(at)', nat, mean(bjd_tot(at))
  ;;bjd_tot(at) = bjd_tot(at) + .04
 ;; print, 'after', mean(bjd_tot(at))

  aat = where(bjd_tot lt 5.60417E4, naat)
 ;; print, 'n_elements(aat)', naat, mean(bjd_tot(aat))
  ;;bjd_tot(aat) = bjd_tot(aat) - .05
;;  print, 'after', mean(bjd_tot(aat))

  
  ;;normalize
  se = where(phase_tot gt 0.47 and phase_tot lt 0.51, secount)
  ph = where(phase_tot lt 0.6 and phase_tot gt 0.55, phcount)
  plot_corrnorm =  mean(flux_tot(se),/nan)
  ;;plot_corrnorm =  0.057605;;0.05815;0.05761 ;
  print, 'plot_corrnorm', plot_corrnorm
  flux_tot = flux_tot / plot_corrnorm
  err_tot = err_tot/plot_corrnorm

  ;;plot
 ;; p1 = errorplot(bjd_tot, flux_tot, err_tot, $
 ;;                xtitle = 'BJD', ytitle = 'Relative Flux', $
 ;;                yrange = [0.99, 1.005])
  
  ;setup for running fitting code
  infile = strcompress('/Users/jkrick/irac_warm/pcrs_planets/' + planetname + '/'+ 'fitting_input_phot_ch'+$
                       chname+'_'+string(apradius)+'_' + pmapname + '.sav',/remove_all)
  outfile = strcompress('/Users/jkrick/irac_warm/pcrs_planets/' + planetname + '/'+ 'fitting_output_phot_ch'+$
                        chname+'_'+string(apradius)+'_' + pmapname + '.sav',/remove_all)

  nonlin_c=[0.523357, -0.74367, 0.801920, -0.316680]

   ;;I don't really need the rest so make them up for now
  nbr_ind = findgen(n_elements(bjd_tot))
  gw = 1.0
  time_tot = bjd_tot
     
  save, bjd_tot, flux_tot, nbr_ind, gw, err_tot, time_tot, nonlin_c, filename=infile

  ;;---------------------------------------------------------------
  ;;call Nikole's code
  ;;---------------------------------------------------------------
  
  ;;driver_trans_sec_phase_ramp_jk_nkl, 2, 'flat', 0, 0, outfile, infile
  driver_trans_sec_phase_ramp_jk_nkl, 2, 'cowan', 0, 0, outfile, infile

  restore, outfile
  finalp = p
  finalperror = perror


  pfit = plot(bjd_tot, trans,  overplot = p1, color = 'cyan', thick = 2)

  ;;now what about error bars on phase amplitude and phase shift
  if keyword_set(error_bars) then begin
     print, '------------------------'
     print, 'starting uncertainty testing'
     ;;finalp = finalp - finalperror ;;minus 1 sigma for all of the parameters
     ;;now try something more sophisticated where I vary some to all
     ;;of them, and make a distribution of amplitudes and phase shifts.
     ;;can have a random distribution of a random choice of each
     ;;parameter between it's 1 sigma values

     ;;pad the errors on some parameters - theses are just too small
     ;;                                    since I fix most of the
     ;;                                    variables.
     ;;using those from fitting wong et al. data where snapshots not available.
     ;;finalperror = [ 0.00000046, 0.,0.13,0.,0.,0.0001728, .00439, 6.78-5, 7.15E-5, 4.02E-5, 4.14E-5, 4.34E-5, 3.8E-5, 0.,0.0,0.]

     
     nruns = 1000
     amparr = fltarr(nruns)
     shiftarr = amparr
     for r = 0, nruns - 1 do begin
        rand = (2* randomu(seed, 16)) - 1 ; these will be the number of sigma for each of 16 parameters
        randomp = finalp + rand*finalperror
        driver_trans_sec_phase_ramp_fixed, 2, 'cowan', 0, 0, outfile, randomp,infile
        restore, outfile
        amparr(r) = ( ph_params(2) - ph_params(0))/2.
        shiftarr(r) = (p(5) - ph_params(1))*(1./p(0))*360.

     endfor
     ;;plot the last iteration
     ;;pfit = plot(bjd_tot, trans,  overplot = p1, color = 'red', thick = 2)

     ;;what are the distributions
     plothist, amparr, xhist, yhist, /noplot, bin = 1E-6
     ph1 = plot(xhist, yhist, xtitle = 'Amplitude pf phase curve', ytitle = 'number')
     ;;fit with a gaussian?
     start = [1.0E-3,1E-4, 100.]
     noise = fltarr(n_elements(yhist))
     noise[*] = 1                                              ;equally weight the values
     result= MPFITFUN('mygauss',xhist,yhist, noise, start) ;/quiet   ; fit a gaussian to the histogram 
     ph1 = plot(xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.), overplot = ph1)
     print, 'amplitude gauss results', result(0), result(1), result(2)
     print, 'amplitude - 2sigma', .000964 - 2*(result(1))
     print, 'shiftarr', shiftarr
     
     plothist, shiftarr, xhist, yhist, /noplot, bin = 2E-1
     ph2 = plot(xhist, yhist, xtitle = 'Phase Curve Shift(degrees)', ytitle = 'number')
          ;;fit with a gaussian?
     start = [13.,5., 100.]
     noise = fltarr(n_elements(yhist))
     noise[*] = 1                                              ;equally weight the values
     result= MPFITFUN('mygauss',xhist,yhist, noise, start) ;/quiet   ; fit a gaussian to the histogram 
     ph1 = plot(xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.), overplot = ph2)
     print, 'phase shift gauss results', result(0), result(1), result(2)

  endif

  if keyword_set(shuffle) then begin
     ;;want to shuffle around the fluxes in
     ;;each snapshot within 1 sigma, then
     ;;refit
     print, 'starting shuffle'
     nruns = 10
     amparr = fltarr(nruns)
     shiftarr = amparr
     colorarr = ['black','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']
     for r = 0, nruns - 1 do begin
        rand = (2* randomu(seed, n_elements(flux_tot))) - 1 ; these will be the number of sigma 
        randflux_tot = flux_tot + rand*err_tot
        save, bjd_tot, randflux_tot, nbr_ind, gw, err_tot, time_tot, nonlin_c, filename=infile
        testplot = plot(bjd_tot, randflux_tot, '1s', color = colorarr(r), overplot = testplot)
        driver_trans_sec_phase_ramp_jk_nkl, 2, 'cowan', 0, 0, outfile, infile
        restore, outfile
        testplot = plot(bjd_tot, trans,  overplot = testplot, color = 'cyan', thick = 2)

        amparr(r) = ( ph_params(2) - ph_params(0))/2.
        shiftarr(r) = (p(5) - ph_params(1))*(1./p(0))*360.
     endfor
     
    ;;what are the distributions
     plothist, amparr, xhist, yhist, /noplot, bin = 1E-6
     ph1 = plot(xhist, yhist, xtitle = 'Amplitude pf phase curve', ytitle = 'number')
     ;;fit with a gaussian?
     start = [8E-4,1E-4, 100.]
     noise = fltarr(n_elements(yhist))
     noise[*] = 1                                              ;equally weight the values
     result= MPFITFUN('mygauss',xhist,yhist, noise, start) ;/quiet   ; fit a gaussian to the histogram 
 ;    ph1 = plot(xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.), overplot = ph1)
     print, 'amplitude gauss results', result(0), result(1), result(2)

     print, 'shiftarr', shiftarr
     
     plothist, shiftarr, xhist, yhist, /noplot, bin = 2E-1
     ph2 = plot(xhist, yhist, xtitle = 'Phase Curve Shift(degrees)', ytitle = 'number')
          ;;fit with a gaussian?
     start = [13.,5., 100.]
     noise = fltarr(n_elements(yhist))
     noise[*] = 1                                              ;equally weight the values
     result= MPFITFUN('mygauss',xhist,yhist, noise, start) ;/quiet   ; fit a gaussian to the histogram 
;     ph1 = plot(xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.), overplot = ph2)
     print, 'phase shift gauss results', result(0), result(1), result(2)


  endif
  


  
  ;;overplot the Wong fit from his fitting routine
  normfactor = 0.9985
  readcol, '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/wong/ch2_data.dat', bjd, flux, format = '(D0, D0)'
  bjd = bjd + 56000 - 0.5

  readcol, '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/wong/ch2.dat', bjdmodel, fluxmodel, format = '(D0,D0)'
  pmodel = plot(bjd, fluxmodel/normfactor, overplot = p1, color = 'yellow', thick = 2)
     

  ;;overplot the fit from the wong et al. dataset
  outfile = '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/fitting_output.sav'
  restore, outfile
  pfit = plot(bjd_tot, trans,  overplot = p1, color = 'red', thick = 2)

end


