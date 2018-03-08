pro run_driver_trans_sec_phase_ramp_test,error_bars = error_bars, wong = wong, pmap = pmap, snapshots=snapshots
  
  COMMON data, bjd_tot, flux_tot, nbr_ind, gw, err_tot, time_tot, nonlin_c

  dirname = '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/'
  normfactor = 0.9985
  colorarr = ['gray', 'gray','gray','gray','gray','burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']
  
  if keyword_set(pmap) then begin
     ;;read in the unbinned data
     restore, dirname + 'WASP-14b_phot_ch2_2.25000_150723_newtime.sav'
     aorname = ['r45426688','r45428224', 'r45428480', 'r45428736','r45428992' ] ; start with stares only
     print, n_elements(planethash[aorname(0),'bmjdarr'])
     bmjd_pmap = [planethash[aorname(1),'bmjdarr'],planethash[aorname(2),'bmjdarr'],planethash[aorname(3),'bmjdarr'],$
                  planethash[aorname(4),'bmjdarr']]
     flux_pmap = [planethash[aorname(1),'corrflux'],planethash[aorname(2),'corrflux'],planethash[aorname(3),'corrflux'],$
                  planethash[aorname(4),'corrflux']]
     err_pmap = [planethash[aorname(1),'corrfluxerr'],planethash[aorname(2),'corrfluxerr'],planethash[aorname(3),'corrfluxerr'],$
                 planethash[aorname(4),'corrfluxerr']]
     ;;normalize
     flux_pmap = flux_pmap/0.0577; 697;80051
     
     ;;is there some way to use the first AOR fluxes from Wong?
     readcol, dirname + '/wong/ch2_data.dat', bjd, flux, format = '(D0, D0)'
     err = fltarr(n_elements(bjd)) + 0.00388070540434
     bjd =  bjd + 56000. - 0.5
     flux = flux / normfactor  ; re-normalize

     bjd_all = [ bjd[0:n_elements(planethash[aorname(0),'bmjdarr'])], bmjd_pmap]
     flux_all = [ flux[0:n_elements(planethash[aorname(0),'bmjdarr'])], flux_pmap]
     fluxerr_all = [err[0:n_elements(planethash[aorname(0),'bmjdarr'])],err_pmap]
;     print, 'n', n_elements(bjd_all), n_elements(flux_all), n_elements(fluxerr_all)

     ;;testing without wong
     bjd_all = bmjd_pmap
     flux_all = flux_pmap
     fluxerr_all = err_pmap
     
     ;;remove bad data
     good = where(finite(bjd_all) gt 0 and finite(flux_all) gt 0 and finite(fluxerr_all) gt 0)
     bjd_all = bjd_all(good)
     flux_all = flux_all(good)
     fluxerr_all =fluxerr_all(good)
;     print, 'n after', n_elements(bjd_all), n_elements(flux_all), n_elements(fluxerr_all)

     ;;bin
     bin_level = 63L
     numberarr = findgen(n_elements(bjd_all))
     h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
     bin_flux = dblarr(n_elements(h))
     bin_fluxerr = bin_flux
     bin_bjd = bin_flux
     c = 0
     for j = 0L, n_elements(h) - 1 do begin
        
        ;;get rid of the bins with no values and low numbers, meaning
        ;;low overlap
        if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
           meanclip, flux_all[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
           bin_flux[c] = meanflux ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           idataerr = fluxerr_all[ri[ri[j]:ri[j+1]-1]]
           bin_fluxerr[c] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
           
           meanbmjdarr = mean( bjd_all[ri[ri[j]:ri[j+1]-1]],/nan)
           bin_bjd[c]= meanbmjdarr
           c = c + 1
        endif else begin
           print, 'not binning'
        endelse
        
     endfor
     bin_flux = bin_flux[0:c-1]
     bin_fluxerr = bin_fluxerr[0:c-1]
     bin_bjd = bin_bjd[0:c-1]

     ;;plot them to check
;     ptest = errorplot(bjd_all, flux_all, fluxerr_all, yrange = [0.97, 1.03])
;     ptest = errorplot(bin_bjd, bin_flux, bin_fluxerr, yrange = [0.99, 1.005])

         ;;use same naming as NKL
     bjd_tot = bin_bjd
     err_tot = bin_fluxerr
     flux_tot = bin_flux

     infile = strcompress('/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/fitting_input_continuous.sav',/remove_all)
     outfile = dirname + 'fitting_output_continuous.sav'

  endif
  

  if keyword_set(wong) then begin
  
     infile = strcompress('/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/fitting_input_phot_wong.sav',/remove_all)
     outfile = dirname + 'fitting_output.sav'
     readcol, dirname + '/wong/ch2_data.dat', bjd, flux, format = '(D0, D0)'
     
     ;;re-normalize
     flux = flux / normfactor
     
     bjd = bjd + 56000 - 0.5
     fluxerr = fltarr(n_elements(flux)) +0.00388070540434  ; from email from Wong
     ;;fluxerr = randomn(seed, n_elements(flux))*1E-3
     ;;bin
     bin_level = 63L
     numberarr = findgen(n_elements(bjd))
     h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
     bin_flux = dblarr(n_elements(h))
     bin_fluxerr = bin_flux
     bin_bjd = bin_flux
     c = 0
     for j = 0L, n_elements(h) - 1 do begin
        
        ;;get rid of the bins with no values and low numbers, meaning
        ;;low overlap
        if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
           meanclip, flux[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
           bin_flux[c] = meanflux ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           idataerr = fluxerr[ri[ri[j]:ri[j+1]-1]]
           bin_fluxerr[c] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
           
           meanbmjdarr = mean( bjd[ri[ri[j]:ri[j+1]-1]],/nan)
           bin_bjd[c]= meanbmjdarr
           c = c + 1
        endif else begin
           print, 'not binning'
        endelse
        
     endfor
     bin_flux = bin_flux[0:c-1]
     bin_fluxerr = bin_fluxerr[0:c-1]
     bin_bjd = bin_bjd[0:c-1]

     ;;use same naming as NKL
     bjd_tot = bin_bjd
     err_tot = bin_fluxerr
     flux_tot = bin_flux

  endif                        ; kwyword set = wong

  
  
  
  ;;plot the data
  p1 = errorplot(bjd_tot, flux_tot, err_tot, xtitle = 'bjd', ytitle = 'flux', yrange = [0.99, 1.005])
  
  
  nonlin_c=[0.523357, -0.74367, 0.801920, -0.316680]
  ;;I don't really need the rest so make them up for now
  nbr_ind = findgen(n_elements(bjd_tot))
  gw = 1.0
  time_tot = bjd_tot
  
  save, bjd_tot, flux_tot, nbr_ind, gw, err_tot, time_tot, nonlin_c, filename=infile
  
  ;;call Nikole's fitting code
;;  driver_trans_sec_phase_ramp_jk, 2, 'flat', 0, 0, outfile
  ;;call it a bunch of times and plot the distribution of answers
  nruns = 1
  amparr = fltarr(nruns)
  shiftarr = amparr
  for nc = 0, nruns - 1 do begin
     driver_trans_sec_phase_ramp_jk_nkl, 2, 'cowan', 0, 0, outfile, infile
     
     restore, outfile
     amparr(nc) = ( ph_params(2) - ph_params(0))/2.
     shiftarr(nc) = (p(5) - ph_params(1))*(1./p(0))*360.
  endfor
;;  print, 'amparr', amparr
;;  plothist, amparr, xhist, yhist, /noplot, bin = 1E-6
;;  ph1 = plot(xhist, yhist, xtitle = 'Amplitude pf phase curve', ytitle = 'number')
;;  plothist, shiftarr, xhist, yhist, /noplot, bin = 1E-1
;;  ph2 = plot(xhist, yhist, xtitle = 'Phase Curve Shift(degrees)', ytitle = 'number')

  
  finalp = p
  finalperror = perror
  
  ;;overplot the results
  pfit = plot(bjd_tot, trans,  overplot = p1, color = 'cyan', thick = 2)

  ;;overplot the Wong fit from his fitting routine
  readcol, dirname + '/wong/ch2.dat', bjdmodel, fluxmodel, format = '(D0,D0)'
  pmodel = plot(bjd, fluxmodel/normfactor, overplot = p1, color = 'green', thick = 2)


  ;;now what about error bars on phase amplitude and phase shift
  if keyword_set(error_bars) then begin
     print, '------------------------'
     print, 'starting fixed fitting'
     ;;finalp = finalp - finalperror ;;minus 1 sigma for all of the parameters
     ;;now try something more sophisticated where I vary some to all
     ;;of them, and make a distribution of amplitudes and phase shifts.
     ;;can have a random distribution of a random choice of each
     ;;parameter between it's 1 sigma values
     nruns = 1000
     amparr = fltarr(nruns)
     shiftarr = amparr
     for r = 0, nruns - 1 do begin
        rand = (2* randomu(seed, 16)) - 1 ; these will be the number of sigma for each of 16 parameters
        randomp = finalp + rand*finalperror
        driver_trans_sec_phase_ramp_fixed, 2, 'cowan', 0, 0, outfile, randomp, infile
        restore, outfile
        amparr(r) = ( ph_params(2) - ph_params(0))/2.
        shiftarr(r) = (p(5) - ph_params(1))*(1./p(0))*360.

     endfor
     ;;plot the last iteration
     pfit = plot(bjd_tot, trans,  overplot = p1, color = 'red', thick = 2)

     ;;what are the distributions
     plothist, amparr, xhist, yhist, /noplot, bin = 1E-6
     ph1 = plot(xhist, yhist, xtitle = 'Amplitude pf phase curve', ytitle = 'number')
     ;;fit with a gaussian?
     start = [7.85E-4,1E-5, 1000.]
     noise = fltarr(n_elements(yhist))
     noise[*] = 1                                              ;equally weight the values
     result= MPFITFUN('mygauss',xhist,yhist, noise, start) ;/quiet   ; fit a gaussian to the histogram 
     ph1 = plot(xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.), overplot = ph1)
     print, 'amplitude gauss results', result(0), result(1), result(2)
     
     plothist, shiftarr, xhist, yhist, /noplot, bin = 2E-1
     ph2 = plot(xhist, yhist, xtitle = 'Phase Curve Shift(degrees)', ytitle = 'number')
          ;;fit with a gaussian?
     start = [16.,1.0, 500.]
     noise = fltarr(n_elements(yhist))
     noise[*] = 1                                              ;equally weight the values
     result= MPFITFUN('mygauss',xhist,yhist, noise, start) ;/quiet   ; fit a gaussian to the histogram 
     ph1 = plot(xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.), overplot = ph2)
     print, 'phase shift gauss results', result(0), result(1), result(2)

  endif
  
end


