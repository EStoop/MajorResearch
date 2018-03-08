pro sample_stare_wasp14_nkl, n_samples, exodata = exodata
  t = systime(1)
  colorarr = ['black','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']
;;want code to read in the staring mode corrected photometry of
;;wasp-14b, then randomly sample it into 43 different 30 min AORs that are similar to
;;the snapshot observations. 

;;then want to fit the "snapshot style" light curves, and make
;;distributions of the resulting best-fit parameters

  ;;setup some known quantities for wasp14
  period = 2.24376543D          ;days
  utmjd_center = 56042.687d0; 2456034.21290D - 2400000.5D
  ;;ecc = 0.0828D
  ;;argperi = 252.11D             ;degrees
  ;;ar_semimaj = 5.99D
  exosystem = 'WASP-14 b'
  verbose = 0

  ;;read in the staring mode photometry with pmap reduction
;;  starefile = '/Users/jkrick/irac_warm/pcrs_planets/wasp-14b/wasp14_stare_150226.txt'
;;  readcol, starefile, starebmjd, starephase, corrflux, corrfluxerr, format = 'D, F, F, F',/nan`
;;     normfactor = 0.05775

  ;;or try reading in the wong et al reduction
  readcol, '/Users/jkrick/irac_warm/pcrs_planets/wasp-14b/wong/ch2_data.dat', bjd, flux, format = '(D0, D0)'
  normfactor = 0.9985
  corrfluxerr = fltarr(n_elements(bjd)) + 0.00388070540434
  starebmjd =  bjd + 56000. - 0.5
  corrflux = flux ;/ normfactor      ; re-normalize

  ;;and convert to phase
  bmjd_dist = starebmjd - utmjd_center ; how many UTC away from the transit center
  starephase =( bmjd_dist / period )- fix(bmjd_dist/period)
  print, 'min/max phase', min(starephase) + .005, max(starephase) - .005
  print, 'min/max bmjd', min(starebmjd), max(starebmjd)
  ;;setup arrays for free parameters
  amparr = fltarr(n_samples)
  shiftarr = amparr

  ;;pick 43 random start phases from -0.5 to 0.5
  fifteen = 0.005

  ;;run this sampling n_samples times
  for ji = 0, n_samples - 1 do begin
     phasecenter = (1.1)*(randomu(seed, 40)) - 0.55 ;1.1 ..55  ;;these will be different every time this is called
     phasecenter = [phasecenter, 0.5, -0.5, 0.0]  ; force transit and eclipse
     plothist, phasecenter, xhist, yhist, /noplot, bin = 0.1
     phtest = plot(xhist, yhist, xtitle = 'random phase choices', overplot = phtest)
     
     print, 'phasecenter', phasecenter
     ;;want to grab 30 min around these phasecenters from the actual data
     for i = 0, n_elements(phasecenter) - 1 do begin
        a = where(starephase lt phasecenter(i) + fifteen and starephase gt phasecenter(i) - fifteen)
        unbincorrflux = corrflux[a]
        unbincorrfluxerr = corrfluxerr[a]
        unbinphase = starephase[a]
        unbinbmjd = starebmjd[a]
;        print, 'unbinphase', unbinphase
;        print, 'unbinbmjd', unbinbmjd
       
        
        ;;remove NaNs
        good = where(finite(unbincorrflux) gt 0, ngood, COMPLEMENT=nbad)
        ;;print, 'number of nan', nbad
        unbincorrflux = unbincorrflux(good)
        unbinbmjd = unbinbmjd(good)
        unbincorrfluxerr = unbincorrfluxerr(good)


        ;;binning
        bin_level = 63L
        n_bins = 15L
        numberarr = findgen(n_elements(unbinbmjd))
;;        h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
        h = histogram(numberarr, OMIN=om, nbins = n_bins, reverse_indices = ri)
;;        print, 'n_elements', n_elements(h), n_elements(numberarr), n_elements(numberarr) / 3.
        bin_flux = dblarr(n_elements(h))
        bin_fluxerr = bin_flux
        bin_bmjd = bin_flux
        c = 0
        for j = 0L, n_elements(h) - 2 do begin
           
           ;;get rid of the bins with no values and low numbers, meaning
           ;;low overlap
;;           print, 'ri[j+1]',ri[j+1]
;;           print, 'ri[j]',ri[j]
           if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
              meanclip, unbincorrflux[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
              bin_flux[c] = meanflux ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
              
              idataerr = unbincorrfluxerr[ri[ri[j]:ri[j+1]-1]]
              bin_fluxerr[c] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
              
              meanbmjdarr = mean( unbinbmjd[ri[ri[j]:ri[j+1]-1]],/nan)
              bin_bmjd[c]= meanbmjdarr
              c = c + 1
           endif else begin
              print, 'not binning'
           endelse
           
        endfor
        bin_flux = bin_flux[0:c-1]
        bin_fluxerr = bin_fluxerr[0:c-1]
        bin_bmjd = bin_bmjd[0:c-1]

;;        print, 'i, bin_bmjd', i,  bin_bmjd
        
        if i eq 0 then begin
           simulcorrflux = bin_flux
           simulcorrfluxerr = bin_fluxerr
           simulbmjd = bin_bmjd
        endif else begin
           simulcorrflux = [simulcorrflux, bin_flux]
           simulcorrfluxerr = [simulcorrfluxerr, bin_fluxerr]
           simulbmjd = [simulbmjd, bin_bmjd]

        endelse
     endfor


     ;;sorting
     sa = sort(simulbmjd)
     flux_tot = simulcorrflux[sa]
     err_tot = simulcorrfluxerr[sa]
     bjd_tot = simulbmjd[sa]

     ;;some errant values too low in bmjd?
     ;;because eccentric, then phase doesn't quite correlate
     ;;with bjd in expected way
     ;;good = where(bjd_tot gt 5.60417E4); and bjd_tot lt 5.60438E4)
     ;;bjd_tot = bjd_tot(good)
     ;;err_tot = err_tot(good)
     ;;flux_tot = flux_tot(good)
 
     ;;need to normalize
     ;;hardcode since it is not agiven that I will get samples during secondary
     flux_tot = flux_tot / normfactor
     err_tot = err_tot / normfactor
     
     ;;test output by plotting
     ;;p1 = errorplot(bjd_tot, flux_tot, err_tot, '1s-', sym_size =0.2, sym_filled= 1, xtitle = 'BMJD', $
     ;;               ytitle = 'corrflux', color = colorarr[ji], errorbar_color = colorarr[ji], overplot = p1)
     
     infile = '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/continuous_subsampled.sav'
     outfile ='/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/fitting_output_continuous_subsampled.sav'
     nonlin_c=[0.523357, -0.74367, 0.801920, -0.316680]
     ;;I don't really need the rest so make them up for now
     nbr_ind = findgen(n_elements(bjd_tot))
     gw = 1.0
     time_tot = bjd_tot
     
     save, bjd_tot, flux_tot, nbr_ind, gw, err_tot, time_tot, nonlin_c, filename=infile
     driver_trans_sec_phase_ramp_jk_nkl, 2, 'cowan', 0, 0, outfile, infile
     restore, outfile
     print, 'ji', ji
     amparr(ji) = ( ph_params(2) - ph_params(0))/2.
     shiftarr(ji) = (p(5) - ph_params(1))*(1./p(0))*360.

     ;;overplot the results
     ;;pfit = plot(bjd_tot, trans,  overplot = p1, color = colorarr[ji], thick = 2)

     ;;wait, 5s
  endfor
  ;;overplot the fit from the wong et al. dataset
  ;;outfile = '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/fitting_output.sav'
  ;;restore, outfile
  ;;pfit = plot(bjd_tot, trans,  overplot = p1, color = 'red', thick = 2)
  
  
  ;;what are the distributions
  plothist, amparr, xhist, yhist, /noplot, bin = 1E-6
  ph1 = plot(xhist, yhist, xtitle = 'Amplitude pf phase curve', ytitle = 'number')
  ;;fit with a gaussian?
  start = [8E-4,1E-4, 100.]
  noise = fltarr(n_elements(yhist))
  noise[*] = 1                                             ;equally weight the values
  result= MPFITFUN('mygauss',xhist,yhist, noise, start)    ;/quiet   ; fit a gaussian to the histogram 
                                ;    ph1 = plot(xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.), overplot = ph1)
  print, 'amplitude gauss results', result(0), result(1), result(2)
  
  
  plothist, shiftarr, xhist, yhist, /noplot, bin = 2E-1
  ph2 = plot(xhist, yhist, xtitle = 'Phase Curve Shift(degrees)', ytitle = 'number')
  ;;fit with a gaussian?
  start = [13.,5., 100.]
  noise = fltarr(n_elements(yhist))
  noise[*] = 1                                             ;equally weight the values
  result= MPFITFUN('mygauss',xhist,yhist, noise, start)    ;/quiet   ; fit a gaussian to the histogram 
;     ph1 = plot(xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.), overplot = ph2)
  print, 'phase shift gauss results', result(0), result(1), result(2)
  
end
