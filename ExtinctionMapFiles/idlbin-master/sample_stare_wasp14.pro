pro sample_stare_wasp14, n_samples, exodata = exodata
  t = systime(1)
;;want code to read in the staring mode corrected photometry of
;;wasp-14b, then randomly sample it into 43 different 30 min AORs that are similar to
;;the snapshot observations. 

;;then want to fit the "snapshot style" light curves, and make
;;distributions of the resulting best-fit parameters

  ;;setup some known quantities for wasp14
  porbit = 2.24376543D          ;days
  mjd_transit = 2456034.21290D - 2400000.5D
  ecc = 0.0828D
  argperi = 252.11D             ;degrees
  ar_semimaj = 5.99D
  exosystem = 'WASP-14 b'
  verbose = 0
  trad = randomu(seed, /double)*2.0D + 1.0D    ;; randomize initial guess
  omrot = randomu(seed,/Double)*0.15D + 0.85D   ;; randomize initial guess
  fstar = 0.0577d ;; initial guess
  albedo = !NULL  ;; set to null to fix the value at zero in the fit
  IF N_ELEMENTS(RP_RSTAR) EQ 0 AND ~KEYWORD_SET(FIT_SYSTEM) THEN rp_rstar = 0.09420d
  IF N_ELEMENTS(INCLINATION) EQ 0 AND ~KEYWORD_SET(FIT_SYSTEM) THEN inclination=84.63

  ;;read in the staring mode photometry
  starefile = '/Users/jkrick/irac_warm/pcrs_planets/wasp-14b/wasp14_stare_150226.txt'
  readcol, starefile, starebmjd, starephase, corrflux, corrfluxerr, format = 'D, F, F, F',/nan

  ;;setup arrays for free parameters
  tradarr = fltarr(n_samples)
  omegaarr = tradarr
  chisquaredarr = tradarr
  fstararr = tradarr
  amparr = tradarr
  phase_maxarr = tradarr
  phase_minarr = tradarr
 ;;pick 43 random start phases from -0.5 to 0.5
 ;;actually don't want to be too close to the edge, 30 min AOR,
 ;;so want random phases to run from -0.495 to 0.495
  edge = 0.45; 0.495
  fifteen = 0.005

  ;;run this sampling n_samples times
  for j = 0, n_samples - 1 do begin
   
     ;phasecenter = (edge*2)*(randomu(seed, 43)) - edge  ;;these will be different every time this is called
     phasecenter = (0.82)*(randomu(seed, 43)) - 0.37  ;;these will be different every time this is called

     ;;want to grab 30 min around these phasecenters from the actual data
     for i = 0, n_elements(phasecenter) - 1 do begin
        a = where(starephase lt phasecenter(i) + fifteen and starephase gt phasecenter(i) - fifteen)
        if i eq 0 then begin
           simulcorrflux = corrflux[a]
           simulcorrfluxerr = corrfluxerr[a]
           simulphase = starephase[a]
           simulbmjd = starebmjd[a]
        endif else begin
           simulcorrflux = [simulcorrflux, corrflux[a]]
           simulcorrfluxerr = [simulcorrfluxerr, corrfluxerr[a]]
           simulphase = [simulphase, starephase[a]]
           simulbmjd = [simulbmjd, starebmjd[a]]

        endelse
     endfor
     
     ;;test output by plotting
     ;;p1 = plot(simulphase, simulcorrflux, '1s', sym_size =0.2,
     ;;sym_filled= 1, xtitle = 'phase', ytitle = 'corrflux')
     
     ;;now do some fitting

     ;;try iteratively
     for iter = 0, 4 do begin
        FIT_PHASE_CURVE,simulbmjd,simulcorrflux,exosystem,fstar,albedo,trad,omrot,tperi,tphase,model,$
                        EXODATA=exodata,LAMBDA_BAND='IRAC4.5',VERBOSE=verbose,PARAM_ERR=param_err,$
                        REDUCED_CHI2=reduced_chi2,ECC=ecc,ARGPERI=argperi,$
                        PORBIT=porbit,MJD_TRANSIT=mjd_transit,RP_RSTAR=rp_rstar,AR_SEMIMAJ=ar_semimaj,$
                        INCLINATION=inclination,TPERI_2=tperi_2,MODEL_2=model_2,FUNC=simulcorrfluxerr,$
                        FIT_SYSTEM=fit_system,TPHASE_2=tphase_2,AMPLITUDE=amplitude,PHASE_MAX=phase_max,$
                        PHASE_MIN=phase_min
     endfor

     ;;want to keep the fitted parameters in an array
     tradarr[j] = trad
     omegaarr[j] = omrot
     fstararr[j] = fstar
     chisquaredarr[j] = reduced_chi2
     amparr[j] = amplitude
     phase_maxarr[j] = phase_max
     phase_minarr[j] = phase_min
     if j mod 100 eq 0 then print, 'testing final values', j, trad, omrot, reduced_chi2, systime(1) - t
  endfor
  savename = strcompress('/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/fit_phase_summary_'+string(n_samples)+'.sav',/remove_all)
  save, tradarr, omegaarr, fstararr, chisquaredarr,amparr, phase_maxarr, phase_minarr, filename=savename

  ;;plot histograms of the fitted parameters
  plothist, tradarr, xhist, yhist, bin = 0.1, /noprint, /noplot
  p2 = barplot(xhist, yhist, xtitle = 'Radiative Timescale', ytitle = 'Number', fill_color = 'blue', $
               title = string(n_samples) + ' realizations')

  plothist, omegaarr, xhist, yhist, bin = 0.01, /noprint, /noplot
  p2 = barplot(xhist, yhist, xtitle = 'Omega rotation', ytitle = 'Number', fill_color = 'blue', $
               title = string(n_samples) + ' realizations')
  
  plothist, chisquaredarr, xhist, yhist, bin = 0.01, /noprint, /noplot
  p2 = barplot(xhist, yhist, xtitle = 'Chi-squared', ytitle = 'Number', fill_color = 'blue', $
               title = string(n_samples) + ' realizations')
  
  plothist, fstararr, xhist, yhist, bin = 5E-6, /noprint, /noplot
  p2 = barplot(xhist, yhist, xtitle = 'F star', ytitle = 'Number', fill_color = 'blue', $
               title = string(n_samples) + ' realizations')
  
   plothist, amparr, xhist, yhist, bin = 5E-5, /noprint, /noplot
   p2 = barplot(xhist, yhist, xtitle = 'Amplitude of Phase Curve', ytitle = 'Number', fill_color = 'blue', $
                title = string(n_samples) + ' realizations')

   plothist, phase_maxarr, xhist, yhist, bin =0.005, /noprint, /noplot
   p2 = barplot(xhist, yhist, xtitle = 'Location of Flux maximum in Phase Units', ytitle = 'Number', fill_color = 'blue', $
                title = string(n_samples) + ' realizations')

   plothist, phase_minarr, xhist, yhist, bin =0.01, /noprint, /noplot
   p2 = barplot(xhist, yhist, xtitle = 'Location of Flux minimum in Phase Units', ytitle = 'Number', fill_color = 'blue', $
                title = string(n_samples) + ' realizations')
  
  ;;check duration of code
  print, 'time check: ',systime(1) - t, 'seconds'
end
