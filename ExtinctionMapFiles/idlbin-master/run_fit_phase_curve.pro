pro run_fit_phase_curve, exodata = exodata
 ;;wrapper to fit phase curves with Jim's model of C&A11

  t = systime(1)

  ;;setup some known quantities for wasp14
  porbit = 2.24376507D          ;days
  mjd_transit = 2456034.21261D - 2400000.5D
  ecc = 0.083D
  argperi = 252.67D             ;degrees
  ar_semimaj = 5.99d
  rp_rstar = 0.09420d
  inclination=84.61d
  exosystem = 'WASP-14 b'
  verbose = 1
  trad =1.d ;2.  ;; initial guess
  omrot =1.d ;0.9  ;; initial guess
  fstar = 0.0577d ;; initial guess


  ;;setup some known quantities for XO3
;  porbit = 3.1915239D           ;days
;  mjd_transit = 2454449.86816D - 2400000.5D
;  ecc = 0.26D
;  argperi = 345.8D             ;degrees
;  ar_semimaj = 5.99d ;XXX
;  rp_rstar = 0.09420d ;XXX
;  inclination=84.79d
;  exosystem = 'XO-3 b'
;  verbose = 1
;  trad =1.d ;2.  ;; initial guess
;  omrot =1.d ;0.9  ;; initial guess
;  fstar = 0.049d ;; initial guess
;  teff_star = 6781.0d

  startvals = [trad, omrot, fstar]
  albedo = !NULL  ;; set to null to fix the value at zero in the fit
  IF N_ELEMENTS(RP_RSTAR) EQ 0 AND ~KEYWORD_SET(FIT_SYSTEM) THEN rp_rstar = 0.09420d
  IF N_ELEMENTS(INCLINATION) EQ 0 AND ~KEYWORD_SET(FIT_SYSTEM) THEN inclination=84.79d

  ;;read in the snapshot mode photometry
  starefile = '/Users/jkrick/irac_warm/pcrs_planets/wasp-14b/wasp14_snap_latenterr_150226.txt'
  starefile = '/Users/jkrick/irac_warm/pcrs_planets/wasp-14b/wasp14_snap_phaseonly_150226.txt'
  readcol, starefile, starebmjd, starephase, corrflux, corrfluxerr, format = 'D, D, D, D',/nan

  ;;in this instance want to not include transit and eclipse
 ; good1 = where(starephase gt -0.47 and starephase lt -0.04)
 ; good2 = where(starephase gt 0.03 and starephase lt 0.45)
 ; goodphase = [starephase(good1), starephase(good2)]
 ; goodcorrflux = [corrflux(good1), corrflux(good2)]
 ; goodcorrfluxerr = [corrfluxerr(good1), corrfluxerr(good2)]
 ; goodbmjd = [starebmjd(good1), starebmjd(good2)]
 
  ;;read in the save file
 ; restore, '/Users/jkrick/external/irac_warm/XO3/XO3_phot_ch2_2.25000_150226.sav'
 ; corrflux = planethash['r46471424', 'corrflux']
 ; corrfluxerr = planethash['r46471424', 'corrfluxerr']
 ; starebmjd =  planethash['r46471424', 'bmjdarr']

  ;;test output by plotting
  ;;p1 = plot(goodphase, goodcorrflux, '1s', sym_size =0.2, sym_filled= 1, xtitle = 'phase', ytitle = 'corrflux')
     

  ;;run the fitting iteratively. start with fixed number of iterations
  for i = 0, 0 do begin
     print, '-------------------------'
     print, 'starting fit', fstar, trad, omrot
     ;;now do some fitting
     FIT_PHASE_CURVE,starebmjd,corrflux,exosystem,fstar,albedo,trad,omrot,tperi,tphase,model,$
                     EXODATA=exodata,LAMBDA_BAND='IRAC4.5',VERBOSE=verbose,PARAM_ERR=param_err,$
                     REDUCED_CHI2=reduced_chi2,ECC=ecc,ARGPERI=argperi,$
                     PORBIT=porbit,MJD_TRANSIT=mjd_transit,$;LIMB_DARKENING=limb_darkening, $
                     RP_RSTAR=rp_rstar,AR_SEMIMAJ=ar_semimaj,$
                     INCLINATION=inclination,TPERI_2=tperi_2,MODEL_2=model_2,FUNC=corrfluxerr,$
                     TPHASE_2=tphase_2,AMPLITUDE=amplitude,PHASE_MAX=phase_max,$
                     PHASE_MIN=phase_min,TRANSIT_DEPTH=transit_depth,ECLIPSE_DEPTH=eclipse_depth ;/FIT_SYSTEM ,


     print, 'final fit amplitude', amplitude; amplitude[0],amplitude[1]
     print, 'phase_max', phase_max[0], phase_max[1]
     print, 'phase_min', phase_min[0], phase_min[1]
     ;;print, 'param err albedo, trad, omrot', param_err
     ;;print, 'other params', rp_rstar, inclination, mjd_transit
     print, 'transit and eclipse depth', transit_depth[0], transit_depth[1], eclipse_depth[0], eclipse_depth[1]
     print, 'trad, omrot, fstar, chi2', trad, omrot, fstar, reduced_chi2
  ;;check duration of code
     print, 'time check: ',systime(1) - t, 'seconds'

  endfor

end
