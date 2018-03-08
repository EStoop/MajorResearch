function function_fit_lightcurve, planetname, phasearr, fluxarr, errarr, savefilename

 ;first pull the info we know from the internet
  exoplanet_data_file = '/Users/jkrick/idlbin/exoplanets.csv'
  exosystem = strmid(planetname, 0, 7) + ' b' ;'HD 209458 b' ;
  if planetname eq 'WASP-52b' then teq_p = 1315
  if planetname eq 'WASP-14b' then exosystem = 'WASP-14 b'
  if planetname eq 'XO3' then exosystem = 'XO-3 b' 
  if planetname eq 'simul_XO3' then exosystem = 'XO-3 b' 

  get_exoplanet_data,EXOSYSTEM=exosystem,MSINI=msini,MSTAR=mstar,TRANSIT_DEPTH=transit_depth,RP_RSTAR=rp_rstar,AR_SEMIMAJ=ar_semimaj,$
                       TEQ_P=teq_p,TEFF_STAR=teff_star,SECONDARY_DEPTH=secondary_depth,SECONDARY_LAMBDA='4.5',$
                       INCLINATION=inclination,MJD_TRANSIT=mjd_transit,P_ORBIT=p_orbit,EXODATA=exodata,RA=ra,DEC=dec,VMAG=vmag,$
                       DISTANCE=distance,ECC=ecc,T14=t14,F36=f36,F45=f45,FP_FSTAR0=fp_fstar0;,/VERBOSE;verbose

  IF N_ELEMENTS(mjd_start) EQ 0 THEN BEGIN
     mjd_start = get_exoplanet_start(p_orbit,mjd_transit,START_HR_BEFORE_TRANSIT=start_hr_before_transit,START_PHASE=start_phase)
  ENDIF 

 
 ;normalize for more understandable plots
; normcorr = mean(fluxarr)
; fluxarr = fluxarr / normcorr
; errarr = errarr / normcorr

 ;actually want the bottom of eclipse to be at 1.0
 ;fake this for now  XXX
 ;;fluxarr = fluxarr + 0.0015

;and sort them into phase-order
; sp = sort(phasearr)
; phasearr = phasearr[sp]
; fluxarr = fluxarr[sp]
; errarr = errarr[sp]

;remove nan's in error array
 print, 'n before nan', n_elements(phasearr)
 bad = where(finite(errarr) lt 1, badcount, complement = good)
 phasearr = phasearr(good)
 fluxarr = fluxarr(good)
 errarr = errarr(good)
 print, 'n after nan', n_elements(phasearr)

  print, 'n before nan', n_elements(phasearr)
 bad = where(finite(fluxarr) lt 1, badcount, complement = good)
 phasearr = phasearr(good)
 fluxarr = fluxarr(good)
 errarr = errarr(good)
 print, 'n after nan', n_elements(phasearr)


;delete this if I want to overplot on the plot_pixphasecorr plot
 ;XXX will need to add the color offsets.
; testplot = errorplot(phasearr, fluxarr, errarr, '1s', yrange = [0.97, 1.005])

;Initial guess, start with those from the literature
 amplitude = .001             ; messing around with amplitude of the phase curve.
 phase_offset = 1.
 linear = 1.0
 params0 = [fp_fstar0, rp_rstar, ar_semimaj, inclination, amplitude, phase_offset] ;, linear]
 print, 'params0', params0
; params0 = [fp_fstar0, .08, ar_semimaj, inclination, amplitude,
; phase_offset]

;;overplot the starting values before fitting to check where we start
;;out
;;  model = mandel_agol(params0, t=phasearr, flux=fluxarr, err=errarr, p_orbit=p_orbit, mjd_start=mjd_start, mjd_transit=mjd_transit ,savefilename = savefilename, /overplot )


 parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(params0))
; ;limit fp_f_star0 to be positive and less than 1.
 parinfo[0].limited[0] = 1
 parinfo[0].limits[0] = 0.0 ;;
 parinfo[0].limited[1] = 1
 parinfo[0].limits[1] = 1.0
;  parinfo[0].fixed = 1
; ;limit rp/rf_star to be positive and less than 1.
 parinfo[1].limited[0] = 1
 parinfo[1].limits[0] = 0.0 ;;
 parinfo[1].limited[1] = 1
 parinfo[1].limits[1] = 1.0
; parinfo[1].fixed = 1  ; just testing
; ;limit ar_semimaj to be positive 
 parinfo[2].limited[0] = 1
 parinfo[2].limits[0] = 0.0 ;;
; limit inclination to be positive 
 parinfo[3].limited[0] = 1
 parinfo[3].limits[0] = 0.0 ;;
;parinfo[4].fixed = 1
;parinfo[5].fixed = 1
;do the fitting
 afargs = {t:phasearr, flux:fluxarr, err:errarr, p_orbit:p_orbit, mjd_start:mjd_start, mjd_transit:mjd_transit, savefilename:savefilename} ;, rp_rstar:rp_rstar ar_semimaj:ar_semimaj,,inclination:inclination}
 ;;print, 'calling mandel_agol'
 pa = mpfit('mandel_agol', params0, FUNCTARGS=afargs, PERROR=spa, BESTNORM=achi, DOF=adof, COVAR = COV, status = status, errmsg = errmsg, parinfo = parinfo,/quiet) 
;, savefilename = savefilename)

  print, 'status', status
  print, errmsg
  print, 'reduced chi squared',  achi / adof, adof
  print, 'perror', spa

  print, 'primary depth', pa(1)^2, 2*pa(1)*spa(1)  ;error propagation for pa(1)^2
  print, 'secondary depth', pa(0), spa(0)

;want to overplot the fitted curve
  params0 = pa  ; just give it the answer, and run with overplot
  ph = findgen(100) / 100. - 0.5   ; give it a nice set of phases for a nice plot
  model = mandel_agol(params0, t=phasearr, flux=fluxarr, err=errarr, p_orbit=p_orbit, mjd_start=mjd_start, mjd_transit=mjd_transit ,savefilename = savefilename, /overplot );, rp_rstar=rp_rstar, ar_semimaj=ar_semimaj,inclination=inclination,


return, 0
end






