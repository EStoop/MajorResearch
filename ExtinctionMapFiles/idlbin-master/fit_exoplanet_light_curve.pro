pro fit_exoplanet_light_curve, planetname, bin_level, apradius, chname, snapshots = snapshots


 ;first pull the info we know from the internet
  exoplanet_data_file = '/Users/jkrick/idlbin/exoplanets.csv'
  exosystem = strmid(planetname, 0, 7) + ' b' ;'HD 209458 b' ;

  get_exoplanet_data,EXOSYSTEM=exosystem,MSINI=msini,MSTAR=mstar,TRANSIT_DEPTH=transit_depth,RP_RSTAR=rp_rstar,AR_SEMIMAJ=ar_semimaj,$
                       TEQ_P=teq_p,TEFF_STAR=teff_star,SECONDARY_DEPTH=secondary_depth,SECONDARY_LAMBDA='4.5',$
                       INCLINATION=inclination,MJD_TRANSIT=mjd_transit,P_ORBIT=p_orbit,EXODATA=exodata,RA=ra,DEC=dec,VMAG=vmag,$
                       DISTANCE=distance,ECC=ecc,T14=t14,F36=f36,F45=f45,FP_FSTAR0=fp_fstar0,/VERBOSE;verbose

  IF N_ELEMENTS(mjd_start) EQ 0 THEN BEGIN
     mjd_start = get_exoplanet_start(p_orbit,mjd_transit,START_HR_BEFORE_TRANSIT=start_hr_before_transit,START_PHASE=start_phase)
  ENDIF 

                                ;then I need final phase and flux for
                                ;the fits, so open up the best
                                ;photometry to date


  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  stareaor = planetinfo[planetname, 'stareaor']
  snapaor = stareaor + 1        ; the first AOR in the list which is snapshots
  basedir = planetinfo[planetname, 'basedir']
  dirname = strcompress(basedir + planetname+ '/'); +'/hybrid_pmap_nn/')


  if keyword_set(snapshots) then begin
     savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
;  filename = '/Users/jkrick/irac_warm/pcrs_planets/wasp-14b/hybrid_pmap_nn/wasp14_phot_ch2_2.50000.sav'
     print, 'inside fit+exopl', savefilename
     restore, savefilename
     for a = snapaor, n_elements(aorname) -1 do begin
                                ; plothist, planethash[aorname(a),'corrflux'], xhist, yhist, bin = 0.0005, /noplot,/nan
                                ; ap = plot(xhist, yhist)
        meanclip, planethash[aorname(a),'corrflux_d'], meancorr, sigmacorr, clipsig = 2.5 ;,/verbose
        meanclip, planethash[aorname(a),'corrfluxerr'], meancorrerr, sigmacorrerr, clipsig = 2.5
                                ; meanerr, planethash[aorname(a),'corrflux'], planethash[aorname(a),'corrfluxerr'], meancorr2, sigmam, sigmad
        meancorrerr = meancorrerr / 7.;sqrt(n_elements(planethash[aorname(a),'corrflux'])) ;
                                ;print, 'compare mean & error', meancorr, meancorrerr, meancorr2, sigmam, sigmad
        if a eq snapaor then phasearr = [median(planethash[aorname(a),'phase'])] else phasearr = [phasearr, median(planethash[aorname(a),'phase'])]
        if a eq snapaor then fluxarr = [meancorr] else fluxarr = [fluxarr, meancorr]
        if a eq snapaor then errarr = [meancorrerr] else errarr = [errarr, meancorrerr]
     endfor
;now sort the arrays since they are taken at random phase
     sp = sort(phasearr)
     phasearr = phasearr[sp]
     fluxarr = fluxarr[sp]
     errarr = errarr[sp]
  endif else begin   
; working with staring mode, not snapshots
     startaor = 1               ; 0th aor is a peakup aor
     for a = startaor, n_elements(aorname) -1 do begin
        ;restore the pixphasecorr_np save file
        filename =strcompress(dirname +'pixphasecorr_ch'+chname+'_'+aorname(a) +string(apradius)+'.sav',/remove_all)
        print, a, ' ', aorname(a), 'restoring', filename
        restore, filename

        ;pull together the whole light curve
        if a eq startaor then phasearr = phase else phasearr = [phasearr, phase]
        if a eq startaor then fluxarr = flux_np else fluxarr = [fluxarr, flux_np] ; want this to be the nearest_neighbor with np light curve
        if a eq startaor then errarr = fluxerr_np else errarr = [errarr, fluxerr_np]  ; and error bars on above
     endfor
; if not snapshots, need to do some binning XXX

  endelse
  
 


 ;normalize for more understandable plots
 normcorr = mean(fluxarr)
 fluxarr = fluxarr / normcorr
 fluxarr = fluxarr + 0.001
 errarr = errarr / normcorr

;what happens if I intentionally offset the minimum of the sin curve
 low = where(phasearr lt -0.05 and phasearr gt -0.3)
; fluxarr(low) = fluxarr(low) - 0.001


 testplot = errorplot(phasearr, fluxarr, errarr, '1s', yrange = [0.991, 1.005])

;Initial guess, start with those from the literature
 amplitude = .001               ; messing around with amplitude of the phase curve.
 phase_offset = 1.
 c1 = -.04
 c2 = .04
 c3 = 1.0
 c4 = 1.0
; params0 = [fp_fstar0, rp_rstar, ar_semimaj, inclination, amplitude, phase_offset]
 params0 = [fp_fstar0, rp_rstar, ar_semimaj, inclination, c1, c2];, c3, c4]
;  params0 = [fp_fstar0, rp_rstar, ar_semimaj, inclination]
 parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(params0))
; ;limit fp_f_star0 to be positive and less than 1.
  parinfo[0].limited[0] = 1
  parinfo[0].limits[0] = 0. ;;
  Parinfo[0].limited[1] = 1
  parinfo[0].limits[1] = 1.0
;  parinfo[0].fixed = 1
; ;limit rp/rf_star to be positive and less than 1.
  parinfo[1].limited[0] = 1
  parinfo[1].limits[0] = 0.0 ;;
  parinfo[1].limited[1] = 1
  parinfo[1].limits[1] = 1.0
; ;limit ar_semimaj to be positive 
  parinfo[2].limited[0] = 1
  parinfo[2].limits[0] = 0.0 ;;
; limit inclination to be positive 
  parinfo[3].limited[0] = 1
  parinfo[3].limits[0] = 0.0 ;;
 ;limit c1 and c2 to be gt -1 and lt 1
;  parinfo[4].limited[0] = 1
;  parinfo[4].limits[0] = -1
;  parinfo[4].limited[1] = 1
;  parinfo[4].limits[1] = 1
;  parinfo[5].limited[0] = 1
;  parinfo[5].limits[0] = -1
;  parinfo[5].limited[1] = 1
;  parinfo[5].limits[1] =1

;do the fitting
  modelfilename = strcompress(dirname + planetname +'_model_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)

  afargs = {t:phasearr, flux:fluxarr, err:errarr, p_orbit:p_orbit, mjd_start:mjd_start, mjd_transit:mjd_transit, savefilename:modelfilename};, rp_rstar:rp_rstar ar_semimaj:ar_semimaj,,inclination:inclination}
  pa = mpfit('mandel_agol', params0, FUNCTARGS=afargs, PERROR=spa, BESTNORM=achi, DOF=adof, COVAR = COV, status = status, errmsg = errmsg, parinfo = parinfo);, savefilename = savefilename)

  print, 'status', status
  print, errmsg
  print, 'reduced chi squared',  achi / adof
  

;want to overplot the fitted curve
  params0 = pa  ; just give it the answer, and run with overplot
  ph = findgen(100) / 100. - 0.5   ; give it a nice set of phases for a nice plot
  model = mandel_agol(params0, t=ph, flux=fluxarr, err=errarr, p_orbit=p_orbit, mjd_start=mjd_start, mjd_transit=mjd_transit ,savefilename = modelfilename, /overplot);, rp_rstar=rp_rstar, ar_semimaj=ar_semimaj,inclination=inclination,
end






