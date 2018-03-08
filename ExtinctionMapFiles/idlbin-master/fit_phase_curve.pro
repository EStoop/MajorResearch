;
;;FIT_PHASE_CURVE.PRO
;;
;SYNTAX:  
;       IDL> FIT_PHASE_CURVE,t,f,exosystem,fstar,albedo,trad,omrot,tperi,tphase,model,/TIDALLY_LOCKED,$
;                            EXODATA=exodata,LAMBDA_BAND=lambda_band,FUNC=func,/VERBOSE,$
;                            PARAM_ERR=param_err,SYSTEM_ERR=system_err,/FIT_SYSTEM,$
;                            REDUCED_CHI2=reduced_chi2,AR_SEMIMAJ=ar_semimaj,ECC=ecc,$
;                            ARGPERI=argperi,RP_RSTAR=rp_rstar,TEFF_STAR=teff_star,$
;                            PORBIT=porbit,INCLINATION=inclination,MJD_TRANSIT=mjd_transit,$
;                            LIMB_DARKENING=limb_darkening,TPERI_2=tperi_2,MODEL_2=model_2,TPHASE_2=tphase_2,$
;                            PHASE_MAX=phase_max,PHASE_MIN=phase_min,T_MAX=t_max,$
;                            T_MIN=t_min,AMPLITUDE=amplitude,/FIXED_FSTAR,/FIXED_TRAD,$
;                            /FIXED_OMROT,TRANSIT_DEPTH=transit_depth,ECLIPSE_DEPTH=eclipse_depth
;                            
;   Fit the phase curve of exoplanet system EXOSYSTEM
;
;   **INPUTS**
;         T  :  Time, in MJD.
;         F  :  Measured Flux of exoplanet system
;  EXOSYSTEM :  Exoplanet system name (string), as listed in exoplanets.org database.
;               Type
;               get_exoplanet_data,EXODATA=exodata,LIST_NAMES='*'
;               to list all the names available.  Or narrow down the list using wild cards
;               (LIST_NAMES = 'WASP*' , LIST_NAMES = 'HD*' etc)
;    FSTAR   :  (can be input as undefined) - Stellar flux, same units as F, initial guess.  Output as fit value.
;    ALBEDO  :  (can be input as undefined) - Albedo, input as a number if you want its value to vary in the fit.  
;               If input as !null, or an undefined variable, it is held fixed at 0.0.  Output as fit value.
;      TRAD  :  (can be input as undefined) - Radiative timescale, initial guess. Output as fit value.
;     OMROT  :  (can be input as undefined) - Rotational angular velocity in units of the periastron orbital angular velocity,
;                            initial guess.  Output as fit value.
;
;   **OUTPUTS**
;      TPERI :  Input time, converted to time relative to periastron, days.
;     TPHASE :  Input time expressed as an orbital "phase", i.e., fraction of orbit since transit
;      MODEL :  Array of model values of relative flux as a function of input time
;
;  **KEYWORDS**
;          FUNC    : Uncertainty in flux values
;     /FIXED_FSTAR : Hold fixed the value of fstar in the fit to the input value.
;      /FIXED_TRAD : Hold fixed the value of trad in the fit to the input value.
;     /FIXED_OMROT : Hold fixed the value of omrot in the fit to the input value.
;  /TIDALLY_LOCKED : Force OMROT to be equal to the instantaneous orbital angular velocity [not currently available]
;         EXODATA  : Exoplanet data structure.  Useful to pass this to avoid having to reread the database file.
;    LAMBDA_BAND = : wavelength band at which the measurement is to be made.  Values can be
;                   'U', 'B', 'V', 'R', 'I', 'J', 'H', 'K', 'L', 'M', 'N', 'O', 'J2MASS', 'H2MASS', 'Ks2MASS',
;                   'IRAC3.6','IRAC4.5', 'IRAC5.7', 'IRAC8.0', 'g', 'r', 'i', 'z'.  Default is 'IRAC4.5'.
;      /FIT_SYSTEM : Fit exoplanet system parameters RP_RSTAR, INCLINATION, and MJD_TRANSIT.
;       /VERBOSE   : Print updates to screen, plot fitting progress.
;       PARAM_ERR  : Vector with formal errors in Albedo, T_RAD, OMEGA_ROT
;      SYSTEM_ERR  : Vector with formal errors in  RP_RSTAR, INCLINATION, and MJD_TRANSIT, if /FIT_SYSTEM is set
;    REDUCED_CHI2  : Reduced Chi^2 in fit
;    TPERI_2       : Periastron time in 1000 equally spaced intervals for smooth model curve
;    TPHASE_2      : Time in 1000 equally spaced intervals for smooth model curve, expressed in terms of orbital phase,
;                    or fraction of orbit since transit.
;    MODEL_2       : Model values at TPERI_2 or TPHASE_2
;;         T_MAX   : Model TPERI of peak amplitude.  May not be an actual time in the returned TPERI_2 array.
;;         T_MIN   : Model TPERI of minimum amplitude.  May not be an actual time in the returned TPERI_2 array.
;;     PHASE_MAX   : [TPHASE_2 of peak amplitude, uncertainty] (two-element vector).  
;;                   May not be an actual value in the TPHASE_2 array.
;;     PHASE_MIN   : [TPHASE_2 of minimum amplitude, uncertainty] (two-element vector). 
;;                    May not be an actual value in the TPHASE_2 array.
;;     AMPLITUDE   : [Peak-to-valley amplitude of phase curve, in units of FSTAR, uncertainty] (two-element vector)
;;   TRANSIT_DEPTH : [Depth of transit in the model, uncertainty] (two-element vector), i.e., the flux at mid-transit divided by the
;;                       unobstructed star + planet flux at that time.
;;   ECLIPSE_DEPTH : [Depth of eclipse in the model, uncertainty] (two-element vector), i.e., the flux at mid-eclipse divided by the
;;                       unobstructed star + planet flux at that time.
;
;  **EXOPLANET SYSTEM PARAMETERS**
;    (Except for LIMB_DARKENING, these do not need to be input because they will be read from the database, BUT
;     they may be input to override values from database)
;    (Parameters indicated by asterisks "*" are fittable if /FIT_SYSTEM keyword is set.
;     In this case, input or database values will be used as a guess.  If neither is available, use default guess.)
;   AR_SEMIMAJ     : Semimajor axis, in units of the stellar radius
;          ECC     : Orbital eccentricity
;     ARGPERI      : Argument of periastron (degrees) (angle between periastron and plane of sky (transit minus 90 deg)
;     RP_RSTAR*    : Planet radius, in units of the stellar radius
;    TEFF_STAR     : Stellar effective temperature (K)  -- NO LONGER USED IN FIT
;       PORBIT     : Orbital period (dy)
;  INCLINATION*    : Orbital inclination (degrees)
;  MJD_TRANSIT*    : MJD of transit (dy)
; LIMB_DARKENING   : The set of stellar limb darkening coefficients to use to compute the transit profile.
;                    See Claret (2000) for a description of the model, and Mandel & Agol (2002) for the
;                    effect on transit profiles.  Leave undefined for no limb darkening.  
;                    NOT IN DATABASE -- MUST BE INPUT.
;



PRO fpc_time,t,mjd_transit,porbit,ecc,argperi,ttransit,tperi,tphase
   degrad = !dpi / 180d 
   tphase = ((t - mjd_transit)/porbit mod 1)
   ttransit = tphase * porbit ;;; number of days since the most recent transit
   theta_transit = !dpi/2. - argperi*degrad   ;; Angle between periastron and transit
   tperi_0 = [0:1:0.01] * porbit ;; grid of values of tperi
   ELLIPTICAL_ORBIT,tperi_0,porbit,ecc,r,theta  ;; Solve for the elliptical orbit
   ;; Predict the value of tperi of transit
   st = SORT(theta)
   tperi_transit = INTERPOL(tperi_0[st],theta[st]*degrad,theta_transit,/NAN)
;; Convert ttransit to tperi
   tperi = ttransit + tperi_transit
;; Condition TPERI values so they are within one orbit
   tperi = porbit * ( tperi / porbit mod 1 )
RETURN
END

FUNCTION PHASE_CURVE_FUNCTION,p,TIME=t,FLUX=flux,FUNC=func,AR_SEMIMAJ=ar_semimaj,ECC=ecc,ARGPERI=argperi, RP_RSTAR=rp_rstar,$
                              TEFF_STAR=teff_star,PORBIT=porbit, LAMBDA_BAND=lambda_band, INCLINATION=inclination,$
                              MODEL=model,A_AU=a_au,LIMB_DARKENING=limb_darkening
  FSTAR = P[0]                              
  ALBEDO=P[1]
  TRAD=P[2]
  OMROT=P[3]  
  IF N_ELEMENTS(P) EQ 7 THEN BEGIN
     rp_rstar=P[4]
     inclination=P[5]
     mjd_transit=P[6]
     fpc_time,t,mjd_transit,porbit,ecc,argperi,ttransit,tperi,tphase
  ENDIF ELSE tperi = t
  ikeep = where(FINITE(flux),nkeep,complement=inan)
  MODEL=MAKE_ARRAY(N_ELEMENTS(flux),/DOUBLE,VALUE=0)
  resid = model
  IF nkeep ne 0 THEN BEGIN      
    ;;; Compute a model light curve for one orbit.                        
     exoplanet_phase_curve,ar_semimaj,ecc,argperi,rp_rstar,teff_star,porbit,lambda_band,albedo,trad,omrot,tperi_model,$
                           fp_fstar,relative_flux,INCLINATION=inclination,A_AU=a_au,LIMB_DARKENING=limb_darkening
     st = SORT(tperi_model)
     nt = N_ELEMENTS(tperi_model)                           
     ;;; triple the orbit in the model, so interpolation works at the boundaries
     tperi_model = [tperi_model[st]-porbit,tperi_model[st],tperi_model[st]+porbit]
     relative_flux = [relative_flux[st],relative_flux[st],relative_flux[st]]
     model[ikeep] = INTERPOL(relative_flux,tperi_model,tperi[ikeep],/NAN) * fstar
     resid[ikeep] = (flux[ikeep]-model[ikeep])/func[ikeep]
  ENDIF
  model[inan] = !values.d_nan
RETURN,resid
END

PRO phase_curve_show_fit,myfunct,p,iter,fnorm,FUNCTARGS=fa,PARINFO=parinfo,QUIET=quiet,DOF=dof,PFORMAT=pformat,UNIT=unit

  IF ~KEYWORD_SET(QUIET) THEN BEGIN
    resid = PHASE_CURVE_FUNCTION(p,TIME=fa.time,FLUX=fa.flux,FUNC=fa.func,AR_SEMIMAJ=fa.ar_semimaj,ECC=fa.ecc,$
                                 ARGPERI=fa.argperi,RP_RSTAR=fa.rp_rstar,TEFF_STAR=fa.teff_star,PORBIT=fa.porbit,$
                                 LAMBDA_BAND=fa.lambda_band,INCLINATION=fa.inclination,MODEL=model,A_AU=fa.a_au,$
                                 LIMB_DARKENING=fa.limb_darkening)
    set_plot,'x'
    st = sort(fa.time)
    CGplot,fa.time[st],fa.flux[st]/p[0],psym=3,xstyle=1,ystyle=1,yrange=[0.98,1.02]
    CGplot,fa.time[st],model[st]/p[0],/OVERPLOT,color='Medium Gray';,psym=3
    ifin = WHERE(FINITE(fa.flux) and FINITE(fa.func),nfin) 
    print,'Iteration ',iter,' CHISQ: ',fnorm/(nfin-N_ELEMENTS(P))
    print,'FSTAR: ',p[0],' ALBEDO: ',p[1],' T_RAD: ',P[2],' OMEGA_ROT: ',P[3]
    IF N_ELEMENTS(P) EQ 7 THEN print,'RP_RSTAR: ',p[4],' inclination: ',p[5],' MJD_TRANSIT: ',P[6]
    cr = ''
    IF iter GT 1 and ITER mod 20 EQ 0 THEN read,'CR to Continue',cr
  ENDIF
END

PRO FIT_PHASE_CURVE,t,f,exosystem,fstar,albedo,trad,omrot,tperi,tphase,model,TIDALLY_LOCKED=tidally_locked,$
                    EXODATA=exodata,LAMBDA_BAND=lambda_band,FUNC=func,VERBOSE=verbose,PARAM_ERR=param_err,$
                    REDUCED_CHI2=reduced_chi2,AR_SEMIMAJ=ar_semimaj,ECC=ecc,ARGPERI=argperi,RP_RSTAR=rp_rstar,$
                    TEFF_STAR=teff_star,PORBIT=porbit,INCLINATION=inclination, MJD_TRANSIT=mjd_transit,$
                    TPERI_2=tperi_2,MODEL_2=model_2,TPHASE_2=tphase_2,FIT_SYSTEM=fit_system,SYSTEM_ERR=system_err,$
                    PHASE_MAX=phase_max,PHASE_MIN=phase_min,T_MAX=t_max,T_MIN=t_min,AMPLITUDE=amplitude,$
                    FIXED_FSTAR=fixed_fstar,FIXED_TRAD=fixed_trad,FIXED_OMROT=fixed_omrot,$
                    TRANSIT_DEPTH=transit_depth,ECLIPSE_DEPTH=eclipse_depth,LIMB_DARKENING=limb_darkening
                    
  IF N_ELEMENTS(LIMB_DARKENING) EQ 0 THEN limb_darkening = 0                    
;
;  (0) Set values of constants
   degrad = !dpi / 180d  
   au = 149597870700d2  ;; cm
   rsun = 6.955d10 ;; cm
   IF N_ELEMENTS(FUNC) EQ 0 THEN func = MAKE_ARRAY(SIZE=SIZE(f),VALUE=1)
   v = KEYWORD_SET(VERBOSE) 
   SIMPLIFIED_LAMBDA = HASH('J','J','H','H','Ks2MASS','KS','K','KP','IRAC3.6','3.6','IRAC4.5','4.5','IRAC5.7','5.8','IRAC8.0','8.0')
      
   IF V THEN PRINT,'(1) Gather data for the exoplanet system'
   get_exoplanet_data,EXOSYSTEM=exosystem,MSINI=msini,MSTAR=mstar,RSTAR=rstar,TRANSIT_DEPTH=transit_depth,RP_RSTAR=rp_rstar0,AR_SEMIMAJ=ar_semimaj0,$
                      TEQ_P=3000.0,TEFF_STAR=teff_star0,SECONDARY_DEPTH=secondary_depth,SECONDARY_LAMBDA=SIMPLIFIED_LAMBDA[lambda_band],OM=argperi0,$
                      INCLINATION=inclination0,MJD_TRANSIT=mjd_transit0,P_ORBIT=porbit0,EXODATA=exodata,RA=ra,DEC=dec,VMAG=vmag,$
                      DISTANCE=distance,ECC=ecc0,VERBOSE=verbose

   ;; If no input values do not override database values               
   IF ~N_ELEMENTS(AR_SEMIMAJ) THEN ar_semimaj=ar_semimaj0
   IF ~N_ELEMENTS(ECC) THEN ecc=ecc0
   IF ~N_ELEMENTS(ARGPERI) THEN argperi=argperi0
   IF ~N_ELEMENTS(RP_RSTAR) THEN rp_rstar=rp_rstar0
   IF ~N_ELEMENTS(TEFF_STAR) THEN teff_star=teff_star0
   IF ~N_ELEMENTS(PORBIT) THEN porbit=porbit0
   IF ~N_ELEMENTS(INCLINATION0) THEN inclination0 = 90.0
   IF ~N_ELEMENTS(INCLINATION) THEN inclination=inclination0
   IF ~N_ELEMENTS(MJD_TRANSIT) THEN mjd_transit=mjd_transit0
   IF ~KEYWORD_SET(FIT_SYSTEM) THEN BEGIN
      IF V THEN PRINT,'(2) Convert to time in days since periastron'
      fpc_time,t,mjd_transit,porbit,ecc,argperi,ttransit,tperi,tphase
      nsyspar = 0
   ENDIF ELSE BEGIN
      tperi = t  ;;; placeholder
      nsyspar = 3
   ENDELSE
   
   IF V THEN PRINT,'(3) Set up parameter initial guess and constraints.'
   parinfo = replicate({value:0.D, fixed:0, limited:[0,0], $
                         limits:[0.D,0], parname:''}, 4+nsyspar)
   PARINFO[0].parname = 'F_STAR'
   IF N_ELEMENTS(fstar) NE 0 THEN PARINFO[0].value = fstar
   PARINFO[0].limited = [1,0]
   IF KEYWORD_SET(FIXED_FSTAR) THEN PARINFO[0].fixed=1
   PARINFO[1].parname = 'ALBEDO'
   IF N_ELEMENTS(ALBEDO) NE 0 THEN BEGIN
      PARINFO[1].limited = [1,1]
      PARINFO[1].limits = [0,1d]
      PARINFO[1].value = albedo
   ENDIF ELSE BEGIN
      PARINFO[1].value = 0d0
      PARINFO[1].FIXED = 1
   ENDELSE
   PARINFO[2].parname = 'T_RAD'
   IF N_ELEMENTS(trad) NE 0 THEN PARINFO[2].value = trad
   PARINFO[2].limited = [1,0]
   IF KEYWORD_SET(FIXED_TRAD) THEN PARINFO[2].fixed=1
   ;PARINFO[2].limits = [0,0.5*PORBIT]  ;; Force radiative timescales to be less than 1/2 orbital period
   PARINFO[3].parname = 'OMEGA_ROT'
   IF N_ELEMENTS(omrot) NE 0 THEN PARINFO[3].value = omrot
   PARINFO[3].limited = [1,0]
   IF KEYWORD_SET(FIXED_OMROT) THEN PARINFO[3].fixed=1
   IF KEYWORD_SET(FIT_SYSTEM) THEN BEGIN
      PARINFO[4].parname='RP_RSTAR'
      PARINFO[4].value = rp_rstar
      PARINFO[4].limited=[1,0]
      
      PARINFO[5].parname='INCLINATION'
      PARINFO[5].value=inclination
      PARINFO[5].limited=[1,1]
      PARINFO[5].limits=[0,180.0]
      
      PARINFO[6].parname='MJD_TRANSIT'
      PARINFO[6].value=mjd_transit
      PARINFO[6].limited=[1,0]
   ENDIF
   IF V THEN PRINT,'(4) Set values of function arguments'
   FUNCTARGS = {AR_SEMIMAJ:ar_semimaj, ECC:ecc, ARGPERI:argperi, RP_RSTAR:rp_rstar, TEFF_STAR:teff_star, MODEL: DBLARR(N_ELEMENTS(t)), $
                PORBIT:porbit, LAMBDA_BAND:lambda_band, INCLINATION:inclination, A_AU:ar_semimaj * rstar * rsun/au, $
               LIMB_DARKENING: limb_darkening, TIME:tperi, FLUX:f, FUNC:func}
   IF V THEN PRINT,'(5) Run the fit'
   PARAM_RESULT = MPFIT('PHASE_CURVE_FUNCTION',FUNCTARGS=functargs,PARINFO=parinfo,AUTODERIVATIVE=1,$
                         QUIET=~V,PERROR=perror,BESTNORM=bestnorm,status=status,errmsg=errmsg,$
                         iterproc='phase_curve_show_fit',XTOL=1d-6,COVAR=covar,NOCATCH=0)
   IF v THEN print,'(6) Done fit.  Status=',status,' ',errmsg
   iFIXED = WHERE(parinfo[*].FIXED EQ 1,nfixed,ncomplement=nfloat,complement=ifloat)
   DOF = N_ELEMENTS(tperi) - N_ELEMENTS(params) - nfixed
   reduced_chi2 = SQRT(BESTNORM/DOF)
   param_err = perror*reduced_chi2
   fstar = param_result[0]
   albedo = param_result[1]
   trad = param_result[2]
   omrot = param_result[3]
   IF KEYWORD_SET(FIT_SYSTEM) THEN BEGIN
      rp_rstar = param_result[4]
      inclination = param_result[5]
      mjd_transit = param_result[6]
      fpc_time,t,mjd_transit,porbit,ecc,argperi,ttransit,tperi,tphase
   ENDIF
   exoplanet_phase_curve,ar_semimaj,ecc,argperi,rp_rstar,teff_star,porbit,lambda_band,albedo,trad,omrot,tperi_2,$
     fp_fstar,model_2,INCLINATION=inclination,A_AU=ar_semimaj * rstar * rsun/au,PHASE_FRACTION=tphase_2,$
     amplitude=amplitude,phase_max=phase_max,phase_min=phase_min,t_max=t_max,t_min=t_min,$
     TRANSIT_DEPTH=transit_depth,ECLIPSE_DEPTH=eclipse_depth,LIMB_DARKENING=limb_darkening
   IF V THEN PRINT,'(7) Compute uncertainties in calculated light curve parameters.'
   ;; GET THE PARTIAL DERIVATIVES and the UNCERTAINTIES in amplitude, phase_max, phase_min, transit_depth, and eclipse_depth
   amplitude_unc = 0d
   phase_max_unc = 0d
   phase_min_unc = 0d
   transit_depth_unc = 0d
   eclipse_depth_unc = 0d
   FOR i = 0,nfloat-1 DO BEGIN ;; Loop over all floating parameters
      dx = param_err[ifloat[i]]
      please_continue=0
      CASE PARINFO[ifloat[i]].parname OF
        'T_RAD': BEGIN
          x = dblarr(3)
          IF TRAD EQ 0 THEN BEGIN
              x[0] = trad+dx
              x[1] = trad
              x[2] = trad+dx*2d0   
          ENDIF ELSE BEGIN
              x[0] = trad-dx > 0
              x[1] = trad
              x[2] = trad+dx
          ENDELSE
          exoplanet_phase_curve,ar_semimaj,ecc,argperi,rp_rstar,teff_star,porbit,lambda_band,albedo,x[2],omrot,tperi_2,$
            fp_fstar,model_2,INCLINATION=inclination,A_AU=ar_semimaj * rstar * rsun/au,PHASE_FRACTION=tphase_2,$
            amplitude=amplitude2,phase_max=phase_max2,phase_min=phase_min2,t_max=t_max,t_min=t_min,$
            TRANSIT_DEPTH=transit_depth2,ECLIPSE_DEPTH=eclipse_depth2,LIMB_DARKENING=limb_darkening
          exoplanet_phase_curve,ar_semimaj,ecc,argperi,rp_rstar,teff_star,porbit,lambda_band,albedo,x[0],omrot,tperi_2,$
            fp_fstar,model_2,INCLINATION=inclination,A_AU=ar_semimaj * rstar * rsun/au,PHASE_FRACTION=tphase_2,$
            amplitude=amplitude0,phase_max=phase_max0,phase_min=phase_min0,t_max=t_max,t_min=t_min,$
            TRANSIT_DEPTH=transit_depth0,ECLIPSE_DEPTH=eclipse_depth0,LIMB_DARKENING=limb_darkening
        END
        'OMEGA_ROT': BEGIN
          x = dblarr(3)
          x[0] = omrot-dx > 0
          x[1] = omrot
          x[2] = omrot+dx
          exoplanet_phase_curve,ar_semimaj,ecc,argperi,rp_rstar,teff_star,porbit,lambda_band,albedo,trad,x[2],tperi_2,$
            fp_fstar,model_2,INCLINATION=inclination,A_AU=ar_semimaj * rstar * rsun/au,PHASE_FRACTION=tphase_2,$
            amplitude=amplitude2,phase_max=phase_max2,phase_min=phase_min2,t_max=t_max,t_min=t_min,$
            TRANSIT_DEPTH=transit_depth2,ECLIPSE_DEPTH=eclipse_depth2,LIMB_DARKENING=limb_darkening
          exoplanet_phase_curve,ar_semimaj,ecc,argperi,rp_rstar,teff_star,porbit,lambda_band,albedo,trad,x[0],tperi_2,$
            fp_fstar,model_2,INCLINATION=inclination,A_AU=ar_semimaj * rstar * rsun/au,PHASE_FRACTION=tphase_2,$
            amplitude=amplitude0,phase_max=phase_max0,phase_min=phase_min0,t_max=t_max,t_min=t_min,$
            TRANSIT_DEPTH=transit_depth0,ECLIPSE_DEPTH=eclipse_depth0,LIMB_DARKENING=limb_darkening
        END
        ELSE: please_continue=1
          ;;; HERE -- add cases
      ENDCASE
      IF please_continue THEN CONTINUE
      sx = SORT(x) 
      damp = (DERIV(x[sx],([amplitude0,amplitude,amplitude2])[sx]))[1]
      dpmx = (DERIV(x[sx],([phase_max0,phase_max,phase_max2])[sx]))[1]
      dpmn = (DERIV(x[sx],([phase_min0,phase_min,phase_min2])[sx]))[1]
      dtran = (DERIV(x[sx],([transit_depth0,transit_depth,transit_depth2])[sx]))[1]
      decl = (DERIV(x[sx],([eclipse_depth0,eclipse_depth,eclipse_depth2])[sx]))[1]
      amplitude_unc += (damp * dx)^2
      phase_max_unc += (dpmx * dx)^2
      phase_min_unc += (dpmn * dx)^2
      transit_depth_unc += (dtran * dx)^2
      eclipse_depth_unc += (decl * dx)^2
   ENDFOR
  amplitude_unc = SQRT(amplitude_unc)
  phase_max_unc = SQRT(phase_max_unc)
  phase_min_unc = SQRT(phase_min_unc)
  transit_depth_unc = SQRT(transit_depth_unc)
  eclipse_depth_unc = SQRT(eclipse_depth_unc) 
  amplitude = [amplitude,amplitude_unc]
  phase_max = [phase_max,phase_max_unc]
  phase_min = [phase_min,phase_min_unc]
  transit_depth = [transit_depth,transit_depth_unc]
  eclipse_depth = [eclipse_depth,eclipse_depth_unc]
RETURN
END