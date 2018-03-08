PRO exoplanet_light_curve,t,rel_flux,EXOSYSTEM=exosystem,MSINI=msini,MSTAR=mstar,TRANSIT_DEPTH=transit_depth,$
    TEQ_P=teq_p,TEFF_STAR=teff_star,INCLINATION=inclination,MJD_TRANSIT=mjd_transit,AR_SEMIMAJ=ar_semimaj,P_ORBIT=p_orbit,$
    NPERIOD=nperiod,MJD_START=mjd_start,SECONDARY_DEPTH=secondary_depth,SECONDARY_LAMBDA=secondary_lambda,RA=ra,DEC=dec,$
    VMAG=vmag,DIST=distance,ECC=ecc,RP_RSTAR=rp_rstar,T14=t14,NT=nt,PLOT=plot,PDF=pdf,VERBOSE=verbose,F36=f36,F45=f45,$
    EXODATA=exodata,TPHASE=tphase,DURATION_HR=duration_hr,LABEL=label,START_HR_BEFORE_TRANSIT=start_hr_before_transit,$
    START_PHASE=start_phase,_extra=extra
;;
;;   Model the light curve of an exoplanet, including transits, occultations, and phase curves.  Either input the parameters of the
;;   system by hand or use an actual planet from exoplanets.org.  Right now we are assuming eccentricity of zero (circular), 
;;   and that all systems are tidally locked and thermally emitting.   
;; 
;  OUTPUT PARAMETERS
;    T - time array, in MJD
;    REL_FLUX - flux of system relative to the flux of the star 
;
;; INPUT KEYWORDS
;;
;;  EXOSYSTEM - the name of the exoplanet system, a string found in the exoplanets.org database (eg., "WASP-46 b").  
;;              If this is given then the other keywords can be passed as outputs.  Otherwise they will all have to be input
;;              to specify the parameters of the orbits.
;;  NPERIOD - the number of orbital periods to return.  Default is 1.  If DURATION_HR is set, then this returns the actual number of periods returned.
;;  DURATION_HR - the duration in hours to return.  Default is the orbital period.  Overrides NPERIOD.
;;  MJD_START - modified julian date at the start of the observation.  Default is MJD_TRANSIT - P/2, if MJD_TRANSIT
;;              exists.  Otherwise default is now.
;;  START_HR_BEFORE_TRANSIT - number of hours prior to transit to start the observation.  Phased to the nearest time subsequent to now.
;;  START_PHASE - phase of orbit at which to start the observation.  Phased to the nearest time subsequent to now.
;;  NT - number of elements in time array, default is 10000
;;  LABEL - extra label for the plot file (in addition to the exoplanet system name)
;;  /PLOT - set this if we want to plot the lightcurve
;;  /PDF  - send the plot to a pdf file, otherwise to screen
;;  /VERBOSE - print some progress to the standard output
;;  /TPHASE - return time in orbital phase 
;;
;; INPUTS if not available via the exoplanet database (or EXOSYSTEM is not given), OUTPUTS if determined within the code
;;
;;  MSINI - Minimum mass of planet, MEarth
;;  MSTAR - Mass of host star, MSun
;;  TRANSIT_DEPTH - Squared Radius of planet, in units of the host star radius squared
;;  RP_RSTAR - radius of planet, in units of star radius (alternate to TRANSIT_DEPTH 
;;             in the event that no transit has been detected, to allow for phase curves without
;;             transits at lower inclination). 
;;  AR_SEMIMAJ - semimajor axis of orbit, in units of the host star radius
;;  TEQ_P - Equilibrium temperature of the planet, Kelvins
;;  TEFF_STAR - Effective temperature of star, Kelvins
;;  SECONDARY_DEPTH - Depth of secondary transit.  If given, then TEQ_P and TEFF_STAR will be ignored.
;;  SECONDARY_LAMBDA - band at which the secondary depth was measured, or is desired.  Values can be 
;;                     'J', 'H', 'KS', 'KP', '3.6', '4.5', '5.8', '8.0'
;;  INCLINATION - orbital inclination (measured from axis of revolution), degrees
;;  MJD_TRANSIT - Modified Julian date of mid-transit (JD-2400000) (days)
;;  P_ORBIT - orbital period (days)
;;  EXODATA - The exoplanet database structure, returned from a previous run
;;  
;; OUTPUT KEYWORDS (output when EXOSYSTEM is input)
;;  RA - Right Ascension string (decimal hr, J2000, epoch 2000, from database)
;;  DEC - Declination string (decimal degrees, J2000, epoch 2000, from database)
;;  VMAG - V magnitude of the host star
;;  DISTANCE - distance of system, parsecs
;;  ECC - orbital eccentricity
;;  T14 - time between first and fourth contacts (total transit time)
;;  F36 - Estimated 3.6 micron flux density (mJy) of star, extrapolated from VMAG and TEFF_STAR
;;  F45 - Estimated 4.5 micron flux density (mJy) of star, extrapolated from VMAG and TEFF_STAR
;;  EXODATA - Passes back the exoplanet database as a structure
;;;
;;;; Systems known to have both transits and secondary eclipses
;WASP-1 b,WASP-18 b,WASP-33 b,XO-3 b,WASP-12 b,CoRoT-1 b,XO-2 b,55 Cnc e,HD 80606 b,WASP-19 b,
;OGLE-TR-113 b,GJ 436 b,WASP-14 b,WASP-24 b,WASP-17 b,XO-1 b,HD 149026 b,TrES-3 b,TrES-4 b,WASP-3 b,
;TrES-1 b,Kepler-12 b,TrES-2 b,KOI-13 b,Kepler-7 b,CoRoT-2 b,HAT-P-7 b,Kepler-6 b,Kepler-17 b,
;Kepler-5 b,HD 189733 b,WASP-2 b,HD 209458 b,HAT-P-8 b,HAT-P-1 b,WASP-4 b,HAT-P-6 b,XO-4 b

exoplanet_data_file = '/Users/jkrick/idlbin/exoplanets.csv'

get_exoplanet_data,EXOSYSTEM=exosystem,MSINI=msini,MSTAR=mstar,TRANSIT_DEPTH=transit_depth,RP_RSTAR=rp_rstar,AR_SEMIMAJ=ar_semimaj,$
                       TEQ_P=teq_p,TEFF_STAR=teff_star,SECONDARY_DEPTH=secondary_depth,SECONDARY_LAMBDA=4.5,$
                       INCLINATION=inclination,MJD_TRANSIT=mjd_transit,P_ORBIT=p_orbit,EXODATA=exodata,RA=ra,DEC=dec,VMAG=vmag,$
                       DISTANCE=distance,ECC=ecc,T14=t14,F36=f36,F45=f45,FP_FSTAR0=fp_fstar0,VERBOSE=verbose

IF N_ELEMENTS(mjd_start) EQ 0 THEN BEGIN
   mjd_start = get_exoplanet_start(p_orbit,mjd_transit,START_HR_BEFORE_TRANSIT=start_hr_before_transit,START_PHASE=start_phase)
ENDIF 

IF N_ELEMENTS(mjd_transit) EQ 0 THEN mjd_transit = mjd_start + p_orbit/4.
IF N_ELEMENTS(nperiod) EQ 0 THEN nperiod = 1.0
IF N_ELEMENTS(duration_hr) NE 0 THEN nperiod = duration_hr/24./p_orbit
IF N_ELEMENTS(nt) EQ 0 THEN nt = 1000L
 
t = DINDGEN(nt) / (nt-1) * nperiod * p_orbit + mjd_start
tdiff = (t-mjd_transit)/P_ORBIT MOD 1   ;;; number of orbits since nearest transit (can be negative)

phase = tdiff
ilow = WHERE(phase LT -0.5,nlow)
IF nlow NE 0 THEN phase[ilow] += 1
ihi = WHERE(phase GT 0.5,nhi)
IF nhi NE 0 THEN phase[ihi] -= 1
 
ineg = WHERE(tdiff LT 0 ,nneg)
IF nneg NE 0 THEN tdiff[ineg]  += 1   ;; "positivize" the number of orbits
phi_orbit = 2*!dPI * tdiff  ;; phi goes from 0 to 2*!pi
z = ar_semimaj * SQRT(1 - COS(phi_orbit)^2 * SIN(inclination/!radeg)^2)

p = rp_rstar
kappa_1 = ACOS( (1-p^2+z^2)/(2*z) )
kappa_2 = ACOS( (p^2+z^2-1)/(2*p*z) )
lambda_2 = (kappa_1 + p^2*kappa_2 - 0.5*SQRT( 4*z^2 - (1 + z^2 - p^2)^2 ))/!dpi

i1 = WHERE( (z GT (1+P)) OR (0.5*!dpi LT phi_orbit AND 1.5*!dpi GE phi_orbit),n1)     ;; Out of transit
i2 = WHERE( z GT (1-P) AND z LE (1+P) AND (phi_orbit GT 1.5*!dpi OR phi_orbit LE 0.5*!dpi),n2)   ;;; Partial transit
i3 = WHERE( z LE (1-P)  AND (phi_orbit GT 1.5*!dpi OR phi_orbit LE 0.5*!dpi),n3)  ;; Full Transit
lambda_star = DBLARR(nt)
IF n1 NE 0 THEN lambda_star[i1] = 0.0
IF n2 NE 0 THEN lambda_star[i2] = lambda_2[i2]
IF n3 NE 0 THEN lambda_star[i3] = p^2

i1 = WHERE( (z GT (1+P)) OR (phi_orbit GT 1.5*!dpi OR phi_orbit LE 0.5*!dpi),n1)   ;; Out of contact
i2 = WHERE( z GT (1-P) AND z LE (1+P) AND (0.5*!dpi LT phi_orbit AND 1.5*!dpi GE phi_orbit),n2)   ;;; Partial occultation
i3 = WHERE( z LE (1-P)  AND (0.5*!dpi LT phi_orbit AND 1.5*!dpi GE phi_orbit),n3)    ;;; full occultaiton
lambda_planet = DBLARR(nt)
IF n1 NE 0 THEN lambda_planet[i1] = 0.0
IF n2 NE 0 THEN lambda_planet[i2] = lambda_2[i2]/P^2
IF n3 NE 0 THEN lambda_planet[i3] = 1.0

alpha = ABS(phi_orbit-!dpi)
phase_curve = (sin(alpha) + (!dpi-alpha)*cos(alpha))/!dpi  ;; Phase curve for a Lambert sphere
rel_flux = (1-lambda_star) + fp_fstar0 * (1-lambda_planet) * phase_curve
IF KEYWORD_SET(TPHASE) THEN BEGIN
   is = sort(phase)
   phase = phase[is]
   rel_flux = rel_flux[is]
   torbit = phase
   t = phase 
ENDIF ELSE torbit = (t-mjd_start)
   
IF KEYWORD_SET(PLOT) OR KEYWORD_SET(PDF) THEN BEGIN
   IF N_ELEMENTS(exosystem) NE 0 THEN BEGIN
      exoplan_string = STRJOIN(STRSPLIT(exosystem,' ',/EXTRACT),'_',/SINGLE) 
      plot_title = exoplan_string
   ENDIF ELSE BEGIN
      exoplan_string = 'exosystem'
      plot_title = ''
   ENDELSE
   IF N_ELEMENTS(LABEL) EQ 0 THEN lab='' ELSE  lab = '_'+label
IF KEYWORD_SET(tphase) then xtitlename = 'Phase' else xtitlename = 'Time (dy)'
   lc = plot(torbit,rel_flux,title=plot_title,xtitle=xtitlename,ytitle='F/F*', thick =4 ,color='Sky Blue')
;   psfile = exoplan_string+lab+'_predicted_lightcurve.eps'
;   ifps,psfile,12,6,/ENCAPSULATED,NOPOST=~KEYWORD_SET(PDF),/COLOR
;      CGPLOT,torbit,rel_flux,title=plot_title,xtitle='Time (dy)',ytitle='F/F*',xstyle=1,ystyle=16,color='Sky Blue',thick=4,_extra=extra;,psym=1
;   ENDPS,psfile,PDF=pdf,NOPOST=~KEYWORD_SET(PDF)
ENDIF

RETURN
END

 

   
