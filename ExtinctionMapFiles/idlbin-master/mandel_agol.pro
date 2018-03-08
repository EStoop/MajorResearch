function mandel_Agol,params0, t = t, flux = flux, err = err, p_orbit =p_orbit, mjd_start =mjd_start, mjd_transit =mjd_transit, savefilename = savefilename, overplot = overplot ;, rp_rstar = rp_rstar , ar_semimaj = ar_semimaj,inclination = inclination,

;; fixing to values from the literature
;p_orbit
;mjd_start  ; XXX do I really need this if using my own phase as input anyway?
;mjd_transit

;; fitting
;fp_fstar0
;ar_semimaj
;inclination
;rp_rstar


;taken from Jim;s exoplanet_light_curve.pro
;nt = 1000.
;nperiod = 1.
;t = DINDGEN(nt) / (nt-1) * nperiod * p_orbit + mjd_start
;tdiff = (t-mjd_transit)/P_ORBIT MOD 1   ;;; number of orbits since nearest transit (can be negative)

;phase = tdiff
;ilow = WHERE(phase LT -0.5,nlow)
;IF nlow NE 0 THEN phase[ilow] += 1
;ihi = WHERE(phase GT 0.5,nhi)
;IF nhi NE 0 THEN phase[ihi] -= 1
phase = t
tdiff = t
nt = n_elements(phase)

ineg = WHERE(tdiff LT 0 ,nneg)
IF nneg NE 0 THEN tdiff[ineg]  += 1   ;; "positivize" the number of orbits
phi_orbit = 2*!dPI * tdiff  ;; phi goes from 0 to 2*!pi
z = params0[2] * SQRT(1 - COS(phi_orbit)^2 * SIN(params0[3]/!radeg)^2)

p = params0[1]; rp_rstar
kappa_1 = ACOS( (1-p^2+z^2)/(2*z) )
kappa_2 = ACOS( (p^2+z^2-1)/(2*p*z) )
lambda_2 = (kappa_1 + p^2*kappa_2 - 0.5*SQRT( 4*z^2 - (1 + z^2 - p^2)^2 ))/!dpi

i1 = WHERE( (z GT (1+P)) OR (0.5*!dpi LT phi_orbit AND 1.5*!dpi GE phi_orbit),n1)     ;; Out of transit
i2 = WHERE( z GT (1-P) AND z LE (1+P) AND (phi_orbit GT 1.5*!dpi OR phi_orbit LE 0.5*!dpi),n2)   ;;; Partial transit
i3 = WHERE( z LE (1-P)  AND (phi_orbit GT 1.5*!dpi OR phi_orbit LE 0.5*!dpi),n3)  ;; Full Transit
lambda_star = DBLARR(nt)
;IF n1 NE 0 THEN lambda_star[i1] = params0[6]*phi_orbit ;0.0
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

;now add in a phase curve
alpha = ABS(phi_orbit-!dpi)

;phase_curve = (sin(alpha) + (!dpi-alpha)*cos(alpha))/!dpi  ;; Phase curve for a Lambert sphere from Seager book 

phase_curve = 1+ params0[0] + params0[4]*cos(2*!dpi*tdiff/p_orbit) + params0[5]*sin(2*!dpi*tdiff/p_orbit);+ params0[6]*cos(4*!dpi*tdiff/p_orbit) + params0[7]*sin(4*!dpi*tdiff/p_orbit)

;testing
;phase_curve = params0[4]*sin(alpha + params0[5])

rel_flux = (1-lambda_star) + params0[0] * (1-lambda_planet) ;* phase_curve
;rel_flux = (1-lambda_star) + fp_fstar0 * (1-lambda_planet) * phase_curve


model = rel_flux
model = (flux - model) / err


if keyword_set(overplot) then begin
   test = plot(phase, rel_flux, color = 'Sky Blue', thick = 4, /overplot )
   openw, outlun, '/Users/jkrick/irac_warm/pcrs_planets/WASP-52b/jk_model.txt',/GET_LUN
   for nm = 0, n_elements(phase) - 1 do printf, outlun, phase(nm), rel_flux(nm)
   free_lun, outlun
   
endif

;change 'phot' to 'model'
;print, 'at end of  mandel', savefilename
save, phase, rel_flux, filename = savefilename
return, model
end
