

pro driver_trans_sec_phase_ramp_fixed, channel, phase_func, ramp_flag, cowan_opt, sfile, finalp, filename

COMMON data, bjd_tot, flux_tot, nbr_ind, gw, err_tot, time_tot, nonlin_c
print, 'running fixed'
restore,filename

;;restore,'/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/fitting_input_phot_ch2_2.25000_150723.sav'
;;flux_tot=flux_tot+0.002d0
;nonlin_c=[0.0225, 0.3828, -0.2748, -0.0522] ;wasp-14b specific
nonlin_c=[0.523357, -0.74367, 0.801920, -0.316680]
Tc=56042.687d
p0_test=[2.24376507d0, 84.63d0*!dtor, 5.98d0, -.0247236, -0.0792322, Tc,  0.09421d0, 0.002115d0, 0.002367d0] ; for inc *!dtor


p0_ramp=[0.00016d0, 0.07d0]
;initial fit from models
p0_phase=[[0.00265d0, 1.0d0, 0.5d0, 0.5d0, 0.0001d0],$ ;flat
          [0.000575d0, 0.69d0, 0.2500d0, 0.2500d0, 0.0004096d0],$ ;cosf
          [0.001184d0, 0.4114d0, 0.88d0, 0.2500d0, 6.962d-05],$ ;lorentz_sym
          [0.001178d0, 0.1654d0, 0.438d0, 1.13d0, 6.403d-05],$ ;lorentz_asym
          [0.0001988d0, 0.6602d0, 0.2500d0, 0.2500d0, 0.0001087d0],$  ;dist_sq
          [0.0003122d0, 0.0002354d0, 1.352d-06, 9.158d-05, 0.0004129d0]]  ;cowan

ph_opt=['flat','cosf','lorentz_sym', 'lorentz_asym', 'dist_sq','cowan']
pind=where(ph_opt eq phase_func)
p0_cowan=[0.d0, 0.d0, 0.d0, 0.d0]

;values for orbital fit
;p(0) -> orbital period (days)
;p(1) -> inc
;p(2) -> a/R_*
;p(3) -> k=e cos (w) (eccentricity * sin (longitude of pericenter))
;p(4) -> h=e sin (w) (eccentricity * sin (longitude of pericenter))
;p(5) -> mid-point of transit
;p(6)  -> (R_p/R_*) (radius of planet in units of radius of star)
;p(7) -> depth of 1st secondary eclipse 
;p(8) -> depth of 1st secondary eclipse
;p(9:13) -> phase curve parameters
;p(14:15) -> ramp parameters
;;pa=replicate({value:0.d0, fixed:0, limited:[0,0], limits:[0.d0, 0.d0], mpmaxstep:0.d0},16)
pa=replicate({value:0.d0, fixed:0, limited:[0,0], limits:[0.d0, 0.d0], mpmaxstep:0.d0},16)
;fix starting values to the output from the last run
pa(*).value=finalp ;[p0_test, p0_phase(*,pind),p0_ramp]

pa([0,1,2,3,4,5, 6, 7, 9, 10,11,12,13,14]).fixed=1 ;, 1, 2, 3, 4
;apparently can't fix all parameters, need to allow one to change
;find one that doesn't change the amplitude and phase shift parameters
pa(15).limits(*)=[0.01d0, 1.0d0]

fa={phase_func:phase_func}

p = mpfit('calc_trans_sec_phase', parinfo=pa, functargs=fa, perror=perror, status=status,  DOF=adof,bestnorm=bestnorm, ERRMSG=ERRMSG, /fastnorm, maxiter = 0)
;print, 'perror', perror
print, 'status', status
print, errmsg
print, 'reduced chi squared',  bestnorm / adof
print, 'pa', pa.fixed
porb=p(0)
ecc=sqrt(p(3)^2+p(4)^2)
om=atan(p(4),p(3))
;cirrange, om,/radians
;omega=om-!dpi
omega=om
a_Rs=p(2)
inc=p(1)

print, [transpose(p), transpose(perror)]

final_answer, p, perror, bestnorm, phase_func, sfile


end






