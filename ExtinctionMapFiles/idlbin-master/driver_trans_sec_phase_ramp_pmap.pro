function calc_trans_sec_phase, p, phase_func=phase_func

COMMON data, bjd_tot, flux_tot, nbr_ind, gw, err_tot, nonlin_c

xt=p(0:8)
xph=p(9:12)

trans=transit_phase_orb_nl(bjd_tot, xt, xph, phase_func, nonlin_c)

flux_mod=trans*p(13)

wf1=flux_tot/flux_mod
res=(wf1-1.d0)/err_tot

return, res

end

pro final_answer, p, perror, phase_func, sfile

COMMON data, bjd_tot, flux_tot, nbr_ind, gw, err_tot, nonlin_c

xt=p(0:8)
xph=p(9:12)

trans=transit_phase_orb_nl(bjd_tot, xt, xph, phase_func, nonlin_c)

flux_mod=trans*p(13)


save, p, perror, bjd_tot, flux_tot, flux_mod, trans, filename=sfile

end


pro driver_trans_sec_phase_ramp_pmap, channel, phase_func, sfile

COMMON data, bjd_tot, flux_tot, nbr_ind, gw, err_tot, nonlin_c

if (channel eq 1) then begin
restore,'ch1_filter_nbr_np_centr_apr6.sav'   ;np apr
nonlin_c=[0.495904,-0.479051,0.490810,-0.201760]
Tc=575.93243d0
p0_test=[3.52474859d0, 86.710d0*!dtor, 8.760d0, 0.d0, 180.d0*!dtor, Tc, sqrt(0.014607d0), 0.00094d0, 0.00094d0]
endif

if (channel eq 2) then begin
restore,'ch2_filter_nbr_nonp_centr_apr3.sav'
nonlin_c=[0.544494,-0.740935,0.771975,-0.295253]
Tc=216.40808d0
p0_test=[3.52474859d0, 86.710d0*!dtor, 8.760d0, 0.d0, 180.d0*!dtor, Tc, sqrt(0.014607d0), 0.00213d0, 0.00213d0]
endif


p0_phase=[0.d0, 0.d0, 0.d0, 0.d0]

;values for orbital fit
;p(0) -> orbital period (days)
;p(1) -> inc = inclination angle (radians)
;p(2) -> a/R_* = semi-major axis divided by R_*
;p(3) -> e (eccentricity)
;p(4) -> omega (radians)
;p(5) -> center of transit
;p(6) -> (R_p/R_*) (radius of planet in units of radius of star)
;p(7) -> depth of 1st secondary eclipse 
;p(8) -> depth of 2nd secondary eclipse
;p(9:12) -> phase curve parameters
;p(13) -> stellar flux
pa=replicate({value:0.d0, fixed:0, limited:[0,0], limits:[0.d0, 0.d0], mpmaxstep:0.d0},14)
;starting values from Pal et al. (2010)
pa(*).value=[p0_test, p0_phase, median(flux_tot)]

pa(*).limited(*)=[1,1]
pa(0).limits(*)=[3.52474d0,3.52475d0]
pa(1).limits(*)=[85.d0,90.d0]*!dtor
pa(2).limits(*)=[8.5d0, 9.d0]
pa(3).limits(*)=[0.0,1.0]
pa(4).limits(*)=[-!dpi, !dpi]
pa(5).limits(*)=Tc+[-1.d0,1.d0]*0.02d0
pa(6).limits(*)=[0.1d0,0.2d0] 
pa([7,8]).limits(*)=[0.0d0,0.005d0]
pa(9:12).limits(*)=[-0.005d0,0.005d0]
pa(13).limits(*)=minmax(flux_tot)

pa([0,1,2,3,4]).fixed=1

fa={phase_func:phase_func} 

p = mpfit('calc_trans_sec_phase', parinfo=pa, functargs=fa, perror=perror, status=status, bestnorm=bestnorm, ERRMSG=ERRMSG, maxiter=25, /fastnorm)

print, status

final_answer, p, perror, phase_func, sfile

end






