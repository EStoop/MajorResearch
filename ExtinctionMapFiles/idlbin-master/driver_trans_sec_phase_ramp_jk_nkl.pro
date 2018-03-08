function ekepler_nl,ms,e
eps=1.d-10
; This routine solves Kepler's equation for E as a function of (e,M)
; as laid out by Marc A. Murison (USNO).  Faster and more accurate 
; than method in Murray & Dermott, because of better starting value 
; for E0 
t34=e^2
t35=e*t34
t33=cos(ms)
E0=ms+(-0.5*t35+e+(t34+1.5*t33*t35)*t33)*sin(ms)
eps3=eps+1.0
while(max(abs(eps3)) gt eps) do  begin
  ;eps1=(e0-e*sin(e0)-ms)/(1.0-e*cos(e0))
  ;t1=-1.0+e*cos(e0)
  ;t2=e*sin(e0)
  ;t3=-e0+t2+ms
  ;eps2=t3/(0.5*t3*t2/t1+t1)
  t1=cos(E0)
  t2=-1.0+e*t1
  t3=sin(E0)
  t4=e*t3
  t5=-E0+t4+ms
  t6=t5/(0.5*t5*t4/t2+t2)
  eps3=t5/((0.5*t3-(1.0/6.0)*t1*t6)*e*t6+t2)
  E0=E0-eps3
endwhile

return,E0
end

pro occultunifast,z,w,muo1
if(abs(w-0.5d0) lt 1.d-3) then w=0.5d0
; This routine computes the lightcurve for occultation
; of a uniform source without microlensing  (Mandel & Agol 2002).
;Input:
;
; rs   radius of the source (set to unity)
; z    impact parameter in units of rs
; w    occulting star size in units of rs
;
;Output:
; muo1 fraction of flux at each b0 for a uniform source
;
; Now, compute pure occultation curve:
nz=n_elements(z)
muo1=dblarr(nz)
; the source is unocculted:
; Table 3, I.
indx=where(z gt 1.d0+w)
if(total(indx) ge 0) then muo1(indx)=1.d0
; the  source is completely occulted:
; Table 3, II.
; (already zero, so do nothing)
; the source is partly occulted and the occulting object crosses the limb:
; Equation (26):
indx=where(z ge abs(1.d0-w) and z le 1.d0+w)
if(total(indx) ge 0) then begin
 zt=z(indx)
 xt=(1.d0-w^2+zt^2)/2.d0/zt
 kap1=acos(xt*double(xt lt 1.d0)+1.d0*double(xt ge 1.d0))
 xt=(w^2+zt^2-1.d0)/2.d0/w/zt
 kap0=acos(xt*double(xt lt 1.d0)+1.d0*double(xt ge 1.d0))
 lambdae=w^2*kap0+kap1
 xt=4.d0*zt^2-(1.d0+zt^2-w^2)^2
 lambdae=(lambdae-0.5d0*sqrt(xt*double(xt ge 0.d0)))/!dpi
 muo1(indx)=1.d0-lambdae
endif
indx=where(z le 1.d0-w)
if(total(indx) ge 0) then muo1(indx)=1.d0-w(indx)^2
return
end

pro occultarb_fn,b,p,c,lc,lcuni,fname
; If you make use of this routine, please cite Mandel & Agol (2002).
;
; Computes a light curve for a limb-darkening model
; which is determined by an arbitrary function in
; the routine fname.  The function iname should
; take as input the mu=cos(theta) values for the
; star and output the specific intensity at these values;
; it should work for all values of 0 <= mu <= 1.
;
; Input:
;   b    vector of impact parameters as a function of time
;        (distance between center of two bodies in units
;         of the radius of the occulted body)
;   c    limb-darkening coefficients
;   p    ratio of radii of the foreground to background body
;
; fname  Name of function that computes I(mu)
; 
; Output:
;  lc   Lightcurve that is normalized to unity outside of transit
;
; Compute the limb-darkening at the edge of the star:
i0=call_function(fname,0d0,c)
; Compute light-curve for a uniform source with a 
; surface brightness equal to that at the edge of the star:
occultunifast,b,p,lcuni

; Now divide star up into nr annuli with radius ranging
; from r_n to r_n+1 in each annulus, and average radius
; of r_(n+1/2).  Then, compute the light curve:

; Calculate maximum depth for a uniform source:
fac=max(abs(1d0-lcuni))
nr=8L
; Compute light curve to a fractional precision of depth of ~10^-3:
dlcmax=1d0
; Only compute the limb-darkened light curve where necessary:
indx=where(lcuni ne 1d0 and lcuni ne 0d0)
nindx=n_elements(indx)
lc1=lcuni[indx]
lc=lcuni
rn  =[0d0,1d0]
imun=[1d0, i0]
while(dlcmax gt fac*.5d-3) do begin
; First solve for points that are evenly spaced in theta:
 dtheta = .5d0*!dpi/double(nr)
 thetan = dindgen(nr+1)*dtheta
 rn     = sin(thetan)
 mun    = cos(thetan)
; Compute radius at center of each annulus:
 rnhalf=sin(.5d0*(thetan[1:nr]+thetan[0:nr-1]))
; Compute mu at the boundary of each annulus:
 mun=sqrt(1d0-rn^2)
; Compute intensities at boundaries:
 imun=call_function(fname,mun,c)
; Compute the flux as a function of impact parameter:
 fofb=lcuni[indx]
; Now compute the light curve for each annulus:
 dimunh=(imun[1L:nr]-imun[0L:nr-1L])/i0
; Compute the unocculted flux:
 f0=1d0-total(rnhalf^2*dimunh)
 for ir=0L,nr-1L do begin
   occultunifast,b[indx]/rnhalf[ir],p/rnhalf[ir],lcnh
   fofb=fofb-rnhalf[ir]^2*dimunh[ir]*lcnh
 endfor
 lc[indx]=fofb/f0
 dlcmax=max(abs(lc[indx]-lc1))
 lc1=lc[indx]
; Double the resolution until convergence:
 nr=nr*2L
endwhile
return
end

function iofmu,mu,c
icalc=1.d0
;four-parameter nonliner limb dk
for j=0,3 do icalc=icalc-c(j)*(1.d0-mu^((j+1.d0)/2.d0))
return,icalc
;return,1d0-c[0]*(1d0-mu)+c[1]*alog((mu+c[2])/(1d0+c[2]))
end

function flat, x, p
phase=x*0.d0+p(0)

return,phase

end

function cowan, x, p

;x is true anomaly, or f+omega, or 
;f+offset

phase=p(0)*cos(x)+p(1)*sin(x)+ $
           p(2)*cos(2.d0*x)+p(3)*sin(2.d0*x)

;phase=p(4)+p(0)*cos(x)

return,phase

end

function cosf, x, p

phase=p(4)+p(0)*cos(x-p(1))

return, phase

end


function cosf_sq, x, p

;phase=1.d0+p(0)*(1.d0+p(2)*cos(x-p(1)))^2
phase=p(0)*cos(x-p(1))^2

return, phase

end

function lorentz_sym, x0, p, porb
;symmetric lorentzian

del=[-1.d0, 0.d0, 1.d0, 2.d0]
peaks=x0*0.d0
for n=0,3 do begin
   x=x0+del(n)*porb
   u = (x-p[1])/p[2]
   peaks=peaks+p[0]/(u^2+1.d0)
endfor


phase=p[4]+peaks

return, phase

end

function lorentz_asym, x0, p, porb
;asymmetric lorentzian

del=[0.d0, -1.d0, 1.d0, 2.d0]
peaks=x0*0.d0
for n=0,3 do begin
   x=x0+del(n)*porb
   u1 = (x - p[1])/p[2]
   step1=x*0.d0
   ind1=where(u1 lt 0.d0, c1)
   if (c1 gt 0.0) then step1(ind1)=1.d0

   u2 = (x - p[1])/p[3]
   step2=x*0.d0
   ind2=where(u2 ge 0.d0, c2)
   if (c2 gt 0.0) then step2(ind2)=1.d0

   peaks=peaks+p[0]*(step1/(u1^2+1.d0)+step2/(u2^2+1.d0))
endfor


phase=p[4]+peaks

return, phase

end


function dist_sq, x, p

;dist=(1d0-e1^2)/(1d0+e1*cos(x-p(1)))
;peak=p(0)*(1.d0/dist^2)
;peak=p(0)*(1.d0+e1*cos(x-p(1)))^2
peak=p(0)*(1.d0+cos(x-p(1)))^2
;peak=peak-min(peak)
phase=p[4]+peak

return, phase

end



function star, t, p

flux=p(0)*t+p(1)*t^2

return,flux

end

function star_sine, t, p

prot=3.95d0
flux=p(0)*sin(2.d0*!dpi/prot*t-p(1))

return,flux

end

function ramp_func, t, p

;ramp=1.d0-exp(-p(0)*t+p(1))
ramp=1.d0-p(0)*exp(-t/p(1))

return,ramp

end

pro transit_phase_orb, xt, x_ph, phase_func, trans, orb_params, ph_params

COMMON data, bjd_tot, flux_tot, nbr_ind, gw, err_tot, time_tot, nonlin_c

;;ON_ERROR, 2

e1=xt(3)
;transit information
f_trans=0.5d0*!dpi-xt(4)
f_trans=f_trans
E_trans=2.d0*atan(sqrt(1d0-e1^2)*tan(.5d0*f_trans),1d0+e1)
M_trans=E_trans-e1*sin(E_trans)
;M_trans=M_trans+2.d0*!dpi
;eclipse information
f_sec=-1.d0*(xt(4)+0.5d0*!dpi)
E_sec=2.d0*atan(sqrt(1d0-e1^2)*tan(.5d0*f_sec),1d0+e1)-2.d0*!dpi
M_sec=E_sec-e1*sin(E_sec)
tp=xt(5)-xt(0)/(2.d0*!dpi)*M_trans
t_sec1=xt(5)-xt(0)/(2.d0*!dpi)*(M_trans-M_sec)
t_sec2=t_sec1+xt(0)
phi_sec=(M_sec-M_trans)/(2.d0*!dpi)

m=2d0*!dpi/xt(0)*(bjd_tot-tp)
ekep=ekepler_nl(m,e1)
f=2.d0*atan(sqrt((1.d0+e1)/(1.d0-e1))*tan(0.5d0*ekep))
nm0=where(m eq 0.d0)
if(total(nm0) ge 0) then f(nm0)=0.d0

radius=xt(2)*(1d0-xt(3)^2)/(1d0+xt(3)*cos(f))

; Now compute sky separation:
z0=radius*sqrt(1d0-(sin(xt(1))*sin(xt(4)+f))^2)

;get transit curve with limb darkening
occultarb_fn,z0,xt(6),nonlin_c,flux_trans,flux_sec,'iofmu'

if (strmid(phase_func,0,4) eq 'cosf') then phase=call_function(phase_func, f, x_ph)
if (strmid(phase_func,0,7) eq 'lorentz') then phase=call_function(phase_func, bjd_tot-tp, x_ph, xt(0)) 
if (phase_func eq 'flat') then phase=call_function(phase_func, f, x_ph)
if (strmid(phase_func,0,4) eq 'dist') then phase=call_function(phase_func, f, x_ph)
if (phase_func eq 'cowan') then phase=cowan(f, x_ph(0:4))
if (phase_func eq 'cowan2') then phase=cowan(f+xt(4)+!dpi, x_ph(0:4))

print, 'tsec1, tsec2, xt(5)', t_sec1, t_sec2, xt(5)

ind_trans=where(abs(bjd_tot-xt(5)) lt 0.25 and flux_trans lt 1.d0, c0)
ind_sec1=where(abs(bjd_tot-t_sec1) lt 0.25 and flux_sec lt 1.d0, c1)
ind_sec2=where(abs(bjd_tot-t_sec2) lt 0.25 and flux_sec lt 1.d0, c2)
;;ind_sec2=where( flux_sec lt 1.d0, c2)
;;ind_sec2=where(abs(bjd_tot-t_sec2) lt 0.25 , c2)
;;print, 'testing c', c0, c1, c2
flux_trans([ind_sec1,ind_sec2])=1.d0
flux_sec(ind_trans)=1.d0


;transit duration, ingress, egress
t14_trans=bjd_tot(ind_trans(c0-1.0))-bjd_tot(ind_trans(0))
ingress=where(z0(ind_trans) ge abs(1.d0-xt(6)) and z0(ind_trans) le 1.d0+xt(6) and bjd_tot(ind_trans)-xt(5) lt 0.d0, ci)
egress=where(z0(ind_trans) ge abs(1.d0-xt(6)) and z0(ind_trans) le 1.d0+xt(6) and bjd_tot(ind_trans)-xt(5) gt 0.d0, ce)
t12_trans=bjd_tot(ind_trans(ingress(ci-1.0)))-bjd_tot(ind_trans(ingress(0)))
t34_trans=bjd_tot(ind_trans(egress(ce-1.0)))-bjd_tot(ind_trans(egress(0)))

;eclipse 1 duration, ingress, egress
t14_sec1=(bjd_tot(ind_sec1(c1-1.0))-bjd_tot(ind_sec1(0)))
ingress=where(z0(ind_sec1) ge abs(1.d0-xt(6)) and z0(ind_sec1) le 1.d0+xt(6) and bjd_tot(ind_sec1)-t_sec1 lt 0.d0, ci)
egress=where(z0(ind_sec1) ge abs(1.d0-xt(6)) and z0(ind_sec1) le 1.d0+xt(6) and bjd_tot(ind_sec1)-t_sec1 gt 0.d0, ce)
t12_sec1=bjd_tot(ind_sec1(ingress(ci-1.0)))-bjd_tot(ind_sec1(ingress(0)))
t34_sec1=bjd_tot(ind_sec1(egress(ce-1.0)))-bjd_tot(ind_sec1(egress(0)))

;eclipse 2 duration, ingress, egress
t14_sec2=(bjd_tot(ind_sec2(c2-1.0))-bjd_tot(ind_sec2(0)))
ingress=where(z0(ind_sec2) ge abs(1.d0-xt(6)) and z0(ind_sec2) le 1.d0+xt(6) and bjd_tot(ind_sec2)-t_sec2 lt 0.d0, ci)
egress=where(z0(ind_sec2) ge abs(1.d0-xt(6)) and z0(ind_sec2) le 1.d0+xt(6) and bjd_tot(ind_sec2)-t_sec2 gt 0.d0, ce)
t12_sec2=bjd_tot(ind_sec2(ingress(ci-1.0)))-bjd_tot(ind_sec2(ingress(0)))
t34_sec2=bjd_tot(ind_sec2(egress(ce-1.0)))-bjd_tot(ind_sec2(egress(0)))


flux_sec_new=flux_sec
flux_sec_new(ind_sec1)=1.d0+(flux_sec(ind_sec1)-1.d0)*xt(7)/xt(6)^2
flux_sec_new(ind_sec2)=1.d0+(flux_sec(ind_sec2)-1.d0)*xt(8)/xt(6)^2

ph1_sec1=phase(ind_sec1(0)-1.0)
ph4_sec1=phase(ind_sec1(c1-1.0)+1.0)
ph0_sec1=(ph1_sec1+ph4_sec1)/2.d0

ph1_sec2=phase(ind_sec2(0)-1.0)
print, 'phase', n_elements(phase), ind_sec2(c2-1.0)+1.0
subvar = ind_sec2(c2-1.0)+1.0
if ind_sec2(c2-1.0)+1.0 ge n_elements(phase) then subvar = n_elements(phase) - 1
ph4_sec2=phase(subvar)
ph0_sec2=(ph1_sec2+ph4_sec2)/2.d0

trans2=flux_sec_new-1.d0
ph_sec=phase
ph_sec(ind_sec1)=(trans2(ind_sec1)+xt(7))*((phase(ind_sec1)-ph0_sec1)/xt(7)+1.d0)-(xt(7)-ph0_sec1)
ph_sec(ind_sec2)=(trans2(ind_sec2)+xt(8))*((phase(ind_sec2)-ph0_sec2)/xt(8)+1.d0)-(xt(8)-ph0_sec2)
;offset=min([ph_sec(ind_sec1), ph_sec(ind_sec2)])
offset=min(ph_sec(ind_sec2))

ph_sec=ph_sec-offset
ind=where(ph_sec lt 0.d0, count)
if (count gt 0.0) then ph_sec(ind)=0.d0

trans=ph_sec+flux_trans

ph_min=min(phase-offset, mind)
if (ph_min) lt 0.d0 then ph_min=0.d0
tph_min=bjd_tot(mind)

ph_max=max(phase-offset, mxind)
tph_max=bjd_tot(mxind)

orb_params=[tp, phi_sec, t_sec1, t_sec2, t14_trans, t12_trans, t34_trans, t14_sec1, t12_sec1, t34_sec1, t14_sec2, t12_sec2, t34_sec2]
ph_params=[ph_min, tph_min, ph_max, tph_max]

;;jk added
;;test plot
;help, bjd_tot
;help, phase8
;print, 'phase', phase[0:10], min(phase), max(phase)
;print, 'ending transit phase orb'
;print, 'shouldnt see this'

end

function calc_trans_sec_phase, p, phase_func=phase_func
COMMON data, bjd_tot, flux_tot, nbr_ind, gw, err_tot, time_tot, nonlin_c
porb=p(0)
ecc=sqrt(p(3)^2+p(4)^2)
om=atan(p(4),p(3))
omega=om
;ecc=0.0830d0
;omega=252.67*!dtor

;a_Rs=p(2)*(porb/(2.d0*!dpi))*sqrt(1.d0-p(1))*sqrt(1.d0-ecc^2)/(1.d0+p(4))
;cosi=(sqrt(p(1))/a_Rs)*(1.d0+p(4))/(1.d0-ecc^2)
;inc=acos(cosi)
a_Rs=p(2)
inc=p(1)

xt=[porb, inc, a_Rs, ecc, omega, p(5:8)] ;parameters relevant to transit
xph=p(9:13)                              ;parameters relevant to phase
xr=p(14:15)                              ;parameters relevant to ramp

;trans=transit_phase_orb_nl(bjd_tot, xt, xph, phase_func, nonlin_c)
transit_phase_orb, xt, xph, phase_func, trans, orb_params, ph_params

flux1=flux_tot/(trans)
;;nl does this to apply gaussian weights
;;w1=total(flux1(nbr_ind)*gw,2)
w1 = 1.d0
wf1=flux1/w1
res=(wf1-1.d0)/err_tot


;;testing different code
;;xph=p(9:12)

;;trans=transit_phase_orb_nl(bjd_tot, xt, xph, phase_func);, nonlin_c)

;;flux_mod=trans;*p(13)

;;wf1=flux_tot/flux_mod
;;res=(wf1-1.d0)/err_tot


return, res

end

pro final_answer, p, perror, bestnorm, phase_func, sfile

COMMON data, bjd_tot, flux_tot, nbr_ind, gw, err, time_tot, nonlin_c

porb=p(0)
ecc=sqrt(p(3)^2+p(4)^2)
om=atan(p(4),p(3))
omega=om
;a_Rs=p(2)*(porb/(2.d0*!dpi))*sqrt(1.d0-p(1))*sqrt(1.d0-ecc^2)/(1.d0+p(4))
;cosi=(sqrt(p(1))/a_Rs)*(1.d0+p(4))/(1.d0-ecc^2)
;inc=acos(cosi)
a_Rs=p(2)
inc=p(1)

xt=[porb, inc, a_Rs, ecc, omega, p(5:8)] 
xph=p(9:13)
xr=p(14:15)

;;trans=transit_phase_orb_nl(bjd_tot, xt, xph, phase_func);, nonlin_c)
transit_phase_orb, xt, xph, phase_func, trans, orb_params, ph_params

plot, bjd_tot-p(5), flux_tot-1.d0, psym=3
oplot, bjd_tot-p(5), trans-1.d0, thick=3

flux1=flux_tot/(trans)
;;w1=total(flux1(nbr_ind)*gw,2)
;;jk making this up
w1 = 1.0
;plot, bjd_tot-xt(5), trans-1.0

save, p, perror, bestnorm, xt, xph, xr, bjd_tot, time_tot, w1, flux_tot, trans, orb_params, ph_params, phase_func, filename=sfile
;ph_params=[ph_min, tph_min, ph_max, tph_max]

print,'ph_params', ph_params
print, 'amplitude',( ph_params(2) - ph_params(0))/2.
print, 'phase shift', (p(5) - ph_params(1))*(1./p(0))*360.
end


pro driver_trans_sec_phase_ramp_jk_nkl, channel, phase_func, ramp_flag, cowan_opt, sfile, filename

COMMON data, bjd_tot, flux_tot, nbr_ind, gw, err_tot, time_tot, nonlin_c
;;restore,'/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/fitting_input_phot_ch2_2.25000_150723.sav'
restore,filename; '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/fitting_input_phot_wong.sav'

;;flux_tot=flux_tot+0.002d0
;nonlin_c=[0.0225, 0.3828, -0.2748, -0.0522] ;wasp-14b specific
nonlin_c=[0.523357, -0.74367, 0.801920, -0.316680]
Tc=56042.687d0;56042.65d0
p0_test=[2.24376507d0, 84.63d0*!dtor, 5.98d0, -.0247236, -0.0792322, Tc,  0.09421d0, 0.002115d0, 0.002367d0] ; for inc *!dtor


p0_ramp=[0.00016d0, 0.07d0]
;initial fit from models
p0_phase=[[0.00265d0, 1.0d0, 0.5d0, 0.5d0, 0.0001d0],$ ;flat
          [0.000575d0, 0.69d0, 0.2500d0, 0.2500d0, 0.0004096d0],$ ;cosf
          [0.001184d0, 0.4114d0, 0.88d0, 0.2500d0, 6.962d-05],$ ;lorentz_sym
          [0.001178d0, 0.1654d0, 0.438d0, 1.13d0, 6.403d-05],$ ;lorentz_asym
          [0.0001988d0, 0.6602d0, 0.2500d0, 0.2500d0, 0.0001087d0],$  ;dist_sq
          [0.00076689572d0, 0.00011690945d0, 2.3513864d-05,2.4894987d-05, 0.00041290000d0]]  ;cowan

ph_opt=['flat','cosf','lorentz_sym', 'lorentz_asym', 'dist_sq','cowan']
pind=where(ph_opt eq phase_func)
p0_cowan=[0.00076689572d0, 0.00011690945d0, 2.3513864d-05,2.4894987d-05];, 0.00041290000d0]; [0.d0, 0.d0, 0.d0, 0.d0]

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
;starting values from Pal et al. (2010)
pa(*).value=[p0_test, p0_phase(*,pind),p0_ramp]
pa(*).limited(*)=[1,1]
pa(0).limits(*)=[2.d0,2.4d0]
pa(1).limits(*)=[70.d0, 90.d0]*!dtor
pa(2).limits(*)=[5.0d0, 7.0d0]
pa(3).limits(*)=[-0.1d0,0.1d0]
pa(4).limits(*)=[-0.1d0, 0.0d0]
pa(5).limits(*)=Tc+[-1.d0,1.d0]* 0.05d0
pa(6).limits(*)=[0.01d0,0.15d0] 
pa([7,8]).limits(*)=[0.001d0,0.007d0]
pa(9).limits(*)=[0.d0,0.011d0]
pa(10).limits(*)=[0.d0,1.8d0] 
pa(11:12).limits(*)=[0.04d0,6.0d0]
pa(13).limits(*)=[-0.0015d0,0.0015d0]
;pa(14).limits(*)=[0.d0, 0.01d0]
pa(14).limits(*)=[-0.01d0, 0.01d0]
pa(15).limits(*)=[0.01d0, 1.0d0]

if (strmid(phase_func,0,5) eq 'cowan') then begin
   pa(9:13).limits(*)=[-0.0015d0,0.0015d0]
   pa(9:12).value=p0_cowan
endif

if (strmid(phase_func,0,4) eq 'cosf') then pa(11:12).fixed=1
if (phase_func eq 'lorentz_sym') then pa(12).fixed=1
if (phase_func eq 'flat') then pa(9:12).fixed=1
;;if (phase_func eq 'flat') then pa(10:13).fixed=1
if (phase_func eq 'dist_sq') then pa(11:12).fixed=1


pa([0, 1,2,3,4,5, 6,7, 8, 14, 15]).fixed=1 ;, 1, 2, 3, 4 ,9,10,11,12,13 ,9,10,11,12,13
;pa(4).fixed=1
;pa(7:8).fixed=1
;pa(1).mpmaxstep=0.0001d0
if (ramp_flag eq 0) then begin
   pa(14).value=0.d0
   pa(14).fixed=1
endif
if (cowan_opt eq 1) then begin
   pa(11:12).value=0.d0
   pa(11:12).fixed=1
endif
if (cowan_opt eq 2) then begin
   pa([9,11]).value=0.d0
   pa([9,11]).fixed=1
endif
if (cowan_opt eq 3) then begin
   pa([9,12]).value=0.d0
   pa([9,12]).fixed=1
endif
if (cowan_opt eq 4) then begin
   pa(9).value=0.d0
   pa(9).fixed=1
endif


fa={phase_func:phase_func}
message, /reset
p = mpfit('calc_trans_sec_phase', parinfo=pa, functargs=fa, perror=perror, status=status,  DOF=adof,bestnorm=bestnorm, ERRMSG=ERRMSG, /fastnorm,/nocatch)
print, 'perror', perror
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

;save, p, perror, xt, bestnorm, filename=sfile

;final_answer, p, perror, data_ch1, nbr_ch1, gw_ch1, data_ch2, nbr_ch2, gw_ch2, E0_ch1, E0_ch2, phase_func, thresh, sfile


end






