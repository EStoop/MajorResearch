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

function ellpic_bulirsch,n,k
; Computes the complete elliptical integral of the third kind using
; the algorithm of Bulirsch (1965):
kc=sqrt(1d0-k*k) &  p=n+1d0
if(min(p) lt 0d0) then print,'Negative p'
m0=1d0 & c=1d0 & p=sqrt(p) & d=1d0/p & e=kc
L1:
f = c & c = d/p+f & g = e/p & d = (f*g+temporary(d))*2d0
p = g + temporary(p) & g = m0 & m0 = kc + temporary(m0)
if(max(abs(1d0-kc/g)) gt 1.d-13) then begin
 kc = 2d0*sqrt(e) & e=kc*m0
 goto,L1
endif
return,.5d0*!dpi*(c*m0+d)/(m0*(m0+p))
end

function ellk,k
; Computes polynomial approximation for the complete elliptic
; integral of the first kind (Hasting's approximation):
m1=1.d0-k^2
a0=1.38629436112d0
a1=0.09666344259d0
a2=0.03590092383d0
a3=0.03742563713d0
a4=0.01451196212d0
b0=0.5d0
b1=0.12498593597d0
b2=0.06880248576d0
b3=0.03328355346d0
b4=0.00441787012d0
ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*alog(m1)
return,ek1-ek2
end

function ellec,k
; Computes polynomial approximation for the complete elliptic
; integral of the second kind (Hasting's approximation):
m1=1.d0-k^2
a1=0.44325141463d0
a2=0.06260601220d0
a3=0.04757383546d0
a4=0.01736506451d0
b1=0.24998368310d0
b2=0.09200180037d0
b3=0.04069697526d0
b4=0.00526449639d0
ee1=1.d0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*alog(1.d0/m1)
return,ee1+ee2
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

pro occultquadfast,z0,u1,u2,p,F,F0
;  This routine computes the lightcurve for occultation of a quadratically 
;  limb-darkened source.  Please cite Mandel & Agol (2002) if you make use 
;  of this routine in your research.  
;  Please report errors or bugs to agol@astro.washington.edu
;
; Input:
;
; z0   impact parameter in units of the source radius, rs (vector)
; p    occulting star/planet size in units of rs (scalar)
; (This routine has only been tested for p > 0.003)
; u1   linear    limb-darkening coefficient (gamma_1 in paper)
; u2   quadratic limb-darkening coefficient (gamma_2 in paper)
;
; Output:
;
; F   fraction of flux at each z0 for a limb-darkened source
; F0  fraction of flux at each z0 for a uniform source
;
; Limb darkening has the form:
;  I(r)=[1-u1*(1-sqrt(1-(r/rs)^2))-u2*(1-sqrt(1-(r/rs)^2))^2]/(1-u1/3-u2/6)/pi
; 
; MA02 equation (1): uniform source
occultunifast,z0,p,F0
lambdae=1.d0-F0
nz=n_elements(z0)  & pp=dblarr(nz)+p & lambdad=dblarr(nz)
etad=dblarr(nz) & F=dblarr(nz)
; substitute z=z0(i) to shorten expressions
z=z0 & etad=0.5d0*p^2*(p^2+2.d0*z^2)
; Table 1, Case 1:
indx=where((z ge 1.d0+p) or (pp lt 0.d0))
if(total(indx) ge 0) then etad(indx)=0.d0
occult=((z le 1.d0+p) and (pp gt 0.d0))
; Table 1, Case 11 (corrected typos):
indx=where((pp ge 1.d0) and (z lt p-1.d0))
if(total(indx) ge 0) then etad(indx)=0.5d0
; Table 1, Cases 2 & 8: (These are the most important cases.)
indx=where((z gt abs(1.d0-p)) and (z lt (1.d0+p)) and (z ne p) and occult)
if(total(indx) ge 0) then begin
 zz=z(indx)
 kap1=acos((1.d0-p^2+zz^2)/2.d0/zz)
 kap0=acos((p^2+zz^2-1.d0)/2.d0/p/zz)
; Equation 7 in MA02, lambda_1:
 a=(zz-p)^2 & b=(zz+p)^2 & q=p^2-zz^2
 k=sqrt((1.d0-a)/4.d0/zz/p)
 en=1.d0/a-1.d0
 lambdad(indx)=(((1.d0-b)*(2.d0*b+a-3.d0)-3.d0*q*(b-2.d0))*ellk(k)+$
;      4.d0*p*zz*(zz^2+7.d0*p^2-4.d0)*ellec(k)-3.d0*q/a*ellpic(en,k))/9.d0/!dpi/sqrt(p*zz)
     4.d0*p*zz*(zz^2+7.d0*p^2-4.d0)*ellec(k)-3.d0*q/a*ellpic_bulirsch(en,k))/9.d0/!dpi/sqrt(p*zz)
; Equation 7 in MA02, eta_1:
 etad(indx)=0.5d0*(kap1+2.d0*etad(indx)*kap0-0.25d0*(1.d0+5.d0*p^2+zz^2)*sqrt((1.d0-a)*(b-1.d0)))/!dpi
endif
; Table 1, Cases 3 & 9, MA02:
case3=((pp gt 0.d0) and (pp lt .5d0) and (z gt p) and (z lt 1.d0-p))
case9=((pp gt 0.d0) and (pp lt 1.d0) and (z gt 0.d0) and (z lt .5d0-abs(p-.5d0)))
indx=where((case3 or case9) and occult)
if(total(indx) ge 0) then begin
 zz=z(indx)
 a=(zz-p)^2 & b=(zz+p)^2 & q=p^2-zz^2
 kinv=2.d0*sqrt(p*zz/(1.d0-a))
; Equation 7 in MA02, lambda_2:
 en=b/a-1.d0
 lambdad(indx)=2.d0/9.d0/!dpi/sqrt(1.d0-a)*((1.d0-5.d0*zz^2+p^2+q^2)*$
    ellk(kinv)+(1.d0-a)*(zz^2+7.d0*p^2-4.d0)*ellec(kinv)-3.d0*q/a*$
;     ellpic(en,kinv))
    ellpic_bulirsch(en,kinv))
endif
; Table 1, Case 4 (corrected typo in equation 7 & include p > 0.5d0):
indx=where((pp gt 0.d0) and (pp lt 1.d0) and (z eq (1.d0-p)) and occult)
if(total(indx) ge 0) then lambdad(indx)=2.d0/3.d0/!dpi*acos(1.d0-2.d0*p)-$
  4.d0/9.d0/!dpi*(3.d0+2.d0*p-8.d0*p^2)*sqrt(p*(1.d0-p))-2.d0/3.d0*(p gt 0.5d0)
; Table 1, Case 5:
indx=where((pp gt 0.d0) and (pp lt 0.5d0) and (z eq p) and occult)
; lambda_4, equation 7 (corrected 2k -> 2p):
if(total(indx) ge 0) then lambdad(indx)=1.d0/3.d0+2.d0/9.d0/!dpi*$
    (4.d0*(2.d0*p^2-1)*ellec(2.d0*p)+(1.d0-4.d0*p^2)*ellk(2.d0*p))
; Table 1, Case 6:
indx=where((pp eq 0.5d0) and (z eq 0.5d0) and occult)
if(total(indx) ge 0) then lambdad(indx)=1.d0/3.d0-4.d0/9.d0/!dpi
; Table 1, Case 7:
indx=where((pp gt 0.5d0) and (z eq p) and occult)
if(total(indx) ge 0) then begin
 kap1=acos(0.5d0/p) & kap0=acos(1.d0-0.5d0/p^2)
; lambda_3, equation 7 (corrected typo 1/2k -> 1/2p):
 lambdad(indx)=1.d0/3.d0+16.d0*p/9.d0/!dpi*(2.d0*p^2-1.d0)*$
    ellec(0.5d0/p)-(1.d0-4.d0*p^2)*(3.d0-8.d0*p^2)/9.d0/!dpi/p*$
    ellk(0.5d0/p)
 etad(indx)=0.5d0*(kap1+2.d0*etad(indx)*kap0-0.25d0*(1.d0+6.d0*p^2)*$
    sqrt(4.d0*p^2-1.d0))/!dpi
endif
; Table 1, Case 10:
indx=where((z eq 0.d0) and (p lt 1.d0) and occult)
if(total(indx) ge 0) then lambdad(indx)=-2.d0/3.d0*(1.d0-p^2)^1.5d0
; Now, using equation (33):
omega=1.d0-u1/3.d0-u2/6.d0
F=1.d0-((1.d0-u1-2.d0*u2)*lambdae+(u1+2.d0*u2)*(lambdad+2.d0/3.d0*(p gt z))+u2*etad)/omega
iflag=where((finite(F) eq 0) or (F lt 0.d0))
if(total(iflag) ge 0) then F(iflag)=0.5d0*(F(iflag-1)+F(iflag+1))
return
end

pro occultarb_fn,b,p,c,lc,fname
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
lc=lcuni
indx=where(lcuni ne 1d0 and lcuni ne 0d0, count)
if (count gt 0.0) then begin
   nindx=n_elements(indx)
   lc1=lcuni[indx]
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
endif 
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

phase=x*0.d0

return,phase

end

function cowan_func, x, p

;x is true anomaly, or f+omega, or 
;f+offset

phase=p(0)*cos(x)+p(1)*sin(x)+ $
           p(2)*cos(2.d0*x)+p(3)*sin(2.d0*x)

return,phase

end



function transit_phase_orb_nl,t,x, x_ph, phase_func;, nonlin_c

  COMMON data, bjd_tot, flux_tot, nbr_ind, gw, err_tot, time_tot, nonlin_c

;new routine based on Eric Agol's transit_orb_quicker routine
;Modified by NKL on 12/1/10 for use with HAT-P-2b CH 1&2
;Further modified by NKL on 01/21/12 to limit round-off erro

; Computes a transit lightcurve, normalized to unity, as a
; function of time, t, from transit center
;
; Input parameters (x) are:
; x(0) = P  (units of day)
; x(1) = inc = inclination angle (radians)
; x(2) = a/R_* = semi-major axis divided by R_*
; x(3) = e (eccentricity)
; x(4) = omega (radians)
; x(5) = center of transit
; x(6) = R_p/R_*
; x(7) = depth of secondary eclipse #1
; x(8) = depth of secondary eclipse #2
; x(9) = light travel time across system (seconds)
; 
; x_ph = phase curve parameters

;correct time to light travel

;tp=x(5)
;f_trans=1.5d0*!dpi-x(4)
f_trans=0.5d0*!dpi-x(4)
e1=x(3)
;tanE=sqrt(1.d0-e1^2)*cos(x(4))/(sin(x(4))+e1)
;E_trans=atan(tanE)
E_trans=2.d0*atan(sqrt(1d0-e1^2)*tan(.5d0*f_trans),1d0+e1)
M_trans=E_trans-e1*sin(E_trans)
tp=x(5)-x(0)/(2.d0*!dpi)*M_trans
;tp=x(5)+x(0)*sqrt(1d0-e1^2)/2d0/!dpi*(e1*sin(f_trans)/(1d0+e1*cos(f_trans))$
;   -2d0/sqrt(1d0-e1^2)*atan(sqrt(1d0-e1^2)*tan(.5d0*f_trans),1d0+e1))

m=2d0*!dpi/x(0)*(t-tp)
ekep=ekepler_nl(m,e1)
f=2.d0*atan(sqrt((1.d0+e1)/(1.d0-e1))*tan(0.5d0*ekep))
nm0=where(m eq 0.d0)
if(total(nm0) ge 0) then f(nm0)=0.d0

;phase correction due to light travel time
;;df=2.d0*!dpi*(x(9)/86400.d0/x(0))*(sin(ekep)+1.0)/2.d0
;f=f-df

radius=x(2)*(1d0-x(3)^2)/(1d0+x(3)*cos(f))

; Now compute sky separation:
z0=radius*sqrt(1d0-(sin(x(1))*sin(x(4)+f))^2)

if min(z0) ge 1.d0 then fit=t*0.d0

;get transit curve with limb darkening
occultarb_fn,z0,x(6),nonlin_c,flux_trans,'iofmu'

;use uniform disk version for secondary eclipses
occultunifast,z0,x(6),flux_sec

phase=call_function(phase_func, f, x_ph) 

ind_trans=where(abs(t-x(5)) gt 0.25)
flux_trans(ind_trans)=1.0
ind_sec=where(abs(t-x(5)) lt 0.25)
flux_sec(ind_sec)=1.0

ind_sec1=where(flux_sec lt 1.d0 and t-tp lt 0.0, c1)
ind_sec2=where(flux_sec lt 1.d0 and t-tp gt 0.0, c2)

flux_sec_new=flux_sec
flux_sec_new(ind_sec1)=1.d0+(flux_sec(ind_sec1)-1.d0)*x(7)/x(6)^2
flux_sec_new(ind_sec2)=1.d0+(flux_sec(ind_sec2)-1.d0)*x(8)/x(6)^2

ph1_sec1=phase(ind_sec1(0)-1.0)
ph4_sec1=phase(ind_sec1(c1-1.0)+1.0)
ph0_sec1=(ph1_sec1+ph4_sec1)/2.d0

ph1_sec2=phase(ind_sec2(0)-1.0)
ph4_sec2=phase(ind_sec2(c2-1.0)+1.0)
ph0_sec2=(ph1_sec2+ph4_sec2)/2.d0

trans2=flux_sec_new-1.d0
ph_sec=phase
ph_sec(ind_sec1)=(trans2(ind_sec1)+x(7))*((phase(ind_sec1)-ph0_sec1)/x(7)+1.d0)-(x(7)-ph0_sec1)
ph_sec(ind_sec2)=(trans2(ind_sec2)+x(8))*((phase(ind_sec2)-ph0_sec2)/x(8)+1.d0)-(x(8)-ph0_sec2)
;offset=min([ph_sec(ind_sec1), ph_sec(ind_sec2)])
offset=min(ph_sec(ind_sec2))

ph_sec=ph_sec-offset
ind=where(ph_sec lt 0.d0, count)
if (count gt 0.0) then ph_sec(ind)=0.d0

flux=ph_sec+flux_trans
;flux=(ph_sec+1.d0)*flux_trans

ph_min=min(phase-offset, mind)
if (ph_min) lt 0.d0 then ph_min=0.d0
tph_min=bjd_tot(mind)

ph_max=max(phase-offset, mxind)
tph_max=bjd_tot(mxind)

;;orb_params=[tp, phi_sec, t_sec1, t_sec2, t14_trans, t12_trans, t34_trans, t14_sec1, t12_sec1, t34_sec1, t14_sec2, t12_sec2, t34_sec2]
ph_params=[ph_min, tph_min, ph_max, tph_max]
print,'ph_params', ph_params
print, 'amplitude',( ph_params(2) - ph_params(0))/2.
print, 'phase shift', (x(5) - ph_params(1))*(1./x(0))*360.

return,flux

end

