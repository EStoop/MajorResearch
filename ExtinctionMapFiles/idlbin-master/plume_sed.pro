
;Get Phot in W/m^2/Hz
cnst = 1.e-3 ;erg/s/cm^2/Hz -> W/m^2/Hz
cnst = 1.e-3*1.d26 ;erg/s/cm^2/Hz -> Jy
dbread,'','plume_phot.dat','DOUBLE',dbphot,db,nl
fnu_fuv = dbphot[0] * cnst
fnu_fuverr = dbphot[1] * cnst
fnu_nuv = dbphot[2] * cnst
fnu_nuverr = dbphot[3] * cnst
fnu_B = dbphot[4] * cnst
fnu_Berr = dbphot[5] * cnst
fnu_V = dbphot[6] * cnst
fnu_Verr = dbphot[7] * cnst

;NUV is UL
nsigul = 3.
fnu_nuvul = nsigul*fnu_nuverr
fnu_nuv = fnu_nuvul
fnu_nuverr = !VALUES.F_NAN
;fnu_nuvul = !VALUES.F_NAN

lam = [1528.,2271.,4400., 5500.]
nu = 2.99e18/lam
f_plume = [fnu_fuv,fnu_nuv,fnu_B,fnu_V] ;Jy
f_plume_err = [fnu_fuverr,fnu_nuverr,fnu_Berr,fnu_Verr] ;Jy
fnu_nuverr = 1d32 ;Make error large for fitting procedure

;Making SED Plot
;---------------------------------------------------------------------------------------
showfont,6,''
xtit=textoidl('\lambda (!6!sA!r!u!9 %!6!n)')
ytit=textoidl('\nu!8f!X_{\nu} (W/m^{2})')


;# Age (yr)      
; 0: 0.000E+00, 1: 1.000E+07, 2: 3.000E+07, 3: 5.000E+07, 4: 1.015E+08
; 5: 2.861E+08, 6: 5.088E+08, 7: 1.015E+09, 8: 3.000E+09  9: 5.000E+09,
;10: 1.000E+10
age=[0.,10.,30.,50.,100.,286.,509.,1015.,3000.,5000.,10000.] ;Myr
specf = 'cb2007_hr_stelib_m62_salp_ssp.1'
dbread,'',specf,'DOUBLE',dbsp,nc
lambda = reform(dbsp[0,*])
nage = nc-1

d2scl = 4.*!pi*5.d25^2 
modnu = 2.99e18/lambda
lam_um = lambda/1.e4 ;microns
flam_mod = dbsp[1:nc-1,*]*3.826e33/d2scl ;ergs/s/cm^2/A
fnu_mod = flam_mod
for i=0, nage-1 do fnu_mod[i,*] = lambda^2*flam_mod[i,*]*(1.e-3*1.e26)/2.99e18 ; 3.3356e4 ;Jy
Av = dindgen(50)/10.
nAv = n_elements(Av)
AlamAv = wd_extcalc(lam_um,tau_lam=tau_lam,Alam=Alam)

nlam = n_elements(lam_um)
vmod = dblarr(nage*nAv)
nuvmod = Vmod
fuvmod = Vmod
Bmod = vmod
COMMON FILTERCURVE,filters
;get_filts

snrfuv = fnu_fuv/fnu_fuverr
snrnuv = fnu_nuv/fnu_nuverr
snrB =  fnu_B/fnu_Berr
snrV = fnu_V/fnu_Verr

for i=0, nage-1 do $
   for j=0, nav-1 do begin
   fuvmod[nav*i + j] = obsv_filter(lam_um,modnu*fnu_mod[i,*]*exp(-AlamAv/1.086*Av[j]),/box,lrange=[.1300,.1700])
   nuvmod[nav*i + j] = obsv_filter(lam_um,modnu*fnu_mod[i,*]*exp(-AlamAv/1.086*Av[j]),/box,lrange=[.1900,.2700])
   Bmod[nav*i + j] = obsv_filter(lam_um,modnu*fnu_mod[i,*]*exp(-AlamAv/1.086*Av[j]),/box,lrange=[.3800,0.4900])
   vmod[nav*i + j] = obsv_filter(lam_um,modnu*fnu_mod[i,*]*exp(-AlamAv/1.086*Av[j]),/box,lrange=[.4950,0.6050])
endfor

snrtot = (snrfuv + snrnuv + snrB + snrV)
wfuv = snrfuv/snrtot
wnuv = snrnuv/snrtot
wB = snrB/snrtot
wV = snrV/snrtot
;fnu_nuv = 3*fnu_nuverr

;wfuv = 0.25
;wnuv = 0.25
;wB = 0.25
;wV = 0.25

tscl = wfuv*fnu_fuv/fuvmod + wnuv*fnu_nuv/nuvmod + wB*fnu_B/Bmod + wV*fnu_V/Vmod 
fdiff = ((fnu_fuv - tscl*fuvmod)/fnu_fuverr)^2 + ((fnu_nuv - tscl*nuvmod)/fnu_nuverr)^2 + $
        ((fnu_B - tscl*Bmod)/fnu_Berr)^2 + ((fnu_V - tscl*Vmod)/fnu_Verr)^2

;Force fit below 3-sigma upper limit of NUV pt
ord = sort(fdiff)
ind_ord = 0.
ind = ord[ind_ord]
while tscl[ind]*nuvmod[ind] GT fnu_nuv AND ind_ord LT nAv*nage do begin
   ind = ord[ind_ord++]
   chisq = fdiff[ind]
endwhile
if ind_ord EQ nAv*nage then chisq = min(fdiff,ind)
;chisq = min(fdiff,ind)
   
rc,ind,nav,r,c
fitfnu = lambda & fitfnu[*] = 0.0
fitfnu = tscl[ind]*fnu_mod[r,*]*exp(-AlamAv/1.086*Av[c])
minAv = Av[c]
minage = age[r]

;GET UV/opt luminosity before/after extinction (difference is IR luminosity)
dl = 16.6 ;Mpc Virgo
;dl = .001 ;Mpc MW
d2fac = 4.*!pi*dl^2*3.086d22^2 ;m^2
fit_llam = 2.99d18/lambda^2*1.d-26*d2fac*(tscl[ind]*fnu_mod[r,*]) ;W/m^2/A
fit_llam_ext = fit_llam*exp(-AlamAv/1.086*Av[c]) ;W/m^2/A
lamUVopt = where(lambda GE 100. and lambda LE 500000.)
L_UVopt =  int_tabulated(lambda[lamUVopt],fit_llam[lamUVopt])/3.826d26 ;Lsun
L_UVopt_ext =  int_tabulated(lambda[lamUVopt],fit_llam_ext[lamUVopt])/3.826d26 ;Lsun
L_IR = L_UVopt-L_UVopt_ext

;Monte Carlo Error Estimation
nmc = 200
seed=[5L,10L,15L,20L]
mcfnu_fuv = randomn(seed[0],nmc,/norm)*fnu_fuverr+fnu_fuv
mcfnu_nuv = randomn(seed[1],nmc,/norm)*fnu_nuverr+fnu_nuv
mcfnu_B = randomn(seed[2],nmc,/norm)*fnu_Berr+fnu_B
mcfnu_V = randomn(seed[3],nmc,/norm)*fnu_Verr+fnu_V

mcminage = dblarr(nmc) 
mcminAv = dblarr(nmc)
   
for mc=0, nmc-1 do begin

   mctscl = wfuv*mcfnu_fuv[mc]/fuvmod + wnuv*mcfnu_nuv[mc]/nuvmod + wB*mcfnu_B[mc]/Bmod + wV*mcfnu_V[mc]/Vmod 
   mcfdiff = ((mcfnu_fuv[mc] - mctscl*fuvmod)/fnu_fuverr)^2 + ((mcfnu_nuv[mc] - mctscl*nuvmod)/fnu_nuverr)^2 + $
             ((mcfnu_B[mc] - mctscl*Bmod)/fnu_Berr)^2 + ((mcfnu_V[mc] - mctscl*Vmod)/fnu_Verr)^2
   
   ;Force fit below 3-sigma upper limit of NUV pt
   ;mcord = sort(mcfdiff)
   ;mcind_ord = 0.
   ;mcind = mcord[mcind_ord]
   ;while mctscl[mcind]*nuvmod[mcind] GT fnu_nuv AND mcind_ord LT nAv*nage do begin
   ;   mcind = ord[mcind_ord++]
   ;   mcchisq = mcfdiff[mcind]
   ;endwhile
   ;if mcind_ord EQ nAv*nage then mcchisq = min(mcfdiff,mcind)
   mcchisq = min(mcfdiff,mcind)

   rc,mcind,nav,mcr,mcc
   mcfitfnu = lambda & mcfitfnu[*] = 0.0
   mcfitfnu = mctscl[mcind]*fnu_mod[mcr,*]*exp(-AlamAv/1.086*Av[mcc])
   mcminAv[mc] = Av[mcc]
   mcminage[mc] = age[mcr]
endfor
Averr = stddev(mcminAv)
ageerr = stddev(mcminage)

set_plot,'ps'
device,/enc,filename='plume-sed.eps'
plot, lam,modnu*fitfnu,/ylog,psym=-1,chars=1.5,xtit=xtit,ytit=ytit,thick=5,syms=2,/xlog,$
      yr=[1.e-16,1.e-13],/nodata,xr=[1.e3,2.e4],/xst,charthick=5,ythick=5,xthick=5
;Blam = 4500.+4500.*[-0.1,-0.1,0.1,0.1]
;polyfill,Blam,[1.e-16,1.e-12,1.e-12,1.e-16],col=0.75*!D.N_COLORS
;irac1lam = 36000.+36000.*[-0.1,-0.1,0.1,0.1]
;polyfill,irac1lam,[1.e-16,1.e-13,1.e-13,1.e-16],col=0.75*!D.N_COLORS
;irac2lam = 45000.+45000.*[-0.1,-0.1,0.1,0.1]
;polyfill,irac2lam,[1.e-16,1.e-13,1.e-13,1.e-16],col=0.75*!D.N_COLORS
 
nsig=3.
plotsym,0,/fill
oplot,lam,nu*f_plume*1.e-26,psym=8,syms=1.
oploterror, lam,nu*f_plume*1.e-26,nsig*nu*f_plume_err*1.e-26,psym=3
plotsym,1
plots,lam[1],nu[1]*fnu_nuvul*1.e-26,psym=8,syms=2.5,thick=3.0
oplot, lambda,modnu*fitfnu*1.e-26,thick=5
;oplot, lambda,modnu*fitfnu_Av0*1.e-29,thick=3,linest=1
;oplot, lambda,modnu*tscl[0]*fnu[0,*]*exp(-AlamAv/1.086*Av)*1.e-29,linest=3,thick=3
;oplot, lambda,modnu*tscl[6]*fnu[6,*]*exp(-AlamAv/1.086*Av)*1.e-29,linest=4,thick=3
;oplot, lambda,modnu*tscl[7]*fnu[7,*]*exp(-AlamAv/1.086*Av)*1.e-29,linest=5,thick=3
agetxt = textoidl('Age = '+strtrim(string(fix(minage)),2)+'\pm'+ $
                  strtrim(string(fix(ageerr)),2)+' Myr')
Avtxt = textoidl('!8A_{V}!X = '+strtrim(string(minAv,format='(D3.1)'),2)+'\pm'+$
                 strtrim(string(Averr,format='(D3.1)'),2)+' mag')
legend,[agetxt,Avtxt],chars=1.25,charthick=4,box=0

device,/close
set_plot,'x'

spawn,'open plume-sed.eps'

;caption:

;The FUV-bright plume fit by a young SED having non-negligible extinction. 
;Assuming a single stellar population, the FUV, NUV, and V-band photometry are best fit by a 100 Myr old SED having an internal extinction of A_V = 1 mag.  The dot-dashed line indicates a zero-age SED while the triple dot-dashed and dashed lines indicate SEDs 500 Myr and 1 Gyr.  The FUV-NUV color is too steep to be fit by a young stellar population alone.  The wavelength coverage of the new broad-band imagery currently being obtained, along with the proposed GALEX bands, are shown.  

END
 


