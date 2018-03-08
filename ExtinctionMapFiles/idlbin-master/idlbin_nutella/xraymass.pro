pro xraymass
;cluster 1 net counts = 171.6  err = 30.5   counts/s =1.73E-3  L(0.3-7, 5Kev) = 1.722E-14 erg/s/cm^2  T=
;cluster 2                       71.6            28.8			   7.21E-4			     7.178E-15
;cluster 3			    9.28		13.9   
;(had to use a smaller area because of bright source on cluster 3)

counts = [171.6, 71.6]  ;in 0.3-7kev range

exptime = 99340.3

countpersec = counts/exptime 
print, countpersec

Fx = [1.722D-14, 7.178D-15]  ;erg/s/cm^2
z=0.9


dl = lumdist(z)  ;mpc
da = dl/(1+z)^2 ;mpc

dl = dl *1D6*3.08D18  ;cm
da = da *1D6*3.08D18  ;cm

Lx = Fx*4*!PI*dl^2

;T calculated from Lx - T relation of maughan06 and Vikhlinin02
B = 3.32
A = 6.24D44
T = 6*10^[(1/B)*alog10(.456*Lx/A)]
print, T

;newly calculated Fx in 0.5 - 2Kev
Fx = [ 6.944D-15,3.094D-15]
Lx = Fx*4*!PI*dl^2
print, 'final Lx in 0.5-2kev', Lx


;use maughan 07 M-Lx relation
omegam = 0.3
omegal = 0.7
E = (omegam*(1+z)^3 + (1-omegam-omegal)*(1+z)^2 + omegal)^0.5
C=1.6D44
Cerr=0.1D44
alpha=7./3.
B=1.5
Berr=0.1
print, E, C, alpha, B

m500=4.D14*10^[(1./B)*alog10(Lx/(C*E^alpha))]

print, 'mass', m500

;error propagation
dLx = [0.6D43, 0.7D43]
dB=0.1
dC=0.1D44

dmlx = 4.D14*10^[(1./B)*alog10((Lx+dLx)/(C*E^alpha))] - m500
dmb = 4.D14*10^[(1./(B+db))*alog10(Lx/(C*E^alpha))] - m500
dmc = 4.D14*10^[(1./B)*alog10(Lx/((C+dc)*E^alpha))] - m500

dmass = sqrt(dmlx^2 + dmb^2 + dmc^2)
print, dmlx, dmb, dmc 
print, 'dmass',dmass

end
