pro getdist,z,dl,t0

; procedure to find the angular
; size and luminosity distances for
; z particular redshift for (0.3,0.7,75) cosmology
; t0 = 14.44 Gyr, true age is 13.67

h0 = 71.
cons = 3d5/h0
dz = z/10000.
omega_m = 0.27
omega_l = 0.73

tot = 0.
tot2 = 0.
tz = 0.
for ct1=0,10000 do begin
 temp = dz/sqrt(omega_m*(1+tz)^3 + omega_l) 
 tot = tot+cons*temp
 tot2 = tot2+temp/(1.+tz)
 tz = tz+dz
endfor

;print,'da=',tot/(1.+z),' dL=',tot*(1.+z)
dl = tot*(1.+z) ; in Mpc
t0 = tot2*3.0857d22/(h0*1d3*3.15d7) ; in years
end
