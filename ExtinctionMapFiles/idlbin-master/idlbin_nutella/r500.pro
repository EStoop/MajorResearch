pro r500

m500=[6.2D13, 3.6D13, 3.6D13]
m500 = 3.6D13
m500err = [1.4D13, 1.4D13, 1.4D13]
print, 'm500', m500
;m500 = 500*rhocrit*4/3 * Pi * r500^3
;r500^3 = (m500 * 3) / (500 * rhocrit* 4 * Pi)

G = 6.672D-8 ;cm^3s^-2g^-1  - km & solar masses
G = G * 1.989D33 ;(g/Msolar)
G = G * (1./1.D5)^3.;(km/cm)^3
G = G * (3.24D-20)^3. ; (Mpc/Km)^3
print, "G,",G

H0 = 70.  ;km s^-1 Mpc^-1
H0 = H0 * 3.24D-20  ; (Mpc/Km)
z = 1.
omegam = 0.27
Hofz = H0^2*[omegam*(1.+z)^3 + (1.-omegam)]
print, "hofz", hofz, hofz^2

rhocrit = (3*Hofz) / (8*!Pi*G)
print, 'rhocrit', rhocrit

r500 = ((m500 * 3.) / (500. * rhocrit* 4. * !Pi)) ^(1./3.)

print, "r500", r500
end
