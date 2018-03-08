pro clustermass
;determine cluster mass from  richness
;sdss hansen et al. method


;G = 6.672E-8 ;(cm^3s^-2g^-1)
G = 6.672D-8*(3.24D-25)^3 * 1.989D33 ;Mpc^3s^-2Msolar^-1
z=1.0
omegam = 0.27
H0 = 72 ;(km/s/Mpc)
H0 = 72*3.24D-20  ;Mpc/s/Mpc
H2 = H0^2*(omegam*(1+z)^3 + 1-omegam)
rhoc = 3*H2/(8*!PI*G)

h=0.72
;-----------
;calculate N200
;number of red galaxies with rest-frame i l>0.4L* in an aperture r

;Ngals = number pf galaxies within 1h Mpc within +- 2sigma of the eso redigeline, brighter than M* + 1
MstarK = -22.9
dm = 44.1  ;distance modulus 
maglimit = mstark +dm +1  ;22.2

;how many galaxies are within 720kpc in the red sequence brighter than irac1=22.2
;take this from cluster.pro
Ngals=[ 7,9,7]
r = 0.156*h*Ngals^0.6  ;Mpc
 print, 'r', r
;N200 is the number of galaxies inside r
N200 = [7,9,7]
r200 = 0.182*h*N200^(0.42) ;  (Mpc)

m200 = 200*rhoc*(4./3.)*!pi*(r200^3)
print, 'sdss', m200

;==================================
;yee et al method
;background corrected number counts
;test values
ncluster =[41,50,40]; [64.,72.,73.]
nback =14.44;48.3
Nnet = ncluster - nback

gamma = 1.8

;angular diameter distance to cluster
;luminosity distance / (1+z)^2
z = 1.0
D = lumdist(z) / (1+z)^2     ;Mpc

;angular size of counting aperture
theta = 1.5    ;Mpc   ;0.013degrees
;theta = 0.013
;angular area of counting aperture
Atheta = !Pi*theta^2

Igamma = 3.78

;integrated cluster LF
kstar = -22.9
alpha = -0.84
phistar = 1.78E-3  ;Mpc^-3

phi = int_lumfunc(-25.,-21.9, phistar, kstar, alpha)
print, 'phi', phi
;===========

Bgc = Nnet*((3-gamma)*(D^(gamma - 3))*(theta^(gamma-1)))/(2*Atheta*Igamma*phi)

deltaBgc = Bgc*(sqrt(Nnet +  1.69*nback)/Nnet)


m200 = 10^(1.62*alog10(Bgc) + 9.86)

print, 'yee', m200, bgc
end
