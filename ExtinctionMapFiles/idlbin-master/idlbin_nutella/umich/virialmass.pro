PRO virialmass
close,/all

G = 6.672E-8                    ;cm^3 s^-2,g^-1
R = 3.96E6*3.0856E18   ;cm
sigma= 1328.  -59.              ;km/s
M = double((3.*!PI*R*((sigma)^2.)))
div = double(2.*G)
M = M*1E10 /div           ;(km to cm)^2
M = M / 1.989E33          ;grams to solar masses
print, "virial mass for sigma = ", sigma," is", M," solar masses"

END
