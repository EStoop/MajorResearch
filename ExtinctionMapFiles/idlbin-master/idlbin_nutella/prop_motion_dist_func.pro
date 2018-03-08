function prop_motion_dist_func, x

;THIS FUNCTION DEFINES THE INTEGRAND IN THE INTEGRAL THAT HAS TO BE COMPUTED
;IN ORDER TO COMPUTE MOST INTERESTING COSMOLOGICAL PARAMETERS.  SEE
;'THE COSMOLOGICAL CONSTANT' BY SEAN CARROLL, EQUATION 42 FOR THE PROPER
;MOTION DISTANCE.

;For most accurate cosmology, uncomment the following two lines
omega_m = 0.3
omega_v = 0.7
H_0 = 70.

;For q = 0.5, uncomment the following two lines (Dressler)
;omega_m = 1.0
;omega_v = 0.0
;H_0 = 50.

H = H_0 * 1000.0 * 100.0 / ((10.0^6.0) * (3.0 * 10.0^18.0)) ;in sec^-1
H_a = H * (omega_m / x^3 + omega_v)^0.5

return, 1.0 / (x^2.0 * H_a)

end
